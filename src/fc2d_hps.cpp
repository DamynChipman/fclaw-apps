#include "fc2d_hps.hpp"
// #include <HPS/fc2d_hps.hpp>
#include "fc2d_hps_options.h"
#include "fc2d_hps_physical_bc.h"
#include "fc2d_hps_fort.h"
#include "fc2d_hps_diagnostics.h"
// #include <Methods/fc2d_hps_methods.hpp>
// #include <Structures/fc2d_hps_vector.hpp>

#include <EllipticForestApp.hpp>
#include <HPSAlgorithm.hpp>
#include <FISHPACK.hpp>

static fc2d_hps_vtable_t s_hps_vt;

using HPSAlgorithm = EllipticForest::HPSAlgorithm<EllipticForest::FISHPACK::FISHPACKFVGrid, EllipticForest::FISHPACK::FISHPACKFVSolver, EllipticForest::FISHPACK::FISHPACKPatch, double>;

/* --------------------- Hps solver (required) ------------------------- */

static
void hps_setup_solver(fclaw2d_global_t *glob)
{
	
    // Get EllipticForest app
    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    app.log("Setting up EllipticForest...");

    // Get glob from user and get options
	fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
	fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);

    // Create root patch
    int nx = clawpatch_opt->mx;
    int ny = clawpatch_opt->my;
    double x_lower = fclaw_opt->ax;
    double x_upper = fclaw_opt->bx;
    double y_lower = fclaw_opt->ay;
    double y_upper = fclaw_opt->by;
    EllipticForest::FISHPACK::FISHPACKFVGrid root_grid(nx, ny, x_lower, x_upper, y_lower, y_upper);
    EllipticForest::FISHPACK::FISHPACKPatch root_patch(root_grid);
    root_patch.level = 0;
    root_patch.isLeaf = true;

    // Create PDE problem to solve
    // TODO: Link to ForestClaw stuff
    // EllipticForest::FISHPACK::FISHPACKProblem pde{};
    // pde.setU([&](double x, double y){
    //     return x + y;
    // });
    // pde.setF([&](double x, double y){
    //     return 0.0;
    // });

    // Create patch solver
    EllipticForest::FISHPACK::FISHPACKFVSolver solver{};

    // Create new HPS algorithm
    // TODO: delete HPS in clean up function (where...?)
    HPSAlgorithm* HPS = new HPSAlgorithm(root_patch, solver);

    // Save HPS into ForestClaw glob
    // TODO: Should I put this somewhere else?
    glob->user = (HPSAlgorithm*) HPS;

    // Call setup stage
    fclaw2d_domain_t* domain = glob->domain;
    p4est_wrap_t* p4est_wrap = (p4est_wrap_t*) domain->pp;
    p4est_t* p4est = p4est_wrap->p4est;
    HPS->setupStage(p4est);

}


static
void hps_rhs(fclaw2d_global_t *glob,
             fclaw2d_patch_t *patch,
             int blockno,
             int patchno)
{

    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    app.log("Setting up RHS...");

    int mx,my,mbc;
    double dx,dy,xlower,ylower;
	fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
								&xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;
	fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);
	FCLAW_ASSERT(mfields == 1);

	/* Compute right hand side */
    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
    FCLAW_ASSERT(hps_vt->fort_rhs != NULL); /* Must be initialized */

	hps_vt->fort_rhs(&blockno,&mbc,&mx,&my,&mfields,
                    &xlower,&ylower,&dx,&dy,rhs);
}

static
void hps_solve(fclaw2d_global_t *glob)
{

    // Get EllipticForest app
    EllipticForest::EllipticForestApp& app = EllipticForest::EllipticForestApp::getInstance();
    app.log("Beginning HPS solve...");

    // Get HPS algorithm from glob
    // TODO: Should I get this from somewhere else?
    HPSAlgorithm* HPS = (HPSAlgorithm*) glob->user;

    // Call build stage
    HPS->buildStage();

    // Call upwards stage
    HPS->upwardsStage([&](EllipticForest::FISHPACK::FISHPACKPatch& leafPatch){
        EllipticForest::FISHPACK::FISHPACKFVGrid& grid = leafPatch.grid();
        fclaw2d_patch_t* fc_patch = &(glob->domain->blocks->patches[leafPatch.leafID]);
        int mfields;
        double* rhs;
        fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);
        leafPatch.vectorF() = EllipticForest::Vector<double>(grid.nPointsX() * grid.nPointsY());
        for (auto i = 0; i < grid.nPointsX(); i++) {
            for (auto j = 0; j < grid.nPointsY(); j++) {
                int idx = j + i*grid.nPointsY();
                int idx_T = i + j*grid.nPointsX();
                leafPatch.vectorF()[idx] = rhs[idx_T];
            }
        }
        return;
    });

    // Call solve stage; provide Dirichlet data via function
    HPS->solveStage([&](EllipticForest::FISHPACK::FISHPACKPatch& rootPatch){
        fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
        fc2d_hps_vtable_t* hps_vt = fc2d_hps_vt();

        EllipticForest::FISHPACK::FISHPACKFVGrid& grid = rootPatch.grid();
        rootPatch.vectorG() = EllipticForest::Vector<double>(2*grid.nPointsX() + 2*grid.nPointsY());

        EllipticForest::Vector<double> gWest(grid.nPointsY());
        EllipticForest::Vector<double> gEast(grid.nPointsY());
        EllipticForest::Vector<double> gSouth(grid.nPointsX());
        EllipticForest::Vector<double> gNorth(grid.nPointsX());

        int dirichletBC = 1;

        for (auto j = 0; j < grid.nPointsY(); j++) {
            double y = grid(1, j);
            double x_lower = grid.xLower();
            double x_upper = grid.xUpper();
            gWest[j] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x_lower, &y);
            gEast[j] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x_upper, &y);
        }
        for (auto i = 0; i < grid.nPointsX(); i++) {
            double x = grid(0, i);
            double y_lower = grid.yLower();
            double y_upper = grid.yUpper();
            gSouth[i] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x, &y_lower);
            gNorth[i] = hps_vt->fort_eval_bc(&dirichletBC, &glob->curr_time, &x, &y_upper);
        }

        rootPatch.vectorG().setSegment(0*grid.nPointsX(), gWest);
        rootPatch.vectorG().setSegment(1*grid.nPointsX(), gEast);
        rootPatch.vectorG().setSegment(2*grid.nPointsX(), gSouth);
        rootPatch.vectorG().setSegment(3*grid.nPointsX(), gNorth);

        return;
    });

    // Copy data to ForestClaw patch
    HPS->quadtree.traversePreOrder([&](EllipticForest::FISHPACK::FISHPACKPatch& patch){
        if (patch.isLeaf) {
            fclaw2d_patch_t* fc_patch = &(glob->domain->blocks->patches[patch.leafID]);

            int mbc;
            int Nx, Ny;
            double x_lower, y_lower, dx, dy;
            double* q;
            int meqn, mfields;
            double* rhs;
            fclaw2d_clawpatch_grid_data(glob, fc_patch, &Nx, &Ny, &mbc, &x_lower, &y_lower, &dx, &dy);
            fclaw2d_clawpatch_soln_data(glob, fc_patch, &q, &meqn);
            fclaw2d_clawpatch_rhs_data(glob, fc_patch, &rhs, &mfields);

            EllipticForest::FISHPACK::FISHPACKFVGrid& grid = patch.grid();
            int nx = grid.nPointsX() + 2*mbc;
            int ny = grid.nPointsY() + 2*mbc;
            for (auto i = 0; i < nx; i++) {
                for (auto j = 0; j < ny; j++) {
                    int idx = j + i*ny;
                    int idx_T = i + j*nx;
                    if (i > mbc-1 && i < nx-mbc && j > mbc-1 && j < ny-mbc) {
                        q[idx_T] = patch.vectorU()[idx];
                        rhs[idx_T] = patch.vectorU()[idx];
                    }
                }
            }
        }
        return;
    });

    return;

}


/* ---------------------------------- Output functions -------------------------------- */

static
void hps_output(fclaw2d_global_t *glob, int iframe)
{
	const fc2d_hps_options_t* hps_opt;
	hps_opt = fc2d_hps_get_options(glob);

	// if (hps_opt->ascii_out != 0)
	// 	fclaw2d_clawpatch_output_ascii(glob,iframe);

	// if (hps_opt->vtk_out != 0)
	// 	fclaw2d_clawpatch_output_vtk(glob,iframe);
}


/* ---------------------------------- Tagging functions ------------------------------- */

int hps_tag4refinement(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int blockno, int patchno,
                       int initflag)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);

    int tag_patch;
    double refine_threshold;

    refine_threshold = fclaw_opt->refine_threshold;

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,this_patch,&rhs,&mfields);

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_tag4refinement != NULL);


    /* Use default fortran tagging routines.  Choose refinement based on criteria
       set in configuration files (clawpatch:refinement-criteria) */
    tag_patch = 0;
    clawpatch_vt->fort_tag4refinement(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs, &refine_threshold,
                                      &initflag, &tag_patch);
    return tag_patch;
}


static
int hps_tag4coarsening(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *fine_patches,
                       int blockno,
                       int patchno,
                       int initflag)
{
    fclaw2d_patch_t *patch0 = &fine_patches[0];
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs[4];
    int mfields;
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_rhs_data(glob,&fine_patches[igrid],&rhs[igrid],&mfields);
    }

    fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_tag4coarsening != NULL);

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    double coarsen_threshold = fclaw_opt->coarsen_threshold;

    int tag_patch = 0;
    clawpatch_vt->fort_tag4coarsening(&mx,&my,&mbc,&mfields,&xlower,&ylower,&dx,&dy,
                                      &blockno, rhs[0],rhs[1],rhs[2],rhs[3],
                                      &coarsen_threshold,&initflag,&tag_patch);
    return tag_patch == 1;
}

/* -------------------------------- Diagnostic functions ------------------------------ */
static
void hps_compute_error(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *patch,
                          int blockno,
                          int patchno,
                          void *user)
{
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;

    fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
    {
        fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
        FCLAW_ASSERT(clawpatch_vt->fort_compute_patch_error != NULL);

        int mx, my, mbc;
        double xlower,ylower,dx,dy;
        fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,&xlower,
                                    &ylower,&dx,&dy);

        double *area = fclaw2d_clawpatch_get_area(glob,patch);  /* Might be null */

        /* Computing solution is stored in the RHS; true solution is stored in soln */
        int mfields;

        /* Computed solution */
        double *rhs;
        fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

        double *err;
        fclaw2d_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);

        /* True solution */
        double *soln;
        fclaw2d_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);
        double t = glob->curr_time;
        clawpatch_vt->fort_compute_patch_error(&blockno, &mx,&my,&mbc,
                                               &mfields,&dx,&dy,
                                               &xlower,&ylower, &t, rhs, err, soln);
        /* Accumulate sums and maximums needed to compute error norms */

        FCLAW_ASSERT(clawpatch_vt->fort_compute_error_norm != NULL);
        clawpatch_vt->fort_compute_error_norm(&blockno, &mx, &my, &mbc, &mfields, 
                                              &dx,&dy, area, err,
                                              error_data->local_error);

    }
}


static
void hps_conservation_check(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *patch,
                            int blockno,
                            int patchno,
                            void *user)
{
    fc2d_hps_error_info_t* error_data = (fc2d_hps_error_info_t*) user;
    int mx, my, mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mfields;
    double *rhs;  /* Solution is stored in the right hand side */ 
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    FCLAW_ASSERT(clawpatch_vt->fort_conservation_check != NULL);


    /* Need a better way to determine which diagnostic to do */
    double* area = fclaw2d_clawpatch_get_area(glob,patch);  
    clawpatch_vt->fort_conservation_check(&mx, &my, &mbc, &mfields, &dx,&dy,
                                          area, rhs, error_data->rhs,
                                          error_data->c_kahan);
    fc2d_hps_options_t *hps_opt = fc2d_hps_get_options(glob);

    int intersects_bc[4];
    fclaw2d_physical_get_bc(glob,blockno,patchno,intersects_bc);

    double t = glob->curr_time;
    int cons_check = 1;

    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt();
    FCLAW_ASSERT(hps_vt->fort_apply_bc != NULL);

    /* Sum up the normal derivative around the boundary */
    hps_vt->fort_apply_bc(&blockno, &mx, &my, &mbc, &mfields, 
                         &xlower, &ylower, &dx,&dy,&t, intersects_bc,
                         NULL,rhs, hps_vt->fort_eval_bc,
                         &cons_check, error_data->boundary);
}


/* ------------------------------ Virtual functions  ---------------------------------- */

static
fc2d_hps_vtable_t* hps_vt_init()
{
	FCLAW_ASSERT(s_hps_vt.is_set == 0);
	return &s_hps_vt;
}

void fc2d_hps_solver_initialize(fclaw2d_global_t* glob)
{
	int claw_version = 4; /* solution data is organized as (i,j,m) */
	fclaw2d_clawpatch_vtable_initialize(glob, claw_version);

	/* Patch : These could be over-written by user specific settings */
	fclaw2d_patch_vtable_t*   patch_vt = fclaw2d_patch_vt(glob);  
	patch_vt->rhs            = hps_rhs;   /* Calls FORTRAN routine */
    patch_vt->initialize     = hps_rhs;   /* Get an initial refinement */
	patch_vt->setup          = NULL;

    /* Tagging functions : Base refinement on the right hand side */
    patch_vt->tag4refinement = hps_tag4refinement;
    patch_vt->tag4coarsening = hps_tag4coarsening;

    /* Clawpatch and ForestClaw : Output functions */
    fclaw2d_vtable_t*   fclaw_vt = fclaw2d_vt(glob);
    fclaw_vt->output_frame = hps_output;

    /* Elliptic specific functions */
    fclaw2d_elliptic_vtable_t *elliptic_vt = fclaw2d_elliptic_vt(glob);
    elliptic_vt->setup = hps_setup_solver;

    /* Solver doesn't do anything so far */
    elliptic_vt->solve = hps_solve;    
    elliptic_vt->apply_bc = fc2d_hps_physical_bc;

    /* BCs : Homogeneous BCs by default */
	fc2d_hps_vtable_t*  hps_vt = hps_vt_init();	
    hps_vt->fort_apply_bc = &FC2D_HPS_FORT_APPLY_BC_DEFAULT;
    hps_vt->fort_eval_bc  = &FC2D_HPS_FORT_EVAL_BC_DEFAULT;

    /* Diagnostics : Error, conservation */
    fclaw2d_clawpatch_vtable_t *clawpatch_vt = fclaw2d_clawpatch_vt(glob);
    clawpatch_vt->compute_error = hps_compute_error;  /* calls user-defined fortran routine */

    /* Conservation check : Compares sum(rhs) with sum of normal fluxes around the boundary
       of the solution.   (uses divergence theorem) */
    clawpatch_vt->conservation_check = hps_conservation_check;        

    /* These are specialized for the elliptic problem */
    fclaw2d_diagnostics_vtable_t *diag_vt = fclaw2d_diagnostics_vt(glob);
    diag_vt->patch_init_diagnostics     = fc2d_hps_diagnostics_initialize;
    diag_vt->patch_reset_diagnostics    = fc2d_hps_diagnostics_reset;
    diag_vt->patch_compute_diagnostics  = fc2d_hps_diagnostics_compute;
    diag_vt->patch_gather_diagnostics   = fc2d_hps_diagnostics_gather;
    diag_vt->patch_finalize_diagnostics = fc2d_hps_diagnostics_finalize;

	hps_vt->is_set = 1;
}


/* ----------------------------- User access to solver functions --------------------------- */

fc2d_hps_vtable_t* fc2d_hps_vt()
{
	FCLAW_ASSERT(s_hps_vt.is_set != 0);
	return &s_hps_vt;
}





