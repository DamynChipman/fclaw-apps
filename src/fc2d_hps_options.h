#ifndef FC2D_HPS_OPTIONS_H
#define FC2D_HPS_OPTIONS_H

#include <fclaw_base.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fclaw2d_patch.h>
#include <fclaw_package.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

struct fclaw2d_global;

typedef struct fc2d_hps_options fc2d_hps_options_t;


typedef enum {
    LAPLACE = 0,      /* Laplacian (no beta) */    
    VARPOISSON,         /* Variable Poisson operator */
    HEAT,               /* Variable Poisson operator */
    USER_OPERATOR
} fc2d_hps_operator_types;

typedef enum {
    FFT = 0,    /* Must use fivepoint */
    FISHPACK,   /* Must use fivepoint */
    DST,        /* Must use fivepoint */
    BICG,       /* Can be used with variable coefficient */
    USER_SOLVER
} fc2d_hps_solver_types;


struct fc2d_hps_options
{
    /* Boundary conditions */
    int *boundary_conditions;
    const char *bc_cond_string;

    /* Output */
    int ascii_out;
    int vtk_out;
    // int mmio_out;

    /* Caching */
    int cache_T;

    /* Timing */
    int time_setup;
    int time_build;
    int time_upwards;
    int time_solve;
    // int time_copy;

    /* Homogeneous vs. non-homogeneous */
    int nonhomogeneous_rhs;

    // Patch solver short-circuit
    int only_patch_solver;

    int operator_type;  /* laplace, varpoisson, etc ... */
    sc_keyvalue_t *kv_operator_type;

    int patch_solver;                /* FFT, DST, FISHPACK, ... */
    sc_keyvalue_t *kv_patch_solver;

    int is_registered;
};

#if 0
fclaw_exit_type_t fc2d_hps_postprocess (fc2d_hps_options_t *
                                               mg_opt);

fclaw_exit_type_t fc2d_hps_check (fc2d_hps_options_t * mg_opt);

void fc2d_hps_reset (fc2d_hps_options_t * mg_opt);
#endif

fc2d_hps_options_t*  fc2d_hps_options_register (fclaw_app_t * app,
                                                const char *configfile);

void fc2d_hps_package_register(fclaw_app_t* app,
                               fc2d_hps_options_t* mg_opt);

fc2d_hps_options_t* fc2d_hps_get_options(struct fclaw2d_global *glob);

void fc2d_hps_options_store (struct fclaw2d_global* glob, 
                             fc2d_hps_options_t* mg_opt);

#ifdef __cplusplus
}
#endif

#endif
