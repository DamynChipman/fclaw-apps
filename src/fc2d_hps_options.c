/*
Copyright (c) 2019-2021 Carsten Burstedde, Donna Calhoun, Damyn Chipman
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "fc2d_hps_options.h"

static int s_hps_options_package_id = -1;

static void*
hps_register (fc2d_hps_options_t* hps_opt, sc_options_t * opt)
{

    /* Array of NumFaces=4 values */
    fclaw_options_add_int_array (opt, 0, "boundary_conditions", 
                                 &hps_opt->bc_cond_string, "1 1 1 1",
                                 &hps_opt->boundary_conditions, 4,
                                 "[hps] Physical boundary condition type [1 1 1 1]");

    sc_options_add_bool (opt, 0, "ascii-out", &hps_opt->ascii_out, 0,
                           "Output ASCII formatted data [F]");

    sc_options_add_bool (opt, 0, "vtk-out", &hps_opt->vtk_out, 0,
                           "Output VTK formatted data [F]");

    // sc_options_add_bool(opt, 0, "mmio_out", &hps_opt->mmio_out, 0,
    //                         "Output matrices in Matrix Market format [F]");

    sc_options_add_bool(opt, 0, "cache_T", &hps_opt->cache_T, 1,
                            "Cache DtN matrix T by level [T]");

    sc_options_add_bool(opt, 0, "time_setup", &hps_opt->time_setup, 0,
                            "Time setup stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_build", &hps_opt->time_build, 0,
                            "Time build stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_upwards", &hps_opt->time_upwards, 0,
                            "Time upwards stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "time_solve", &hps_opt->time_solve, 0,
                            "Time solve stage of HPS method [F]");

    // sc_options_add_bool(opt, 0, "time_copy", &hps_opt->time_copy, 0,
    //                         "Time copy stage of HPS method [F]");

    sc_options_add_bool(opt, 0, "nonhomogeneous_rhs", &hps_opt->nonhomogeneous_rhs, 0,
                            "Flag for non-homogeneous RHS; determines if upwards pass is done [F]");

    sc_options_add_bool(opt, 0, "only_patch_solver", &hps_opt->only_patch_solver, 0,
                            "Option to bypass HPS method and just use the patch solver [F]");


    /* Set operator type (laplace, varpoisson, heat, ...) */
    sc_keyvalue_t *kv_op = hps_opt->kv_operator_type = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv_op, "laplace",  LAPLACE);     /* Uses FFT or BICG */
    sc_keyvalue_set_int (kv_op, "varpoisson", VARPOISSON);   /* Uses BICG */
    sc_keyvalue_set_int (kv_op, "heat",       HEAT);   /* Uses BICG */
    sc_keyvalue_set_int (kv_op, "user_solver",  USER_SOLVER);   /* Uses BICG */
    sc_options_add_keyvalue (opt, 0, "operator-type", &hps_opt->operator_type,
                             "laplace", kv_op, "Set operator type [laplace]");

    /* Set solver type (FFT, DST, BICG, ...) */
    sc_keyvalue_t *kv_s = hps_opt->kv_patch_solver = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv_s, "bicg", BICG);
    sc_keyvalue_set_int (kv_s, "fishpack", FISHPACK);
    sc_keyvalue_set_int (kv_s, "dst", DST);
    sc_keyvalue_set_int (kv_s, "fft",  FFT);     
    sc_keyvalue_set_int (kv_s, "user_solver",  USER_SOLVER);     
    sc_options_add_keyvalue (opt, 0, "patch-solver", &hps_opt->patch_solver,
                             "fishpack", kv_s, "Set patch solver [fishpack]");

    hps_opt->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
hps_postprocess (fc2d_hps_options_t * hps_opt)
{
    fclaw_options_convert_int_array (hps_opt->bc_cond_string, 
                                     &hps_opt->boundary_conditions,4);
    
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
hps_check(fc2d_hps_options_t *hps_opt,
                 fclaw2d_clawpatch_options_t *clawpatch_opt)
{
    return FCLAW_NOEXIT;
}

static
void hps_destroy (fc2d_hps_options_t * hps_opt)
{
    fclaw_options_destroy_array (hps_opt->boundary_conditions);

    FCLAW_ASSERT (hps_opt->kv_patch_solver != NULL);
    sc_keyvalue_destroy (hps_opt->kv_patch_solver);

    FCLAW_ASSERT (hps_opt->kv_operator_type != NULL);
    sc_keyvalue_destroy (hps_opt->kv_operator_type);

}

/* ------------------------------------------------------
   This is boiler plate below here .... don't change 

   Generic calls to options handling;  each calls 
   package specific options call back
   ------------------------------------------------------ */

static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);

    hps_opt = (fc2d_hps_options_t*) package;

    return hps_register(hps_opt,opt);
}


static fclaw_exit_type_t
options_postprocess (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    return hps_postprocess (hps_opt);
}


static fclaw_exit_type_t
options_check (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_hps_options_t *hps_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    clawpatch_opt = (fclaw2d_clawpatch_options_t *)
        fclaw_app_get_attribute(app,"clawpatch",NULL);
    FCLAW_ASSERT(clawpatch_opt->is_registered);

    return hps_check(hps_opt,clawpatch_opt);    
}

static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    hps_opt = (fc2d_hps_options_t*) package;
    FCLAW_ASSERT (hps_opt->is_registered);

    hps_destroy (hps_opt);

    FCLAW_FREE (hps_opt);
}

static const fclaw_app_options_vtable_t hps_options_vtable = {
    options_register,
    options_postprocess,
    options_check,
    options_destroy,
};

/* ----------------------------------------------------------
   Public interface to clawpack options
   ---------------------------------------------------------- */
fc2d_hps_options_t*  fc2d_hps_options_register (fclaw_app_t * app,
                                                              const char *configfile)
{
    fc2d_hps_options_t *hps_opt;

    FCLAW_ASSERT (app != NULL);

    hps_opt = FCLAW_ALLOC (fc2d_hps_options_t, 1);
    fclaw_app_options_register (app, "hps", configfile,
                                &hps_options_vtable, hps_opt);
    
    fclaw_app_set_attribute(app,"hps",hps_opt);
    return hps_opt;
}

fc2d_hps_options_t* fc2d_hps_get_options(fclaw2d_global_t *glob)
{
    int id = s_hps_options_package_id;
    return (fc2d_hps_options_t*) fclaw_package_get_options(glob,id);
}

void fc2d_hps_options_store (fclaw2d_global_t* glob, fc2d_hps_options_t* hps_opt)
{
    int id = fclaw_package_container_add_pkg(glob,hps_opt);
    s_hps_options_package_id = id;
}
