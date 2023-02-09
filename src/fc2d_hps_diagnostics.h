#ifndef FC2D_HPS_DIAGNOSTICS_H
#define FC2D_HPS_DIAGNOSTICS_H

#include <fclaw2d_include_all.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_diagnostics.h>

#include "fc2d_hps.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct {
    double* local_error;  /* meqn x 3 array of errors on a patch */
    double* global_error; /* meqn x 3 array of errors after gather */
    double *mass0;  /* Mass at initial time */
    double *mass;
    double area;
    double *rhs;       /* Sum of rhs hand side */
    double *boundary;  /* sum around boundary */
    double *c_kahan;  
} fc2d_hps_error_info_t;

/* --------------------------- Problem dependent functions -----------------------------*/

void fc2d_hps_diagnostics_initialize(fclaw2d_global_t *glob, void **acc_patch);


void fc2d_hps_diagnostics_reset(fclaw2d_global_t *glob, void* patch_acc);

void fc2d_hps_diagnostics_compute(fclaw2d_global_t* glob,
                                           void* patch_acc);

void fc2d_hps_diagnostics_gather(fclaw2d_global_t *glob, void* patch_acc,
                               int init_flag);

void fc2d_hps_diagnostics_finalize(fclaw2d_global_t *glob, void** patch_acc);

void fc2d_hps_compute_diagnostics(fclaw2d_domain_t *domain,
                                fclaw2d_patch_t *patch,
                                int blockno,
                                int patchno,
                                void* user);


#ifdef __cplusplus
}
#endif

#endif
