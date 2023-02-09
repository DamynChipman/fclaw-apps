#ifndef FC2D_HPS_PHYSICAL_BC_H
#define FC2D_HPS_PHYSICAL_BC_H

#include <fclaw2d_elliptic_solver.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_physical_bc.h>

#include "fc2d_hps.hpp"
#include "fc2d_hps_options.h"

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_global;
struct fclaw2d_domain;
struct fclaw2d_patch;


typedef struct fc2d_hps_time_info
{
    double t;
} fc2d_hps_time_info_t;


/* Inhomogeneous boundary conditions */
void fc2d_hps_physical_bc(struct fclaw2d_global *glob);


#ifdef __cplusplus
}
#endif

#endif
