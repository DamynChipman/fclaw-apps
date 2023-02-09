#ifndef FC2D_HPS_FORT_H
#define FC2D_HPS_FORT_H

#include <fclaw_base.h>         /* Needed for FCLAW_F77_FUNC */
#include "fc2d_hps.hpp"           /* Needed for fc2d_hps_fort_eval_bc_t def */

#ifdef __cplusplus
extern "C"
{
#endif


/* - -------------------------------- BC functions ------------------------------------ */

#define FC2D_HPS_FORT_APPLY_BC_DEFAULT FCLAW_F77_FUNC(fc2d_hps_fort_apply_bc_default, \
                                                       FC2D_HPS_FORT_APPLY_BC_DEFAULT)

void FC2D_HPS_FORT_APPLY_BC_DEFAULT(const int* blockno, const  int* mx, const  int* my, 
                                    const  int* mbc, const  int* mfields, 
                                    const double* xlower, const double* ylower,
                                    const double* dx, const double* dy, 
                                    const double* t, 
                                    int intersects_bc[], int mthbc[], 
                                    double rhs[], fc2d_hps_fort_eval_bc_t g_bc, 
                                    int *cons_check, double flux_sum[]);


#define FC2D_HPS_FORT_EVAL_BC_DEFAULT FCLAW_F77_FUNC(fc2d_hps_fort_eval_bc_default, \
                                                FC2D_HPS_FORT_EVAL_BC_DEFAULT)

double FC2D_HPS_FORT_EVAL_BC_DEFAULT(const int* iface, 
                                     const double* t, 
                                     const double* x, const double* y);


/* Specialized to do conservation check as well */
#define FC2D_HPS_FORT_APPLY_BC FCLAW_F77_FUNC(fc2d_hps_fort_apply_bc, \
                                                       FC2D_HPS_FORT_APPLY_BC)

void FC2D_HPS_FORT_APPLY_BC(const int* blockno, const  int* mx, const  int* my, 
                                    const  int* mbc, const  int* mfields, 
                                    const double* xlower, const double* ylower,
                                    const double* dx, const double* dy, 
                                    const double* t, 
                                    int intersects_bc[], int mthbc[], 
                                    double rhs[], fc2d_hps_fort_eval_bc_t g_bc, 
                                    int *cons_check, double flux_sum[]);


#define FC2D_HPS_FORT_EVAL_BC FCLAW_F77_FUNC(fc2d_hps_fort_eval_bc, \
                                                FC2D_HPS_FORT_EVAL_BC)

double FC2D_HPS_FORT_EVAL_BC(const int* iface, 
                                     const double* t, 
                                     const double* x, const double* y);


#ifdef __cplusplus
}
#endif

#endif

