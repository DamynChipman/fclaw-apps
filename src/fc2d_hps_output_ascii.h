#ifndef FC2D_HPS_OUTPUT_H
#define FC2D_HPS_OUTPUT_H

#include "fc2d_hps.hpp"
#include "fc2d_hps_options.h"

#ifdef __cplusplus
extern "C"
{
#endif

void fc2d_hps_time_header_ascii(fclaw2d_global_t* glob, int iframe);

void cb_hps_output_ascii(fclaw2d_domain_t * domain,
                         fclaw2d_patch_t * patch,
                         int blockno, int patchno,
                         void *user);


#ifdef __cplusplus
}
#endif


#endif /* FC2D_HPS_OUTPUT_H */
