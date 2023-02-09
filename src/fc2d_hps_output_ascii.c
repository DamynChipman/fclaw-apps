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


#include "fc2d_hps_output_ascii.h"

void fc2d_hps_time_header_ascii(fclaw2d_global_t* glob, int iframe)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = 
    fclaw2d_clawpatch_get_options(glob);
    char matname1[20];
    sprintf(matname1,"fort.q%04d",iframe);

    FILE *f1 = fopen(matname1,"w");
    fclose(f1);

    char matname2[20];
    sprintf(matname2,"fort.t%04d",iframe);

    double time = glob->curr_time;

    int ngrids = glob->domain->global_num_patches;

    int mfields = clawpatch_opt->rhs_fields;  
    int maux = clawpatch_opt->maux;

    int mf;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    if (fclaw_opt->compute_error)
        mf = mfields + 2;  /* Print out error and true solution */
    else
        mf = mfields;  /* Only print out computed solution */


    FILE *f2 = fopen(matname2,"w");
    fprintf(f2,"%12.6f %23s\n%5d %30s\n%5d %30s\n%5d %30s\n%5d %30s\n",time,"time",
            mf,"mfields",ngrids,"ngrids",maux,"num_aux",2,"num_dim");
    fclose(f2);

}


void cb_hps_output_ascii(fclaw2d_domain_t * domain,
                         fclaw2d_patch_t * patch,
                         int blockno, int patchno,
                         void *user)
    {
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;
    int iframe = *((int *) s->user);

    /* Get info not readily available to user */
    int global_num, local_num;
    int level;
    fclaw2d_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num,&local_num, &level);
    
    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *rhs;
    int mfields;
    fclaw2d_clawpatch_rhs_data(glob,patch,&rhs,&mfields);

    double *err;
    fclaw2d_clawpatch_elliptic_error_data(glob,patch,&err,&mfields);

    double *soln;
    fclaw2d_clawpatch_elliptic_soln_data(glob,patch,&soln,&mfields);

    char fname[BUFSIZ];
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    snprintf (fname, BUFSIZ, "%s.q%04d", fclaw_opt->prefix, iframe);

    fc2d_hps_vtable_t*  hps_vt = fc2d_hps_vt(); 
    FCLAW_ASSERT(hps_vt->fort_output != NULL);

    if (fclaw_opt->compute_error)
    {        
        hps_vt->fort_output(fname,&mx,&my,&mfields,&mbc,
                            &xlower,&ylower,&dx,&dy,rhs,
                            soln, err, &global_num, &level,&blockno,
                            &glob->mpirank);
    }
    else
    {
        fclaw2d_clawpatch_vtable_t* clawpatch_vt = fclaw2d_clawpatch_vt(glob);
        clawpatch_vt->fort_output_ascii(fname,&mx,&my,&mfields,&mbc,
                                        &xlower,&ylower,&dx,&dy,rhs,
                                        &global_num,&level,&blockno,
                                        &glob->mpirank);        
    }

}
