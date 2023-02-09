!! This is a specialized routine for elliptic problems. For these problems, 
!! we tag based on the right hand side

subroutine fc2d_hps_fort_tag4refinement(mx,my,mbc, & 
    mfields, xlower,ylower,dx,dy,blockno, & 
    rhs, tag_threshold, init_flag, tag_patch)
    IMPLICIT NONE

    INTEGER :: mx,my, mbc, mfields, tag_patch, init_flag
    INTEGER :: blockno
    DOUBLE PRECISION :: xlower, ylower, dx, dy
    DOUBLE PRECISION :: tag_threshold
    DOUBLE PRECISION :: rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    INTEGER :: i,j, mq
    DOUBLE PRECISION :: qmin, qmax

    INTEGER :: exceeds_th, fclaw2d_clawpatch_exceeds_threshold
    INTEGER :: ii,jj
    DOUBLE PRECISION :: xc,yc,quad(-1:1,-1:1), qval
    logical :: is_ghost

    !! # Assume that we won't refine      
    tag_patch = 0

    !! # Default : Refinement based only on first variable in system.  
    !! # Users can modify this by creating a local copy of this routine
    !! # and the corresponding tag4coarsening routine.

    is_ghost = .false.  !! We don't try to tag ghost cells

    mq = 1    !! Tag on the first field
    qmin = rhs(1,1,mq)
    qmax = rhs(1,1,mq)
    do j = 1,my
        do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(qmin,rhs(i,j,mq))
            qmax = max(qmax,rhs(i,j,mq))
            qval = rhs(i,j,mq)
            do jj = -1,1
                do ii = -1,1
                    quad(ii,jj) = rhs(i+ii,j+jj,mq)
                end do
            end do
            exceeds_th = fclaw2d_clawpatch_exceeds_threshold( & 
                        blockno, qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                        tag_threshold, init_flag, is_ghost)
            !! # -1 : Not conclusive (possibly ghost cell); don't tag for refinement
            !! # 0  : Does not pass threshold (don't tag for refinement)      
            !! # 1  : Passes threshold (tag for refinement)
            if (exceeds_th .gt. 0) then
                tag_patch = 1
                return
            endif
        enddo
    enddo

end subroutine fc2d_hps_fort_tag4refinement
