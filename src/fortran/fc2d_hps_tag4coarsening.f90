subroutine fc2d_hps_fort_tag4coarsening(mx,my,mbc,mfields, & 
           xlower,ylower,dx,dy, blockno, rhs0, rhs1, rhs2, rhs3, & 
           coarsen_threshold, initflag, tag_patch)
    implicit none

    integer :: mx,my, mbc, mfields, tag_patch, initflag
    integer :: blockno
    double precision :: xlower(0:3), ylower(0:3), dx, dy
    double precision :: coarsen_threshold
    double precision :: rhs0(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision :: rhs1(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision :: rhs2(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision :: rhs3(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    integer :: mq
    double precision :: qmin, qmax

    !! # Don't coarsen when initializing the mesh
    if (initflag .ne. 0) then
        tag_patch = 0
        return
    endif

    !! # Assume that we will coarsen a family unless we find a grid
    !! # that doesn't pass the coarsening test.
    tag_patch = 1
    mq = 1
    qmin = rhs0(1,1,mq)
    qmax = rhs0(1,1,mq)

    call fc2d_hps_test_refine(blockno,mx,my,mbc,mfields, & 
            mq,rhs0,qmin,qmax, dx,dy,xlower(0), ylower(0), & 
            coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fc2d_hps_test_refine(blockno,mx,my,mbc,mfields, & 
           mq,rhs1,qmin,qmax,dx,dy,xlower(1), ylower(1), & 
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fc2d_hps_test_refine(blockno,mx,my,mbc,mfields, & 
           mq,rhs2,qmin,qmax,dx,dy,xlower(2), ylower(2), & 
           coarsen_threshold,initflag, tag_patch)
    if (tag_patch == 0) return

    call fc2d_hps_test_refine(blockno,mx,my,mbc,mfields, & 
           mq,rhs3,qmin,qmax,dx,dy,xlower(3), ylower(3), & 
           coarsen_threshold,initflag, tag_patch)

end subroutine fc2d_hps_fort_tag4coarsening

subroutine fc2d_hps_test_refine(blockno,mx,my,mbc, & 
    mfields,mq,rhs, qmin,qmax,dx,dy,xlower,ylower, & 
    coarsen_threshold,init_flag,tag_patch)
    
    implicit none
    integer :: mx,my,mbc,mfields,mq,tag_patch, init_flag, blockno
    double precision :: coarsen_threshold
    double precision :: qmin,qmax, dx, dy, xlower, ylower
    double precision :: rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    double precision :: xc,yc,quad(-1:1,-1:1),qval

    integer i,j, ii, jj

    integer :: exceeds_th, fclaw2d_clawpatch_exceeds_threshold
    logical(kind=4) :: is_ghost

    is_ghost = .false.
    do i = 1,mx
        do j = 1,my
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(rhs(i,j,mq),qmin)
            qmax = max(rhs(i,j,mq),qmax)
            qval = rhs(i,j,mq)
            do ii = -1,1               
                do jj = -1,1
                    quad(ii,jj) = rhs(i+ii,j+jj,mq)
                end do
            end do
            exceeds_th = fclaw2d_clawpatch_exceeds_threshold( & 
                  blockno, qval,qmin,qmax,quad, dx,dy,xc,yc, & 
                  coarsen_threshold, init_flag, is_ghost)
            
            !! # -1 : Not conclusive (possibly ghost cell) (do not tag for coarsening)
            !! # 0  : Does not pass threshold (tag for coarsening)      
            !! # 1  : Passes threshold (do not tag for coarsening)
            !! # Note : exceeds_th = -1 leads to over-refining, so it is 
            !! # ignored here.  Logic of regridding (coarsening then 
            !! # refining) isn't clear.
            if (exceeds_th .gt. 0) then
               tag_patch = 0
               return
            endif
        enddo
    enddo

end subroutine fc2d_hps_test_refine
