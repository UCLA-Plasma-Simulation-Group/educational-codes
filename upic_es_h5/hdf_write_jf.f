      module hdf_write_jf
      
      use globals, only: LINEAR,QUADRATIC
      use pinit2d_jf
! FST
      use hdf5io_class  ! for HDF5 I/O from Weiming
! FST

      implicit none
      private
      public :: writef, write_jf
      public :: EX,EY,JX,JY,DEN,VDOTEX,VDOTEY,PHSL_XX,PHSL_XY,PHSL_YX,PHSL_YY
      public :: ESPOYNT_X,ESPOYNT_Y,FIN_VS_INIT_ENE,FIN_VS_INIT_VX,FIN_VS_INIT_VY
      public :: VDOTEX_PART,VDOTEY_PART,VDOTEX_INT,VDOTEY_INT,FVX_XY_LABEL,FVY_XY_LABEL
      public :: VDOTEX_FOLLOW_PART,VDOTEY_FOLLOW_PART,ESPOYNT_INT_X,ESPOYNT_INT_Y
      public :: EX_ENE_INT,EY_ENE_INT,DIV_ESPOYNTINT,VXVY,PH_VX_X,PH_VY_X,PH_VX_Y,PH_VY_Y
      public :: DIVESPOYNT,BIN_EX_LABEL,BIN_EY_LABEL,POT,EZ,BX,BY,BZ
      PUBLIC :: THROUGH_VX_VX,THROUGH_VX_VY,THROUGH_VY_VX,THROUGH_VY_VY,DEXDT,DEYDT
      public :: GRAD_PHI_X,GRAD_PHI_Y,KE_LABEL,THROUGH_DENE
		
		save
      integer,parameter :: EX = 1,EY = 2,JX =3,JY = 4,DEN = 5,VDOTEX = 6,VDOTEY = 7
      integer,parameter :: PHSL_XX = 8,PHSL_XY=9,PHSL_YX=10,PHSL_YY=11
      integer,parameter :: ESPOYNT_X=12,ESPOYNT_Y=13,FIN_VS_INIT_ENE=14,FIN_VS_INIT_VX=15
      integer,parameter :: FIN_VS_INIT_VY=16,VDOTEX_PART=17,VDOTEY_PART=18
      integer,parameter :: VDOTEX_INT=19,VDOTEY_INT=20,FVX_XY_LABEL=21,FVY_XY_LABEL=22
      integer,parameter :: VDOTEX_FOLLOW_PART = 23,VDOTEY_FOLLOW_PART=24
      integer,parameter :: ESPOYNT_INT_X=25,ESPOYNT_INT_Y=26,EX_ENE_INT=27,EY_ENE_INT=28
      integer,parameter :: DIV_ESPOYNTINT=29,VXVY=30,PH_VX_X=31,PH_VY_X=32,PH_VX_Y=33,PH_VY_Y=34
      integer,parameter :: DIVESPOYNT=35,BIN_EX_LABEL=36,BIN_EY_LABEL=37,POT=38,EZ=39
      integer,parameter :: BX=40,BY=41,BZ=42
      integer,parameter :: THROUGH_VX_VX=43,THROUGH_VX_VY=44,THROUGH_VY_VX=45,THROUGH_VY_VY=46
      integer,parameter :: DEXDT=47,DEYDT=48,GRAD_PHI_X=49,GRAD_PHI_Y=50,KE_LABEL=51,THROUGH_DENE=52
      interface writef
!         module procedure ipwrite2_hdf
	 module procedure ipwrite2_h5
!         module procedure ipwrite2
      end interface

     interface write_jf
 !			module procedure write_jf2d_hdf
 !			module procedure write_jf3d_hdf
 			module procedure write_1core_2d_h5
 			module procedure write_1core_3d_h5
     end interface
      
      contains

! FST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_1core_2d_h5(f,nx,ny,fname,label_code,time,it,dvx,dvy,xoff,yoff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine writes real data f to a file
! f = input data to be written, modified on node 0
! nx/ny = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! name = file name (used only if nren < 0)
! dvx and dvy are for changing the space between points on the x and y axis, ie, 
!		dvx = xaxis(2)-xaxis(1)
! xoff and yoff are the positions of the xaxis(1) and yaxis(1)
! input: f, nx, ny, fname
      implicit none
      integer nx,ny,it
      real f
      character*(*) fname
      dimension f(nx,ny,1)
      real time, dvx,dvy,xoff,yoff
      integer label_code
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld

! lstat = length of status array
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      character(len=6) :: nlabel
      character(len=100) :: name
! HDF stuff
!
!     Function declaration.
!
! 
!**** Variable declaration *******************************************
!

! file handles for HDF
! array informations needed for HDF
      integer dims(2)
      integer offset
      integer ledge(2)
      
      character*70 title
      character*70 axistitles(3),axisunits(3)
      character*70 timeunits
      
! arrays needed for HDF-write      
      real*8 xaxis(nx),yaxis(ny)
      real*8 ftemp(nx,ny)
      integer ix,iy
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double

! local variables for HDF5
      real axismin(3),axismax(3)
      type (hdf5file) h5_header
      
! DEBUG
      write(*,*)'in pwrite_1core_h5'
      write(*,*)'(pwrite_h5) nproc= ',nproc
      write(*,*)'(pwrite_h5) nx = ',nx
      write(*,*)'(pwrite_h5) ny = ',ny
      write(*,*)'(pwrite_h5) xoff,yoff =',xoff,yoff
      write(*,*)'(pwrite_h5) dvx, dvy =', dvx, dvy
! DEBUG
! define HDF5 "chunk" here
!     start(1)=0
!     start(2)=0
!     stride(1)=1
!     stride(2)=1
!     ledge(1)=nx
!     ledge(2)=ny
!
      dims(1)= nx 
      dims(2)= ny 
      offset = 0
      
	call do_labels(label_code,title,axistitles(1),axisunits(1),axistitles(2),axisunits(2),&
        & axistitles(3),axisunits(3))
	! write (nlabel,'(i6)') it + 100000
	write (nlabel,'(i6.6)') it 
	name = trim(fname)//'_'//trim(adjustl(nlabel))//'.h5'
	timeunits = '\omega_p^{-1}'
	axismin(1) = xoff
	axismin(2) = yoff
	axismax(1) = xoff + (nx-1) * dvx
	axismax(2) = yoff + (ny-1) * dvy


	 		
	call h5_header%new(filename=name, timeunits=timeunits,ty='grid', &
	& n=it, t=time, dt = 1.0, axisname=axistitles,axislabel=axisunits, &
	& axismin=axismin,axismax=axismax, dataname =  title, units = axisunits(3), &
	& label = axistitles(3), rank=2)
	
	do i = 1,nx
	    do j = 1,ny
	        ftemp(i,j) = f(i,j,1)
            end do
	end do
	call pwfield_serial(h5_header,ftemp(:,:), dims, ierr)

	if(ierr .ne. 0) then
	    write(*,*) 'HDF5 write failed'
	    stop
	endif
!

      return
      end subroutine write_1core_2d_h5

!
! FST
! modified by JF to write a single array without collecting from other nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_1core_3d_h5(f,nx,ny,nz,fname,label_code,time,it,&
      & dvx,dvy,dvz,xoff,yoff,zoff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine writes real data f to a file
! f = input data to be written, modified on node 0
! nx/ny = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! iunit = fortran unit number
! nrec = current record number for write (if negative, open file with
! recl=-nren)
! name = file name (used only if nren < 0)
! dvx and dvy are for changing the space between points on the x and y axis, ie, 
!		dvx = xaxis(2)-xaxis(1)
! xoff and yoff are the positions of the xaxis(1) and yaxis(1)
! input: f, nx, ny, fname
      implicit none
      integer nx,ny,nz,it
      real f
      character*(*) fname
      dimension f(nx,ny,nz)
      real time, dvx,dvy,dvz,xoff,yoff,zoff
      integer label_code
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      character(len=6) :: nlabel
      character(len=100) :: name

! 
!**** Variable declaration *******************************************
!

! array informations needed for HDF
      integer start(3), edges(3), stride(3), dims(3)
      integer ledge(3)
      integer offset(2)
      
      character*70 title
      character*70 axistitles(3),axisunits(3)
      character*70 timeunits
      
! arrays needed for HDF-write      
      integer ix,iy,iz
      real axismin(3),axismax(3)
      type (hdf5file) h5_header
      
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double
      
! DEBUG
!     write(*,*)'in pwrite_hdf'
!     write(*,*)'(pwrite_hdf) nproc= ',nproc
!     write(*,*)'(pwrite_hdf) nx = ',nx
!     write(*,*)'(pwrite_hdf) ny = ',ny
!     write(*,*)'(pwrite_hdf)'
! DEBUG
! define HDF "chunk" here
      start(1)=0
      start(2)=0
      start(3)=0
      stride(1)=1
      stride(2)=1
      stride(3)=1
      ledge(1)=nx
      ledge(2)=ny
      ledge(3)=nz
! for another implemenation, not used
!      gedge(1)=nx
!      gedge(2)=ny*nproc
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=(ledge(2)-start(2))/stride(2)
      dims(3)=(ledge(3)-start(3))/stride(3)

      offset(1) = 0
      offset(2) = 0 

      
      ! write (nlabel,'(i6)') it + 100000
        write (nlabel,'(i6.6)') it 
        name = trim(fname)//'_'//trim(adjustl(nlabel))//'.h5'
	call do_labels(label_code,title,axistitles(1),axisunits(1),axistitles(2),axisunits(2),&
        & axistitles(3),axisunits(3))
	timeunits = '\omega_p^{-1}'
	axismin(1) = xoff
	axismin(2) = yoff
	axismin(3) = zoff
	axismax(1) = xoff + (nx-1) * dvx
	axismax(2) = yoff + (ny-1) * dvy
	axismax(3) = zoff + (nz-1) * dvz


	 		
	call h5_header%new(filename=name, timeunits=timeunits,ty='grid', &
	& n=it, t=time, dt = 1.0, axisname=axistitles,axislabel=axisunits, &
	& axismin=axismin,axismax=axismax, dataname =  title, units = 'a.u.', &
	& label = title , rank=3)
	call pwfield(h5_header,f(:,:,:), dims, dims, offset,ierr)

	if(ierr .ne. 0) then
	    write(*,*) 'HDF5 write failed'
	    stop
	endif


      return
      end subroutine write_1core_3d_h5

! FST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine PWRITE2_H5(f,nx,nxv,kyp,name,label_code,iter,time,iorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine collects distributed real data f and writes to a file
! f = input data to be written, modified on node 0
! nx/kyp = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! name = file name (used only if nren < 0)
! input: f, nx, kyp, nxv, fname
! output: nrec
      implicit none
      integer nx,kyp, nxv,iter
      real f
      character*(*) name
      dimension f(nxv,kyp,1)
      real time
      integer iorder, label_code
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k,l,temp
      integer ierr,ny, nstart
      integer per
      dimension per(2)
      

! 
!**** Variable declaration *******************************************
!

! array informations needed for HDF5
      integer start(2), edges(2), stride(2), dims(2), dims_local(2)
      integer offset
      integer ledge(2)
      type(hdf5file) h5_header
      
      character*70 title
      character*70 xtitle,xunits
      character*70 ytitle,yunits
      character*70 ztitle,zunits
      character*70 axistitles(3), axisunits(3)
      character*70 timeunits
      
! arrays needed for HDF-write      
!      real*4 xaxis(nx),yaxis((kyp-3)*nproc)
!      real*4 ftemp(nx,kyp-3)
!      real*4 ftemprcv(nxv,kyp)
      real,dimension(:), allocatable :: xaxis,yaxis
      real,dimension(:,:), allocatable :: ftemp,ftemprcv
			
      integer ix,iy
      real axismin(3),axismax(3)
!      
! temporary attributes, needed for Ricardo's IDL routines
!
      real*8 time_double
      
      allocate(xaxis(nx),yaxis((kyp-3)*nproc))
      allocate(ftemp(nx,kyp-3))
!     allocate(ftemprcv(nxv,kyp))
      

!
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)

! DEBUG
!     write(*,*)'in pwrite_h5'
!     write(*,*)'(pwrite_hdf) nproc,idproc= ',nvp,idproc
!     write(*,*)'(pwrite_hdf) nx = ',nx
!     write(*,*)'(pwrite_hdf) nxv = ',nxv
!     write(*,*)'(pwrite_hdf) kyp = ',kyp,ny
!     write(*,*)'(pwrite_hdf) iorder = ',iorder
! DEBUG
      ny = kyp-3
! define HDF "chunk" here
      start(1)=0
      start(2)=0
      stride(1)=1
      stride(2)=1
      ledge(1)=nx
      ledge(2)=ny
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=((ledge(2)-start(2))/stride(2))*nproc
      dims_local(1) = dims(1)
      dims_local(2) = ((ledge(2) - start(2))/stride(2))

      axismin(1) = 0
      axismin(2) = 0
      axismin(3) = 0
      axismax(1) = real(nx)
      axismax(2) = real(ny)
      axismax(3) = 0
      timeunits = '\omega_p^{-1}'


            call do_labels(label_code,title,axistitles(1), axisunits(1),&
 	    & axistitles(2), axisunits(2), axistitles(3), axisunits(3))
  	    call h5_header%new(filename = name,timeunits=timeunits,&
            & ty='grid',n=iter,t=time,&
            & dt = 1.0,&
  	    & axisname = axistitles, axislabel = axistitles, &
            & axisunits = axisunits,axismin = axismin,axismax = axismax, &
            & dataname = title,units =axisunits(3),&
            & label = axistitles(3),rank = 2)
!           call h5_header%new(filename=name,rank=2)

            if (iorder .eq. QUADRATIC) then
            do ix=1,nx
                do iy=1,ny
                    ftemp(ix,iy)=f(ix+1,iy+1,1)
                enddo
            enddo            
            else
	    ! FST --> this does not seem possible, since ftemp does not go to kyp
            do ix=1,nx
                do iy=1,ny
                    ftemp(ix,iy)=f(ix,iy,1)
                enddo
            enddo                     
            endif 

! new HDF5 write
            offset = idproc*ny
! DEBUG
!           write(*,*) 'before h5 call'
!           write(*,*) 'global dims=',dims
!           write(*,*) 'local dims=',dims_local
!           write(*,*) 'offset=',offset
! DEBUG
            call pwfield(h5_header,ftemp(:,:),dims,&
            & dims_local,offset,ierr)
	    if (ierr .ne. 0) then
	        write(*,*)'HDF5 write failed'
		stop
            endif
         
         
      
      deallocate(xaxis,yaxis)
      deallocate(ftemp)
!     deallocate(ftemprcv)

      
      return
      end subroutine PWRITE2_H5
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HDF5
! this is the wrapper function for pwrite2_hdf
! HDF5
!! END OF PWRITE_H5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine ipwrite2_h5(f,nxv,nypmx,it,time,label_code,name,inorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 2d data f 
! and it writes to a file.  Only some of the data can be written.
! name = file name (used only if new file is to be opened)
! it = integer to append to file name
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
         implicit none
         real, dimension(:,:,:) :: f
         integer :: nxv, nypmx
         character(len=*), optional :: name
         integer :: it, label_code
         integer, optional :: inorder
         real :: time
! local data
         integer :: mx, order,i,temp
         character(len=6) :: nlabel
         character(len=60) :: fname
         character(len=10), save :: sname = ':ipwrite2:'

! unpack arguments (not necessary here since BEPS feeds it)
!         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
!	  nblok = 1	!This is only appropriate for when there is no shared memory
!         mx = 2*this%nd1; kyp = this%nd2p; nxv = size(f,1)
         fname = ' '
         if (present(name)) then
            fname = name
!            if (present(it)) then
				!write (nlabel,'(i6)') it + 100000
				write (nlabel,'(i6.6)') it 
				fname = trim(name)//'_'//trim(adjustl(nlabel))//'.h5'
!            endif
         endif
         order = QUADRATIC
         if (present(inorder)) then
            if ((inorder >= 1) .and. (inorder <= 2)) then
               order = inorder
            endif
         endif

         if (order==QUADRATIC) then
 !           if ((mx+2) <= nxv) mx = mx + 1
            call PWRITE2_H5(f,nxv-4,nxv,nypmx,trim(fname),label_code,it,time,QUADRATIC)
         else
 !           if ((mx+1) <= nxv) mx = mx + 1
            call PWRITE2_H5(f,nxv-1,nxv,nypmx,trim(fname),label_code,it,time,LINEAR)
         endif
         end subroutine ipwrite2_h5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!takes the label_code and generates titles  and units for all the different diagnostics
      subroutine do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
      	implicit none
      	integer :: label_code
      	character(len=*) :: title,xtitle,xunits,ytitle,yunits,ztitle,zunits
      	
      	select case(label_code)
      	
      	case(EX)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(EY)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(EZ)
      		title = 'Ez'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dz!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(BX)
      		title = 'Bx'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'B!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(BY)
      		title = 'By'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'B!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(BZ)
      		title = 'Bz'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'B!Dz!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(JX)
      		title = 'jx'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dx!N'
      		zunits = 'env!Dth!N'
      	case(JY)
      		title = 'jy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dy!N'
      		zunits = 'env!Dth!N'
      	case(VDOTEX)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(DEN)
      		title = 'Electron Density'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Density'
      		zunits = 'n!D0!N'
      	case(PHSL_XX)
      		title = 'Phase Space Slice - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_XY)
      		title = 'Phase Space Slice - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YX)
      		title = 'Phase Space Slice - vy vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YY)
      		title = 'Phase Space Slice - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(ESPOYNT_X)
      		title = 'ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_Y)
      		title = 'ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(FIN_VS_INIT_ENE)
      		title = 'Initial vs. Final Energy'
      		xtitle = 'Initial Ene'
      		xunits = 'mv!Dth!N!U2!N'
      		ytitle = 'Final Ene'
      		yunits = 'mv!Dth!N!U2!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VX)
      		title = 'Initial vs. Final vx'
      		xtitle = 'Initial v!Dx!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dx!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VY)
      		title = 'Initial vs. Final vy'
      		xtitle = 'Initial v!Dy!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dy!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_PART)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_PART)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEX_INT)
      		title = 'Time Integrated vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_INT)
      		title = 'Time Integrated vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(FVX_XY_LABEL)
      		title = 'Phase Space vx vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N'
      		zunits = '1/v!Dth!N'
      	case(FVY_XY_LABEL)
      		title = 'Phase Space vy vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N'
      		zunits = '1/v!Dth!N'
      	case(VXVY)
      		title = 'Phase Space vx vs vy'
      		xtitle = 'v!Dx!N'
      		xunits = '1/v!Dth!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_FOLLOW_PART)
      		title = 'Summed by following particle vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_FOLLOW_PART)
      		title = 'Summed by following particle vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(ESPOYNT_INT_X)
      		title = 'Time Integrated ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_INT_Y)
      		title = 'Time Integrated ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(EX_ENE_INT)
      		title = 'x Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ex!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(EY_ENE_INT)
      		title = 'y Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ey!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(DIV_ESPOYNTINT)
      		title = 'Time Integrated div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(DIVESPOYNT)
      		title = 'div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(PH_VX_X)
      		title = 'Phase Space - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_X)
      		title = 'Phase Space - vy vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VX_Y)
      		title = 'Phase Space - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_Y)
      		title = 'Phase Space - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(BIN_EX_LABEL)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(BIN_EY_LABEL)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(POT)
      		title = 'Potential'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Potential'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N!U2!N'
      	case(THROUGH_VX_VX)
      		title = 'Initial vx vs. Final vx'
      		xtitle = 'Initial v!Dx!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dx!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(THROUGH_VX_VY)
      		title = 'Initial vx vs. Final vy'
      		xtitle = 'Initial v!Dx!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dy!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(THROUGH_VY_VX)
      		title = 'Initial vy vs. Final vx'
      		xtitle = 'Initial v!Dy!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dx!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(THROUGH_VY_VY)
      		title = 'Initial vy vs. Final vy'
      		xtitle = 'Initial v!Dy!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dy!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(THROUGH_DENE)
      		title = 'Initial vy vs. Final Energy in Wave Frame'
      		xtitle = 'Initial Ene'
      		xunits = 'mv!Dth!U2!N'
      		ytitle = 'Final Ene'
      		yunits = 'mv!Dth!U2!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(DEXDT)
      		title = 'dEx / dt'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'dE!Dx!N / dt'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(DEYDT)
      		title = 'dEy / dt'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'dE!Dy!N / dt'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(GRAD_PHI_X)
      		title = 'grad phi x comp'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'grad phi!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(GRAD_PHI_Y)
      		title = 'grad phi y comp'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'grad phi!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(KE_LABEL)
      		title = 'Kinetic Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Kin. Ene.'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'

      	end select
      end subroutine do_labels


		end module hdf_write_jf
