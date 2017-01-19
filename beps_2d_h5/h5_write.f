      module hdf_write_jf
      
      use globals, only: LINEAR,QUADRATIC
      use pinit2d_jf
      use hdf5
      implicit none
      private
      public :: writef, write_jf
      public :: EX,EY,JX,JY,DEN,VDOTEX,VDOTEY,PHSL_XX,PHSL_XY,PHSL_YX,PHSL_YY
      public :: ESPOYNT_X,ESPOYNT_Y,FIN_VS_INIT_ENE,FIN_VS_INIT_VX,FIN_VS_INIT_VY
      public :: VDOTEX_PART,VDOTEY_PART,VDOTEX_INT,VDOTEY_INT,FVX_XY_LABEL,FVY_XY_LABEL
      public :: VDOTEX_FOLLOW_PART,VDOTEY_FOLLOW_PART,ESPOYNT_INT_X,ESPOYNT_INT_Y
      public :: EX_ENE_INT,EY_ENE_INT,DIV_ESPOYNTINT,VXVY,PH_VX_X,PH_VY_X,PH_VX_Y,PH_VY_Y
      public :: DIVESPOYNT,BIN_EX_LABEL,BIN_EY_LABEL,POT
		
		save
      integer,parameter :: EX = 1,EY = 2,JX =3,JY = 4,DEN = 5,VDOTEX = 6,VDOTEY = 7
      integer,parameter :: PHSL_XX = 8,PHSL_XY=9,PHSL_YX=10,PHSL_YY=11
      integer,parameter :: ESPOYNT_X=12,ESPOYNT_Y=13,FIN_VS_INIT_ENE=14,FIN_VS_INIT_VX=15
      integer,parameter :: FIN_VS_INIT_VY=16,VDOTEX_PART=17,VDOTEY_PART=18
      integer,parameter :: VDOTEX_INT=19,VDOTEY_INT=20,FVX_XY_LABEL=21,FVY_XY_LABEL=22
      integer,parameter :: VDOTEX_FOLLOW_PART = 23,VDOTEY_FOLLOW_PART=24
      integer,parameter :: ESPOYNT_INT_X=25,ESPOYNT_INT_Y=26,EX_ENE_INT=27,EY_ENE_INT=28
      integer,parameter :: DIV_ESPOYNTINT=29,VXVY=30,PH_VX_X=31,PH_VY_X=32,PH_VX_Y=33,PH_VY_Y=34
      integer,parameter :: DIVESPOYNT=35,BIN_EX_LABEL=36,BIN_EY_LABEL=37,POT=38
      interface writef
         module procedure ipwrite2_hdf
!         module procedure ipwrite2
      end interface

		interface write_jf
			module procedure write_jf2d_hdf
			module procedure write_jf3d_hdf
		end interface
      
			interface add_h5_attribute
				module procedure add_h5_attribute_str
				module procedure add_h5_attribute_str_v1
				module procedure add_h5_attribute_logical
				module procedure add_h5_attribute_logical_v1
				module procedure add_h5_attribute_int
			
				module procedure add_h5_attribute_double
				module procedure add_h5_attribute_v1_double
			
			end interface

      
      contains
!
! FST
! modified by JF to write a single array without collecting from other nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_jf2d_hdf(f,nx,ny,fname,label_code,time,it,dvx,dvy,xoff,yoff)
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
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      dimension istatus(lstat)
      character(len=5) :: nlabel
      character(len=100) :: name
! HDF stuff
!
!     Function declaration.
!
      integer sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      integer sfcreate
      external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      external sfcreate
      integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
      integer DFNT_FLOAT64
      parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
      parameter(DFNT_FLOAT32 = 5)
      parameter(DFNT_FLOAT64 = 6)
      integer SD_UNLIMITED
      parameter(SD_UNLIMITED = 0)

! 
!**** Variable declaration *******************************************
!

! file handles for HDF
      integer sd_id, sds_id, sds_index, status
      integer dimid1,dimid2,dimid3
! array informations needed for HDF
      integer start(2), edges(2), stride(2), dims(2)
      integer ledge(2)
      
      character*40 title
      character*40 xtitle,xunits
      character*40 ytitle,yunits
      character*40 ztitle,zunits
      
! arrays needed for HDF-write      
      real*8 xaxis(nx),yaxis(ny)
      real*8 ftemp(nx,ny)
      integer ix,iy
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
      stride(1)=1
      stride(2)=1
      ledge(1)=nx
      ledge(2)=ny
! for another implemenation, not used
!      gedge(1)=nx
!      gedge(2)=ny*nproc
!
      dims(1)=(ledge(1)-start(1))/stride(1)
      dims(2)=(ledge(2)-start(2))/stride(2)
      
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
			call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
			write (nlabel,'(i5)') it + 10000
			name = trim(fname)//'_'//trim(adjustl(nlabel))//'.hdf'

      if (idproc.eq.0) then
         sd_id = sfstart(trim(name), DFACC_CREATE)
! new hdf write
         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 2, dims)
! HDF
         status = sfwdata(sds_id,start,stride,ledge,f)
         ! now write HDF header needed for Ricardo's package
         time_double=time
         status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
         status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,it)

         do ix=1,nx
             xaxis(ix)=(ix-1)*dvx + xoff
         enddo
         dimid1 = sfdimid(sds_id,0)
         status = sfsdmname(dimid1, 'x1')
         status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
	 			 status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
         do iy=1,(ny)
             yaxis(iy)=(iy-1)*dvy + yoff
         enddo
         
         if ((label_code .eq. PHSL_XX) .OR. (label_code .eq. PHSL_XY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvxmax / (real(ny)/2.)
				 endif
         if ((label_code .eq. PHSL_YX) .OR. (label_code .eq. PHSL_YY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvymax / (real(ny)/2.)
         endif
         
         dimid2 = sfdimid(sds_id,1)
         status = sfsdmname(dimid2, 'x2')
         status = sfsdscale(dimid2, (ny), DFNT_FLOAT64, yaxis)
	 			 status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')
	 		
				 status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
	 		
! HDF equivalent of the "close" statement
         status=sfendacc(sds_id)
         status=sfend(sd_id)
!

      endif
      return
      end subroutine write_jf2d_hdf

!
! FST
! modified by JF to write a single array without collecting from other nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
      subroutine write_jf3d_hdf(f,nx,ny,nz,fname,label_code,time,it,&
      	&dvx,dvy,dvz,xoff,yoff,zoff)
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
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k
      integer ierr
      dimension istatus(lstat)
      character(len=5) :: nlabel
      character(len=100) :: name
! HDF stuff
!
!     Function declaration.
!
      integer sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      integer sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      integer sfcreate
      external sfstart, sfselect, sfwdata, sfendacc, sfend, sfsdtstr
      external sfdimid,sfsdmname,sfsdscale,sfsdmstr, sfsnatt
      external sfcreate
      integer DFACC_CREATE,DFNT_INT32,DFNT_FLOAT32
      integer DFNT_FLOAT64
      parameter(DFACC_CREATE = 4, DFNT_INT32 = 24)
      parameter(DFNT_FLOAT32 = 5)
      parameter(DFNT_FLOAT64 = 6)
      integer SD_UNLIMITED
      parameter(SD_UNLIMITED = 0)

! 
!**** Variable declaration *******************************************
!

! file handles for HDF
      integer sd_id, sds_id, sds_index, status
      integer dimid1,dimid2,dimid3
! array informations needed for HDF
      integer start(3), edges(3), stride(3), dims(3)
      integer ledge(3)
      
      character*40 title
      character*40 xtitle,xunits
      character*40 ytitle,yunits
      character*40 ztitle,zunits
      
! arrays needed for HDF-write      
      real*4 xaxis(nx),yaxis(ny),zaxis(nz)
      real*4 ftemp(nx,ny)
      integer ix,iy,iz
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
      
! determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
			call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
			write (nlabel,'(i5)') it + 10000
			name = trim(fname)//'_'//trim(adjustl(nlabel))//'.hdf'

      if (idproc.eq.0) then
         sd_id = sfstart(trim(name), DFACC_CREATE)
! new hdf write
!         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT32, 3, dims)
         sds_id = sfcreate(sd_id,trim(title),DFNT_FLOAT64, 3, dims)
! HDF
         status = sfwdata(sds_id,start,stride,ledge,f)
         ! now write HDF header needed for Ricardo's package
         time_double=time
         status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
         status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,it)

         do ix=1,nx
             xaxis(ix)=(ix-1)*dvx + xoff
         enddo
         dimid1 = sfdimid(sds_id,0)
         status = sfsdmname(dimid1, 'x1')
         status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
	 			 status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
         do iy=1,(ny)
             yaxis(iy)=(iy-1)*dvy + yoff
         enddo
         do iz=1,(nz)
             zaxis(iz)=(iz-1)*dvz + zoff
         enddo
         
         if ((label_code .eq. PHSL_XX) .OR. (label_code .eq. PHSL_XY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvxmax / (real(ny)/2.)
				 endif
         if ((label_code .eq. PHSL_YX) .OR. (label_code .eq. PHSL_YY)) then
         	yaxis = (yaxis - real(ny)/2.) * fvymax / (real(ny)/2.)
         endif
         
         dimid2 = sfdimid(sds_id,1)
         status = sfsdmname(dimid2, 'x2')
         status = sfsdscale(dimid2, (ny), DFNT_FLOAT64, yaxis)
	 			 status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')

         dimid3 = sfdimid(sds_id,2)
         status = sfsdmname(dimid3, 'x3')
         status = sfsdscale(dimid3, (nz), DFNT_FLOAT64, zaxis)
	 			 status = sfsdmstr(dimid3, trim(ztitle),trim(zunits),'F5.2')
	 		
!				 status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
	 		
! HDF equivalent of the "close" statement
         status=sfendacc(sds_id)
         status=sfend(sd_id)
!

      endif
      return
      end subroutine write_jf3d_hdf

      subroutine PWRITE2_HDF(f,nx,nxv,kyp,name,label_code,iter,time,iorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
! this subroutine collects distributed real data f and writes to a file
! f = input data to be written, modified on node 0
! nx/kyp = length of data f in x/y on each processor to write
! nxv = first dimension of data array f, must be >= nx
! name = file name (used only if nren < 0)
! input: f, nx, kyp, nxv, fname
      implicit none
      integer nx,kyp, nxv,iter
      real f
      character*(*) name
      dimension f(nxv,kyp,1)
      real time
      integer iorder, label_code
! common block for parallel processing
!      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
!      parameter(lstat=10)
! nproc = number of real or virtual processors obtained
! lgrp = current communicator
! lworld = MPI_COMM_WORLD communicator
!      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
      integer istatus, lrec, nvp, idproc, np, ioff, id, nrec0, i, j, k,l,temp
      integer nproc
      integer ierr,ny
      integer per
      dimension per(2)

			integer, parameter :: izero = ichar('0')

! 
!**** Variable declaration *******************************************
!

      real,dimension(2) :: xaxis_range, yaxis_range
  integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID, fileID

! array informations needed for HDF
 !     integer start(2), edges(2), stride(2)!, dims(2)
		  integer(hsize_t), dimension(1) :: dims

!      integer ledge(2)
      
      character*70 title
      character*70 xtitle,xunits
      character*70 ytitle,yunits
      character*70 ztitle,zunits
      
!			real,dimension(:), allocatable :: xaxis,yaxis
			real,dimension(:,:), allocatable :: ftemp!,ftemprcv
			
      integer ix,iy

      real*8 time_double
      
!      allocate(xaxis(nx),yaxis((kyp-3)*nproc))
      allocate(ftemp(nx,kyp-3))
!      allocate(ftemprcv(nxv,kyp))
      
			nproc = 1
			ny = kyp-3
			! define HDF "chunk" here
!      start(1)=0
!      start(2)=0
!      stride(1)=1
!      stride(2)=1
!      ledge(1)=nx
!      ledge(2)=ny

!      dims(1)=(ledge(1)-start(1))/stride(1)
!      dims(2)=((ledge(2)-start(2))/stride(2))*nproc

			call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)

      xaxis_range(1) = 0; xaxis_range(2) = nx - 1
      yaxis_range(1) = 0; yaxis_range(2) = ny - 1
      
			call create_h5_dataset(fileID, datasetID, rootID, dataspaceID, plistID, name, nx,ny, &
				& label_code,iter,time, xaxis_range,yaxis_range)

! no special diagnostic node
!			 if (nvp.eq.nproc) then
!					np = nvp
!					ioff = 1
! special diagnostic node present
!			 else
!					np = nvp - nproc
!					ioff = 0
!					id = 1
!					call MPI_RECV(ftemprcv,nx*kyp,mreal,id,99,lworld,istatus,ierr)
!			 endif
! first write data for node 0
!			 nrec0 = nrec
!			 start(2)=(nrec-1)*(ny)
			if (iorder .eq. QUADRATIC) then
				do ix=1,nx
					do iy=1,ny
						ftemp(ix,iy)=f(ix+1,iy+1,1)
					enddo
				enddo            
			else
				do ix=1,nx
					do iy=1,kyp
						ftemp(ix,iy)=f(ix,iy,1)
					enddo
				enddo                     
			endif 

			call add_h5_dataset_2d( fileID, 'data', ftemp, units = zunits, long_name = ztitle)

      call close_hdf5_file( fileID )

! then write data from remaining nodes
!			 do i = 2, np
!					id = i - ioff
! calculate the starting point for             
!					start(2)=(nrec-1)*(ny)
!					call MPI_RECV(ftemprcv,nxv*kyp,mreal,id,99,lworld,istatus,ierr)
!					do ix=1,nx
!						do iy=1,ny
!							ftemp(ix,iy) = ftemprcv(ix+1,iy+1)
!						enddo
!					enddo
! old binary write            
!           write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,kyp)
! 
! new hdf write
!					status = sfwdata(sds_id,start,stride,ledge,ftemp)
!           write(*,*)'status (sfwdata_) = ', status
!					nrec = nrec + 1
!			 enddo
!			 per(1) = 1
!			 per(2) = nvp
			 
			 ! now write HDF header needed for Ricardo's package
!			 time_double=time
!			 status = sfsnatt(sds_id,'TIME',DFNT_FLOAT64,1,time_double)
!			 status = sfsnatt(sds_id,'ITER',DFNT_INT32,1,iter)

!			 do ix=1,nx
!					 xaxis(ix)=(ix-1)
!			 enddo
!			 dimid1 = sfdimid(sds_id,0)
!			 status = sfsdmname(dimid1, 'x1')
!			 status = sfsdscale(dimid1, nx, DFNT_FLOAT64, xaxis)
!			status = sfsdmstr(dimid1, trim(xtitle),trim(xunits),'F5.2')         
!			 do iy=1,(ny)*nproc
!					 yaxis(iy)=(iy-1)
!			 enddo
!			 dimid2 = sfdimid(sds_id,1)
!			 status = sfsdmname(dimid2, 'x2')
!			 status = sfsdscale(dimid2, (ny)*nproc, DFNT_FLOAT64, yaxis)
!			status = sfsdmstr(dimid2, trim(ytitle),trim(yunits),'F5.2')
		
!			status = sfsdtstr(sds_id,trim(ztitle),trim(zunits),'F5.2','cartesian')
		
! HDF equivalent of the "close" statement
!			 status=sfendacc(sds_id)
!			 status=sfend(sd_id)
!

      
!      deallocate(xaxis,yaxis)
      deallocate(ftemp)
 !     deallocate(ftemprcv)

      
      return
      end subroutine pwrite2_hdf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine ipwrite2_hdf(this,f,iunit,nrec,name,it,time,inorder)
         subroutine ipwrite2_hdf(f,nxv,nypmx,it,time,label_code,name,inorder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 2d data f 
! and it writes to a file.  Only some of the data can be written.
! this = ufield2d descriptor of data
! f = input data to be written
! name = file name (used only if new file is to be opened)
! it = integer to append to file name
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
         implicit none
!         type (ufield2d), intent(in) :: this
         real, dimension(:,:,:) :: f
         integer :: nxv, nypmx
         character(len=*), optional :: name
         integer :: it, label_code
         integer, optional :: inorder
         real :: time
! local data
         integer :: mx, order,nblok,i,temp
         character(len=5) :: nlabel
         character(len=60) :: fname
         character(len=10), save :: sname = ':ipwrite2:'
! check for errors
!         if (monitor > 0) then
!         if (monitor==2) call werrfl(class//sname//' started')
!         if (checkinit(sname,this) /= 0) return
!         if (this%layout /= XLOCAL) then
!            erstr = ' invalid layout'
!            UFIELD2D_ERR = 2; EXCEPTION = EXCEPTION + 1
!            call ehandler(EDEFAULT,class//sname//erstr); return
!         endif
!         endif
! unpack arguments
!         nxv = size(f,1); nypmx = size(f,2); nblok = size(f,3)
			nblok = 1	!This is only appropriate for when there is no shared memory
!         mx = 2*this%nd1; kyp = this%nd2p; nxv = size(f,1)
         fname = ' '
         if (present(name)) then
            fname = name
!            if (present(it)) then
				write (nlabel,'(i5)') it + 10000
!				fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
				fname = trim(name)//'_'//trim(adjustl(nlabel))//'.h5'
!            endif
         endif
         order = QUADRATIC
         if (present(inorder)) then
            if ((inorder >= 1) .and. (inorder <= 2)) then
               order = inorder
            endif
         endif
!         if (.not.present(time)) then
!             time = 0.0
!         endif
! write guard cells if present

         if (order==QUADRATIC) then
 !           if ((mx+2) <= nxv) mx = mx + 1
            call PWRITE2_HDF(f,nxv-4,nxv,nypmx,trim(fname),label_code,it,time,QUADRATIC)
!            call PWRITE2_HDF(f,mx,kyp,nxv,iunit,lrec,&
!     &      trim(fname),time,QUADRATIC)
         else
 !           if ((mx+1) <= nxv) mx = mx + 1
            call PWRITE2_HDF(f,nxv-1,nxv,nypmx,trim(fname),label_code,it,time,LINEAR)
!            call PWRITE2_HDF(f,mx,kyp,nxv,iunit,lrec,&
!     &      trim(fname),time,LINEAR)
         endif
!         if (monitor==2) call werrfl(class//sname//' complete')
         end subroutine ipwrite2_hdf
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine create_h5_dataset(fileID, datasetID, rootID, dataspaceID, plistID, name,nx,ny, &
		& label_code,iter,time,xaxis_range,yaxis_range)

	  integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID, fileID
	  integer :: label_code, nx, ny, iter
	  real :: time
		character*(*) name
		real,dimension(:) :: xaxis_range,yaxis_range
	  
	  real,dimension(2) :: box_min, box_max
		integer(hsize_t), dimension(1) :: dims
	  integer(KIND=4) :: ierr
	  integer(KIND=4) :: iter_i4
	  real :: dt
		integer, parameter :: izero = ichar('0')
!		real, dimension(2) :: axis_range
		character*70 title
		character*70 xtitle,xunits
		character*70 ytitle,yunits
		character*70 ztitle,zunits
		
		box_min(1) = 0;box_min(2) = 0;
		box_max(1) = nx-1;box_max(2) = ny-1;
	
		call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
		call h5pcreate_f(H5P_FILE_ACCESS_F, plistID, ierr)
		! setup alignment 
		! although this is not parallel I/O specific it is only recommended for this type of I/O
!		if ( h5_tune%alignment(1) /= -1 ) then
!			call h5pset_alignment_f( plistID, h5_tune%alignment(1), h5_tune%alignment(2), ierr )
!		endif

		 ! Create the file
		 call h5fcreate_f(name, H5F_ACC_TRUNC_F, fileID, ierr, access_prp = plistID)
	
		 ! Close the property list
		 call h5pclose_f(plistID, ierr)
	
		! this is required for older versions of hdf5 that don't allow setting attributes
		! for the file, only the root group
		call h5gopen_f( fileID, '/', rootID, ierr )
	 
		! add name property
		call add_h5_attribute( rootID, 'NAME', name ) 
	
		! set a default write property list for axis values
		plistID = H5P_DEFAULT_F
	
		call add_h5_attribute( rootID, 'TYPE', 'grid' ) 
		
		! add axis information
		call h5gcreate_f( rootID, 'AXIS', axisGroupID, ierr) 
		
		dims(1) = 2
		call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
	
	
		!xaxis
		call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+1), H5T_NATIVE_DOUBLE, dataspaceID, &
							datasetID, ierr )		
		call add_h5_attribute( datasetID, 'TYPE', 'linear' ) 
		call add_h5_attribute( datasetID, 'UNITS', xunits ) 
		call add_h5_attribute( datasetID, 'NAME', xtitle ) 
		call add_h5_attribute( datasetID, 'LONG_NAME', xtitle ) 
!		axis_range(1) = 0
!		axis_range(2) = nx-1
		call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, xaxis_range, dims, ierr, xfer_prp = plistID )
		call h5dclose_f( datasetID, ierr )
		!yaxis
		call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+2), H5T_NATIVE_DOUBLE, dataspaceID, &
							datasetID, ierr )
		call add_h5_attribute( datasetID, 'TYPE', 'linear' ) 
		call add_h5_attribute( datasetID, 'UNITS', yunits ) 
		call add_h5_attribute( datasetID, 'NAME', ytitle ) 
		call add_h5_attribute( datasetID, 'LONG_NAME', ytitle ) 
!		axis_range(1) = 0
!		axis_range(2) = ny-1
		call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, yaxis_range, dims, ierr, xfer_prp = plistID )
		call h5dclose_f( datasetID, ierr )
	  
		! simulation information
		! possibly move to a separate group
		if (iter /= 0) then
			call add_h5_attribute( rootID, 'DT', time / real(iter) )
		else
			call add_h5_attribute( rootID, 'DT', 0.2 )
		endif
			
		call add_h5_attribute( rootID, 'TIME UNITS', '!Nt!Mw!N!Dp!N' ) 
		call add_h5_attribute( rootID, 'XMIN', box_min ) 
		call add_h5_attribute( rootID, 'XMAX', box_max ) 
		call add_h5_attribute( rootID, 'PERIODIC', .TRUE. ) 
!		call add_h5_atribute( rootID, 'MOVE C', diagFile%move_c ) 
		call add_h5_attribute( rootID, 'TIME', time ) 
		iter_i4 = iter
		call add_h5_attribute( rootID, 'ITER', iter_i4 ) 

		call h5gclose_f( axisGroupID, ierr ) 
	  call h5gclose_f( rootID, ierr )

	end subroutine create_h5_dataset


!From osiris
!---------------------------------------------------------------------------------------------------
subroutine add_h5_dataset_2d( parentID, name, dataset, units, long_name, tag  )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 2
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  
  real, dimension(:,:), intent(in) :: dataset

  character( len = * ), intent(in), optional :: units, long_name, tag

  
  integer(hsize_t), dimension(rank) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer(KIND=4) :: ierr
  
  dims(1) = size( dataset, 1 )
  dims(2) = size( dataset, 2 )
  
  ! create the dataset
  call h5screate_simple_f(rank, dims, dataspaceID, ierr ) 
  call h5dcreate_f( parentID, name, H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr )
  
  ! add optional attributes
  if ( present(units) ) then
	call add_h5_attribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_attribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_attribute( datasetID, 'TAG', long_name ) 
  endif
  
  ! write the data
  call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, dataset, dims, ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  call h5dclose_f( datasetID, ierr )
  
end subroutine add_h5_dataset_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_double( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real, intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(KIND=4) :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_v1_double( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  real, dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(KIND=4) :: ierr
  
  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_v1_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_logical( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
!  integer :: bool
  integer(KIND=4) :: ierr,bool
  
  dims(1) = 1
  if ( attribute ) then
    bool = 1
  else
    bool = 0
  endif
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_logical
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_logical_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer :: i
  integer(KIND=4) :: ierr
  integer(hsize_t), dimension(1) :: dims
  integer(KIND=4), dimension(:), allocatable :: bool
  
  dims(1) = size(attribute)
  allocate( bool(dims(1)) )
  do i = 1, dims(1)
	if ( attribute(i) ) then
	  bool(i) = 1
	else
	  bool(i) = 0
	endif
  enddo
  
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

  deallocate(bool)
  
end subroutine add_h5_attribute_logical_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_str( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: size
  
  integer(KIND=4) :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  
  size = len(attribute)
  call h5tset_size_f(typeID, size, ierr)
  
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5tclose_f( typeID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_str
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_str_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: maxlen
  integer :: i
  integer(KIND=4) :: ierr
  
  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  
  maxlen = 0
  do i = 1, size(attribute)-1
    if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
  enddo
  
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  call h5tset_size_f(typeID, maxlen, ierr)
 
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5tclose_f( typeID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_str_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_attribute_int( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  integer(KIND=4), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(KIND=4) :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_attribute_int
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine close_hdf5_file( file_id )
!---------------------------------------------------

   implicit none
   
   ! dummy variables
   
   integer(HID_T),    intent(in) :: file_id
   
   ! local variables
   
   integer(KIND=4) :: ierr
	   
   ! executable statements
   call h5fclose_f(file_id, ierr)

end subroutine close_hdf5_file
!---------------------------------------------------

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
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
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
      		ytitle = 'v!Dx!N'
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


      	end select
      end subroutine do_labels


		end module hdf_write_jf
