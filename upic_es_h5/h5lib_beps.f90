! hdf5io module for QuickPIC Open Source 1.0
!> @author
!> Weiming An
! DESCRIPTION:
!> HDF5 I/O class originally for QuickPIC, extended for UPIC-EMMA
! update: 04/18/2016
! 
!
! update:
!
! 1/9/2017  --> added serial attribute and serial write functions for diagnostics such as phase space
!               (FST)
!               
!               new functions
!               pwfield_serial
!               wrattr_file_serial
!               
!

      module hdf5io_class

!     use parallel_pipe_class
      use HDF5
      use mpi
         
      implicit none
      
      private
      
      public :: hdf5file, pwfield ,pwfield_serial, pwpart
      
      type hdf5file
         
         private
         
         character(len=100) :: filename, timeunits, dataname, units, label
         character(len=100) :: ty
         integer :: n, rank
         real :: t, dt
         character(len=100), dimension(3) :: axisname, axislabel, axisunits
         real, dimension(3) :: axismax, axismin

         
         contains
         
         generic :: new => init_hdf5file
         procedure, private :: init_hdf5file
         
      end type 

      interface add_h5_atribute
         module procedure add_h5_atribute_str
         module procedure add_h5_atribute_str_v1
         module procedure add_h5_atribute_single
         module procedure add_h5_atribute_v1_single
         module procedure add_h5_atribute_int
         module procedure add_h5_atribute_v1_int
      end interface
      
      interface pwfield
        module procedure pwfield_3d
        module procedure pwfield_2d
      end interface

      interface pwfield_serial
        module procedure pwfield_3d_serial
        module procedure pwfield_2d_serial
      end interface

      
      interface pwpart
        module procedure pwpart_2d
      end interface

      
      contains
      
       subroutine init_hdf5file(this,filename,timeunits,ty,n,t,dt,axisname,&
      &axislabel,axisunits,axismax,axismin,dataname,units,label,rank)

         implicit none
         
         class(hdf5file), intent(inout) :: this
         character(len=*), intent(in), optional :: filename, timeunits
         character(len=*), intent(in), optional :: ty, dataname, units, label
         integer, intent(in), optional :: n, rank
         real, intent(in), optional :: t, dt
         character(len=*), dimension(3), intent(in), optional :: axisname, &
         &axislabel, axisunits
         real, dimension(3), intent(in), optional :: axismax, axismin
         
         if (present(filename)) then
            this%filename = filename
         else
            this%filename = 'file.h5'
         end if

         if (present(dataname)) then
            this%dataname = dataname
         else
            this%dataname = 'Data'
         end if

         if (present(units)) then
            this%units = units
         else
            this%units = ''
         end if

         if (present(label)) then
            this%label = label
         else
            this%label = 'Data'
         end if

         if (present(timeunits)) then
            this%timeunits = timeunits
         else
            this%timeunits = 'a.u.'
         end if

         if (present(ty)) then
            this%ty = ty
         else
            this%ty = 'grid'
         end if
         
         if (present(n)) then
            this%n = n
         else
            this%n = 1
         end if

         if (present(rank)) then
            this%rank = rank
         else
            this%rank = 2
         end if

         if (present(t)) then
            this%t = t
         else
            this%t = 1.0
         end if

         if (present(dt)) then
            this%dt = dt
         else
            this%dt = 1.0
         end if

         if (present(axisname)) then
            this%axisname = axisname
         else
            this%axisname = (/'x1','x2','x3'/)
         end if

         if (present(axislabel)) then
            this%axislabel = axislabel
         else
            this%axislabel = (/'x1','x2','x3'/)
         end if

         if (present(axisunits)) then
            this%axisunits = axisunits
         else
            this%axisunits = (/'a.u.','a.u.','a.u.'/)
         end if

         if (present(axismax)) then
            this%axismax = axismax
         else
            this%axismax = (/1.0,1.0,1.0/)
         end if

         if (present(axismin)) then
            this%axismin = axismin
         else
            this%axismin = (/0.0,0.0,0.0/)
         end if

      end subroutine init_hdf5file
!
      subroutine add_h5_atribute_str( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         character( len = * ), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, typeID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer(size_t) :: size
          
         integer :: ierr
          
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
        
      end subroutine add_h5_atribute_str        
!
      subroutine add_h5_atribute_str_v1( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         character( len = * ), dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, typeID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer(size_t) :: maxlen
         integer :: i, ierr
          
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
        
      end subroutine add_h5_atribute_str_v1
!
      subroutine add_h5_atribute_single( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         real, intent(in) :: attribute
         
         integer(hid_t) :: dataspaceID, attrID
         integer(hid_t) :: d_float 
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         d_float = detect_precision()
         dims(1) = 1
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, d_float, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, d_float, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_single
!
      subroutine add_h5_atribute_v1_single( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         real, dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hid_t) :: d_float
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         d_float = detect_precision()
         dims(1) = size(attribute)
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, d_float, dataspaceID, attrID, ierr )
         call h5awrite_f( attrID, d_float, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_v1_single
!
      subroutine add_h5_atribute_int( objID, name, attribute )

         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         integer, intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         dims(1) = 1
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID,&
         &attrID, ierr )
         call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_int
!
      subroutine add_h5_atribute_v1_int( objID, name, attribute )
        
         implicit none
          
         integer(hid_t), intent(in) :: objID
         character( len = * ), intent(in) :: name
         integer, dimension(:), intent(in) :: attribute
          
         integer(hid_t) :: dataspaceID, attrID
         integer(hsize_t), dimension(1) :: dims
         integer :: ierr
          
         dims(1) = size(attribute)
         call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
         call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID,&
         &attrID, ierr )
         call h5awrite_f( attrID, H5T_NATIVE_INTEGER, attribute, dims, ierr)
         call h5aclose_f( attrID, ierr )
         call h5sclose_f( dataspaceID, ierr )
        
      end subroutine add_h5_atribute_v1_int
!      
      subroutine wrattr_file(this,file_id,xferID)
      
         implicit none
         
         class(hdf5file), intent(in) :: this
         integer(hid_t), intent(in) :: file_id
! local data
          integer(hid_t) :: rootID, aid, dspace_id, dset_id, treal, xferID
          integer(hsize_t), dimension(1) :: dims
          integer, parameter :: zero = ichar('0')          
          integer :: i, ierr
          
          treal = detect_precision()
          
          call h5gopen_f(file_id, '/', rootID, ierr)

          call add_h5_atribute(rootID, 'NAME', this%dataname) 
          call add_h5_atribute(rootID, 'TYPE', this%ty) 
          call add_h5_atribute(rootID, 'TIME', this%t)
          call add_h5_atribute(rootID, 'ITER', this%n) 
          call add_h5_atribute(rootID, 'DT', this%dt)
          call add_h5_atribute(rootID, 'TIME UNITS', this%timeunits) 
          
          if (this%ty == 'grid') then
             call h5gcreate_f(rootID, 'AXIS', aid, ierr) 

             dims(1) = 2 
             call h5screate_simple_f(1, dims, dspace_id, ierr ) 
             do i = 1, this%rank
                call h5dcreate_f(aid, 'AXIS'//char(zero+i), treal, dspace_id,&
                &dset_id, ierr )

                call add_h5_atribute(dset_id, 'TYPE', 'linear') 
                call add_h5_atribute(dset_id, 'UNITS', this%axisunits(i)) 
                call add_h5_atribute(dset_id, 'NAME', this%axisname(i)) 
                call add_h5_atribute(dset_id, 'LONG_NAME', this%axislabel(i)) 

                call h5dwrite_f(dset_id, treal, (/this%axismin(i),this%axismax(i)/),&
                &dims, ierr, xfer_prp=xferID)

                call h5dclose_f(dset_id, ierr)
             enddo
    
             call h5sclose_f(dspace_id, ierr) 
             call h5gclose_f(aid, ierr) 
          endif
             call h5gclose_f(rootID, ierr)

      end subroutine wrattr_file

      subroutine wrattr_file_serial(this,file_id)
      
         implicit none
         
         class(hdf5file), intent(in) :: this
         integer(hid_t), intent(in) :: file_id
! local data
          integer(hid_t) :: rootID, aid, dspace_id, dset_id, treal, xferID
          integer(hsize_t), dimension(1) :: dims
          integer, parameter :: zero = ichar('0')          
          integer :: i, ierr
          
          treal = detect_precision()
          

          call add_h5_atribute(file_id, 'NAME', this%dataname) 
          call add_h5_atribute(file_id, 'TYPE', this%ty) 
          call add_h5_atribute(file_id, 'TIME', this%t)
          call add_h5_atribute(file_id, 'ITER', this%n) 
          call add_h5_atribute(file_id, 'DT', this%dt)
          call add_h5_atribute(file_id, 'TIME UNITS', this%timeunits) 
          
          if (this%ty == 'grid') then
             call h5gcreate_f(file_id, 'AXIS', aid, ierr) 

             dims(1) = 2 
             call h5screate_simple_f(1, dims, dspace_id, ierr ) 
             do i = 1, this%rank
                call h5dcreate_f(aid, 'AXIS'//char(zero+i), treal, dspace_id,&
                &dset_id, ierr )

                call add_h5_atribute(dset_id, 'TYPE', 'linear') 
                call add_h5_atribute(dset_id, 'UNITS', this%axisunits(i)) 
                call add_h5_atribute(dset_id, 'NAME', this%axisname(i)) 
                call add_h5_atribute(dset_id, 'LONG_NAME', this%axislabel(i)) 

                call h5dwrite_f(dset_id, treal, (/this%axismin(i),this%axismax(i)/),&
                &dims, ierr)

                call h5dclose_f(dset_id, ierr)
             enddo
    
             call h5sclose_f(dspace_id, ierr) 
             call h5gclose_f(aid, ierr) 
          endif
             ! call h5gclose_f(rootID, ierr)

      end subroutine wrattr_file_serial
!
      subroutine wrattr_dataset(this,dset_id,unit,name)
      
         implicit none
         
         class(hdf5file), intent(in) :: this
         integer(hid_t), intent(in) :: dset_id
         character(len=*), intent(in), optional :: unit,name

         if (present(unit)) then
            call add_h5_atribute(dset_id, 'UNITS', unit)
         else
            call add_h5_atribute(dset_id, 'UNITS', this%units)
         endif 
          
         if (present(name)) then
            call add_h5_atribute(dset_id, 'LONG_NAME', name) 
         else
            call add_h5_atribute(dset_id, 'LONG_NAME', this%label) 
         endif 
  
      end subroutine wrattr_dataset
!
      subroutine pwfield_3d(file,fd,gs,ls,noff,ierr)

         implicit none
        
         class(hdf5file), intent(in) :: file
         real, dimension(:,:,:), intent(in) :: fd
         integer, dimension(3), intent(in) :: gs, ls
         integer, dimension(2), intent(in) :: noff
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id
         integer(hsize_t), dimension(3) :: start
         integer(hsize_t), dimension(3) :: gsize, lsize
         integer(hsize_t), dimension(2) :: lnoff
         integer :: info
                  
         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff
         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, MPI_COMM_WORLD, info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
      
         call h5fcreate_f(file%filename,H5F_ACC_TRUNC_F,file_id,ierr,&
         &access_prp=flplID) 
         call wrattr_file(file,file_id,xferID)

         call h5screate_simple_f(3, gsize, dspace_id, ierr)
         call h5screate_simple_f(3, lsize, memspaceID, ierr )
         call h5gopen_f(file_id, '/', rootID, ierr)
         call h5dcreate_f(rootID, file%dataname, treal, dspace_id, dset_id,&
         &ierr, dcplID)

         start(1) = 0
         start(2) = lnoff(1)
         start(3) = lnoff(2)
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2),1:lsize(3)),&
         &lsize, ierr, memspaceID, dspace_id, xfer_prp=xferID)

         call wrattr_dataset(file,dset_id)

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_3d
!
      subroutine pwfield_2d(file,fd,gs,ls,noff,ierr)

         implicit none
        
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: fd
         integer, dimension(2), intent(in) :: gs, ls
         integer, intent(in) :: noff
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id
         integer(hsize_t), dimension(2) :: start
         integer(hsize_t), dimension(2) :: gsize, lsize
         integer(hsize_t) :: lnoff
         integer :: info
                  
         ierr = 0
         gsize = gs
         lsize = ls
         lnoff = noff
!	 
!
!
	
         call h5open_f(ierr)
         treal = detect_precision()
         call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
         call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)         
         call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
         info = MPI_INFO_NULL
         call h5pset_fapl_mpio_f(flplID, MPI_COMM_WORLD, info, ierr)
         call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
 ! H5F_ACC_TRUNC_F is equivalent to the "replace" argument in regular Fortran open     
 ! as opposed to H5f_ACC_EXCL
         call h5fcreate_f(file%filename,H5F_ACC_TRUNC_F,file_id,ierr,&
         &access_prp=flplID) 
         call wrattr_file(file,file_id,xferID)

         call h5screate_simple_f(2, gsize, dspace_id, ierr)
         call h5screate_simple_f(2, lsize, memspaceID, ierr )
         call h5gopen_f(file_id, '/', rootID, ierr)
         call h5dcreate_f(rootID, file%dataname, treal, dspace_id, dset_id,&
         &ierr, dcplID)

         start(1) = 0
         start(2) = lnoff
   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, lsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:lsize(1),1:lsize(2)),lsize, ierr,&
         &memspaceID, dspace_id, xfer_prp=xferID)

         call wrattr_dataset(file,dset_id)

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5pclose_f(xferID, ierr)
         call h5pclose_f(dcplID, ierr)
         call h5pclose_f(flplID, ierr)
         call h5gclose_f(rootID, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_2d
!
!
!
      subroutine pwfield_3d_serial(file,fd,gs,ierr)

         implicit none
        
         class(hdf5file), intent(in) :: file
         real, dimension(:,:,:), intent(in) :: fd
         integer, dimension(3), intent(in) :: gs
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id
         integer(hsize_t), dimension(3) :: start
         integer(hsize_t), dimension(3) :: gsize
         integer :: info
                  
         ierr = 0
         gsize = gs
         call h5open_f(ierr)
         treal = detect_precision()
         info = MPI_INFO_NULL
      
         call h5fcreate_f(file%filename,H5F_ACC_TRUNC_F,file_id,ierr)
         call wrattr_file_serial(file,file_id)

         call h5screate_simple_f(3, gsize, dspace_id, ierr)
         ! call h5screate_simple_f(3, lsize, memspaceID, ierr )
         call h5dcreate_f(file_id, file%dataname, treal, dspace_id, dset_id,&
         &ierr)
	 ! call h5dget_space_f(dset_id,file_id,ierr)

   
         call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, gsize,&
         &ierr)

         call h5dwrite_f(dset_id, treal, fd(1:gsize(1),1:gsize(2),1:gsize(3)),&
         &gsize, ierr, memspaceID, dspace_id)

         call wrattr_dataset(file,dset_id)

         call h5sclose_f(memspaceID, ierr)
         call h5sclose_f(dspace_id, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_3d_serial
!
!
!
      subroutine pwfield_2d_serial(file,fd,gs,ierr)

         implicit none
        
         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: fd
         integer, dimension(2), intent(in) :: gs
         integer, intent(inout) :: ierr
! local data
         integer(hid_t) :: treal,flplID, xferID, dcplID, memspaceID 
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id
         integer(hsize_t), dimension(2) :: start
         integer(hsize_t), dimension(2) :: gsize, lsize
         integer(hsize_t) :: lnoff
         integer :: info
                  
         ierr = 0
         gsize = gs
!	 
	
         call h5open_f(ierr)
         treal = detect_precision()
! H5F_ACC_TRUNC_F is equivalent to the "replace" argument in regular Fortran open     
! as opposed to H5f_ACC_EXCL
         call h5fcreate_f(file%filename,H5F_ACC_TRUNC_F,file_id,ierr)
         call wrattr_file_serial(file,file_id)

         call h5screate_simple_f(2, gsize, dspace_id, ierr)
         call h5dcreate_f(file_id, file%dataname, treal, dspace_id, dset_id,&
         &ierr)

         start(1) = 0
         start(2) = 0
   

         call h5dwrite_f(dset_id, treal, fd(1:gsize(1),1:gsize(2)),gsize, ierr)

         call wrattr_dataset(file,dset_id)

         call h5sclose_f(dspace_id, ierr)
         call h5dclose_f(dset_id, ierr)
         call h5fclose_f(file_id, ierr)
         call h5close_f(ierr)
               
      end subroutine pwfield_2d_serial
!
!
    subroutine pwpart_2d(file,part,npp,dspl,delta,ierr)
     
         implicit none

         class(hdf5file), intent(in) :: file
         real, dimension(:,:), intent(in) :: part
         real, dimension(2), intent(in) :: delta
         integer, intent(in) :: npp,dspl
         integer, intent(inout) :: ierr
! local data
         integer :: tnpp, tp, color, pgrp, pid, pnvp, i, j
         integer(hsize_t), dimension(1) :: ldim
         integer, dimension(:), pointer :: np
         integer, dimension(:,:), pointer:: dims
         real, dimension(:), pointer :: buff
         integer(hsize_t), dimension(1) :: start,maxdim
         integer(hid_t) :: treal
         integer(hid_t) :: flplID, xferID, memspaceID, aid
         integer(hid_t) :: file_id, rootID, dset_id, dspace_id, aspace_id
         integer :: info
         integer, dimension(10) :: istat
      
         ierr = 0
         ldim(1) = 1
         call h5open_f(ierr)
         treal = detect_precision()
         
         tnpp = int(npp/dspl)
         tp = 0
         call MPI_ALLREDUCE(tnpp,tp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

         if (tp == 0) then
            call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
            call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
            info = MPI_INFO_NULL
            call h5pset_fapl_mpio_f(flplID, MPI_COMM_WORLD, info, ierr)
            call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)
            call h5fcreate_f(file%filename, H5F_ACC_TRUNC_F, file_id, ierr,&
            &access_prp=flplID) 
            call wrattr_file(file,file_id,xferID)
            call h5gopen_f(file_id, '/', rootID, ierr)
            call h5screate_simple_f(1, ldim, aspace_id, ierr)
            call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
            &aid, ierr )
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
            call h5aclose_f(aid, ierr)
            call h5sclose_f(aspace_id, ierr)
            call h5pclose_f(xferID, ierr)
            call h5pclose_f(flplID, ierr)
            call h5gclose_f(rootID, ierr)
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
            return
         else 
            if (tnpp > 0) then 
               color = 1
            else
               color = MPI_UNDEFINED
            endif
            call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, 0, pgrp, ierr )

            if (tnpp > 0) then
               call MPI_COMM_RANK(pgrp, pid, ierr)
               call MPI_COMM_SIZE(pgrp, pnvp, ierr)
               allocate(np(pnvp), dims(2,pnvp), stat = ierr)
               call MPI_ALLGATHER(tnpp, 1, MPI_INTEGER, np, 1, MPI_INTEGER,&
               &pgrp, ierr)
               dims(1, 1) = 1
               dims(2, 1) = np(1) 
               do i = 2, pnvp
                  dims(1,i) = dims(2,i-1) + 1
                  dims(2,i) = dims(1,i) + np(i) - 1
               enddo
               allocate(buff(tnpp), stat = ierr)

               call h5pcreate_f(H5P_FILE_ACCESS_F, flplID, ierr)         
               call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)  
               info = MPI_INFO_NULL
               call h5pset_fapl_mpio_f(flplID, pgrp, info, ierr)
               call h5pset_dxpl_mpio_f(xferID, H5FD_MPIO_COLLECTIVE_F, ierr)    
               call h5fcreate_f(file%filename, H5F_ACC_TRUNC_F, file_id, ierr,&
               &access_prp=flplID) 
               call wrattr_file(file,file_id,xferID)
               call h5gopen_f(file_id, '/', rootID, ierr)
               call h5screate_simple_f(1, ldim, aspace_id, ierr)
               call h5acreate_f(rootID, 'tp', H5T_NATIVE_INTEGER, aspace_id,&
               &aid, ierr )
               call h5awrite_f(aid, H5T_NATIVE_INTEGER, tp, ldim, ierr)
               call h5aclose_f(aid, ierr)
               call h5sclose_f(aspace_id, ierr)

               do i = 1, 2
                  buff(1:tnpp) = part(i,1:(1+tnpp*dspl-1):dspl)*delta(i)
                  ldim(1) = tp
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  call h5dcreate_f(rootID, 'x'//char(iachar('0')+i), treal,&
                  &dspace_id, dset_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr)
                  start = dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                  &ldim, ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
                  call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                  &name='x_'//char(iachar('0')+i))
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               do i = 1, 3
                  buff(1:tnpp) = part((i+2),1:(1+tnpp*dspl-1):dspl) 
                  ldim(1) = tp
                  call h5screate_simple_f(1, ldim, dspace_id, ierr)
                  call h5dcreate_f(rootID, 'p'//char(iachar('0')+i), treal,&
                  &dspace_id, dset_id, ierr)
                  ldim(1) = tnpp
                  call h5screate_simple_f(1, ldim, memspaceID, ierr)
                  start = dims(1,pid+1) - 1
                  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,start,&
                  &ldim, ierr)
                  call h5dwrite_f(dset_id, treal, buff, ldim, ierr, memspaceID,&
                  &dspace_id, xfer_prp=xferID)
                  call wrattr_dataset(file,dset_id,unit='c/\omega_p',&
                  &name='p_'//char(iachar('0')+i))
                  call h5sclose_f(memspaceID, ierr)
                  call h5sclose_f(dspace_id, ierr)
                  call h5dclose_f(dset_id, ierr)
               enddo

               call h5pclose_f(xferID, ierr)
               call h5pclose_f(flplID, ierr)
               call h5gclose_f(rootID, ierr)
               call h5fclose_f(file_id, ierr)
               deallocate(np,dims,buff)
            endif            
            if (pgrp /= MPI_COMM_NULL) then
               call MPI_COMM_FREE(pgrp, ierr)
            endif
            call h5close_f(ierr)
         endif
         
      end subroutine pwpart_2d


!
      function detect_precision()
         integer(hid_t) :: detect_precision
! local data
         real :: small

         small = 1.0e-12
         small = 1.0 + small
         if (small>1.0) then 
            detect_precision = H5T_NATIVE_DOUBLE 
         else
            detect_precision = H5T_NATIVE_REAL
         endif       

      end function detect_precision     
!      
      end module hdf5io_class

