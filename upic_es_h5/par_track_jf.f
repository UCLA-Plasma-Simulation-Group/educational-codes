! Particle tracking module
! Jay Fahlen
! Based on code from particle tracking in OSIRIS
module par_track_jf
use pinit2d_jf
include 'mpif.h'
implicit none

!-------------------------------------------------------------------------------
! Particle tracks diagnostics and particle track class, taken from OSIRIS by Jay Fahlen
!-------------------------------------------------------------------------------
integer, parameter :: p_max_tracks = 16384
integer,parameter :: p_max_filename_len = 200
integer,parameter :: p_max_tagstring_len = 200
integer, parameter, public :: file_id_tem  =  10
integer, parameter, public :: file_id_rst  =  40

type t_track
  
  integer :: npoints = 0                ! number of points in memory
  integer :: savedpoints = 0            ! number of points saved to disk
  integer, dimension(2) :: tag = 0      ! tag of particle being tracked  
  integer :: part_idx = -1              ! index of particle being tracked (-1) means
                                        ! that particle is not on local node
  
  real, dimension(:,:), pointer :: data => null()  ! track points data
  integer, dimension(:), pointer :: n => null()               ! track points iterations

end type t_track


type t_track_set

  ! maximum number of points to hold
  integer :: maxpoints = -1
  
  ! number of iterations between track points
  integer :: niter = 1			

  ! file holding tags of particles to follow 
  character(len = p_max_filename_len) :: file_tags = 'input_tags'
  
  ! filename to write tracks to
  character(len = p_max_filename_len) :: file_write = ''
  
  ! total number of tracks 
  integer :: ntracks = 0                   
  ! actual tracks
  type( t_track ), dimension(:), pointer :: tracks  => NULL()

  ! list of tags present / missing
  ! each item is the index of the track on the tracks array
  integer, dimension(:), pointer :: present, missing 

  ! list sizes (for simplicity)
  integer :: npresent, nmissing

end type t_track_set

private
integer :: xdim = 2,vdim = 2

public :: assign_tags,t_track,t_track_set, setup, cleanup, missing_particles
public :: new_present_particles,add_track_data,write_tracks,create_file
public :: change_index,update_indexes

!       HDF Function declaration.
integer,external :: sfstart, sfcreate, sfwdata, sfsdtstr, sfsdmstr
integer,external :: sfdimid, sfsdmname, sfsdscale, sfsblsz 
integer,external :: sfendacc, sfend, sfsnatt, sfscompress

!       HDF Constant declaration
integer, parameter :: DFACC_CREATE = 4
integer, parameter :: DFACC_WRITE = 2
integer, parameter :: DFNT_CHAR = 4
integer, parameter :: DFNT_FLOAT32 = 5
integer, parameter :: DFNT_FLOAT64 = 6
integer, parameter :: DFNT_INT32 = 24
integer, parameter :: COMP_CODE_DEFLATE = 4
integer, parameter :: SD_UNLIMITED = 0
integer, parameter :: SD_FILL = 0

! string to id restart data
character(len=*), parameter :: p_track_rst_id = "track rst data - 0x0001"

interface missing_particles
  module procedure missing_particles_list
end interface 

interface new_present_particles
  module procedure new_present_particles_bounds
  module procedure new_present_particles_single
end interface

interface add_track_data
  module procedure add_track_set_data
end interface

interface write_tracks
  module procedure write_tracks
end interface

interface create_file
  module procedure create_single_track_file
  module procedure create_track_set_file
end interface

interface setup
  module procedure setup_single_track
  module procedure setup_track_set_send
  module procedure setup_track_set
end interface

interface cleanup
  module procedure cleanup_single_track
  module procedure cleanup_track_set
end interface

interface save_track_data
  module procedure save_track_data_int
  module procedure save_track_data_double
end interface

interface update_indexes
  module procedure update_indexes_single_track
  module procedure update_indexes_track_set
end interface

interface change_index
  module procedure change_index_track_set
end interface

interface restart_write
  module procedure restart_write_tracks
  module procedure restart_write_single_track
end interface

interface restart_read
  module procedure restart_read_tracks
  module procedure restart_read_single_track
end interface

interface heapsort
!	module procedure heapsort_inplace
	module procedure heapsort_index
end interface


! this should only be called internally
private :: restart_read

 contains

!-------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-------------------------------------------------------------------------------
subroutine update_indexes_single_track( track, new_idx )

   implicit none

   type( t_track ), intent(inout) :: track
   integer, dimension(:), intent(in) :: new_idx

   track%part_idx = new_idx(track%part_idx)

end subroutine update_indexes_single_track

!-------------------------------------------------------------------------------
!  Creates HDF5 group and datasets to save track data
!-------------------------------------------------------------------------------
subroutine create_single_track_file( self, fileID, ndump )
!-------------------------------------------------------------------------------
   
   use hdf5
   
   implicit none
   
   type( t_track ), intent(in) :: self
   integer(hid_t), intent(in) :: fileID
   integer, intent(in) :: ndump
   
   ! local variables
   integer(hsize_t), dimension(1) :: dims, maxdims
   integer(hid_t) :: groupID, dataspaceID, propertyID, datasetID, attrID
   integer :: i, ierr

   
   character( len = p_max_tagstring_len ) :: tagstring
   character( len = p_max_tagstring_len) :: temp1,temp2

	 write(temp1,*) self%tag(1)
	 write(temp2,*) self%tag(2)
	 temp1 = adjustl(temp1)
	 temp2 = adjustl(temp2)

   ! save the track
   tagstring = trim(temp1)//'-'//trim(temp2)
      
   ! create group
   call h5gcreate_f( fileID, tagstring, groupID, ierr) 
   if ( ierr/= 0 ) then 
	 write(*,*)'Unable to create group.'
!	 call abort_program(p_err_diagfile)
   endif
   
   ! add tag as an attribute
   dims(1) = 2
   call h5screate_simple_f( 1, dims, dataspaceID, ierr)
   call h5acreate_f( groupID, 'tag', H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
   call h5awrite_f( attrID, H5T_NATIVE_INTEGER, self%tag, dims, ierr)
   call h5aclose_f( attrID, ierr )
   call h5sclose_f(dataspaceID, ierr)
   
   ! create unlimited dataspace
   dims(1) = 0
   maxdims(1) = H5S_UNLIMITED_F
   call h5screate_simple_f( 1, dims, dataspaceID, ierr, maxdims )
   
   ! set chunk size
   call h5pcreate_f( H5P_DATASET_CREATE_F, propertyID, ierr )
   dims(1) = ndump
   call h5pset_chunk_f( propertyID, 1, dims, ierr )

   ! create datasets
   call h5dcreate_f( groupID, 'n', H5T_NATIVE_INTEGER, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )

   call h5dcreate_f( groupID, 't', H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )
   
   do i = 1, xdim
	 call h5dcreate_f( groupID, 'x'//char(iachar('0')+i), H5T_NATIVE_DOUBLE, dataspaceID, &
	                   datasetID, ierr, propertyID )
     call h5dclose_f( datasetID, ierr )
   enddo

   do i = 1, vdim
	 call h5dcreate_f( groupID, 'p'//char(iachar('0')+i), H5T_NATIVE_DOUBLE, dataspaceID, &
	                   datasetID, ierr, propertyID )
     call h5dclose_f( datasetID, ierr )
   enddo

   call h5dcreate_f( groupID, 'q', H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )

   call h5dcreate_f( groupID, 'ene', H5T_NATIVE_DOUBLE, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )
   
   ! close group
   call h5gclose_f( groupID, ierr )
   
end subroutine create_single_track_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Adds integer track data to file
!-------------------------------------------------------------------------------
subroutine save_track_data_int( fileID, datasetName, newData, savedPoints, npoints )
   
   use hdf5
  
   implicit none

   integer(hid_t), intent(in) :: fileID
   character( len = * ), intent(in) :: datasetName
   
   integer, dimension(:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: datasetID, filespaceID, memspaceID
   integer(hsize_t), dimension(1) :: start, count, size
   integer :: ierr
   
   ! open the dataset
   call h5dopen_f( fileID, datasetName, datasetID, ierr )	  
   
   ! ensure the dataset is big enough for data
   size(1) = savedpoints + npoints
!print*,"i",size(1),fileID
   call h5dextend_f(datasetID, size, ierr)
   
   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = savedPoints
   count(1) = nPoints
   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(1, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, H5T_NATIVE_INTEGER, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        
   call h5dclose_f( datasetID, ierr )
   
end subroutine save_track_data_int
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Adds double track data to file
!-------------------------------------------------------------------------------
subroutine save_track_data_double( fileID, datasetName, newData, savedPoints, npoints )
  
   use hdf5

   implicit none

   integer(hid_t), intent(in) :: fileID
   character( len = * ), intent(in) :: datasetName
   
   real, dimension(:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: datasetID, filespaceID, memspaceID
   integer(hsize_t), dimension(1) :: start, count, size
   integer :: ierr
   
!   print*,newData
   
   ! open the dataset
   call h5dopen_f( fileID, datasetName, datasetID, ierr )	  
   
   ! ensure the dataset is big enough for data
   size(1) = savedpoints + npoints
!print*,"d",size(1)
   call h5dextend_f(datasetID, size, ierr)
   
   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = savedPoints
   count(1) = nPoints
   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(1, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, H5T_NATIVE_DOUBLE, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        
   call h5dclose_f( datasetID, ierr )
   
end subroutine save_track_data_double
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! sort track according to iteration (needed after gathering from all nodes)
!-------------------------------------------------------------------------------
subroutine sort_track( track )
   
   implicit none

   type( t_track ), intent(inout) :: track
   
   integer, dimension(:), allocatable :: idx, temp_int
   real, dimension(:), allocatable :: temp_double
   integer :: i, j, ierr
   
   ! if just 1 point there's no need to sort
   if ( track%npoints > 1 ) then
      allocate( idx(track%npoints), stat = ierr )
	  if ( ierr /= 0 ) then
		 write(*,*)'error allocating data for track object sort.'
!		 call abort_program(p_err_alloc)
	  endif
	  	  
	  ! get sorted index
	  call heapsort( track%n(1:track%npoints), idx )
	  
	  ! reorder iterations
	  allocate( temp_int(track%npoints), stat = ierr )
	  if ( ierr /= 0 ) then
		 write(*,*)'error allocating data for track object sort.'
!		 call abort_program(p_err_alloc)
	  endif
	  do i = 1, track%npoints
	    temp_int(i) = track%n(idx(i))
	  enddo
	  do i = 1, track%npoints
	    track%n(i) = temp_int(i)
	  enddo
      deallocate( temp_int )	
      	  
	  ! reorder remaining data
	  allocate( temp_double(track%npoints), stat = ierr )
	  if ( ierr /= 0 ) then
		 write(*,*)'error allocating data for track object sort.'
!		 call abort_program(p_err_alloc)
	  endif
	  do j = 1, 3 + xdim + vdim
		 do i = 1, track%npoints
		   temp_double(i) = track%data(j,idx(i))
		 enddo
		 do i = 1, track%npoints
		   track%data(j,i) = temp_double(i)
		 enddo
	  enddo
      deallocate( temp_double )	  
	  
	  deallocate( idx )
   endif
   
end subroutine sort_track

!-------------------------------------------------------------------------------
! Pack single track data into comm buffer
!-------------------------------------------------------------------------------
subroutine pack_data_single_track( track, buffer, bufsize, position )
! This function def has changed cause used to have no_co as parameter, but got rid of it   
   implicit none
   
   type( t_track ), intent(inout) :: track
  
  integer, dimension(:), intent(inout) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position
	integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
	common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld


   integer :: i, ierr

   ! pack number of points
   call mpi_pack( track%npoints, 1, MPI_INTEGER, &
                  buffer, bufsize, position, &
			      lworld, ierr )
	
   if ( track%npoints > 0 ) then
	  ! add iterations
	  call mpi_pack( track%n(1:track%npoints), track%npoints, MPI_INTEGER, &
					 buffer, bufsize, position, &
					 lworld, ierr )
	  
	  ! add datasets
	  do i = 1, 3 + xdim + vdim
		 call mpi_pack( track%data(i, 1:track%npoints), track%npoints, MPI_DOUBLE_PRECISION, &
						buffer, bufsize, position, &
						lworld, ierr )
	  enddo
      
      ! free data points
      track%npoints = 0
      
   endif   

end subroutine pack_data_single_track

!-------------------------------------------------------------------------------
! Unpack single track data from comm buffer
!-------------------------------------------------------------------------------
subroutine unpack_data_single_track( track, buffer, bufsize, position )
   
   implicit none
   
   type( t_track ), intent(inout) :: track
   integer, intent(inout), dimension(:) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position

   integer :: npoints, i, ierr
   real, dimension(:), allocatable :: quant_buffer
   integer, dimension(:), allocatable :: n_buffer

   ! unpack number of points
   call mpi_unpack( buffer, bufsize, position, &
                    npoints, 1, MPI_INTEGER, & 
			        mpi_comm_world, ierr )
	
   if ( npoints > 0 ) then
	  allocate( n_buffer( npoints ) )
	  ! unpack iterations
	  call mpi_unpack( buffer, bufsize, position, &
	                   n_buffer, npoints, MPI_INTEGER, &
					   mpi_comm_world, ierr )
      track%n(track%npoints+1: track%npoints+npoints) = n_buffer
      deallocate( n_buffer )
	  
	  ! add datasets
	  allocate( quant_buffer( npoints ) )
	  do i = 1, 3 + xdim + vdim
		 call mpi_unpack( buffer, bufsize, position, &
		                  quant_buffer, npoints, MPI_DOUBLE_PRECISION, &
						  mpi_comm_world, ierr )
		 track%data(i, track%npoints+1: track%npoints+npoints) = quant_buffer
	  enddo
	  deallocate( quant_buffer )
	  
	  track%npoints = track%npoints + npoints
   endif   

end subroutine unpack_data_single_track

!-------------------------------------------------------------------------------
! Save points currently in memory to disk and discards them
!-------------------------------------------------------------------------------
subroutine write_single_track( track, fileID )
!-------------------------------------------------------------------------------

   use hdf5
   
   implicit none

   type( t_track ), intent(inout) :: track
   integer(hid_t), intent(in) :: fileID
   
   
   integer :: i
   character( len = p_max_tagstring_len ) :: tagstring
   character( len = p_max_tagstring_len) :: temp1,temp2


!print*,"track%npoints",track%npoints

   if ( track%npoints > 0 ) then
   
     write(temp1,*) track%tag(1)
     write(temp2,*) track%tag(2)
     temp1 = adjustl(temp1)
     temp2 = adjustl(temp2)
     
	  tagstring = '/'//trim(temp1)//'-'//trim(temp2)

!print*,trim(tagstring)

      ! write iteration
      
!      print*,track%n(1:track%npoints),track%data(1,1:track%npoints),track%data(2,1:track%npoints),&
!print*,"v",track%data(6,1:track%npoints),track%data(7,1:track%npoints)
      call save_track_data( fileID, trim(tagstring)//'/n', track%n(1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/t', track%data(1,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/q', track%data(2,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/ene', track%data(3,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      do i = 1, xdim
		 call save_track_data( fileID, trim(tagstring)//'/x'//char(iachar('0')+i), &
		                       track%data(3+i,1:track%npoints), &
							   track%savedpoints, track%npoints )     
      enddo

      do i = 1, vdim
		 call save_track_data( fileID, trim(tagstring)//'/p'//char(iachar('0')+i), &
		                       track%data(3+xdim+i,1:track%npoints), &
							   track%savedpoints, track%npoints )     
      enddo
      
      
      ! update saved points data
      track%savedpoints = track%savedpoints + track%npoints      
      track%npoints = 0
   endif
   
end subroutine write_single_track
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Add current iteration point to track
!-------------------------------------------------------------------------------
subroutine add_single_track_data( track, part, n, t ) 
! this function def changed too, to have part array not species object
!-------------------------------------------------------------------------------
   
   implicit none
   
   type( t_track ), intent(inout) :: track
!   type( t_species ), intent(in) :: species
 	 real,dimension(:,:,:) :: part
	 integer, intent(in) :: n
   real, intent(in) :: t
   
   real, dimension(3 + xdim + vdim) :: point
   real :: u2
      
   ! debug
   ! check that we are looking at the correct tag
!   if ((track%tag(1) /= species%tag( 1, track%part_idx )) .or. &
!       (track%tag(2) /= species%tag( 2, track%part_idx )) ) then
!     
!     print *, mpi_node(), ' Track is pointing to invalid index, tag ', track%tag
!     call abort_program()   
!       
!   endif
   
   track%npoints = track%npoints + 1 
   track%n( track%npoints ) = n
   
   point(1) = t
!   point(2) = species%q( track%part_idx )
   point(2) = 1.
  
   ! energy must be calculated here
   u2 = part( 3, track%part_idx, 1 )**2 + &
        part( 4, track%part_idx, 1 )**2  
!   point(3) = u2 / ( sqrt(1.0d0 + u2) + 1.0d0) 
! This is not always a relativistic code, so for now just use nonrelativistic energy
	 point(3) = 0.5*u2
   
   point(4 : 5 ) = part( 1:2, track%part_idx,1 )
   point(6 : 7 ) = part( 3:4, track%part_idx, 1 )
      
   track%data( : , track%npoints ) = point
   
end subroutine add_single_track_data
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Cleanup single track object
!-------------------------------------------------------------------------------
subroutine cleanup_single_track( track ) 
!-------------------------------------------------------------------------------
   implicit none
   
   type( t_track ), intent(inout) :: track
   
   integer :: ierr
   
   track%savedpoints = 0
   track%npoints = 0
   track%tag = 0
   deallocate( track%data, track%n, stat = ierr )
   if ( ierr/=0 ) then 
      write(*,*)'error deallocating data for track object.'
!	  call abort_program(p_err_dealloc)
   endif
   
end subroutine cleanup_single_track
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Setup single track object
!-------------------------------------------------------------------------------
subroutine setup_single_track( track, max_points, tag ) 
!-------------------------------------------------------------------------------
   implicit none
   
   type( t_track ), intent(inout) :: track
   integer, intent(in) :: max_points
   integer, dimension(:) :: tag
   
   integer :: ierr, data_size
   
   allocate( track%n(max_points), stat = ierr )
   if ( ierr/=0 ) then 
      write(*,*)'error allocating data for track object.'
!	  call abort_program(p_err_dealloc)
   endif

   !           t +  q + ene +   x     +   p     
   data_size = 1 +  1 +  1  + xdim + vdim

   allocate( track%data(data_size, max_points), stat = ierr )
   if ( ierr/=0 ) then 
      write(*,*)'error allocating data for track object.'
!	  call abort_program(p_err_alloc)
   endif
   
   track%savedpoints = 0
   track%npoints = 0
   track%tag(1) = tag(1)
   track%tag(2) = tag(2)
   track%part_idx = -1

end subroutine setup_single_track


!-------------------------------------------------------------------------------

!*******************************************************************************
!******************************* Track set routines ****************************
!*******************************************************************************

!-------------------------------------------------------------------------------
! Change particle pointer if particle has moved inside particle buffer (usually
! because another one was deleted, sort is handled below)
!-------------------------------------------------------------------------------
subroutine change_index_track_set( track_set, old_idx, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
  
   integer, intent(in) :: old_idx, new_idx

   integer :: j
   
   do j = 1, track_set%npresent
     if ( track_set%tracks(track_set%present(j))%part_idx == old_idx ) then
       track_set%tracks(track_set%present(j))%part_idx = new_idx
       exit
     endif
   enddo

end subroutine change_index_track_set

!-------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-------------------------------------------------------------------------------
subroutine update_indexes_track_set( track_set, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, dimension(:), intent(in) :: new_idx

   integer :: i
      
   do i = 1, track_set%npresent
      call update_indexes( track_set%tracks( track_set%present(i) ), new_idx )
   enddo

end subroutine update_indexes_track_set

!-------------------------------------------------------------------------------
! Particles have left the local node, update present/missing lists
!-------------------------------------------------------------------------------
subroutine missing_particles_list( track_set, idx, n_idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(:), intent(in) :: idx
  integer, intent(in) :: n_idx
  
  integer :: i, j

  ! loop through all the missing particles
  do i = 1, n_idx
	 
	 ! loop through all missing tags 
	 j = 1
	 do
		if (j > track_set%npresent) exit

		if (idx(i) == track_set%tracks(track_set%present(j))%part_idx) then  
		   ! debug
		   !print *, mpi_node(), ' particle has left node, tag = ', &
		   !                       track_set%tracks(track_set%present(j))%tag 
		   
		   ! delete particle index
		   ! (not necessary for production but useful for debugging)
		   track_set%tracks(track_set%present(j))%part_idx = -1

		   ! add track to missing list
		   track_set%nmissing = track_set%nmissing + 1
		   track_set%missing( track_set%nmissing ) = track_set%present(j)

		   ! remove track from present list
		   track_set%present(j) = track_set%present(track_set%npresent) 
		   track_set%npresent = track_set%npresent - 1
		   		   
		   exit
		else 
		   j = j + 1
		endif
	 enddo
	 
	 ! if no tags are present we can end the search
	 if (track_set%npresent == 0) exit
  enddo

end subroutine missing_particles_list

!-------------------------------------------------------------------------------
! New particle is present in local nodes, check if its tag is in the
! missing list and if so store their index and move it to the present list
!-------------------------------------------------------------------------------
subroutine new_present_particles_single( track_set, tag, idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(2), intent(in) :: tag
  integer, intent(in) :: idx
  
  integer :: j
    
  if ( track_set%nmissing > 0 ) then
    j = 1
    do
      if (j > track_set%nmissing) exit
      if ( track_set%tracks(track_set%missing(j))%tag(1) == tag(1) .and. &
				 track_set%tracks(track_set%missing(j))%tag(2) == tag(2) ) then
				 
				 ! debug
				 !print *, mpi_node(), ' particle has entered node node, tag = ', &
				 !                       tag 
		
				 
				 ! store particle index
				 track_set%tracks( track_set%missing(j))%part_idx = idx

				 ! add tag to present list
				 track_set%npresent = track_set%npresent + 1
				 track_set%present( track_set%npresent ) = track_set%missing(j)
		
				 ! remove tag from missing list
				 track_set%missing(j) = track_set%missing(track_set%nmissing) 
				 track_set%nmissing = track_set%nmissing - 1
				 exit
      else
         j = j + 1
      endif
    enddo
  endif
  
end subroutine new_present_particles_single

!-------------------------------------------------------------------------------
! New particles are present in local nodes, check if their tags are in the
! missing list and if so store their index and move them to the present list
!-------------------------------------------------------------------------------
subroutine new_present_particles_bounds( track_set, part, idx0, idx1 )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
!  type( t_species ), intent(in) :: spec
	real,dimension(:,:,:) :: part
  integer, intent(in) :: idx0, idx1
  
  integer, dimension(2) :: tag,int_tag
  real :: real_tag
  equivalence (real_tag,int_tag)
  integer :: i, j

  ! loop through all the new particles
  do i = idx0, idx1
	 
	 ! loop through all missing tags 
	 j = 1
	 do
		if (j > track_set%nmissing) exit

		tag(1) = track_set%tracks(track_set%missing(j))%tag(1)
		tag(2) = track_set%tracks(track_set%missing(j))%tag(2)
		
		real_tag = part(5,i,1)
				
		if (( int_tag(1) == tag(1)) .and. ( int_tag(2) == tag(2))) then  
		   ! store particle index
		   track_set%tracks( track_set%missing(j))%part_idx = i

		   ! add tag to present list
		   track_set%npresent = track_set%npresent + 1
		   track_set%present( track_set%npresent ) = track_set%missing(j)

		   ! remove tag from missing list
		   track_set%missing(j) = track_set%missing(track_set%nmissing) 
		   track_set%nmissing = track_set%nmissing - 1
		   		   
		   exit
		else 
		   j = j + 1
		endif
	 enddo
	 
	 ! if no tags are missing we can end the search
	 if (track_set%nmissing == 0) exit
  enddo

end subroutine new_present_particles_bounds

!-------------------------------------------------------------------------------
! Add current values to track data
!-------------------------------------------------------------------------------
subroutine add_track_set_data( track_set, part, n, t )

  implicit none
  
  type( t_track_set ), intent(inout) :: track_set
!  type( t_species ), intent(in) :: spec
	real,dimension(:,:,:) :: part
  integer, intent(in) :: n
  real, intent(in) :: t
  
  integer :: i

  do i = 1, track_set%npresent
    call add_single_track_data( track_set%tracks( track_set%present(i) ), &
                                part, n, t )
  enddo

end subroutine add_track_set_data

!-------------------------------------------------------------------------------
! Gather single track data from all nodes and sort it
!-------------------------------------------------------------------------------
subroutine gather_data( track_set, nvp, idproc )

   implicit none
  
   type( t_track_set ), intent(inout) :: track_set
!   type( t_node_conf ), intent(in) :: no_co
	 integer :: nvp,idproc
	 
!   integer( p_byte ), dimension(:), allocatable :: buffer
   integer, dimension(:), allocatable :: buffer
   integer :: bufsize, position, totalpoints
   integer :: pack_size_int
   integer :: pack_size_double
   
   
   integer, dimension(mpi_status_size):: stat
   integer :: handle
	integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
	common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
   
   integer :: i, j, ierr

   if ( nvp > 1 ) then  
      
      ! get pack sizes
      call mpi_pack_size( 1, MPI_INTEGER, lworld, pack_size_int, ierr)
      call mpi_pack_size( 1, MPI_DOUBLE_PRECISION, lworld, pack_size_double, ierr)
      
      if ( idproc == 0 ) then
		 ! allocate comm bufer to the max size
		 ! npoints, iteration data, remaining data
				bufsize = track_set%ntracks * ( pack_size_int * (1 + track_set%maxpoints ) + &
		                pack_size_double * track_set%maxpoints * (3 + xdim + vdim) )
		 		 
		 		allocate( buffer( bufsize ), stat = ierr )
		 		if ( ierr /= 0 ) then
					write(*,*)'error allocating data for track object communication.'
!			call abort_program(p_err_alloc)
		 		endif

        do j = 1, nvp-1
           ! print *, ' 1 - sending ping to ', j
 !           call send_ping( no_co, tid(no_co,j), 0 )
            
           ! print *, ' 1 - posted receive from ', j
           call MPI_RECV( buffer, bufsize, MPI_PACKED, j, 0, &
                          lworld, stat, ierr )
                           
           position = 0
           ! print *, ' 1 - unpacking data from ', j
           do i = 1, track_set%ntracks
             call unpack_data_single_track( track_set%tracks(i), buffer, bufsize, position )
           enddo
           ! print *, ' 1 - finished with data from node ',j
         enddo
     
				 ! sort the data
				 do i = 1, track_set%ntracks
						call sort_track( track_set%tracks(i) )
				 enddo
      else
         
         ! wait for ping
         ! print *, my_aid( no_co ) , ' - waiting for ping from node 1'
!         call recv_ping( no_co, tid(no_co, 1), 0 )


				 totalpoints = track_set%tracks(1)%npoints
				 do i = 2, track_set%ntracks
					 totalpoints = totalpoints + track_set%tracks(i)%npoints
				 enddo
				 
				 ! allocate comm bufer to the size of points in this node
				 bufsize = track_set%ntracks * pack_size_int + &
									 totalpoints*(pack_size_int + (3 + xdim + vdim)*pack_size_double )
						 
				 allocate( buffer( bufsize ), stat = ierr )
				 if ( ierr /= 0 ) then
					write(*,*)'error allocating data for track object communication.'
		!			call abort_program(p_err_alloc)
				 endif

         position = 0
         do i = 1, track_set%ntracks
           call pack_data_single_track( track_set%tracks(i), &
                                        buffer, bufsize, position )
         enddo
         
         ! wait for ping
         ! print *, my_aid( no_co ) , ' - waiting for ping from node 1'
         !call recv_ping( no_co, tid(no_co, 1), 0 )
        
         ! send data
         !print *, my_aid( no_co ) , ' - sending data to node 0'
         call MPI_ISEND( buffer, bufsize, MPI_PACKED, 0, 0, &
                         lworld, handle, ierr )
         call MPI_WAIT( handle, stat, ierr )
         
      endif
      
      deallocate( buffer )
     
   endif
   
end subroutine gather_data

!-------------------------------------------------------------------------------
! Write tracks data to file
!-------------------------------------------------------------------------------
subroutine write_tracks( track_set, nvp, idproc )

   use hdf5

   implicit none
   
   type( t_track_set ), intent(inout) :: track_set
!   type( t_node_conf ), intent(in) :: no_co
	 integer :: nvp,idproc
	 
   integer(hid_t) :: fileID
   integer :: i, ierr
   
   if (track_set%ntracks > 0 ) then
		 if ( idproc == 0 ) then 
			! open file
			
			call h5fopen_f(trim(track_set%file_write), H5F_ACC_RDWR_F, fileID, ierr)
		 endif
		 
		 ! gather data from all nodes
		 call gather_data( track_set, nvp, idproc )
		 
		 if ( idproc == 0 ) then 
				! loop through all tracks
				do i = 1, track_set%ntracks
					call write_single_track( track_set%tracks(i), fileID )	 
				enddo
		
				! close file
				call h5fclose_f(fileID, ierr)
		 endif

   endif
   
end subroutine write_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! create tracks file
!-------------------------------------------------------------------------------
subroutine create_track_set_file( track_set, sp_name, ndump, nx, ny,&
                                  dt, periodic, move_c, idproc)
    
   use hdf5

   implicit none
   
   ! parameters
   type( t_track_set ), intent(inout) :: track_set
   character( len = * ), intent(in) :: sp_name
   integer, intent(in) :: ndump,nx,ny, idproc
   
   real, intent(in) :: dt
   logical, dimension(:), intent(in) :: periodic, move_c
   
   ! local vars
   integer(hid_t) :: fileID, dataspaceID, typeID, attrID, groupID
   integer(hsize_t), dimension(1) :: dims
   integer :: i, ierr
   character(len = 128)   :: path
   integer, dimension(xdim) :: bool
   real, dimension(2,2) :: x_bnd
   
   integer, parameter :: p_max_txt_len = 16 
   character( len = p_max_txt_len ), dimension(4 + xdim + xdim ) :: txt 
      
   if (track_set%ntracks > 0 ) then
   
	 ! create tracks data file
	 path  = './DIAG/TRACKS/'
	 track_set%file_write = trim(path) // 'elec-tracks.h5'
  
	 if ( idproc == 0 ) then 
		
		! create directory if directory is missing
!		call mkdir_f_( trim(path), ierr ) !this is made in diag_jf.f

		! create file erasing any previous one
		! note that this routine is only called for n = 0
		call h5fcreate_f(trim(track_set%file_write), H5F_ACC_TRUNC_F, fileID, ierr)
		
		! set global attributes
		call h5gopen_f( fileID, '/', groupID, ierr )
		
		! create a scalar dataspace
		dims(1) = 1
		call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
		
		! write type
		call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
		call h5tset_size_f(typeID, len('tracks'), ierr)
		call h5acreate_f( groupID, 'TYPE', typeID, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, typeID, 'tracks', dims, ierr)
		call h5aclose_f( attrID, ierr )
		call h5tclose_f( typeID, ierr )
		
		! write species name
		call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
		call h5tset_size_f(typeID, len(sp_name), ierr)
		call h5acreate_f( groupID, 'NAME', typeID, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, typeID, sp_name, dims, ierr)
		call h5aclose_f( attrID, ierr )
		call h5tclose_f( typeID, ierr )
		
		! write number of tracks	 
		call h5acreate_f( groupID, 'NTRACKS', H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_INTEGER, track_set%ntracks, dims, ierr)
		call h5aclose_f( attrID, ierr )

		! write number of iterations between dumps	 
		call h5acreate_f( groupID, 'NDUMP', H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_INTEGER, ndump, dims, ierr)
		call h5aclose_f( attrID, ierr )

		! write timestep size	 
		call h5acreate_f( groupID, 'DT', H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, dt, dims, ierr)
		call h5aclose_f( attrID, ierr )

		
		! write box size information
		call h5sclose_f( dataspaceID, ierr )
		dims(1) = xdim
		call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
    
    x_bnd(1,1) = 0.
    x_bnd(1,2) = 0.
    x_bnd(2,1) = nx
    x_bnd(2,2) = ny
		call h5acreate_f( groupID, 'XMIN', H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, x_bnd(1, :), dims, ierr)
		call h5aclose_f( attrID, ierr )
		call h5acreate_f( groupID, 'XMAX', H5T_NATIVE_DOUBLE, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_DOUBLE, x_bnd(2, :), dims, ierr)
		call h5aclose_f( attrID, ierr )
        
        
		! write periodic/moving window information
		bool(:) = 1
		call h5acreate_f( groupID, 'PERIODIC', H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
		call h5aclose_f( attrID, ierr )

		bool(:) = 0
		call h5acreate_f( groupID, 'MOVE C', H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
		call h5aclose_f( attrID, ierr )
		
		! write available quantities
		call h5sclose_f( dataspaceID, ierr )
		dims(1) = 4 + xdim + xdim
		call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
   
		txt(1) = 'n'
		txt(2) = 't'
		do i = 1, xdim
		  txt(2+i) = 'x'//(char(iachar('0')+i))
		enddo
		do i = 1, xdim
		  txt(2+xdim+i) = 'p'//(char(iachar('0')+i))
		enddo
		txt(2+xdim+xdim+1) = 'q'
		txt(2+xdim+xdim+2) = 'ene'
		call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
		call h5tset_size_f(typeID, p_max_txt_len, ierr)
		call h5acreate_f( groupID, 'QUANTS', typeID, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, typeID, txt, dims, ierr)
		call h5aclose_f( attrID, ierr )
		call h5tclose_f( typeID, ierr )
   
		txt(1) = ''
		txt(2) = '1/\omega_p'
		do i = 1, xdim
		  txt(2+i) = 'k_D'
		enddo
		do i = 1, vdim
		  txt(2+xdim+i) = 'm_e v_th'
		enddo
		txt(2+xdim+vdim+1) = 'e'
		txt(2+xdim+vdim+2) = 'm_e v_th^2'
		call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
		call h5tset_size_f(typeID, p_max_txt_len, ierr)
		call h5acreate_f( groupID, 'UNITS', typeID, dataspaceID, attrID, ierr )
		call h5awrite_f( attrID, typeID, txt, dims, ierr)
		call h5aclose_f( attrID, ierr )
		call h5tclose_f( typeID, ierr )
		
		call h5sclose_f( dataspaceID, ierr )
		
		! close root group
		call h5gclose_f( groupID, ierr )
		
		! loop through all tracks
		do i = 1, track_set%ntracks
		  call create_file( track_set%tracks(i), fileID, ndump )	 
		enddo
		
		! close file
		call h5fclose_f(fileID, ierr)
	 endif

   endif
   
end subroutine create_track_set_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! setup track set object 
!-------------------------------------------------------------------------------
subroutine setup_track_set( track_set, ndump, restart )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, intent(in) :: ndump ! number of iterations between each write
   logical, intent(in) :: restart
   
   character(len = 256) :: linebuffer
   integer :: i, tag1,tag2, ierr
   
   if ( restart ) then
   
      call restart_read( track_set, file_id_rst )
   
   else if ( ndump > 0 ) then
   
	  ! read tags
	  open( file_id_tem, file = track_set%file_tags, status = 'old', &
			access = 'sequential', form = 'formatted', iostat = ierr )
	  if ( ierr /= 0 ) then
		 write(*,*)'Unable to open particle tags file'
		 write(*,*)'"',trim(track_set%file_tags),'"'
!		 call abort_program()
	  endif
	  
	  ! skip comments until first valid value
	  do 
		 read( unit = file_id_tem, fmt= '(A)') linebuffer
		 if ( linebuffer(1:1) /= '!' ) then
			backspace file_id_tem
			exit
		 endif
	  enddo
	  
	  read( unit = file_id_tem, fmt= '(i)') track_set%ntracks
	  
	  ! allocate track objects
	  allocate( track_set%tracks(track_set%ntracks), &
				track_set%present(track_set%ntracks), &
				track_set%missing(track_set%ntracks), &
				stat = ierr )
	  if ( ierr /= 0 ) then
		write(*,*)'Unable to allocate track objects'
!		call abort_program( p_err_alloc )
	  endif
	  
	  ! set number of present/missing tracks
	  track_set%npresent = 0
	  track_set%nmissing = track_set%ntracks
	  
	  ! skip comments until first valid value
	  do 
		 read( unit = file_id_tem, fmt= '(A)') linebuffer
		 if ( linebuffer(1:1) /= '!' ) then
			backspace file_id_tem
			exit
		 endif
	  enddo
      
	  track_set%maxpoints = ndump/track_set%niter + 1
   
	  ! read tags and setup individual track objects   
	  do i = 1, track_set%ntracks
		 read( unit = file_id_tem, fmt= '(i,i)') tag1, tag2
		 
!		 print*,tag1,tag2
		 call setup( track_set%tracks(i), track_set%maxpoints, (/tag1,tag2/) )
		 
		 ! store all tracks in the missing list
		 track_set%missing(i) = i
	  enddo
	  
	  ! close tags file
	  close( file_id_tem )
	  
   endif
      
end subroutine setup_track_set

!-------------------------------------------------------------------------------
! setup track set object 
! taken from setup_track_set but modified so that only node 0 reads the tracks file
! and then broadcasts to other nodes. Jay Fahlen
!-------------------------------------------------------------------------------
subroutine setup_track_set_send( track_set, ndump, restart, idproc, nvp )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, intent(in) :: ndump ! number of iterations between each write
   logical, intent(in) :: restart
   integer :: idproc, nvp
   
   character(len = 256) :: linebuffer
   integer :: i, tag1,tag2, ierr
   integer,dimension(:),allocatable :: tagsarray
	 integer :: nproc, lgrp, mreal, mint, mcplx, mdouble
	 integer :: lworld
   integer, parameter :: lstat = 10
   integer,dimension(lstat) :: istatus
	 common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
	 integer :: temp_ntracks	 

   
   if ( restart ) then
   
      call restart_read( track_set, file_id_rst )
   
   else if ( ndump > 0 ) then

! determine the rank of the calling process in the communicator
!		call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
!		call MPI_COMM_SIZE(lworld,nvp,ierr)
   
   	if (idproc == 0) then
   	
			! read tags
			open( file_id_tem, file = track_set%file_tags, status = 'old', &
				access = 'sequential', form = 'formatted', iostat = ierr )
			if ( ierr /= 0 ) then
			 write(*,*)'Unable to open particle tags file'
			 write(*,*)'"',trim(track_set%file_tags),'"'
	!		 call abort_program()
			endif
			
			! skip comments until first valid value
			do 
			 read( unit = file_id_tem, fmt= '(A)') linebuffer
			 if ( linebuffer(1:1) /= '!' ) then
				backspace file_id_tem
				exit
			 endif
			enddo
			
			read( unit = file_id_tem, fmt= '(i)') track_set%ntracks
			
			allocate(tagsarray(2*track_set%ntracks))

			! allocate track objects
			allocate( track_set%tracks(track_set%ntracks), &
					track_set%present(track_set%ntracks), &
					track_set%missing(track_set%ntracks), &
					stat = ierr )
			if ( ierr /= 0 ) then
			write(*,*)'Unable to allocate track objects'
	!		call abort_program( p_err_alloc )
			endif
			
			! set number of present/missing tracks
			track_set%npresent = 0
			track_set%nmissing = track_set%ntracks
			
			! skip comments until first valid value
			do 
			 read( unit = file_id_tem, fmt= '(A)') linebuffer
			 if ( linebuffer(1:1) /= '!' ) then
				backspace file_id_tem
				exit
			 endif
			enddo
				
			track_set%maxpoints = ndump/track_set%niter + 1
		 
			! read tags and setup individual track objects   
			do i = 1, track_set%ntracks
			 read( unit = file_id_tem, fmt= '(i,i)') tag1, tag2
			 
			 tagsarray(2*i-1) = tag1
			 tagsarray(2*i) = tag2
			enddo
			
			if (nvp > 1) then
				do i = 1, nvp-1
					call mpi_send(track_set%ntracks,1,mint,i,0,lworld,ierr)
					call mpi_send(tagsarray,2*track_set%ntracks,mint,i,1,lworld,ierr)
				enddo
			endif

			do i = 1, track_set%ntracks
	!		 print*,tag1,tag2
			 call setup( track_set%tracks(i), track_set%maxpoints, tagsarray(2*i-1:2*i) )
			 
			 ! store all tracks in the missing list
			 track_set%missing(i) = i
			enddo
			
			deallocate(tagsarray)
			! close tags file
			close( file_id_tem )

			
		else !do other nodes
			call mpi_recv(temp_ntracks,1,mint,0,0,lworld,istatus,ierr)
			track_set%ntracks = temp_ntracks
			! set number of present/missing tracks
			track_set%npresent = 0
			track_set%nmissing = track_set%ntracks
			track_set%maxpoints = ndump/track_set%niter + 1
			
			allocate(tagsarray(2*track_set%ntracks))

			! allocate track objects
			allocate( track_set%tracks(track_set%ntracks), &
					track_set%present(track_set%ntracks), &
					track_set%missing(track_set%ntracks), &
					stat = ierr )
			if ( ierr /= 0 ) then
			write(*,*)'Unable to allocate track objects'
	!		call abort_program( p_err_alloc )
			endif
			
			call mpi_recv(tagsarray,2*track_set%ntracks,mint,0,1,lworld,istatus,ierr)
			
			do i = 1, track_set%ntracks
			 call setup( track_set%tracks(i), track_set%maxpoints, tagsarray(2*i-1:2*i) )
			 track_set%missing(i) = i
			enddo
			
			deallocate(tagsarray)
		endif
	  
   endif
      
end subroutine setup_track_set_send

!-------------------------------------------------------------------------------
! destroy track set object 
!-------------------------------------------------------------------------------
subroutine cleanup_track_set( track_set )
  
   implicit none

   type( t_track_set ), intent(inout) :: track_set

   integer :: i, ierr
      
   do i = 1, track_set%ntracks
     call cleanup( track_set%tracks(i) )
   enddo

   if ( track_set%ntracks > 0 ) then 
	  deallocate(track_set%tracks, &
				 track_set%present, &
				 track_set%missing, &
				 stat = ierr )
	  if ( ierr /= 0 ) then
		write(*,*)'Unable to deallocate track objects'
!		call abort_program( p_err_dealloc )
	  endif
   endif
   
end subroutine cleanup_track_set

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! write track restart information
!-------------------------------------------------------------------------------
subroutine restart_write_tracks( track_set, file_id )
  
   implicit none

   type( t_track_set ), intent(in) :: track_set
   integer, intent(in) :: file_id
   
   integer :: i, ierr 
   
   write(unit=file_id, iostat=ierr) p_track_rst_id
   if ( ierr/=0 ) then 
 write(*,*)'error writing restart data for track object.'
!	 call abort_program(p_err_rstwrt)
   endif

   write(unit=file_id, iostat=ierr) track_set%ntracks
   
   if ( track_set%ntracks > 0 ) then
      write(unit=file_id, iostat=ierr) track_set%maxpoints, &
             track_set%npresent, track_set%present, &
             track_set%nmissing, track_set%missing, &
             track_set%file_write
      
      do i = 1, track_set%ntracks
        call restart_write( track_set%tracks(i), file_id )
      enddo
      
   endif
   
   
end subroutine restart_write_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read track restart information
!-------------------------------------------------------------------------------
subroutine restart_read_tracks( track_set, file_id )
  
   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, intent(in) :: file_id
   
   character(len=len(p_track_rst_id)) :: rst_id
   integer :: i, ierr 
   
   read(unit=file_id, iostat=ierr) rst_id
   if ( ierr/=0 ) then 
	 write(*,*)'error reading restart data for track object.'
!	 call abort_program(p_err_rstrd)
   endif
 
   ! check if restart file is compatible
   if ( rst_id /= p_track_rst_id) then
	 write(*,*)'Corrupted restart file, or restart file '
	 write(*,*)'from incompatible binary (tracks)'
!	 call abort_program(p_err_rstrd)
   endif

   read(unit=file_id, iostat=ierr) track_set%ntracks
   
   if ( track_set%ntracks > 0 ) then
      ! allocate track objects
	  allocate( track_set%tracks(track_set%ntracks), &
				track_set%present(track_set%ntracks), &
				track_set%missing(track_set%ntracks), &
				stat = ierr )
	  if ( ierr /= 0 ) then
		write(*,*)'Unable to allocate track objects'
!		call abort_program( p_err_alloc )
	  endif
      
      ! read tracking lists
      read(unit=file_id, iostat=ierr) track_set%maxpoints, &
             track_set%npresent, track_set%present, &
             track_set%nmissing, track_set%missing, &
             track_set%file_write
      
      ! read individual tracks
      do i = 1, track_set%ntracks
        call restart_read( track_set%tracks(i), track_set%maxpoints, file_id )
      enddo
      
   endif
   
end subroutine restart_read_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! write track restart information
!-------------------------------------------------------------------------------
subroutine restart_write_single_track( track, file_id )
  
   implicit none

   type( t_track ), intent(in) :: track
   integer, intent(in) :: file_id
   
   integer :: ierr 
   
   write(unit=file_id, iostat=ierr) track%npoints, track%savedpoints, track%tag, &
     track%part_idx
   
   if ( track%npoints > 0 ) then
      write(unit=file_id, iostat=ierr) track%data( :, 1:track%npoints )
      write(unit=file_id, iostat=ierr) track%n(1:track%npoints )
   endif
   
   
end subroutine restart_write_single_track
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read track restart information
!-------------------------------------------------------------------------------
subroutine restart_read_single_track( track, max_points, file_id )
  
   implicit none

   type( t_track ), intent(inout) :: track
   integer, intent(in) :: max_points, file_id
   
   integer :: ierr 
   
   ! allocate required memory
   call setup_single_track( track, max_points, (/0,0/) )
   
   ! read information from file
   read(unit=file_id, iostat=ierr) track%npoints, track%savedpoints, track%tag, &
     track%part_idx
      
   if ( track%npoints > 0 ) then
      read(unit=file_id, iostat=ierr) track%data( :, 1:track%npoints )
      read(unit=file_id, iostat=ierr) track%n( 1:track%npoints )
   endif
   
   
end subroutine restart_read_single_track
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! index integer array using heapsort (store indexes of sorted array in and index array)
! this needs to be fixed to use 1: based arrays
! Taken from OSIRIS
!---------------------------------------------------------------------------------------------------
subroutine heapsort_index( list, idx ) 
  
  implicit none
  
  integer, dimension(0:), intent(in) :: list
  integer, dimension(0:), intent(inout) :: idx
    
  integer :: t, t_idx
  integer :: n, i, parent, child
  
  n = size(idx) 

  if ( n<2 ) then
    idx(0) = 1
    return
  endif
  
  i = n/2 
  
  do i = 0, n-1 
    idx(i) = i
  enddo
  
  do
    if ( i > 0 ) then
      i = i-1
      t_idx = idx(i)
      t = list(idx(i))
    else
      n = n-1
      if ( n == 0 ) exit
      t_idx = idx(n)
      t = list(idx(n))
      idx(n) = idx(0)
    endif
  
    parent = i
    child = i*2+1
    
    do 
      if (child >= n) exit
      if (child +1 < n ) then
         if (list(idx(child)) < list(idx(child+1))) child = child+1
      endif

      if ( t < list(idx(child)) ) then
        idx(parent) = idx(child)
        parent = child
        child = parent*2 + 1
      else
        exit
      endif
    enddo
    
    idx(parent) = t_idx
  enddo

  do i = 0, size(idx) -1 
    idx(i) = idx(i) + 1
  enddo
  
 
end subroutine heapsort_index
!---------------------------------------------------------------------------------------------------


subroutine assign_tags(part,npp,idproc)
	implicit none
	real,dimension(:,:,:) :: part
	integer,dimension(:) :: npp
	integer :: idproc
	
	integer :: i
	integer,dimension(2) :: tags
	real :: temp_real
	equivalence (temp_real,tags)

	tags(1) = idproc+1
	do i = 1, npp(1)
		tags(2) = i
		part(5,i,1) = temp_real
	enddo
end subroutine assign_tags
		

		
end module par_track_jf