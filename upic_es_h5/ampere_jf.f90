!-----------------------------------------------------------------------
!
      module ampere_jf
      use p0d, only: plsum
      implicit none
      private
      public :: ampere
      
      contains
      
      subroutine ampere(E,part,dt,np,npp,idimp,nblok)
      	integer :: np,idimp,nblok
      	integer,dimension(nblok) :: npp
      	real,dimension(2) :: E
      	real,dimension(:,:,:) :: part
      	real :: dt
         integer :: nproc, lgrp, mreal, mint, mcplx, mdouble
         integer :: lworld,idproc,nvp
         common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworl&
        &d
         !cur(1) is curx, cur(2) = cury
         real,dimension(2) :: cur,tempcur
         real :: Exold, Eyold
         integer :: i,ierr
         integer, parameter :: lstat = 10
         integer,dimension(lstat) :: istatus

! determine the rank of the calling process in the communicator
   	   call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
	      call MPI_COMM_SIZE(lworld,nvp,ierr)

			Exold = E(1)
			Eyold = E(2)
			
         cur(1) = 0.
         cur(2) = 0.
         tempcur = 0.
         
         !add each node's current up
         do i = 1,npp(1)
         	cur(1) = cur(1) + part(3,i,1)
         	cur(2) = cur(2) + part(4,i,1)
         enddo
         ! if node 0, then get current from other nodes and add to cur
!         if (idproc == 0) then
!         	do i = 1, nvp-1
!         		call MPI_RECV(tempcur,2,mreal,i,99,lworld,istatus,ierr)
!!         		cur(1) = cur(1) + tempcur(1)
!         		cur(2) = cur(2) + tempcur(2)
!         	enddo
!         else
!         	!All other nodes send their current
!         	call MPI_SEND(cur,2,mreal,0,99,lworld,ierr)
!         endif
!         print*,"cur",cur(1),cur(2)
			
			!Use Viktor's function that does log(nvp) sends
			call plsum(cur)
			cur = cur / real(np)
!			print*,cur(1),cur(2)
         E(1) = Exold + cur(1)*dt
         E(2) = Eyold + cur(2)*dt
         
         !Now send back to all nodes
!         if (idproc == 0) then
!         	do i = 1,nvp-1
!         		call MPI_SEND(E,2,mreal,i,99,lworld,ierr)
!         	enddo
!         else
!         	call MPI_RECV(E,2,mreal,i,99,lworld,istatus,ierr)
!         endif

      end subroutine ampere
         
      end module ampere_jf
         


      	
      	