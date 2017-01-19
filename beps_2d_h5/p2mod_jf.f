!-----------------------------------------------------------------------
! This module adds some functions written by Jay Fahlen to fix problems
! with the original p2mod.f file by Viktor Decyk.  Feb. 2008
		module p2d_jf
		use par_track_jf
		use pinit2d_jf
		use p0d
		use p2d, only: pwtimer, plsum, plmax, plscan, plbcast, writebf, &
			&readbf, wtimer, wrdata, rddata, dcomp, pcguard, pcguardp, &
			&paguard,paguardp,pncguardp, pnaguardp, pfmove, repart, fnoff,&
			&dblsin, dblcos, hafdbl, zdbl, pmoves, LINEAR, QUADRATIC, &
			&STANDARD, LOOKAHEAD, VECTOR
		implicit none
		private
		public :: pmove,write_jf
		public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
		public :: PPINIT, PPID, PPEXIT, HARTBEAT
		public :: pwtimer, plsum, plmax, plscan, plbcast, writebf, readbf
		public :: wtimer, wrdata, rddata
		public :: dcomp, pcguard, pcguardp, paguard, paguardp
		public :: pncguardp, pnaguardp, pfmove, repart, fnoff
		public :: dblsin, dblcos, hafdbl, zdbl, pmoves, sortp

		interface write_jf
			module procedure write_jf3d
			module procedure write_jf1d
		end interface
		
		interface pmove
			module procedure pmove_fix_buf_track
!			module procedure pmove_fix_buf
!			module procedure pmove_tracks
		end interface
		
		interface sortp
			module procedure ipsortp2y_track
			module procedure ipdsortp2y_track
		end interface
		
		interface plsum
			module procedure ipsum3
		end interface

		contains
		
		subroutine ipsum3(f)
! perform global sum of 3d real array
			implicit none
			real, dimension(:,:,:) :: f
! local data
			integer :: nxyp, nblok
!			real, dimension(size(f,1),size(f,2),size(f,3)) :: g
			real,dimension(:,:,:),allocatable :: g
!			real,dimension(:,:,:),allocatable,save :: g
			
!			if(.not. allocated(g)) then
!				allocate( g(size(f,1),size(f,2),size(f,3)))
!			endif
			allocate( g(size(f,1),size(f,2),size(f,3)))
			
			nxyp = size(f,1)*size(f,2)*size(f,3); nblok = 1
			call PSUM(f,g,nxyp,nblok)
		end subroutine ipsum3

		
		subroutine write_jf1d(f,nx,iunit,nrec,name,order)
! write a single 1 dim array to file, it does not collect from other nodes
         implicit none
         integer :: nx, iunit, nrec,j,idproc,ierr
         real, dimension(:) :: f
         character(len=*), optional :: name
         integer, optional :: order
	      integer :: nproc, lgrp, lstat=10, mreal, mint, mcplx, mdouble,&
	      	lworld
	      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
         endif
         nxv = size(f,1)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
	      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
!   	   call MPI_COMM_SIZE(lworld,nvp,ierr)
! node 0 receives messages from other nodes
			if (idproc.eq.0) then

				if (nrec.lt.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='replace')
					nrec = 1
	! open old file
				else if (nrec.eq.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='old')
				endif
				write (unit=iunit,rec=nrec) (f(j),j=1,nx)
				nrec = nrec + 1
			endif
		end subroutine write_jf1d
		
		subroutine write_jf3d(f,nx,kyp,iunit,nrec,name,order)
! write a single 3 dim array to file, it does not collect from other nodes
         implicit none
         integer :: nx, kyp, iunit, nrec,k,l,j,idproc,ierr
         real, dimension(:,:,:) :: f
         character(len=*), optional :: name
         integer, optional :: order
	      integer :: nproc, lgrp, lstat=10, mreal, mint, mcplx, mdouble,&
	      	lworld
	      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
! local data
         integer :: lrec, nxv, kypmx, nblok, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nx*kyp*lrec
         endif
         nxv = size(f,1); kypmx = size(f,2); nblok = size(f,3)
! this segment is used for mpi computers
! determine the rank of the calling process in the communicator
	      call MPI_COMM_RANK(lworld,idproc,ierr)
! determine the size of the group associated with a communicator
!   	   call MPI_COMM_SIZE(lworld,nvp,ierr)
! node 0 receives messages from other nodes
			if (idproc.eq.0) then
	
				if (nrec.lt.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='replace')
					nrec = 1
	! open old file
				else if (nrec.eq.0) then
					open(unit=iunit,file=name,form='unformatted',access='direct'&
						,recl=lrec,status='old')
				endif
				write (unit=iunit,rec=nrec) (((f(j,k,l),j=1,nx),k=1,kyp),l=1,&
					nblok)
				nrec = nrec + 1
			endif
		end subroutine write_jf3d

		subroutine pmove_fix_buf_track(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vt,&
     &ierr,tracks)
! particle manager: moves particles to appropriate processor
! non-uniform 1d partition boundaries in 2d code
! modified by Jay Fahlen to put buffers in heap rather than stack 
! because OSX limits stack size
			implicit none
         integer :: ny, kstrt, nvp, nbmax, vt, ierr
         real, dimension(2) :: tmove
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp
         type (t_track_set),optional :: tracks

! local data
         integer, dimension(1+vt*(size(part,2)-1),size(part,3)) :: maskp
!         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufl
!         real, dimension(size(part,1),nbmax,size(part,3)) :: sbufr
!         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufl
!         real, dimension(size(part,1),nbmax,size(part,3)) :: rbufr
!         integer, dimension(size(edges,1),size(edges,2)) :: jsl, jsr
!         integer, dimension(size(edges,1),size(edges,2)) :: jss
!         integer, dimension(2*nbmax,size(part,3)) :: ihole
         integer, dimension(7) :: info
         real, dimension(:,:,:),allocatable,save :: sbufl
         real, dimension(:,:,:),allocatable,save :: sbufr
         real, dimension(:,:,:),allocatable,save :: rbufl
         real, dimension(:,:,:),allocatable,save :: rbufr
         integer, dimension(:,:),allocatable,save :: jsl, jsr
         integer, dimension(:,:),allocatable,save :: jss
         integer, dimension(:,:),allocatable,save :: ihole
         integer :: idimp, npmax, nblok, idps, ntmax
         real :: tf, th
         double precision :: dtime

         if (.not. allocated(sbufl)) then
         	allocate(sbufl(size(part,1),nbmax,size(part,3)))
         	allocate(sbufr(size(part,1),nbmax,size(part,3)))
         	allocate(rbufl(size(part,1),nbmax,size(part,3)))
         	allocate(rbufr(size(part,1),nbmax,size(part,3)))
         	allocate(jsl(size(edges,1),size(edges,2)))
         	allocate(jsr(size(edges,1),size(edges,2)))
         	allocate(jss(size(edges,1),size(edges,2)))
         	allocate(ihole(2*nbmax,size(part,3)))
         endif
         	
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         idps = size(edges,1)
         ntmax = 2*nbmax
! initialize timer
         call wtimer(tf,dtime,-1) 
				 if (.not. present(tracks) ) then
						if (vt.eq.1) then
							call WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
			 &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,maskp&
			 &,info)
					 else
							call WPMOVE2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js&
			 &r,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info)
					 endif
				 else 
	   			if (vt.eq.1) then
            call WPXMOV2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,&
            &jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,&
            &maskp,info)
          else
            call WPMOVE2_track(part,edges,npp,sbufr,sbufl,rbufr,rbufl,&
            	&ihole,jsr,jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,&
            	&nbmax,ntmax,info,tracks)
          endif
         endif
				 	
! record time
         call wtimer(tf,dtime)
         tmove(1) = tmove(1) + tf
         tmove(2) = tmove(2) + th
         ierr = info(1)
      end subroutine pmove_fix_buf_track
      
      subroutine WPMOVE2_track(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,&
      	&jsl,jss,th,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info,tracks)
! wrapper function for particle manager
! info(6) = total number of particles on entry
! info(7) = total number of particles on exit
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
      type(t_track_set) :: tracks
      integer, dimension(:), allocatable :: temp_ihole_idx
      
!local data
      integer j, l, npr, nps, nter
      real tf
      double precision dtime
      integer ibflg, iwork
      dimension ibflg(2), iwork(2)
      do 10 j = 1, 7
      info(j) = 0
   10 continue
      th = 0.0
! debugging section: count total number of particles before move
      npr = 0
      do 20 l = 1, nblok
      npr = npr + npp(l)
   20 continue
! find outgoing particles
   30 nter = info(4)
      call PWTIMERA(-1,tf,dtime)
      call PMOVEH2(part,edges,npp,ihole,jss,idimp,npmax,nblok,idps,ntmax)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
      
! ihole has list of indexes of particles that are leaving this processor
! use it to update the missing particle lists for tracking

!			print*,jss(1,1)
!			print*,ihole

			allocate(temp_ihole_idx(jss(1,1)))
			temp_ihole_idx = ihole(1:jss(1,1),1)
!			print*,"a",jss(1,1)
			call missing_particles(tracks, temp_ihole_idx, jss(1,1), part)
!			print*,"b",jss(1,1)
			deallocate(temp_ihole_idx)
      
! send outgoing particles
!      print*,"before PMOVES2_track"
      call PMOVES2_track(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl,&
      	&jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info,tracks)
      	
!      print*,"after PMOVES2_track"
! particle overflow error
      if (info(1).gt.0) return
! buffer overflowed and more particles remain to be checked
      if (info(4).gt.nter) go to 30
! iteration overflow
      if (info(1).lt.0) go to 50
! debugging section: count total number of particles after move
      nps = 0
      do 40 l = 1, nblok
      nps = nps + npp(l)
   40 continue
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(6) = ibflg(1)
      info(7) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
! information
   50 nter = info(4)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end subroutine WPMOVE2_track
      
!-----------------------------------------------------------------------
      subroutine PMOVES2_track(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,&
      	&jsl,jss,ny,kstrt,nvp,idimp,npmax,nblok,idps,nbmax,ntmax,info,tracks)
! this subroutine moves particles into appropriate spatial regions
! periodic boundary conditions
! part(1,n,l) = position x of particle n in partition l
! part(2,n,l) = position y of particle n in partition l
! part(3,n,l) = velocity vx of particle n in partition l
! part(4,n,l) = velocity vy of particle n in partition l
! edges(1,l) = lower boundary of particle partition l
! edges(2,l) = upper boundary of particle partition l
! npp(l) = number of particles in partition l
! sbufl = buffer for particles being sent to lower processor
! sbufr = buffer for particles being sent to upper processor
! rbufl = buffer for particles being received from lower processor
! rbufr = buffer for particles being received from upper processor
! ihole = location of holes left in particle arrays
! jsl(idps,l) = number of particles going down in particle partition l
! jsr(idps,l) = number of particles going up in particle partition l
! jss(idps,l) = scratch array for particle partition l
! ny = system length in y direction
! kstrt = starting data block number
! nvp = number of real or virtual processors
! idimp = size of phase space = 4, 5 if particle tracking is on
! npmax = maximum number of particles in each partition
! nblok = number of particle partitions.
! idps = number of partition boundaries
! nbmax =  size of buffers for passing particles between processors
! ntmax =  size of hole array for particles leaving processors
! info = status information
! info(1) = ierr = (0,N) = (no,yes) error condition exists
! info(2) = maximum number of particles per processor
! info(3) = minimum number of particles per processor
! info(4) = maximum number of buffer overflows
! info(5) = maximum number of particle passes required
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, kstrt, nvp, idimp, npmax, nblok, idps, nbmax, ntmax
      dimension part(idimp,npmax,nblok)
      dimension edges(idps,nblok), npp(nblok)
      dimension sbufl(idimp,nbmax,nblok), sbufr(idimp,nbmax,nblok)
      dimension rbufl(idimp,nbmax,nblok), rbufr(idimp,nbmax,nblok)
      dimension jsl(idps,nblok), jsr(idps,nblok), jss(idps,nblok)
      dimension ihole(ntmax,nblok)
      dimension info(7)
! common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
! lstat = length of status array
      parameter(lstat=10)
! lgrp = current communicator
! mint = default datatype for integers
! mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      type(t_track_set) :: tracks
! local data
! iy = partitioned co-ordinate
      integer iy
      parameter(iy=2)
      integer ierr, l, ks, iter, nps, npt, kb, kl, kr, j, j1, j2, i
      integer nbsize, nter, mter, itermax
      integer msid, istatus
      integer ibflg, iwork
      real any, yt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      integer,dimension(2) :: int_tag
      real :: real_tag
      equivalence(real_tag,int_tag)
      
!      print*,"PMOVES2_track, idimp",idimp
      
      any = float(ny)
      ks = kstrt - 2
      nbsize = idimp*nbmax
      iter = 2
      nter = info(4)
      itermax = 2000
      mter = 0
! buffer outgoing particles
      do 50 l = 1, nblok
      kb = l + ks
      jsl(1,l) = 0
      jsr(1,l) = 0
! load particle buffers
      do 30 j = 1, jss(1,l)
      yt = part(iy,ihole(j,l),l)
! particles going down
      if (yt.lt.edges(1,l)) then
         if (kb.eq.0) yt = yt + any
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            do 10 i = 1, idimp
            sbufl(i,jsl(1,l),l) = part(i,ihole(j,l),l)
   10       continue
            sbufl(iy,jsl(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
!           go to 40
         endif
! particles going up
      else
         if ((kb+1).eq.nvp) yt = yt - any
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            do 20 i = 1, idimp
            sbufr(i,jsr(1,l),l) = part(i,ihole(j,l),l)
   20       continue
            sbufr(iy,jsr(1,l),l) = yt
            ihole(jsl(1,l)+jsr(1,l),l) = ihole(j,l)
         else
            jss(2,l) = 1
!           go to 40
         endif
      endif
   30 continue
   40 jss(1,l) = jsl(1,l) + jsr(1,l)
   50 continue
! check for full buffer condition
      nps = 0
      do 60 l = 1, nblok
      nps = max0(nps,jss(2,l))
   60 continue
      ibflg(3) = nps
! copy particle buffers
   70 iter = iter + 2
      mter = mter + 1
      do 120 l = 1, nblok
! get particles from below and above
      kr = l + ks + 2
      if (kr.gt.nvp) kr = kr - nvp
      kl = l + ks
      if (kl.lt.1) kl = kl + nvp
! this segment is used for shared memory computers
!     jsl(2,l) = jsr(1,kl)
!     do 90 j = 1, jsl(2,l)
!     do 80 i = 1, idimp
!     rbufl(i,j,l) = sbufr(i,j,kl)
!  80 continue
!  90 continue
!     jsr(2,l) = jsl(1,kr)
!     do 110 j = 1, jsr(2,l)
!     do 100 i = 1, idimp
!     rbufr(i,j,l) = sbufl(i,j,kr)
! 100 continue
! 110 continue
! this segment is used for mpi computers
! post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
! send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,l),mreal,kr-1,iter-1,lgrp,msid(3),ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,l),mreal,kl-1,iter,lgrp,msid(4),ierr)
! wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,l) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,l) = nps/idimp
  120 continue
! check if particles must be passed further
      nps = 0
      do 150 l = 1, nblok
! check if any particles coming from above belong here
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 130 j = 1, jsr(2,l)
      if (rbufr(iy,j,l).lt.edges(1,l)) jsl(1,l) = jsl(1,l) + 1
      if (rbufr(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
  130 continue
      if (jsr(1,l).ne.0) write (2,*) 'Info: particles returning up'
! check if any particles coming from below belong here
      do 140 j = 1, jsl(2,l)
      if (rbufl(iy,j,l).ge.edges(2,l)) jsr(1,l) = jsr(1,l) + 1
      if (rbufl(iy,j,l).lt.edges(1,l)) jss(2,l) = jss(2,l) + 1
  140 continue
      if (jss(2,l).ne.0) write (2,*) 'Info: particles returning down'
      jsl(1,l) = jsl(1,l) + jss(2,l)
      nps = max0(nps,jsl(1,l)+jsr(1,l))
  150 continue
      ibflg(2) = nps
! make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 260
! remove particles which do not belong here
      do 250 l = 1, nblok
      kb = l + ks
! first check particles coming from above
      jsl(1,l) = 0
      jsr(1,l) = 0
      jss(2,l) = 0
      do 190 j = 1, jsr(2,l)
      yt = rbufr(iy,j,l)
! particles going down
      if (yt.lt.edges(1,l)) then
         jsl(1,l) = jsl(1,l) + 1
         if (kb.eq.0) yt = yt + any
         rbufr(iy,j,l) = yt
         do 160 i = 1, idimp
         sbufl(i,jsl(1,l),l) = rbufr(i,j,l)
  160    continue
! particles going up, should not happen
      elseif (yt.ge.edges(2,l)) then
         jsr(1,l) = jsr(1,l) + 1
         if ((kb+1).eq.nvp) yt = yt - any
         rbufr(iy,j,l) = yt
         do 170 i = 1, idimp
         sbufr(i,jsr(1,l),l) = rbufr(i,j,l)
  170    continue
! particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 180 i = 1, idimp
         rbufr(i,jss(2,l),l) = rbufr(i,j,l)
  180    continue
      endif
  190 continue
      jsr(2,l) = jss(2,l)
! next check particles coming from below
      jss(2,l) = 0
      do 230 j = 1, jsl(2,l)
      yt = rbufl(iy,j,l)
! particles going up
      if (yt.ge.edges(2,l)) then
         if (jsr(1,l).lt.nbmax) then
            jsr(1,l) = jsr(1,l) + 1
            if ((kb+1).eq.nvp) yt = yt - any
            rbufl(iy,j,l) = yt
            do 200 i = 1, idimp
            sbufr(i,jsr(1,l),l) = rbufl(i,j,l)
  200       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
! particles going down, should not happen
      elseif (yt.lt.edges(1,l)) then
         if (jsl(1,l).lt.nbmax) then
            jsl(1,l) = jsl(1,l) + 1
            if (kb.eq.0) yt = yt + any
            rbufl(iy,j,l) = yt
            do 210 i = 1, idimp
            sbufl(i,jsl(1,l),l) = rbufl(i,j,l)
  210       continue
         else
            jss(2,l) = 2*npmax
            go to 240
         endif
! particles staying here
      else
         jss(2,l) = jss(2,l) + 1
         do 220 i = 1, idimp
         rbufl(i,jss(2,l),l) = rbufl(i,j,l)
  220    continue
      endif
  230 continue
  240 jsl(2,l) = jss(2,l)
  250 continue
! check if move would overflow particle array
  260 nps = 0
      npt = npmax
      do 270 l = 1, nblok
      jss(2,l) = npp(l) + jsl(2,l) + jsr(2,l) - jss(1,l)  
      nps = max0(nps,jss(2,l))
      npt = min0(npt,jss(2,l))
  270 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      do 360 l = 1, nblok
! distribute incoming particles from buffers
! distribute particles coming from below into holes
      jss(2,l) = min0(jss(1,l),jsl(2,l))
      do 290 j = 1, jss(2,l)
      do 280 i = 1, idimp
      
      part(i,ihole(j,l),l) = rbufl(i,j,l)

  280 continue
  
			real_tag = part(addtag_loc,ihole(j,l),l)
			if (int_tag(1) < 0) then
				call new_present_particles(tracks, int_tag(2), ihole(j,l))
			endif
!      call new_present_particles(tracks, int_tag, ihole(j,l))


  290 continue
      if (jss(1,l).gt.jsl(2,l)) then
         jss(2,l) = min0(jss(1,l)-jsl(2,l),jsr(2,l))
      else
         jss(2,l) = jsl(2,l) - jss(1,l)
      endif
      do 320 j = 1, jss(2,l)
! no more particles coming from below
! distribute particles coming from above into holes
      if (jss(1,l).gt.jsl(2,l)) then
         do 300 i = 1, idimp

         part(i,ihole(j+jsl(2,l),l),l) = rbufr(i,j,l)

  300    continue
  			 
  			 real_tag = part(addtag_loc,ihole(j+jsl(2,l),l),l)
				if (int_tag(1) < 0) then
					call new_present_particles(tracks, int_tag(2), ihole(j+jsl(2,l),l))
				endif
!      	 call new_present_particles(tracks, int_tag, ihole(j+jsl(2,l),l))

      else
! no more holes
! distribute remaining particles from below into bottom
         do 310 i = 1, idimp

         part(i,j+npp(l),l) = rbufl(i,j+jss(1,l),l)

  310    continue
  			 
  			 real_tag = part(addtag_loc,j+npp(l),l)
				if (int_tag(1) < 0) then
					call new_present_particles(tracks, int_tag(2), j+npp(l))
				endif

!      	 call new_present_particles(tracks, int_tag, j+npp(l))

      endif
  320 continue
      if (jss(1,l).le.jsl(2,l)) then
         npp(l) = npp(l) + (jsl(2,l) - jss(1,l))
         jss(1,l) = jsl(2,l)
      endif
      jss(2,l) = jss(1,l) - (jsl(2,l) + jsr(2,l))
      if (jss(2,l).gt.0) then
         jss(1,l) = (jsl(2,l) + jsr(2,l))
         jsr(2,l) = jss(2,l)
      else
         jss(1,l) = jss(1,l) - jsl(2,l)
         jsr(2,l) = -jss(2,l)
      endif
      do 350 j = 1, jsr(2,l)
! holes left over
! fill up remaining holes in particle array with particles from bottom
      if (jss(2,l).gt.0) then
         j1 = npp(l) - j + 1
         j2 = jss(1,l) + jss(2,l) - j + 1
         if (j1.gt.ihole(j2,l)) then
! move particle only if it is below current hole
            do 330 i = 1, idimp

            part(i,ihole(j2,l),l) = part(i,j1,l)

  330       continue
						
						real_tag = part(addtag_loc,ihole(j2,l),l)
						if (int_tag(1) < 0) then
							call change_index(tracks, int_tag(2), ihole(j2,l))
						endif

!						call change_index( tracks, j1, ihole(j2,l) )

         endif
      else
! no more holes
! distribute remaining particles from above into bottom
         do 340 i = 1, idimp
         part(i,j+npp(l),l) = rbufr(i,j+jss(1,l),l)
  340    continue
						real_tag = part(addtag_loc,j+npp(l),l)
						if (int_tag(1) < 0) then
							call new_present_particles(tracks, int_tag(2), j+npp(l))
						endif
!		      	call new_present_particles(tracks, int_tag, j+npp(l) )

      endif
  350 continue
      if (jss(2,l).gt.0) then
         npp(l) = npp(l) - jsr(2,l)
      else
         npp(l) = npp(l) + jsr(2,l)
      endif
      jss(1,l) = 0
  360 continue
! check if any particles have to be passed further
      info(5) = max0(info(5),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 70
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
! check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(4) = nter
      endif
      return
      end subroutine PMOVES2_track

!These two sort functions are normally in ppush2mod.f and ppush2lib.f
!Jay Fahlen modified these to handle particle tracking. They're in p2mod_jf.f
!because I haven't had a chance to move them into a new ppush2mod_jf.f file
      subroutine ipsortp2y_track(part,pt,ip,npp,noff,nyp,npic,tsort,tpush,i&
     &norder,tracks)
! sort particles by y grid using memory-conserving bin sort
! sort when current accumulated time, tsort(2), is greater than ideal
! accumulated time, tsort(4), plus last sorting time, tsort(5)
         implicit none
         integer :: inorder
         real :: tpush
         real, dimension(8) :: tsort
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: pt
         integer, dimension(:,:), pointer :: ip, npic
         integer, dimension(:), pointer :: npp, noff, nyp
         type(t_track_set) :: tracks
! local data
         integer :: idimp, npmax, nblok, nypm1, order
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); nypm1 = size(npic,1)
         order = QUADRATIC
!         if (present(inorder)) order = inorder
         order = inorder
! check if sortime exceeded
         tsort(1) = tsort(1) + 1.0
         if (tsort(1) >= min(tsort(7),tsort(8))) then
            tsort(8) = tsort(1)
            tsort(1) = 0.0
         endif
! accumulate push times
         if (tsort(1)==1.0) then
            tsort(2) = tpush
            tsort(3) = tpush
            tsort(4) = tpush
         else if (tsort(1) > 1.0) then
            tsort(2) = tsort(2) + tpush
!           tsort(2) = tsort(1)*tpush
            tsort(4) = tsort(4) + tsort(3)
         endif
! sort
!        if ((tsort(1)==0.0) .or. (tsort(2) > (tsort(4)+tsort(5)))) then
! disable autosort
         if (tsort(1)==0.0) then
! initialize timer
            call wtimer(tsort(5),dtime,-1)
						if (order==LINEAR) then
							 call PSORTP2YL(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,n&
		 &blok,nypm1)
						else
							 call PSORTP2Y(part,pt,ip,npic,npp,noff,nyp,idimp,npmax,nb&
		 &lok,nypm1)
						endif
						
			   		call update_indexes(tracks, part, npp)
						
! record time
            call wtimer(tsort(5),dtime)
            tsort(6) = tsort(6) + tsort(5)
!           tsort(5) = .5*tsort(5)
            if (tsort(1) > 0.0) tsort(8) = tsort(1)
            tsort(1) = 0.0
         endif
      end subroutine ipsortp2y_track
      
			 subroutine ipdsortp2y_track(parta,partb,npp,noff,nyp,npic,tsort,tpush&
	 &,inorder,tracks)
! sort particles by y grid using optimized bin sort
! sort when current accumulated time, tsort(2), is greater than ideal
! accumulated time, tsort(4), plus last sorting time, tsort(5)
			 implicit none
			 integer, optional :: inorder
			 real :: tpush
			 real, dimension(8) :: tsort
			 real, dimension(:,:,:), pointer :: parta, partb
			 integer, dimension(:), pointer :: npp, noff, nyp
			 integer, dimension(:,:), pointer :: npic
			 real, dimension(:,:,:), pointer :: part
			 type(t_track_set) :: tracks
! local data
			 integer :: idimp, npmax, nblok, nypm1, order
			 double precision :: dtime
			 idimp = size(parta,1); npmax = size(parta,2)
			 nblok = size(parta,3); nypm1 = size(npic,1)
			 order = QUADRATIC
			 if (present(inorder)) order = inorder
! check if sortime exceeded
			 tsort(1) = tsort(1) + 1.0
			 if (tsort(1) >= min(tsort(7),tsort(8))) then
					tsort(8) = tsort(1)
					tsort(1) = 0.0
			 endif
! accumulate push times
			 if (tsort(1)==1.0) then
					tsort(2) = tpush
					tsort(3) = tpush
					tsort(4) = tpush
			 else if (tsort(1) > 1.0) then
					tsort(2) = tsort(2) + tpush
!           tsort(2) = tsort(1)*tpush
					tsort(4) = tsort(4) + tsort(3)
			 endif
! sort
!        if ((tsort(1)==0.0) .or. (tsort(2) > (tsort(4)+tsort(5)))) then
! disable autosort
			 if (tsort(1)==0.0) then
! initialize timer
					call wtimer(tsort(5),dtime,-1)
					if (order==LINEAR) then
						 call PDSORTP2YL(parta,partb,npic,npp,noff,nyp,idimp,npmax&
	 &,nblok,nypm1)
					else
						 call PDSORTP2Y(parta,partb,npic,npp,noff,nyp,idimp,npmax,&
	 &nblok,nypm1)
					endif
					part => parta
					parta => partb
					partb => part

					call update_indexes(tracks, parta, npp)

! record time
					call wtimer(tsort(5),dtime)
					tsort(6) = tsort(6) + tsort(5)
!           tsort(5) = .5*tsort(5)
					if (tsort(1) > 0.0) tsort(8) = tsort(1)
					tsort(1) = 0.
			 endif
			 end subroutine ipdsortp2y_track
      
      end module p2d_jf
