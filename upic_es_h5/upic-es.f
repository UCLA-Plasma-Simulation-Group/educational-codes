!-----------------------------------------------------------------------
! * * * periodic 2d electrostatic particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves electrons and ions, with periodic electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms.  the only diagnostic is particle and field energy.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory computers
! update: may 30, 2006
      program pbeps2
      use pinit2d
!     use mprbpush2d
!     use mppush2d
      use prbpush2d
      use ppush2d
!     use mpfft2d
      use pfft2d
      use pfield2d
      use pfield2d_jf
      use pdiag2d
!     use mp2d
!      use p2d
			use pbpush2d, only: djpost
			use ppush2_jf
      use p2d_jf
			use par_track_jf
			use p0d, only:get_funit
      use mp0d, only: mpinit
      use diag_jf
      use pinit2d_jf
      use ampere_jf
      use ext_driver_jf
      use hdf_write_jf
      use hdf5
      implicit none
! idps = number of partition boundaries
! idimp = dimension of phase space = 4
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    2, idimp =   4, mshare =   0
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, vect = 0, ipbc = 1
      integer :: npxy, npxyb, np, npxyi, npxybi, npi
      integer :: nx, ny, nxh, nyh, nyv, nxe, nxeh
      integer :: nloop, nvp, nblok, npmax, npimax = 0, kyp, kxp, nypmx
      integer :: kyb, kxb, kxyb, kbmin, kblok, jbmin, jblok
      integer :: ngds, nxyh, nxhy, nx1, nypm1, nbmax
      integer :: idproc, id0, kstrt, itime, isign, isc, irc, ierr
      integer :: it, modesy2p
      integer :: ntasks
      integer :: i_ens
      real,dimension(2) :: E0xy
      real :: zero = 0.0, we = 0.0, wke = 0.0, wki = 0.0
      real :: tpush = 0.0, tdpost = 0.0
      real :: tpushi = 0.0, tdposti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, qi0, etx
      real :: vtxi, vtyi, vtdxi, vtdyi
      double precision :: dtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:), pointer :: fxye
      complex, dimension(:,:,:), pointer :: qt
      complex, dimension(:,:,:,:), pointer :: fxyt
      complex, dimension(:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer  :: edges
      integer, dimension(:), pointer :: nyp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:), pointer :: sfield
      complex, dimension(:,:,:), pointer :: sfieldt, pott
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
! for diagnostics JF
			real, dimension(:,:,:), pointer :: fvxx, fvyx, fvxy, fvyy,fvxvy
			real, dimension(:,:,:), pointer :: fvxxslice,fvyxslice,fvxyslice,fvyyslice
			real,dimension(:,:,:),pointer :: fv_xy
			integer :: i,ntfieldfunitx,ntfieldfunity
			real,dimension(:,:,:,:), pointer :: cu ! For current diagnostic
			real,dimension(:,:,:,:), pointer :: jE, jE_sum
			real,dimension(:,:,:), pointer :: kE
			complex,dimension(:,:,:,:),pointer :: cut,tempt !For current diagnostic
			complex,dimension(:,:,:),pointer :: phit
			real,dimension(:,:,:,:), pointer :: field_temp
			complex,dimension(:,:,:,:),pointer :: field_tempt !For field diagnostic			
			real :: field_ene
			real,dimension(:,:,:,:), pointer :: grad_phi,grad_phi_dt ! For field diagnostic
			real,dimension(:,:,:,:), pointer :: grad_phi2
			real,dimension(:,:,:,:), pointer :: dEdt, grad_phi_int ! For field diagnostic
			real,dimension(:,:,:,:), pointer :: dEdt_dt
			complex,dimension(:,:,:,:),pointer :: grad_phit !For field diagnostic			
			complex,dimension(:,:,:,:),pointer :: grad_phit2 !For field diagnostic			
			real,dimension(:,:,:),pointer :: div_ESPoynt,phi_old,div_ESPoynt_int !Needed for field and ESPoynt_int diag
			real,dimension(:,:,:),pointer :: div_ESPoynt_sum
			real,dimension(:,:,:,:), pointer :: vdotE_inst,vdotE_sum
			real,dimension(1) :: tot_dEne, tot_mag_dEne
			integer,dimension(1) :: tot_part_count
																			!tempt is like temp space
			real,dimension(:,:,:), pointer :: Ex_bin, Ey_bin
			real :: tdjpost = 0.0, wm = 0.0, tdcjpost = 0.0
			type (t_track_set) :: tracks
			! junk is an array used to pass the boundary conditions and moving window info
			! to the track create_file subroutine.  Since this code is always periodic
			! and doesn't have moving window, they are junk arrays.
		  logical, dimension(2) :: junk
      real, dimension(:,:,:,:), pointer :: ESPoynt,ESPoynt_tempE,ESPoynt_temp,ESPoynt_int
      real,dimension(3,4) :: Pflow
      real,dimension(3) :: totjE, totU, totPflow
      real,dimension(3) :: sub_first2_U=(/0.0,0.0,0.0/),sub_first_jE=(/0.0,0.0,0.0/),sub_first_div=(/0.0,0.0,0.0/)
      real,dimension(3,4) ::sub_first_P
      real,dimension(:),allocatable :: div_P_sumover_x,jE_sumover_x,Py_sumover_x,U_sumover_x,U_sumover_x_fromE
			character(len=5) :: tempname
			integer(KIND=8) :: check_div, nptot, temp8
			integer(KIND=4) :: h5error


      real :: invdt, num_par_cell,inv_num_par_cell,temp_we
      real :: write_tracks_time = 0., write_dtime = 0., temp_time=0.
      real :: add_tracks_time = 0., add_dtime = 0.
      real,dimension(1) :: curx, cury
      
!Debugging Poynting vector variables
!      real :: poynt_tot, ES_ene,ES_ene_h, vdote_tot,vdote_tot_part=0.,kin_ene=0.,kin_ene_tot=0.
 !     real :: kin_ene_cent=0.,kin_ene_0=0.,kin_ene_cent_0=0.,vdote_0=0.,ES_ene_old=0.,ES_ene_h_old
  !    real :: temp_we,tot_ene,tot_ene_old
   !   real :: ene_120 = 0.,poynt_120=0., vdote_120=0.,poynt=0.,allthree=0.,top,left,bot,right,wtotold=0.
    !  real :: kinene=0., kinene_old=0.,kinene_old2=0., tot_par=0.,wkeold=0.,div_poynt
     ! real :: line_on_cell,line_between_cells,line_ave_3,line_weight_3
      integer :: icnt,jcnt
      real,dimension(3) :: div_tot
		  
		  !for debugging
		  real :: real_tag
		  integer, dimension(2) :: int_tag
		  equivalence (real_tag,int_tag)
		  
!
		real, dimension(:,:), pointer :: wt
      integer, dimension(2) :: ktime
! for line diagnostics
		real,dimension(:,:), pointer :: lines
! wtot = total energy
      real, dimension(4) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, totfft = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
      real, dimension(2) :: tpart = 0.0, tparti = 0.0
      real, dimension(2) :: time
      real, dimension(8) :: tsort, tsorti
! msg = heartbeat array
      double precision, dimension(8) :: msg
      character(len=10) :: cdrun
      character(len=100) :: fname
      character(len=12) :: label
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
! nvp = number of real or virtual processors
! initialize for parallel processing
      call PPINIT(idproc,id0,nvp)

      call h5open_f(h5error)

      kstrt = idproc + 1
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      if (id0==0) then
! read namelist
				open(unit=8,file='pinput2',form='formatted',status='old')
				read (8,pinput2)
				read (8,pinput2_jf)
! text output file
				fname = 'energies.'//cdrun
				open(unit=76,file=trim(fname),form='formatted',status='replace')
				fname = 'poutput2.'//cdrun
				open(unit=18,file=trim(fname),form='formatted',status='replace')
				write (18,pinput2,iostat=irc)
				if (irc /= 0) write (18,*) 'pinput2 namelist not written'
				write (18,pinput2_jf,iostat=irc)
				if (irc /= 0) write (18,*) 'pinput2_jf namelist not written'
! open initial diagnostic metafile
				fname = 'pdiag2.init.'//cdrun
				open(unit=19,file=trim(fname),form='formatted',status='replace')
      endif
! broadcast namelist to other nodes
      call sendnml()
      call sendnml_jf()
! np = total number of electrons in simulation
      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
      nptot = npx; nptot = nptot * npy
! check to make sure that the number of particles is evenly divisible by the 
! number of processors.  This is required by the vdistr and fdistr functions below
      check_div = nptot / nvp
      if (check_div*nvp .ne. nptot) then
         write (2,*) 'np must be integer multiple of nvp'
         call PPEXIT
         stop
       endif
! npi = total number of ions in simulation
      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = ny/2
      nyv = ny + 2; nxe = nx + 4
! kyp = number of complex grids in each field partition in y direction
! nypmx = maximum size of particle partition, including guard cells.
      kyp = (ny - 1)/nvp + 1; nypmx = kyp + 3
! ngds = number of guard cells
      ngds = 3*((idps - 1)/2 + 1)
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871
         nxe = nx + 2; ; nypmx = kyp + 1
         ngds = (idps - 1)/2 + 1
      endif
      nxeh = nxe/2
! check if too many processors
      if (nvp > ny) then
         write (2,*) 'Too many processors requested, ny, nvp=', ny, nvp
         call PPEXIT
         stop
      endif      
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      if (id0==0) write (18,*) ntasks+1, ' processors used'
      if (dopt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! nblok = number of particle partitions
      nblok = 1 + mshare*(nvp - 1)
!      npmax = (np/nvp)*1.05 + 7000		!Original
			!Changed by JF to have percentage increase as a input parameter
			!also to prevent integer overflow for large number of particles
			npmax = npx/nvp
			npmax = (npmax*npy + npxyb/nvp)*(1. + npfac) + 7000	
!      if (movion==1) npimax = (npi/nvp)*1.05 + 7000
      if (movion==1) then
      	npimax = npxi/nvp
      	npimax = (npimax*npyi + npxybi/nvp)*(1.0 + npfac) + 7000
      endif
! kxp = number of complex grids in each field partition in x direction
      kxp = (nxh - 1)/nvp + 1
! kyb = number of processors in y
! kxb = number of processors in x
      kyb = ny/kyp; kxb = nxh/kxp
! kxyb = maximum(kxb,kyb)
      kxyb = max(kxb,kyb)
! kblok = number of field partitions in y direction
      kbmin = 1 + (1 - mshare)*(kxyb/kxb - 1)
      kblok = 1 + mshare*(ny/kyp - 1)
! jblok = number of field partitions in x direction
      jbmin = 1 + (1 - mshare)*(kxyb/kyb - 1)
      jblok = 1 + mshare*(nxh/kxp - 1)
! nxyh = maximum(nx,ny)/2
      nxyh = max(nx,ny)/2
! nxhy = maximum(nx/2,ny)
      nxhy = max(nxh,ny)
! dimensions for index and sorting arrays
      nx1 = nx + 1; nypm1 = kyp + 1
! nbmax = size of buffer for passing particles between processors
! Jay Fahlen rewrote to avoid integer overflow with npxy
      nbmax = 1 + (2*(npxy*vty + npxyb*vtdy) + 1.4*npxyb*abs(vdy))*dt/ny
      temp8 = nptot / ny
      nbmax = temp8
      nbmax = 1 + int(real(2*nbmax*vty)*dt)+int(real(2*npxyb*vtdy/ny)*dt)
      nbmax = nbmax + int(1.4*real(npxyb*abs(vdy)/ny)*dt)
!      nbmax = nbmax / 2
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
      endif
      do icnt=1, 3
       do jcnt=1,4
       	sub_first_P(icnt,jcnt) = 0.
       enddo
      enddo
      E0xy=0.
      if (nt_bin_E > 0) then
      	allocate(Ex_bin(bin_E_xrange(2)-bin_E_xrange(1),bin_E_yrange(2)-bin_E_yrange(1),bin_E))
      	allocate(Ey_bin(bin_E_xrange(2)-bin_E_xrange(1),bin_E_yrange(2)-bin_E_yrange(1),bin_E))
      endif
      	
      if (ntlines .ne. 0) then
				do i = 1, 5
					if (linepos(i) .ne. -1) then
						nlines = nlines + 1
					endif
				enddo
      	allocate(nlinerec(nlines))
      	allocate(lines(nlines,nx))
      	lines=0.
      	nlinerec=0.
      endif
      num_par_cell = real(npxy)/real(nx)/real(ny)
      inv_num_par_cell = 1. / num_par_cell
			if ( (nvdotE_part .ne. 0) .OR. (nvdotE_follow_part .ne. 0)) then
				set_part_array_vdotE = 1
			endif
			if (nt_write_U_sumover_x_fromE .ne. 0) then
				allocate(U_sumover_x_fromE(ny))
				U_sumover_x_fromE = 0.
			endif
			tot_part_count(1) = 0
			tot_dEne(1) = 0.
			tot_mag_dEne(1) = 0.
			if ((ntESPoynt_int .ne. 0) .or. (nt_write_ene_int .ne. 0) .or. (nt_write_Py_vs_y .ne. 0) &
			& .or. (nt_write_jE_onlytracked_sumover_x .ne. 0)) then
				allocate(ESPoynt_int(2,nxe,nypmx*kbmin,kblok),ESPoynt_temp(2,nxe,nypmx*kbmin,kblok))
				allocate(ESPoynt(2,nxe,nypmx*kbmin,kblok))
				allocate(Py_sumover_x(ny),div_P_sumover_x(ny),jE_sumover_x(ny),U_sumover_x(ny))
				ESPoynt_int = 0.
				ESPoynt = 0.
				ESPoynt_temp = 0.
				jE_sumover_x = 0.
				Py_sumover_x = 0.
				U_sumover_x = 0.
			endif
      if ((ntj .ne. 0) .OR. (ntESPoynt .ne. 0) .OR. (set_part_array_vdotE .ne. 0) &
      	& .OR. (ntESPoynt_int .ne. 0) .or. (ntvdotE .ne. 0) .or. (nt_write_ene_int .ne. 0) &
      	& .or. (nt_write_Py_vs_y .ne. 0)) then
      	allocate(cu(2,nxe,nypmx*kbmin,kblok))
	      allocate(cut(2,nyv,kxp,jblok),tempt(2,nyv,kxp,jblok))
	      allocate(phit(nyv,kxp,jblok))
      endif
      if ((ntESPoynt_int .ne. 0) .or. (ntvdotE .ne. 0) .or. (nt_write_ene_int .ne. 0) &
      	& .or. (nt_write_Py_vs_y .ne. 0)) then
      	allocate(jE(2,nxe,nypmx*kbmin,kblok))
      	allocate(kE(nxe,nypmx*kbmin,kblok))
      	allocate(jE_sum(2,nxe,nypmx*kbmin,kblok))
      	allocate(div_ESPoynt(nxe,nypmx*kbmin,kblok))
      	allocate(div_ESPoynt_sum(nxe,nypmx*kbmin,kblok))
      	allocate(div_ESPoynt_int(nxe,nypmx*kbmin,kblok))
				jE_sum = 0.
				div_ESPoynt_sum = 0.
      	div_ESPoynt_int = 0.
      	div_ESPoynt = 0.
      endif
      if (ntphsl_x .ne. 0) then
      	do i = 1, 5		!Only five for now, change later if you want
      		if (phsl_x_pos(i) .ne. -1) then
      			num_phsl_x = num_phsl_x + 1
      		endif
      	enddo
			allocate(fvxxslice(nx,2*nphbx+1,num_phsl_x))
			allocate(fvyxslice(nx,2*nphby+1,num_phsl_x))
      endif
      if (ntphsl_y .ne. 0) then
      	do i = 1, 5		!Only five for now, change later if you want
      		if (phsl_y_pos(i) .ne. -1) then
      			num_phsl_y = num_phsl_y + 1
      		endif
      	enddo
			allocate(fvxyslice(ny,2*nphbx+1,num_phsl_y))
			allocate(fvyyslice(ny,2*nphby+1,num_phsl_y))
      endif
!Obsolete with hdf use, but required for certain diagnostics
      !Added by JF to write unformatted file with diag info to ease
      !reading these parameters in IDL
      call write_jf_diag_file(id0,cdrun)
      !End of JF addition
      call mkdir_structure(idproc)

      if ((nphxx .ne. 0) .or. (nphxy .ne. 0) .or. (ntphsl_x .ne. 0)) then
      	nphbx = nphbx*2+1
      	if ((nphxx .ne. 0) .or. (ntphsl_x .ne. 0)) then
	      	allocate(fvxx(nx,nphbx,1))
	      endif
      	if (nphxy .ne. 0) then
	      	allocate(fvxy(ny,nphbx,1))
	      endif
      endif

      if ((nphyx .ne. 0) .or. (nphyy .ne. 0) .or. (ntphsl_x .ne. 0)) then
      	nphby = nphby*2+1
      	if (nphyx .ne. 0 .or. (ntphsl_x .ne. 0)) then
	      	allocate(fvyx(nx,nphby,1))
	      endif
      	if (nphyy .ne. 0) then
	      	allocate(fvyy(ny,nphby,1))
	      endif
	      
      endif
      
			if (nt_phase_vx_vy .ne. 0 .OR. nt_vx_vy_speed .ne. 0) then
      	if ((nphxx .eq. 0) .AND. (nphxy .eq. 0) .AND. (ntphsl_x .eq. 0)) then
					nphbx = nphbx*2+1
				endif
				if ((nphyx .eq. 0) .AND. (nphyy .eq. 0) .AND. (ntphsl_x .eq. 0)) then
					nphby = nphby*2+1
				endif
				allocate(fvxvy(nphbx,nphby,1))
			endif
      
      if (ntESPoynt .ne. 0) then
      	allocate(ESPoynt_tempE(2,nxe,nypmx*kbmin,kblok))
      endif
      if ((nttrack .ne. 0) .or. (ntraw .ne. 0)) then
      	addtag = 1
      else
      	nt_dump_track = 0
      endif
      if (ntphxy .ne. 0) then
      	if (nfvxy_xy_xrange .eq. 0) then
      		nfvxy_xy_xrange = nx
      	endif
      	if (nfvxy_xy_yrange .eq. 0) then
      		nfvxy_xy_yrange = ny
      	endif
      	nphb_xy = 2*nphb_xy+1
      	allocate(fv_xy(nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy))
      endif
      
      ! Since rise, flat, fall, and timerise are in units of wavelength
      ! or period, convert them to actual distance and time here.
      rise = real(nx) / real(wavemode) * rise
      flat = real(nx) / real(wavemode) * flat
      fall = real(nx) / real(wavemode) * fall
      
! diagnostic information needed by diagnostic nodes
! velocity diagnostics
      if (ntv > 0) then
         allocate(fv(2*nmv+2,2,nblok),fvm(2,2,nblok))
         if (id0==0) then
            fname = 'fv2.'//cdrun
            open(unit=10,file=trim(fname),form='formatted',status='unkno&
     &wn')
! write captions
            write (10,*) 'it vdx vdy vtx vty'
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,2,nblok),fvmi(2,2,nblok))
            if (id0==0) then
               fname = 'fvi2.'//cdrun
               open(unit=20,file=trim(fname),form='formatted',status='un&
     &known')
! write captions
                write (20,*) 'it vdxi vdyi vtxi vtyi'
            endif
         endif
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0) .or. (ntfield > 0) .or. (ntj > 0) &
      	& .or. (ntvdotE > 0) .or. (ntESPoynt .ne. 0) .or. (set_part_array_vdotE .ne. 0) &
      	& .or. (ntESPoynt_int .ne. 0) .or. (ntden .ne. 0) .or. (nt_write_grad_phi .ne. 0)) then
         allocate(sfield(nxe,nypmx*kbmin,kblok))
      endif
! potential diagnostics
! Modified by JF to create pott arrays if ntESPoynt is on too
      if ((ntp > 0) .OR. (ntESPoynt .ne. 0) .OR. (ntESPoynt_int .ne. 0)) then
         if (modesxp > nxh) modesxp = nxh
         if (modesyp > nyh) modesyp = nyh
         modesy2p = 2*modesyp - 1
         allocate(pott(modesy2p,min(modesxp,kxp),jblok))
         if (id0==0) then
            write (19,ppot2d,iostat=irc)
            if (irc /= 0) write (18,*) 'ppot2d namelist not written'
         endif
      endif
! energy diagnostics
      if (ntw > 0) allocate(wt((nloop-1)/ntw+1,4))
! open restart files
      if (nustrt /= 1) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='old&
     &')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='old&
     &')
         endif
      else if (ntr > 0) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='unk&
     &nown')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='unk&
     &nown')
         endif
      endif
! open graphics device
      call GROPEN
      call SETNPLT(nplot,irc)
      call STPALIT(idpal)
!
! diagnostic nodes have special processing
      if (idproc < 0) call diag2nodes
!
! switch for determining size of part array
! Three different options: add tags, keep initial velocites, vdotE kept with particles
! 1) Allow any three cases alone
! 2) add tags and keep initial velocities
! 3) add tags and vdotE kept with particles
! 4) do not allow all three together or vdotE with keep initial velocities (too much memory
!			and cases too complicated to be worthwhile)

! These always remain the same, it's for part(>4,n,l) that change
	! part(1,n,l) = position x of particle n in partition l
	! part(2,n,l) = position y of particle n in partition l
	! part(3,n,l) = velocity vx of particle n in partition l
	! part(4,n,l) = velocity vy of particle n in partition l
	
		if ( (nvdotE_part .ne. 0) .AND. (nvdotE_follow_part .ne. 0)) then
			print*,"Can't have both nvdotE_part and nvdotE_follow_part nonzero at once!."
			print*,"nvdotE_part will over ride nvdotE_follow_part."
			nvdotE_follow_part = 0
		endif
		if ( (nvdotE_part .ne. 0) .OR. (nvdotE_follow_part .ne. 0)) then
			set_part_array_vdotE = 1
		endif
		if ( (nvdotE_int .ne. 0) .or. (ntESPoynt_int .ne. 0) .or. (ntfield_ene_int .ne. 0) .or. &
			&	(nt_write_ene_int .ne. 0) .or. (nt_write_Py_vs_y .ne. 0) .or. (nt_dEdt .ne. 0)) then	
			n_jE_u_poynt = 1;
			allocate(grad_phi(2,nxe,nypmx*kbmin,kblok))
			allocate(grad_phi_dt(2,nxe,nypmx*kbmin,kblok))
			allocate(dEdt(2,nxe,nypmx*kbmin,kblok))  !Needed for du/dt, but not integrated du/dt
			allocate(grad_phi_int(2,nxe,nypmx*kbmin,kblok))
			allocate(dEdt_dt(2,nxe,nypmx*kbmin,kblok))
			allocate(grad_phit(2,nyv,kxp,jblok))
			grad_phi = 0.
			dEdt_dt = 0.
			grad_phi_dt = 0.
			dEdt = 0.
			grad_phi_int = 0.
			grad_phit = 0.
		endif

		if (nt_write_grad_phi .ne. 0) then
			allocate(grad_phi2(2,nxe,nypmx*kbmin,kblok))
			allocate(grad_phit2(2,nyv,kxp,jblok))
			grad_phi2 = 0.
			grad_phit2 = 0.
		endif

		if (nt_write_field_ene_in_driver_region .ne. 0) then
			allocate(field_tempt(2,nyv,kxp,jblok))
			allocate(field_temp(2,nxe,nypmx*kbmin,kblok))
			if (idproc == 0) then
				fname = 'field_ene_in_driver_region.'//cdrun
				open(unit=77,file=trim(fname),form='formatted',status='replace')
			endif
		endif

		if (nt_write_ene_int .ne. 0) then
			if (idproc == 0) then
				fname = 'poynting_flow_enes.'//cdrun
				open(unit=78,file=trim(fname),form='formatted',status='replace')
				write (78,*) "Time bot_totU bot_totjE bot_P_b bot_P_l bot_P_t bot_P_r bot_P_tot div_bot bot_tot bot_tot_div ", &
									&	      "mid_totU mid_totjE mid_P_b mid_P_l mid_P_t mid_P_r mid_P_tot div_mid mid_tot mid_tot_div " , &
									&				"top_totU top_totjE top_P_b top_P_l top_P_t top_P_r top_P_tot div_top top_tot top_tot_div "

			endif
		endif
		
		if ((nt_write_jE_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/jE_sumover_x_file'
				nt_write_jE_sumover_x_funit = get_funit(20)
				open(unit=nt_write_jE_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_jE_sumover_x_funit) ny,nt_write_jE_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_jE_onlytracked_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/jE_onlytracked_sumover_x_file'
				nt_write_jE_onlytracked_sumover_x_funit = get_funit(20)
				open(unit=nt_write_jE_onlytracked_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_jE_onlytracked_sumover_x_funit) ny,nt_write_jE_onlytracked_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_kE_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/kE_sumover_x_file'
				nt_write_kE_sumover_x_funit = get_funit(20)
				open(unit=nt_write_kE_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_kE_sumover_x_funit) ny,nt_write_kE_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_U_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/U_sumover_x_file'
				nt_write_U_sumover_x_funit = get_funit(20)
				open(unit=nt_write_U_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_U_sumover_x_funit) ny,nt_write_U_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_U_sumover_x_fromE .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/U_sumover_x_fromE'
				nt_write_U_sumover_x_fromE_funit = get_funit(20)
				open(unit=nt_write_U_sumover_x_fromE_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_U_sumover_x_fromE_funit) ny,nt_write_U_sumover_x_fromE,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_Ux_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/Ux_sumover_x_file'
				nt_write_Ux_sumover_x_funit = get_funit(20)
				open(unit=nt_write_Ux_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_Ux_sumover_x_funit) ny,nt_write_Ux_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_div_P_vs_y .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/div_P_sumover_x_file'
				nt_write_div_P_vs_y_funit = get_funit(20)
				open(unit=nt_write_div_P_vs_y_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_div_P_vs_y_funit) ny,nt_write_div_P_vs_y,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		if ((nt_write_Py_sumover_x .ne. 0)) then
			if (idproc==0) then
				fname = './DIAG/Py_sumover_x_file'
				nt_write_Py_sumover_x_funit = get_funit(20)
				open(unit=nt_write_Py_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
				write (nt_write_Py_sumover_x_funit) ny,nt_write_Py_sumover_x,dt,(timerise+timeflat+timefall),tend
			endif
		endif
		
    if ( (addtag == 0) .AND. (keep_init_vel == 0) .AND. (set_part_array_vdotE == 0) ) then
    	allocate(part(idimp,npmax,nblok))
    endif
    if ( (addtag == 1) .AND. (keep_init_vel == 0) .AND. (set_part_array_vdotE == 0) ) then
    	allocate(part(idimp+1,npmax,nblok))
    	addtag_loc = 5
    endif
    if ( (addtag == 0) .AND. (keep_init_vel .ne. 0 .OR. nt_through_wave .ne. 0) .AND. &
    &(set_part_array_vdotE == 0) ) then
    	if (keep_init_vel .ne. 0 .AND. nt_through_wave .ne. 0) then
    		print*,"Cannot have keep_init_vel and nt_through_wave nonzero simultaneously. Exiting."
    		stop
    	endif
    	allocate(part(idimp+2,npmax,nblok))
    	kivx_loc = 5
    	kivy_loc = 6
    endif
    if ( (addtag == 0) .AND. (keep_init_vel == 0) .AND. (set_part_array_vdotE .ne. 0) ) then
    	allocate(part(idimp+2,npmax,nblok))
     	allocate(vdotE_inst(2,nxe,nypmx*kbmin,kblok))
     	allocate(vdotE_sum(2,nxe,nypmx*kbmin,kblok))
     	vdotE_inst = 0.
     	vdotE_sum = 0.
    	vEx_loc = 5
    	vEy_loc = 6
    endif
    if ( (addtag == 1) .AND. (keep_init_vel .ne. 0) .AND. (set_part_array_vdotE == 0) ) then
    	allocate(part(idimp+3,npmax,nblok))
    	addtag_loc = 5
    	kivx_loc = 6
    	kivy_loc = 7
    endif
    if ( (addtag == 1) .AND. (keep_init_vel == 0) .AND. (set_part_array_vdotE .ne. 0) ) then
    	allocate(part(idimp+3,npmax,nblok))
     	allocate(vdotE_inst(2,nxe,nypmx*kbmin,kblok))
     	allocate(vdotE_sum(2,nxe,nypmx*kbmin,kblok))
     	vdotE_inst = 0.
     	vdotE_sum = 0.
    	addtag_loc = 5
    	vEx_loc = 6
    	vEy_loc = 7
    endif
    if ( (keep_init_vel .ne. 0) .AND. (set_part_array_vdotE .ne. 0) ) then
    	print*,"Cannot have keep_init_vel and nvdotE_part nonzero at the same time. Exiting."
    	stop
    endif
    	

! in real space, qe(j+1,k,l) = charge density at grid point (j,kk)
! in real space, fxye(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk)
! in other words, fxye are the convolutions of the electric field
! over the particle shape, where kk = k + noff(l) - 1
      allocate(qe(nxe,nypmx*kbmin,kblok),fxye(2,nxe,nypmx*kbmin,kblok))
! qt(k,j,l) = complex charge density for fourier mode jj-1,k-1
! fxyt(1,k,j,l) = x component of force/charge for fourier mode jj-1,k-1
! fxyt(2,k,j,l) = y component of force/charge for fourier mode jj-1,k-1
! where jj = j + kxp*(l - 1)
      allocate(qt(nyv,kxp,jblok),fxyt(2,nyv,kxp,jblok))
! ffc = form factor array for poisson solver
      allocate(ffc(nyh,kxp,jblok))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
! edges(1,l) = lower boundary of particle partition l
! edges(2,l) = upper boundary of particle partition l
      allocate(edges(idps,nblok))
! nyp(l) = number of primary gridpoints in particle partition l.
! noff(l) = lowermost global gridpoint in particle partition l.
      allocate(nyp(nblok),noff(nblok))
! npp(l) = number of particles in partition l
! nps(l) = starting address of particles in partition l
      allocate(npp(nblok),nps(nblok))
! sorting arrays
      allocate(pt(max(npmax,npimax),nblok))
      allocate(ip(max(npmax,npimax),nblok),npic(nypm1,nblok))

      tsort = 0.0
      tsort(1) = -1.0
      tsort(7) = real(sortime); tsort(8) = tsort(7)
      if (movion==1) then
         tsorti = 0.0
         tsorti(1) = -1.0
         tsorti(7) = real(sortimi); tsorti(8) = tsorti(7)
      endif   

! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      itime = 0
      qbme = qme
      affp = float(nx)*float(ny)/float(np)
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
      endif
! diagnostics
! velocity diagnostics
      if (ntv > 0) then
         fv(1,:,:) = 8.*max(vtx,vty)
         if (movion==1) fvi(1,:,:) = 8.*max(vtxi,vtyi)
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0) .or. (ntESPoynt .ne. 0) .OR. (ntESPoynt_int .ne. 0) .OR. &
      	& (nt_write_field_ene_in_driver_region .ne. 0) .or. (ntden .ne. 0) &
      	& .or. (nt_write_grad_phi .ne. 0)) then
         allocate(sfieldt(nyv,kxp,jblok))
      endif
! calculate partition variables
      call dcomp(edges,nyp,noff,ny,kstrt,nvp,inorder)
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,kstrt)  
! allocate background charge density
      if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,nblok),nppi(nblok))
         nullify(parti2)
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
!        if (npxy > 0) call distr(part,edges,npp,nps,vtx,vty,vx0,vy0,npx&
!    &,npy,nx,ny,ipbc)
         if (nptot > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npx,npy,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            do i_ens=1,nensemble
                call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)
                nps = 1
                npp = 0
            end do
            call vdistr(part,npp,nps,vtx,vty,vx0,vy0,npx,npy,kstrt,nvp)

         endif
         
! beam electrons
         nps = npp + 1
!        if (npxyb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vdx,vdy,&
!    &npxb,npyb,nx,ny,ipbc)
         if (npxyb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,npxb,npyb,nx,ny,kstrt,nvp,ipbc,ndprof,nsrand)
            call vdistr(part,npp,nps,vtdx,vtdy,vdx,vdy,npxb,npyb,kstrt,n&
     &vp)
         endif
! move electrons into appropriate spatial regions
				 if (nttrack .eq. 0) then
	         call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
	       else
	       	 call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr,tracks)
	       endif
!         call pmove_jf(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! initialize ions
         if (movion==1) then
            nps = 1
            nppi = 0
! background ions
!           if (npxyi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vxi&
!    &0,vyi0,npxi,npyi,nx,ny,ipbc)
            if (npxyi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxi,npyi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vxi0,vyi0,npxi,npyi,&
     &kstrt,nvp)
            endif
! beam ions
            nps = nppi + 1
!           if (npxybi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi,&
!    &vdxi,vdyi,npxbi,npybi,nx,ny,ipbc)
            if (npxybi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,npxbi,npybi,nx,ny,kstrt,nvp,ipbc,ndprofi,nsrandi)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vdxi,vdyi,npxbi,np&
     &ybi,kstrt,nvp)
            endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ie&
     &rr)
!            call pmove_jf(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ie&
!     &rr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! initialize background charge density
				if (movion==0) then
            qi0 = -qme/affp
            call sguard(qi,nyp,zero,nx,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! debug
           call sguard(qi,nyp,qi0,nx,inorder)
!           qi = 0.

! freeze the ions now
         else if ((movion==1).and.(itime==ionoff)) then
            allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
            call sguard(qi,nyp,zero,nx,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! restart
      else
! determine most recent restart file
         if (id0==0) then
            read (16,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            read (17,iostat=ierr) ktime(2)
            if (ierr /= 0) ktime(2) = -1
            if (ktime(1) > ktime(2)) then
               ktime(2) = 16
            else
               ktime(1) = ktime(2)
               ktime(2) = 17
            endif
         endif
         call plbcast(ktime)
         itime = ktime(1)
         if (itime < 0) go to 400
! read restart file
         it = ktime(2)
         call rddata(part,npp,it,ierr)
         if (ierr /= 0) go to 400
         if (movion==1) then
            call rddata(parti,nppi,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (movion==0) then
            call rddata(qi,nvp,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (ntw > 0) then
            call rddata(wt,1,it,ierr)
            if (ierr /= 0) go to 400
            call plbcast(wt)
         endif
         if (ntp > 0) then
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               irc = 0
               fname = 'ppotk2.'//cdrun
               call writebf(pott,modesxp,modesy2p,kxp,11,irc,trim(fname))
            endif
            call plbcast(ktime)
            nprec = ktime(1)
            if (nprec< 0) go to 400
         endif
         if (id0==0) then
            read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            rewind it
         endif
         call plbcast(ktime)
         irc = ktime(1)
         if (irc==itime) go to 490
! handle error
  400    if (id0==0) write (18,*) 'Restart Error'
         go to 3000
      endif
     
		if (addtag .ne. 0) then
			call assign_tags(part,npp,idproc)
		endif
		if ( (raw_dump_tagged .ne. 0) .and. (nttrack .eq. 0) ) then
			call setup_tagged_particles_send(part, idproc, nvp )
		endif

		if (nttrack .ne. 0) then
			! store particle tracking data
			tracks%niter = nttrack
			!the name of the tag file is hard coded in for now
			!tracks%file_tags = file_tags

	!      	call setup(tracks, nt_dump_track, .false.) ! No restart for now!
			call setup(tracks, nt_dump_track, .false., part, idproc, nvp) ! No restart for now!

			if (track_no_h5 .eq. 0) then
				!init hdf 5
				call h5open_f(ierr)
				call create_file(tracks,'elec',nt_dump_track,nx,ny,dt,junk,junk,idproc)
			else
				call create_file_no_h5(tracks,'elec',nt_dump_track,nx,ny,dt,tend,junk,junk,idproc)
			endif
		endif
		
		if (keep_init_vel .ne. 0) then
			call store_init_vel(part,npp)
		endif
		if (nt_through_wave > 0) then
			call init_through_wave(part,npp)
		endif

! GET RID OF THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********************************************************************************
! This is for initializing two particles that are tracked to watch the collisions.
!		do i = 1, npp(1)
!			real_tag = part(5,i,1)
!			if ((int_tag(1) == 1) .AND. (int_tag(2) == 1)) then
!				part(1,i,1) = 16.
!				part(2,i,1) = 33.
!				part(3,i,1) = .707
!				part(4,i,1) = .707
!			endif
!			if ((int_tag(1) == 1) .AND. (int_tag(2) == 2)) then
!				part(1,i,1) = 48.
!				part(2,i,1) = 33.
!				part(3,i,1) = .707
!				part(4,i,1) = .707
!			endif
!		enddo


!		part(1,1,1) = 15.
!		part(2,1,1) = 15.
!		part(3,1,1) = 5.
!		part(4,1,1) = 0.



!*********************************************************************************
! record time
  490 call pwtimer(time,dtime)
! send initial CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      call HARTBEAT(msg,2)
      if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= itime) go to 2000  		
! send time step to diagnostic nodes
      msg(1) = itime
      call HARTBEAT(msg,1)
!      if (id0==0) write (18,991) itime
      write (label,991) itime
      call LOGNAME(label)
! initialize charge density to background
      call sguard(qe,nyp,zero,nx,inorder)
! deposit charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
            sfield = -qe
! add guard cells for density in x
            call aguard(sfield,nyp,nx,inorder)
! add guard cells for density in y
            call paguard(sfield,kstrt,nvp,nx,kyp,ngds)
! transform density to fourier space
            isign = -1
            call fft(sfield,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,&
     &inorder)
! calculate smoothing in fourier space
            call pois(qt,sfieldt,ffc,nx,ny,kstrt)
! transform density to real space
            isign = 1
            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)
! display density
            call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle,nx&
     &,ny,irc,inorder)
            if (irc==1) go to 2000
         endif
      endif

! JF's density diagnostic
      if (ntden > 0) then
         it = itime/ntden
         if (itime==ntden*it) then
            sfield = -qe
! add guard cells for density in x
            call aguard(sfield,nyp,nx,inorder)
! add guard cells for density in y
            call paguard(sfield,kstrt,nvp,nx,kyp,ngds)
! transform density to fourier space
            isign = -1
!            call fft(sfield,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,&
!          & inorder)
! calculate smoothing in fourier space
!            call pois(qt,sfieldt,ffc,nx,ny,kstrt)
! transform density to real space
            isign = 1
!            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
!     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)

						fname = './DIAG/Den/'//'ne.'//cdrun
						call writef(sfield,nxe,nypmx,itime,itime*dt,DEN,trim(fname),inorder)
            
         endif
      endif
! add ion density
      if (movion==1) then
         call dpost(parti,qe,qmi,nppi,noff,tdposti,inorder,dopt)
! Ion Density Diagnostics
         if(itime==ntden*it .AND. ntden > 0) then
           sfield = 0.0
           call dpost(parti,sfield,qmi,nppi,noff,tdposti,inorder,dopt)
!           sfield = qi
! add guard cells for density in x
            call aguard(sfield,nyp,nx,inorder)
! add guard cells for density in y
            call paguard(sfield,kstrt,nvp,nx,kyp,ngds)
! transform density to fourier space
            isign = -1
            call fft(sfield,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,&
     &inorder)
! calculate smoothing in fourier space
            call pois(qt,sfieldt,ffc,nx,ny,kstrt)
! transform density to real space
            isign = 1
            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)

						fname = './DIAG/IDen/'//'ni.'//cdrun
						call writef(sfield,nxe,nypmx,itime,itime*dt,DEN,trim(fname),inorder)  
		  endif
! Ion Density Diagnostics         
      else
         qe = qe + qi
      endif
! add guard cells for density in x
      call aguard(qe,nyp,nx,inorder)
! add guard cells for density in y
      call paguard(qe,kstrt,nvp,nx,kyp,ngds)
! freeze the ions
      if ((movion==1).and.(itime==ionoff)) then
         allocate(qi(nxe,nypmx*kbmin,kblok))
! initialize ion charge density to zero
         call sguard(qi,nyp,zero,nx,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
         deallocate(parti,nppi)
         movion = 0
      endif
! velocity diagnostic
      if (ntv > 0) then
         it = itime/ntv
         if (itime==ntv*it) then
! calculate electron distribution function and moments
            call vdist(part,fv,fvm,npp,nmv,2)
            call plsum(fv(:,:,1))
            fv(1,:,:) = 8.*max(vtx,vty)
! display velocity distributions
            call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)
            if (irc==1) go to 2000
! print out velocity moments
            if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
            if (movion==1) then
! calculate ion distribution function and moments
               call vdist(parti,fvi,fvmi,nppi,nmv,2)
               call plsum(fvi(:,:,1))
               fvi(1,:,:) = 8.*max(vtxi,vtyi)
! display velocity distributions
               call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)
               if (irc==1) go to 2000
! print out velocity moments
               if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
            endif
         endif
      endif
! phase space diagnostic
      if (nts > 0) then
         it = itime/nts
         if (itime==nts*it) then
            isc = 999
! plot electrons x versus y
            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&
     &1,2,irc)
            if (irc==1) go to 2000
! plot electrons vx versus x
            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&
     &3,1,irc)
            if (irc==1) go to 2000
! plot electrons vy versus y
            call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,ny,&
     &4,2,irc)
            if (irc==1) go to 2000
            if (movion==1) then
! plot ions x versus y
               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&
     &1,2,irc)
               if (irc==1) go to 2000
! plot ions vx versus x
               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&
     &3,1,irc)
               if (irc==1) go to 2000
! plot ions vy versus y
               call grasp(parti,nppi,' ION PHASE SPACE',itime,isc,nx,ny,&
     &4,2,irc)
               if (irc==1) go to 2000
            endif
         endif
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,qt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
! potential diagnostic

      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! calculate potential in fourier space
            call pois(qt,sfieldt,ffc,we,nx,ny,kstrt)
! store selected fourier modes
            call gtmodes(sfieldt,pott,nx,ny,modesxp,modesyp,kstrt)
! transform potential to real space
            isign = 1
            call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt&
     &,kyp,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvp,kyp,inorder)
            call cguard(sfield,nyp,nx,inorder)
! display potential
            call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle,nx&
     &,ny,irc,inorder)
     
					fname = './DIAG/pot/pot.'//cdrun
					call writef(sfield,nxe,nypmx,itime,itime*dt,POT,trim(fname),inorder)
     
            if (irc==1) go to 2000

!Added by JF to write plain electron density   

					! fname = './DIAG/den/'//'pden.'//cdrun
					! call writef(sfield,nxe,nypmx,itime,itime*dt,DEN,trim(fname),inorder)

!           g(:,1) = sfield(nxh,:,1)
!           call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorder)
!           if (irc==1) go to 2000
! write diagnostic output
            if (nprec==0) then
               nprec = -1
               fname = 'ppotk2.'//cdrun
               call writebf(pott,modesxp,modesy2p,kxp,11,nprec,trim(fname))
            else
               call writebf(pott,modesxp,modesy2p,kxp,11,nprec)
            endif
         endif
      endif

! calculate force/charge in fourier space
      call pois(qt,fxyt,ffc,we,nx,ny,kstrt)

! transform force/charge to real space
      isign = 1
      call fft(fxye,fxyt,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!     call fftn(fxye,fxyt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inord&
!    &er)

		if (turn_off_self_con .ne. 0) then
			fxye = 0.
		endif

! call external traveling wave driver
		if (driver_select == 1) then
			call plane_wave(fxye,real(itime)*dt,nx,nxe,nvp,idproc)
		else if (driver_select == 2) then
			call plane_finite_wavelen(fxye,real(itime)*dt,nx,nxe,nvp,idproc)
		else if (driver_select == 3) then
			call gauss_tran_finite_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 4) then
			call angle_wave(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 5) then
			call cross_wave(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 6) then
			call rect_tran_finite_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 7) then
			call gauss_tran_per_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 8) then
			call doub_gauss_tran_finite_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 9) then
			call gauss_tran_per_wavelen_circle_wavefronts(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 10) then
			call gauss_tran_per_wavelen_pois(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 11) then
			call rect_tran_per_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 12) then
			call supergauss_tran_per_wavelen(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		else if (driver_select == 13) then
			call plane_finite_wavelen_supergauss(fxye,real(itime)*dt,nx,nxe,nvp,idproc)
		else if (driver_select == 14) then
			call plane_wave_in_y(fxye,real(itime)*dt,nx,nxe,ny,nypmx,nvp,idproc)
		endif

! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxye,kstrt,nvp,kyp,inorder)
      call cguard(fxye,nyp,nx,inorder)
! external pump
      if ((itpon > 0).and.(itime >= itpon)) then
         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))
         fxye(1,:,:,:) = fxye(1,:,:,:) + etx
      endif
! calculate 0-mode efield from ampere's law
			if (ampere_k0 .ne. 0) then
	      call ampere(E0xy,part,dt,np,npp,idimp,nblok)
	    endif

!			curx(1) = 0.
!			cury(1) = 0.
!			do i=1,npp(1)
!				curx(1) = curx(1) + 0.5*part(3,i,1)*part(3,i,1)
!				cury(1) = cury(1) + 0.5*part(4,i,1)*part(4,i,1)
!			enddo
!			call plsum(curx)
!			call plsum(cury)
!			if (idproc == 0) then
!				print*,itime,",",curx(1),",",cury(1)
!			endif
				
!		print*,E0xy(1),E0xy(2)
	  	fxye(1,:,:,:)=fxye(1,:,:,:)+E0xy(1)
	  	fxye(2,:,:,:)=fxye(2,:,:,:)+E0xy(2)
	  	
! efield diagnostic add by JF
			if ((dump_start < 0.) .OR. &
				&((real(itime)*dt >= dump_start) .and. (real(itime)*dt <= dump_end))) then
			if (ntfield > 0) then
				it = itime/ntfield
				if (itime==ntfield*it) then
					fname = './DIAG/Ex/'//'Ex-'//cdrun
					sfield = fxye(1,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,EX,trim(fname),inorder)
					fname = './DIAG/Ey/'//'Ey-'//cdrun
					sfield = fxye(2,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,EY,trim(fname),inorder)
				endif
			endif
			endif

			if ((nt_write_U_sumover_x_fromE .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_U_sumover_x_fromE
				if (itime == nt_write_U_sumover_x_fromE*it) then
					U_sumover_x_fromE = 0.
					call get_2D_sumover_x(0.5*fxye**2,U_sumover_x_fromE,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_U_sumover_x_fromE_funit) U_sumover_x_fromE
					endif
				endif
			endif
			
			! write the efield in bins so that it can be read as a surface plot
			! This diagnostic is a profligate waste of memory and disk space, but can be useful
			! sometimes
			if (nt_bin_E > 0) then
      	it = itime/nt_bin_E
      	if (itime == nt_bin_E*it) then
      		Ex_bin = 0.
      		Ey_bin = 0.
      		call bin_Ex_and_Ey(fxye,Ex_bin,Ey_bin,nx,ny,kyp,itime,itime*dt,nvp,idproc,inorder)
      	endif
      endif
				
! current diagnostic added by JF
      if (ntj > 0) then
      	it = itime/ntj
      	if (itime == ntj*it) then
					!Use Viktor's current deposit
					! the 0. is dth=0. so it doesn't do any push
					! the 1 is ipbc=1, which means periodic 2d boundary condition
					cu = 0.
					call gcjpost(part,fxye,npp,noff,cu,qme,qbme,dt/2.,tdcjpost,inorder)
					! add guard cells for current in x direction
					call aguard(cu,nyp,nx,inorder)
					! add guard cells for current in y direction
					call paguard(cu,kstrt,nvp,nx,kyp,ngds)
					! transform current to fourier space
					isign = -1
					call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
					! smooth current
					call sbpois(cut,tempt,ffc,nx,ny,kstrt)
					! transform current to back
					isign = 1
					call fft(cu,tempt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
					! copy guard cells for current in y direction
					call pcguard(cu,kstrt,nvp,kyp,inorder)
					! copy guard cells for current in x direction
					call cguard(cu,nyp,nx,inorder)

					fname = './DIAG/jx/'//'pcurx.'//cdrun
					!      			call j_diag_deposit(part,sfield,npp,1,ny,nvp,idproc)
					sfield = cu(1,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,JX,trim(fname),inorder)
					!do jy
					fname = './DIAG/jy/'//'pcury.'//cdrun
					!      			call j_diag_deposit(part,sfield,npp,2,ny,nvp,idproc)
					sfield = cu(2,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,JY,trim(fname),inorder)
      	endif
      endif

! nvdotE_part diagnostic added by JF
      if (nvdotE_follow_part > 0) then
      	it = itime/nvdotE_follow_part
      	if (itime == nvdotE_follow_part*it) then
      		!Use Viktor's current deposit
      		! the 0. is dth=0. so it doesn't do any push
					! the 1 is ipbc=1, which means periodic 2d boundary condition
					cu = 0.
					call djpost_jf(part,cu,npp,noff,qme,0.,tdjpost,nx,ny,1,inorder,djopt)
					! add guard cells for current in x direction
					call aguard(cu,nyp,nx,inorder)
					! add guard cells for current in y direction
					call paguard(cu,kstrt,nvp,nx,kyp,ngds)
					! transform current to fourier space
					isign = -1
					call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
					! smooth current
					call sbpois(cut,tempt,ffc,nx,ny,kstrt)
					! transform current to back
					isign = 1
					call fft(cu,tempt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
					! copy guard cells for current in y direction
					call pcguard(cu,kstrt,nvp,kyp,inorder)
					! copy guard cells for current in x direction
					call cguard(cu,nyp,nx,inorder)

					! do x component
					fname = './DIAG/vdotEx_follow_part/vdotEx'
					sfield = cu(1,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,VDOTEX_FOLLOW_PART,trim(fname),inorder)
					!do y component
					fname = './DIAG/vdotEy_follow_part/vdotEy'
					sfield = cu(2,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,VDOTEY_FOLLOW_PART,trim(fname),inorder)
      	endif
      endif

! Dump raw particle data
		! if ((dump_start < 0.) .OR. &
			! &((real(itime)*dt >= dump_start) .and. (real(itime)*dt <= dump_end))) then
		! if ( ntraw > 0) then
			! it = itime / ntraw
			! if (itime == ntraw*it) then
				! call write_raw(part,nx,ny,npp,nvp,itime*dt,itime,itime/ntraw,idproc,tracks)
			! endif
		! endif
		! endif

! Add track data to tracks
		if ( nttrack > 0) then
			it = itime / nttrack
			if (itime == nttrack*it) then
			 temp_time = 0.
       call PWTIMERA(-1,temp_time,add_dtime)
       call add_track_data( tracks, part, itime, real(itime)*dt )
       call PWTIMERA(1,temp_time,add_dtime)
       add_tracks_time = add_tracks_time + temp_time
			endif
		endif
! Write track data to disk

		if ( nt_dump_track > 0) then
			it = itime / nt_dump_track
			if (itime == nt_dump_track*it) then
				temp_time = 0.
				call PWTIMERA(-1,temp_time,write_dtime)
				if (track_no_h5 .eq. 0) then
					call write_tracks( tracks, nvp, idproc )
				else
					call write_tracks_no_h5( tracks, nvp, idproc )
				endif
					
				call PWTIMERA(1,temp_time,write_dtime)
				write_tracks_time = write_tracks_time + temp_time
			endif
		endif
		
! Phasespace of f(vx,x,y) and f(vy,x,y)
		if (ntphxy > 0) then
			it = itime / ntphxy
			if (itime == ntphxy*it) then
				fv_xy=0.
				call phase_space_vxy_vs_x_and_y(part,fv_xy,1,nx,ny,np,npp,itime,itime*dt,idproc)
				fv_xy=0.
				call phase_space_vxy_vs_x_and_y(part,fv_xy,2,nx,ny,np,npp,itime,itime*dt,idproc)
			endif
		endif
				
! Summed over transverse dim phase space diagnostic added by JF
		if (nphxx > 0) then
			it = itime / nphxx
			if (itime==nphxx*it) then
				fvxx = 0.
				call phase_space_vxy_vs_x(part,fvxx,1,nx,npp,real(itime)*dt,itime,nblok,idproc)
			    if (movion==1) then
				fvxx = 0.
			        call phase_space_vxy_vs_x_ion(parti,fvxx,1,nx,nppi,real(itime)*dt,itime,nblok,idproc)
			    endif
			endif
		endif
		if (nphyx > 0) then
			it = itime / nphyx
			if (itime==nphyx*it) then
				fvyx = 0.
				call phase_space_vxy_vs_x(part,fvyx,2,nx,npp,real(itime)*dt,itime,nblok,idproc)
			endif
		endif
		if (nphxy > 0) then
			it = itime / nphxy
			if (itime==nphxy*it) then
				fvxy = 0.
				call phase_space_vxy_vs_y(part,fvxy,1,ny,npp,real(itime)*dt,itime,nblok,idproc)
			endif
		endif
		if (nphyy > 0) then
			it = itime / nphyy
			if (itime==nphyy*it) then
				fvyy = 0.
				call phase_space_vxy_vs_y(part,fvyy,2,ny,npp,real(itime)*dt,itime,nblok,idproc)
			endif
		endif
		
		if (nt_phase_vx_vy > 0) then
			it = itime/nt_phase_vx_vy
			if (itime == nt_phase_vx_vy*it) then
				fvxvy = 0.
				call phase_vx_vy(part,fvxvy,np,npp,real(itime)*dt,itime,idimp,nblok,idproc)
			endif
		endif
		
		if (nt_vx_vy_speed > 0) then
			it = itime/nt_vx_vy_speed
			if (itime == nt_vx_vy_speed*it) then
				fvxvy = 0.
				call phase_vx_vy_x_y_movie_at_each_dump(part,fvxvy,np,npp,real(itime)*dt,itime,idimp,nblok,idproc)
			endif
		endif

! phase space slices diagnostic added by JF
		if (ntphsl_x > 0) then
			it = itime / ntphsl_x
			if (itime==ntphsl_x*it) then
				fvxxslice = 0.
				call phase_slices_x(part,fvxxslice,1,nx,ny,np,npp,itime*dt,itime,idimp,nblok,nvp,idproc)
				fvyxslice = 0.
				call phase_slices_x(part,fvyxslice,2,nx,ny,np,npp,itime*dt,itime,idimp,nblok,nvp,idproc)
			endif
		endif
		if (ntphsl_y > 0) then
			it = itime / ntphsl_y
			if (itime==ntphsl_y*it) then
				fvxyslice = 0.
				call phase_slices_y(part,fvxyslice,1,nx,ny,np,npp,itime*dt,itime,idimp,nblok,nvp,idproc)
				fvyyslice = 0.
				call phase_slices_y(part,fvyyslice,2,nx,ny,np,npp,itime*dt,itime,idimp,nblok,nvp,idproc)
			endif
		endif
! field lines diagnostic added by JF
		if (ntlines > 0) then
			it = itime / ntlines
			if (itime == ntlines*it) then
				call write_Eoft_aty(lines,fxye,nx,ny,nvp,idproc)
			endif
		endif
		
		if (keep_init_vel > 0) then
			it = itime / keep_init_vel
			if (itime == keep_init_vel*it) then			
				call calc_ene_dist(part,npp,itime*dt,itime,idproc)
				call calc_vx_dist(part,npp,itime*dt,itime,idproc)
				call calc_vy_dist(part,npp,itime*dt,itime,idproc)
			endif
		endif

		if (nt_through_wave > 0) then
			it = itime / nt_through_wave
			if (itime == nt_through_wave*it) then
				call calc_through_wave_dist(part,nx,npp,itime*dt,it,idproc, tot_dEne,tot_part_count, tot_mag_dEne)
			endif
		endif
		
			

		
		if (nt_write_field_ene_in_driver_region > 0) then
			it = itime/nt_write_field_ene_in_driver_region
			if (itime == nt_write_field_ene_in_driver_region*it) then
				
				call pois(qt,sfieldt,ffc,temp_we,nx,ny,kstrt)
				call ipgradf2(sfieldt,field_tempt,nx,ny,kstrt)
				isign = 1
				call fft(field_temp,field_tempt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				call get_ene_in_driver_region(field_temp,field_ene,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				if (idproc == 0) then
					write (77,*) real(itime)*dt,field_ene*num_par_cell
				endif
			endif
		endif

!**************************************************************************************************************
! nvdotE_part diagnostic added by JF
      if (nvdotE_part > 0) then

! if nvdotE_part is not 0, then we need to deposit the vdotE info every time step
! so we can add it into the vdotE_sum array.  Then we will write both of these
! at the interval specified by nvdotE_part.
! use cu array as a temp to hold vdotE for the current time step

				! cu will hold jE(t-dt)
				cu = 0.
				call djpost_jf(part,cu,npp,noff,qme,0.,tdjpost,nx,ny,1,inorder,djopt)
				! add guard cells for current in x direction
				call aguard(cu,nyp,nx,inorder)
				! add guard cells for current in y direction
				call paguard(cu,kstrt,nvp,nx,kyp,ngds)

				vdotE_sum = vdotE_sum + cu

      	it = itime/nvdotE_part
      	if (itime == nvdotE_part*it) then
      		!Use Viktor's current deposit
      		! the 0. is dth=0. so it doesn't do any push
	
					! do x component
					fname = './DIAG/vdotEx_part/vdotEx'
					sfield = cu(1,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,VDOTEX_PART,trim(fname),inorder)
					!do y component
					fname = './DIAG/vdotEy_part/vdotEy'
					sfield = cu(2,:,:,:)
					call writef(sfield,nxe,nypmx,itime,itime*dt,VDOTEY_PART,trim(fname),inorder)
      	endif				
      	endif

!**************************************************************************************************************
		if (n_jE_u_poynt .ne. 0) then

!------------------------			
				!This diag is to provide dE/dt time centered to one timestep back
				!1) Copy old grad_phi to grad_phi_dt and old dEdt to dEdt_dt
				!2) Get new grad_phi
				!3) put dE/dt in dEdt, this is centered at half timestep back
				!4) Now, need dE/dt centered at one timestep back, so that is 0.5*(dEdt + dEdt_dt)
			
				!Store old grad phi
				dEdt_dt = dEdt !This is now at t-3/2
				grad_phi_dt = grad_phi				!This is at t-1
				call pois(qt,sfieldt,ffc,temp_we,nx,ny,kstrt)
				call ipgradf2(sfieldt,grad_phit,nx,ny,kstrt)
				isign = 1
				call fft(grad_phi,grad_phit,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
				call pcguard(grad_phi,kstrt,nvp,kyp,inorder)
				call cguard(grad_phi,nyp,nx,inorder)
				grad_phi = -1.0 * grad_phi
				
			if (nt_dEdt .ne. 0) then
				dEdt = (grad_phi - grad_phi_dt) / dt
				it = itime/nt_dEdt
				if (itime == nt_dEdt*it+1) then
					fname = './DIAG/dExdt/'//'dExdt.'//cdrun
					sfield(:,:,1) = (dEdt(1,:,:,1) + dEdt_dt(1,:,:,1)) * 0.5
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,DEXDT,trim(fname),inorder)
					!do y
					fname = './DIAG/dEydt/'//'dEydt.'//cdrun
					sfield(:,:,1) = (dEdt(2,:,:,1) + dEdt_dt(2,:,:,1)) * 0.5
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,DEYDT,trim(fname),inorder)
				endif

				if (nt_write_grad_phi .ne. 0) then
					it = itime/nt_write_grad_phi
					if (itime == nt_write_grad_phi*it+1) then
						fname = './DIAG/grad_phi_x/'//'grad_phi_x.'//cdrun
						sfield(:,:,1) = (grad_phi_dt(1,:,:,1) + grad_phi_dt(1,:,:,1)) * 0.5
						call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,GRAD_PHI_X,trim(fname),inorder)
						!do y
						fname = './DIAG/grad_phi_y/'//'grad_phi_y.'//cdrun
						sfield(:,:,1) = (grad_phi_dt(2,:,:,1) + grad_phi_dt(2,:,:,1)) * 0.5
						call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,GRAD_PHI_Y,trim(fname),inorder)
					endif
				endif


			endif
!------------------------			

		!***** Do file writes and other manipulations to j.E and P here so that time
		!			 centering is correct.

			if ((nt_write_Py_sumover_x .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_Py_sumover_x
				if (itime == nt_write_Py_sumover_x*it) then
					Py_sumover_x = 0.
					sfield = ESPoynt(2,:,:,:)
					call get_1D_sumover_x(sfield,Py_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_Py_sumover_x_funit) Py_sumover_x
					endif
				endif
			endif
			if ((nt_write_jE_sumover_x .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_jE_sumover_x
				if (itime == nt_write_jE_sumover_x*it) then
					jE_sumover_x = 0.
					call get_2D_sumover_x(jE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_jE_sumover_x_funit) jE_sumover_x
					endif
				endif
			endif
			!just use jE_sumover_x as temp array
			if ((nt_write_kE_sumover_x .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_kE_sumover_x
				if (itime == nt_write_kE_sumover_x*it) then
					jE_sumover_x = 0.
					call get_1D_sumover_x(kE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_kE_sumover_x_funit) jE_sumover_x
					endif
				endif
			endif
			if ((nt_write_U_sumover_x .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_U_sumover_x
				if (itime == nt_write_U_sumover_x*it) then
					U_sumover_x = 0.
					call get_2D_sumover_x(0.5*grad_phi**2,U_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_U_sumover_x_funit) U_sumover_x
					endif
				endif
			endif
			if ((nt_write_Ux_sumover_x .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_Ux_sumover_x
				if (itime == nt_write_Ux_sumover_x*it) then
					U_sumover_x = 0.
					call get_2D_sumover_x_xcomp(0.5*grad_phi**2,U_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_Ux_sumover_x_funit) U_sumover_x
					endif
				endif
			endif

!			Pflow = 0.
!			totjE = 0.
!			totU = 0.
!			div_tot = 0.
!			call calc_ene_three_regions(grad_phi_int,ESPoynt_int,jE_sum,div_ESPoynt_sum,totU,Pflow,totjE,div_tot,num_par_cell,dt,&
!				&	nx,nxe,ny,nypmx,nvp,idproc,inorder)

!			if ((itime == int((timerise+timeflat+timefall)/dt)+2) .and. (nvdotE_int .ne. 0)) then		!+2 cause this diagnostic is back one time step
!				! do x component
!				fname = './DIAG/init_cond_vdotE_int/vdotEx'
!				sfield = jE_sum(1,:,:,:)
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEX_INT,trim(fname),inorder)
!				!do y component
!				fname = './DIAG/init_cond_vdotE_int/vdotEy'
!				sfield = jE_sum(2,:,:,:)
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEY_INT,trim(fname),inorder)
!			endif
			
!			if (nvdotE_int .ne. 0) then
!				it = itime/nvdotE_int
!				if (itime == nvdotE_int*it+1) then
!					! do x component
!					fname = './DIAG/vdotEx_int/vdotEx'
!					sfield = jE_sum(1,:,:,:)
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEX_INT,trim(fname),inorder)
!					!do y component
!					fname = './DIAG/vdotEy_int/vdotEy'
!					sfield = jE_sum(2,:,:,:)
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEY_INT,trim(fname),inorder)
!				endif
!			endif
      if (ntvdotE > 0) then      
      	it = itime/ntvdotE
      	if (itime == ntvdotE*it+1) then
					fname = './DIAG/vdotEx/'//'vdotEx.'//cdrun
					sfield(:,:,1) = jE(1,:,:,1)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEX,trim(fname),inorder)
					!do y
					fname = './DIAG/vdotEy/'//'vdotEy.'//cdrun
					sfield(:,:,1) = jE(2,:,:,1)
!					call vdotE(part,fxye,sfield,npp,2,ny,nvp,idproc)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,VDOTEY,trim(fname),inorder)
      	endif
      endif
      if (nt_kE > 0) then      
      	it = itime/nt_kE
      	if (itime == nt_kE*it+1) then
					fname = './DIAG/kE/'//'kE.'//cdrun
					sfield(:,:,1) = kE(:,:,1)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,KE_LABEL,trim(fname),inorder)
      	endif
      endif

!Write the Poynting flux at the first time step after the driver turns off so that 
!subtraction of initial condition can be done
!			if ((itime == int((timerise+timeflat+timefall)/dt)+2) .and. (ntESPoynt_int > 0)) then		!+2 cause this diagnostic is back one time step
!				fname = './DIAG/init_cond_ESPoynt_int/ESPoynt_int_x'
!				sfield = ESPoynt_int(1,:,:,:)
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_INT_X,trim(fname),inorder)
			
!				fname = './DIAG/init_cond_ESPoynt_int/ESPoynt_int_y'
!				sfield = ESPoynt_int(2,:,:,:)
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_INT_Y,trim(fname),inorder)						
!			endif

!			if (ntESPoynt_int > 0) then
!				it = itime/ntESPoynt_int
!				if (itime == ntESPoynt_int*it+1) then  !This is +1 because this diagnostic is 
!					fname = './DIAG/ESPoynt_int_x/ESPoynt_int_x'
!					sfield = ESPoynt_int(1,:,:,:)
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_INT_X,trim(fname),inorder)
				
!					fname = './DIAG/ESPoynt_int_y/ESPoynt_int_y'
!					sfield = ESPoynt_int(2,:,:,:)
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_INT_Y,trim(fname),inorder)							
!				endif
!			endif

			if (nt_write_jE_onlytracked_sumover_x > 0) then
				it = itime / nt_write_jE_onlytracked_sumover_x
				if (itime == nt_write_jE_onlytracked_sumover_x*it) then
					
					jE = 0.
					call pgcjepost2_onlytrack_jf(part,fxye,npp,noff,jE,qme,qbme,dt,tdcjpost,inorder)
	
					jE_sumover_x = 0.
					call get_2D_sumover_x(jE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_jE_onlytracked_sumover_x_funit) jE_sumover_x
					endif
				endif
			endif


			if (ntESPoynt > 0) then
				it = itime/ntESPoynt
				if (itime == ntESPoynt*it+1) then
					fname = './DIAG/ESPoynt_x/ESPoynt_x'
					sfield = ESPoynt(1,:,:,:)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_X,trim(fname),inorder)
				
					fname = './DIAG/ESPoynt_y/ESPoynt_y'
					sfield = ESPoynt(2,:,:,:)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,ESPOYNT_Y,trim(fname),inorder)
				endif
			endif

!For this diag must subtract the first two for initial condition			
!			if (((itime == int((timerise+timeflat+timefall)/dt)+2) .or. (itime == int((timerise+timeflat+timefall)/dt)+3)) .and. (ntfield_ene_int > 0)) then
!				fname = './DIAG/init_cond_u_int/u_x_int'
!				sfield = grad_phi_int(1,:,:,:)*num_par_cell
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,EX_ENE_INT,trim(fname),inorder)
!				fname = './DIAG/init_cond_u_int/u_y_int'
!				sfield = grad_phi_int(2,:,:,:)*num_par_cell
!				call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,EY_ENE_INT,trim(fname),inorder)
!			endif

!			if (ntfield_ene_int > 0) then
!				it = itime/ntfield_ene_int
!				if (itime==ntfield_ene_int*it+1) then
!					fname = './DIAG/u_x_int/u_x_int'
!					sfield = grad_phi_int(1,:,:,:)*num_par_cell
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,EX_ENE_INT,trim(fname),inorder)
!					fname = './DIAG/u_y_int/u_y_int'
!					sfield = grad_phi_int(2,:,:,:)*num_par_cell
!					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,EY_ENE_INT,trim(fname),inorder)
!				endif
!			endif

			!Get div P
!			isign = -1
!			call fft(ESPoynt,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!			call ipdivf2(cut,phit,nx,ny,kstrt)				
!			isign = 1
!			call fft(div_ESPoynt,phit,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
!			call pcguard(div_ESPoynt,kstrt,nvp,kyp,inorder)
!			call cguard(div_ESPoynt,nyp,nx,inorder)

			call div_finite_difference(ESPoynt,div_ESPoynt,nx,ny/nvp)

			!Time integrate it
	!		div_ESPoynt_sum = div_ESPoynt_sum + div_ESPoynt*dt

			if ((nt_write_div_P_vs_y .ne. 0) .and. (itime*dt > (timerise+timeflat+timefall))) then
				it = itime/nt_write_div_P_vs_y
				if (itime == nt_write_div_P_vs_y*it) then
					div_P_sumover_x = 0.
					call get_1D_sumover_x(div_ESPoynt,div_P_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_div_P_vs_y_funit) div_P_sumover_x
					endif
				endif
			endif
			if (nt_div_ESPoynt > 0) then
				it = itime/nt_div_ESPoynt
				if (itime == nt_div_ESPoynt*it+1) then  !This is +1 because this diagnostic is 
					fname = './DIAG/div_ESPoynt/div_ESPoynt'
					sfield = div_ESPoynt(:,:,:)
					call writef(sfield,nxe,nypmx,itime-1,itime*dt-dt,DIVESPOYNT,trim(fname),inorder)
				endif
			endif


		!-----Find P and jE, get P through cuperp*phi
			!Get cu and jE at current time step
		
			cu = 0.
			jE = 0.
			kE = 0.
			!This deposits both the current and j.E at the current timestep by pushing particles
			! half timestep, does not modify part array though 
			! should use dt not dt/2 cause inside function it does dt/2
			call gcjekejpost_jf(part,fxye,npp,noff,cu,jE,kE,qme,qbme,dt,tdcjpost,inorder)!tdcjpost is timer info
!			call gcjejpost_jf(part,fxye,npp,noff,cu,jE,qme,qbme,dt/2.,tdcjpost,inorder)!tdcjpost is timer info
			call aguard(cu,nyp,nx,inorder)
			call paguard(cu,kstrt,nvp,nx,kyp,ngds)
			call aguard(jE,nyp,nx,inorder)
			call paguard(jE,kstrt,nvp,nx,kyp,ngds)
			call aguard(kE,nyp,nx,inorder)
			call paguard(kE,kstrt,nvp,nx,kyp,ngds)
			!Now get cuperp
			isign = -1
			call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
			call cuperp_jf(cut,nx,ny,kstrt)
			call pois(qt,sfieldt,ffc,temp_we,nx,ny,kstrt)
			isign = 1
			!Now get j_perp in real space and put into cu
			call fft(cu,cut,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
			!Now put phi in real space into sfield
			call fft(sfield,sfieldt,isign,mixup,sct,tfft,indx,indy,kstrt,kyp,inorder)
						
			!Now get poynt and put into ESPoynt
			ESPoynt(1,:,:,:) = cu(1,:,:,:)*sfield
			ESPoynt(2,:,:,:) = cu(2,:,:,:)*sfield

			call pcguard(ESPoynt,kstrt,nvp,kyp,inorder)
			call cguard(ESPoynt,nyp,nx,inorder)
			call pcguard(cu,kstrt,nvp,kyp,inorder)
			call cguard(cu,nyp,nx,inorder)
			
			!Now ESPoynt holds the ES Poynting vector at current timestep
			!jE holds j.E at current timestep

!			ESPoynt_int = ESPoynt_int + ESPoynt*dt
!			jE_sum = jE_sum + jE*dt
			
		endif

! -----------------------------------OLD---------------------------------------------------------------
! this section removed to old_poynt.f			

!-----------------------------------------------END OLD-------------------------------------------------

! particle push and charge density update
      wke = 0.
! push electrons
      if (relativity==1) then
         call rpush(part,fxye,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,ipbc,i&
     &norder,popt)
      else
				if (set_part_array_vdotE .eq. 0) then
        	call push(part,fxye,npp,noff,qbme,dt,wke,tpush,nx,ny,ipbc,inorder,popt)
			  else
        	call push_jf(part,fxye,npp,noff,qbme,dt,wke,tpush,nx,ny,ipbc,inorder,popt)
			  endif

      endif
! move electrons into appropriate spatial regions
			if (nttrack .eq. 0) then
	      call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
			else
	      call pmove(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr,tracks)
	    endif

!      call pmove_jf(part,edges,npp,tmove,ny,kstrt,nvp,nbmax,vect,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! push ions
      if (movion==1) then
         wki = 0.
         if (relativity==1) then
            call rpush(parti,fxye,nppi,noff,qbmi,dt,ci,wki,tpushi,nx,ny,&
     &ipbc,inorder,popt)
         else
            call push(parti,fxye,nppi,noff,qbmi,dt,wki,tpushi,nx,ny,ipbc&
     &,inorder,popt)
         endif
         wki = wki*rmass
! move ions into appropriate spatial regions
         call pmove(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr)
!         call pmove_jf(parti,edges,nppi,tmovi,ny,kstrt,nvp,nbmax,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
      endif
! accumulate timings
      tpush = tpush + tdpost
      totpush = totpush + tpush
      tpart = tpart + tmove
      tmove = 0.0
      if (movion==1) then
         tpushi = tpushi + tdposti
         totpushi = totpushi + tpushi
         tparti = tparti + tmovi
         tmovi = 0.0
      endif
      totfft = totfft + tfft
      tfft = 0.0
! sort electrons

      if (sortime > 0) then
      	if (nttrack .eq. 0) then
         if (mod(itime,sortime)==0) then
        	call sortp(part,pt,ip,npp,noff,nyp,npic,tsort(6),inorder)
        endif
				else
				
				!Memory conserving sortp, the other is commented out below
				!Viktor's sortp has changed, so different parameters here
         if (mod(itime,sortime)==0) then
	        call sortp(part,pt,ip,npp,noff,nyp,npic,tsort,tpush,inorder,tracks)
!          call sortp(part,part2,npp,noff,nyp,npic,tsort,tpush,inorder,tracks)
				 endif
	      endif

!         call sortp(part,part2,npp,noff,nyp,npic,tsort,tpush,inorder)
! debug
!        if (tsort(1)==0.) then
!           if (id0==0) write (18,*) 'tsort(8) = ', tsort(8)
!        endif
      endif
! sort ions
      if ((movion==1) .and. (sortimi > 0)) then
         if (mod(itime,sortimi)==0) then
	         call sortp(parti,pt,ip,nppi,noff,nyp,npic,tsorti(6),inorder&
  	   &)
  	   	 endif
!        call sortp(parti,parti2,nppi,noff,nyp,npic,tsorti,tpushi,inorde&
!    &r)
      endif

! energy diagnostic
      if (ntw > 0) then
         it = itime/ntw
         if (itime==ntw*it) then
            wtot(1) = we
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = we + wke + wki
            call plsum(wtot)
! send energy values to diagnostic node
            msg(1:4) = wtot
            call HARTBEAT(msg,4)
            if (id0==0) then
            	write (18,992) wtot(1), wtot(2), wtot(4)
            	write (76,*) "field, kinetic, total ene",wtot(1), wtot(2), wtot(4)
            endif
            wt(it+1,:) = wtot
         endif
      endif
      itime = itime + 1
! restart file
      if (ntr > 0) then
         it = itime/ntr
         if (itime==ntr*it) then
            it = 16 + mod(it-1,2)
            if (id0==0) write (it) itime
            call wrdata(part,npp,it)
            if (movion==1) call wrdata(parti,nppi,it)
            if (movion==0) call wrdata(qi,nvp,it)
            if (ntw > 0) call wrdata(wt,1,it)
            if ((ntp > 0) .and. (id0==0)) write (it) nprec
            if (id0==0) then
               write (it) itime
               end file it
               rewind it
            endif
         endif
      endif
!**************************************************************************************************************

!			tot_ene_old = tot_ene
!			tot_ene = ES_ene_h + kin_ene
!			if (idproc == 0) then
!				print*,itime,vdote_tot,(ES_ene*num_par_cell),vdote_tot+(ES_ene*num_par_cell)
!				print*,itime,vdote_tot,(ES_ene*num_par_cell),(poynt_tot*num_par_cell), &
!					&vdote_tot+(ES_ene*num_par_cell)+(poynt_tot*num_par_cell),wtot
				
!				print*,"d vdote + d ESene - ",
!			endif

!**************************************************************************************************************
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! send QUIT message to diagnostic nodes
      call plsum(tot_dEne)
      call plsum(tot_mag_dEne)
      call plsum(tot_part_count)
      tot_dEne(1) = tot_dEne(1) / real(tot_part_count(1))
      tot_mag_dEne(1) = tot_mag_dEne(1) / real(tot_part_count(1))

      msg = -1.
      call HARTBEAT(msg,1)
! energy diagnostic
      if (ntw > 0) then
         it = (itime - 1)/ntw + 1
         call displayw(wt,dt*real(ntw),it,irc)
! check error return code
         if (irc==1) go to 3000
      endif
      call pwtimer(time,dtime)
! send main CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      msg(3) = totpush; msg(4) = tsort(6); msg(5:6) = tpart
      msg(7:8) = totfft
      call HARTBEAT(msg,8)
      if (id0==0) then
         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'
         write (18,*) 'total electron push time=', totpush, 'sec'
         write (18,*) 'electron sort/move time=', tsort(6), tpart
         tpart(1) = totpush + tsort(6) + tpart(1)
         write (18,*) 'total electron time=', tpart(1), 'sec'
         write (18,*) 'total add_tracks_time=',add_tracks_time
         write (18,*) 'total write_tracks_time=',write_tracks_time
      endif
      if (movion==1) then
         msg(1) = totpushi; msg(2) = tsorti(6); msg(3:4) = tparti
         call HARTBEAT(msg,4)
         if (id0==0) then
            write (18,*) 'total ion push time=', totpushi, 'sec'
            write (18,*) 'ion sort/move time=', tsorti(6), tparti
            tparti(1) = totpushi + tsorti(6) + tparti(1)
            write (18,*) 'total ion time=', tparti(1), 'sec'
         endif
      endif
      if (id0==0) then
         write (18,*) 'total fft time=', totfft, 'sec'
         if (nt_through_wave .ne. 0) write (18,*) 'tot_dEne=',tot_dEne
         if (nt_through_wave .ne. 0) write (18,*) 'tot_mag_dEne=',tot_mag_dEne
! write final diagnostic metafile
         close(unit=19)
         fname = 'pdiag2.final.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
! potential diagnostics
         if (ntp > 0) then
            nprec = nprec - 1
            write (19,ppot2d,iostat=irc)
         endif
			call h5close_f(h5error)

! done
         write (18,*) '* * * q.e.d. * * *'
      endif
! close graphics device

      if (nttrack .ne. 0) then
				call cleanup(tracks)
				if (track_no_h5 .eq. 0) then
					call h5close_f(ierr)
				else
					call close_file_no_h5(tracks)
				endif
			endif

 3000 call PGRCLOSE
      call MP_END
      call PPEXIT
      stop
!
      contains
!
         subroutine diag2nodes
         implicit none
! diagnostic nodes have special processing
  991    format (' T = ',i7)
  992    format (' field, kinetic, total energies = ',3e14.7)
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),nblok))
            allocate(npp(nblok))
         endif
         if (movion==0) allocate(qi(nxe,nypmx*kbmin,kblok))
! restart
         if (nustrt /= 1) then
! determine most recent restart file
            if (id0==0) then
               read (16,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               read (17,iostat=ierr) ktime(2)
               if (ierr /= 0) ktime(2) = -1
               if (ktime(1) > ktime(2)) then
                  ktime(2) = 16
               else
                  ktime(1) = ktime(2)
                  ktime(2) = 17
               endif
            endif
            call plbcast(ktime)
            itime = ktime(1)
            if (itime < 0) go to 30
! read restart file
            it = ktime(2)
            call rddata(part,npp,it,ierr)
            if (ierr /= 0) go to 30
            if (movion==1) then
               call rddata(part,npp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (movion==0) then
               call rddata(qi,nvp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (ntw > 0) then
               call rddata(wt,1,it,ierr)
               if (ierr /= 0) go to 30
               call plbcast(wt)
            endif
            if (ntp > 0) then
               if (id0==0) then
                  read (it,iostat=ierr) ktime(1)
                  if (ierr /= 0) ktime(1) = -1
                  irc = 0
                  fname = 'ppotk2.'//cdrun
                  call writebf(pott,modesxp,modesy2p,kxp,11,irc,trim(fname)&
     &)
               endif
               call plbcast(ktime)
               nprec = ktime(1)
               if (nprec< 0) go to 30
            endif
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
               rewind it
            endif
            call plbcast(ktime)
            irc = ktime(1)
            if (irc==itime) go to 40
! handle error
   30       if (id0==0) write (18,*) 'Restart Error'
            call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         endif
! get initial CPU Time
   40    call HARTBEAT(msg,2)
         time(1) = msg(1); time(2) = msg(2)
         if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
         endif
! get time step
   10    call HARTBEAT(msg,1)
         it = msg(1)
         if (it < 0) then
! energy diagnostic
            if (ntw > 0) then
               it = (itime - 1)/ntw + 1
               call displayw(wt,dt*real(ntw),it,irc)
! check error return code
               if (irc==1) go to 20
            endif
! get main CPU Time
            call HARTBEAT(msg,8)
            time(1) = msg(1); time(2) = msg(2)
            totpush = msg(3); tsort(6) = msg(4); tpart = msg(5:6)
            totfft = msg(7:8)
            if (id0==0) then
               write (18,*) 'main max/min real time=', time(1), time(2),&
     &'sec'
               write (18,*) 'total electron push time=', totpush, 'sec'
               write (18,*) 'electron sort/move time=', tsort(6), tpart
               tpart(1) = totpush + tsort(6) + tpart(1)
               write (18,*) 'total electron time=', tpart(1), 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,4)
               totpushi = msg(1); tsorti(6) = msg(2); tparti = msg(3:4)
               if (id0==0) then
                  write (18,*) 'total ion push time=', totpushi, 'sec'
                  write (18,*) 'ion sort/move time=', tsorti(6), tparti
                  tparti(1) = totpushi + tsorti(6) + tparti(1)
                  write (18,*) 'total ion time=', tparti(1), 'sec'
               endif
            endif
            if (id0==0) then
               write (18,*) 'total fft time=', totfft, 'sec'
! write final diagnostic metafile
               close(unit=19)
               fname = 'pdiag2.final.'//cdrun
               open(unit=19,file=trim(fname),form='formatted',status='re&
     &place')
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  write (19,ppot2d,iostat=irc)
               endif
! done
               write (18,*) '* * * q.e.d. * * *'
            endif
   20       call PGRCLOSE
            call MP_END
            call PPEXIT
            stop
         else
            itime = it
         endif
         if (id0==0) write (18,991) itime
         write (label,991) itime
         call LOGNAME(label)
! density diagnostic
         if (ntd > 0) then
            it = itime/ntd
            if (itime==ntd*it) then
! display density
               call displays(sfield,nvp,' E DENSITY',itime,999,2,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 10
            endif
         endif
! velocity diagnostic
         if (ntv > 0) then
            it = itime/ntv
            if (itime==ntv*it) then
! display velocity distributions
               call displayfv(fv,fvm,' ELECTRON',itime,nmv,2,2,irc)
               if (irc==1) go to 10
! print out velocity moments
               if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
               if (movion==1) then
! display velocity distributions
                  call displayfv(fvi,fvmi,' ION',itime,nmv,2,2,irc)
                  if (irc==1) go to 10
! print out velocity moments
                  if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
              endif
            endif
         endif
! phase space diagnostic
         if (nts > 0) then
            it = itime/nts
            if (itime==nts*it) then
               isc = 999
! plot electrons x versus y
               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&
     &ny,1,2,irc)
               if (irc==1) go to 10
! plot electrons vx versus x
               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&
     &ny,3,1,irc)
               if (irc==1) go to 10
! plot electrons vy versus y
               call grasp(part,npp,' ELECTRON PHASE SPACE',itime,isc,nx,&
     &ny,4,2,irc)
               if (irc==1) go to 10
               if (movion==1) then
! plot ions x versus x
                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&
     &,1,2,irc)
                  if (irc==1) go to 10
! plot ions vx versus x
                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&
     &,3,1,irc)
                  if (irc==1) go to 10
! plot ions vy versus y
                  call grasp(part,npp,' ION PHASE SPACE',itime,isc,nx,ny&
     &,4,2,irc)
                  if (irc==1) go to 10
               endif
            endif
         endif
! potential diagnostic
         if (ntp > 0) then
            it = itime/ntp
            if (itime==ntp*it) then
! display potential
               call displays(sfield,nvp,' POTENTIAL',itime,999,0,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 10
!              call displays(g,nvp,' POT(NXH)',itime,999,0,ny,irc,inorde&
!    &r)
!              if (irc==1) go to 10
! write diagnostic output
               if (nprec==0) then
                  nprec = -1
                  fname = 'ppotk2.'//cdrun
                  call writebf(pott,modesxp,modesy2p,kxp,11,nprec,trim(fnam&
     &e))
               else
                  call writebf(pott,modesxp,modesy2p,kxp,11,nprec)
               endif
            endif
         endif
! efield diagnostic add by JF
!			if (ntfield > 0) then
!				 it = itime/ntfield
!				 if (itime==ntfield*it) then
!						if (nfieldxrec==0) then
!							 fname = 'Ex-'//cdrun
!							 sfield = fxye(1,:,:,:)
!							 call writebf(sfield,nx,kyp,50,nfieldxrec,trim(fname))
!							 fname = 'Ey-'//cdrun
!							 sfield = fxye(2,:,:,:)
!							 call writebf(sfield,nx,kyp,51,nfieldyrec,trim(fname))
!						else
!							 sfield = fxye(1,:,:,:)
!							 call writebf(sfield,nx,kyp,50,nfieldxrec)
!							 sfield = fxye(2,:,:,:)
!							 call writebf(sfield,nx,kyp,51,nfieldyrec)
!						endif
!				 endif
!			endif


! energy diagnostic
      if (ntw > 0) then
         it = itime/ntw
         if (itime==ntw*it) then
! get energy values
            call HARTBEAT(msg,4)
            wtot = msg(1:4)
            if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
         endif
      endif
! restart file
         itime = itime + 1
         if (ntr > 0) then
            it = itime/ntr
            if (itime==ntr*it) then
               it = 16 + mod(it-1,2)
               if (id0==0) write (it) itime
               call wrdata(part,npp,it)
               if (movion==1) call wrdata(part,npp,it)
               if (movion==0) call wrdata(qi,nvp,it)
               if (ntw > 0) call wrdata(wt,1,it)
               if ((ntp > 0) .and. (id0==0)) write (it) nprec
               if (id0==0) then
                  write (it) itime
                  end file it
                  rewind it
               endif
            endif
         endif
         go to 10
         end subroutine
!
      end program pbeps2
