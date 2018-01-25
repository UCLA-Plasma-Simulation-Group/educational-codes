!-----------------------------------------------------------------------
!
      module pinit2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pinit2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: june 24, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR,&
     & PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D, VACUUM_2D,     &
     &VACUUM_3D, NEUMANN_2D
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D
      public :: VACUUM_2D, VACUUM_3D, NEUMANN_2D
      public :: idrun, idrun0, indx, indy, npx, npy, npxb, npyb, inorder
      public :: popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, ntd, nta, ntv, nts, ntm, nte
      public :: ndw, ndp, ndd, nda, ndv, nds, ndm, nde
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0 
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve, relativity, omx, omy, omz, ci, ax, ay
      public :: ndc, movion, npxi, npyi, npxbi, npybi
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, v0, w0
      public :: sortime, sortimi, nplot, idpal, ndstyle, sntasks
      public :: itpon, ionoff, nsrand, ndprof, nsrandi, ndprofi
      public :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy
      public :: ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi
      public :: modesxd, modesyd, modesxp, modesyp, modesxa, modesya
      public :: modesxe, modesye
      public :: imbalance, mpimon
      public :: pinput2, sendnml, distr, pvdistr, ldistr, fdistr, vdistr
      public :: vfdistr, vvdistr, fedges
      public :: t0, ceng, indian, rlprec, inprec, pden2d, ndrec, fdname
      public :: ppot2d, nprec, fpname, pvpot2d, narec, faname, pem2d
      public :: nerec, fename
      public :: waterbag,supergauss4,supergauss5
      public :: nensemble
      public :: b0
!
! Namelist Input
      save
! idrun/idrun0 = run identifier for current/old run
      integer :: idrun = 0, idrun0 = 0
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! npx/npy = initial number of particles distributed in x/y direction
!     integer :: indx =   5, indy =   6, npx =      96, npy =     192
      integer :: indx =   6, indy =   7, npx =     384, npy =     768
!     integer :: indx =   7, indy =   8, npx =    1280, npy =    2560
! npxb/npyb = initial number of particles in beam in x/y direction
!     integer :: npxb =   0, npyb =   0
      integer :: npxb =  32, npyb =  64
!     integer :: npxb = 128, npyb = 256
!     integer :: npxb = 384, npyb = 768
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
      integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
      integer :: djopt = STANDARD
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
      integer :: nustrt = 1, ntr = 0
! ntw, ndw = number of time steps between energy diagnostic
! ntp, ndp = number of time steps between potential diagnostic
! ntd, ndd = number of time steps between ion density diagnostic
! nta, nda = number of time steps between vector potential diagnostic
! ntv, ndv = number of time steps between velocity-space diagnostic
! nts, nds = number of time steps between phase space diagnostic
      integer :: ntw = 1, ntp = 0, ntd = 0, nta = 0, ntv = 0, nts = 0
      integer :: ndw = 1, ndp = 0, ndd = 0, nda = 0, ndv = 0, nds = 0
! ntm, ndm = number of time steps between momentum diagnostic
! nte, nde = number of time steps between electromagnetic diagnostic
      integer :: ntm = 0, nte = 0
      integer :: ndm = 0, nde = 0
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend =  65.000, dt = 0.2000000e+00
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! psolve = type of poisson solver = (1,2,3)
      integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! ci = reciprical of velocity of light
      real :: ci = 1.0
! ax/ay = half-width of particle in x/y direction
!     real :: ax = .816497, ay = .816497
!     real :: ax = .866025, ay = .866025
      real :: ax = .912871, ay = .912871
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
! movion = (0,1) = (no,yes) move the ions
! npxi/npyi = initial number of ions distributed in x/y/z direction
      integer :: movion = 0, npxi =  384, npyi =  768
! npxbi/npybi = initial number of ions in beam in x/y/z direction
      integer :: npxbi =   32, npybi =   64
! qmi = charge on ion, in units of 3
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
! v0 = external pump strength, in units vos/vthermal
! w0 = external pump frequency, in units of wpe
      real :: v0 = 0.0, w0 = 0.0
! sortime = number of time steps between electron sorting
! sortimi = number of time steps between ion sorting
      integer :: sortime = 50, sortimi = 250
! nplot = maximum number of plots per page
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
! ndstyle = (1,2,3) = display (color map,contour plot,both)
      integer :: nplot = 4, idpal = 1, ndstyle = 1
! sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
      integer :: sntasks = -1
! itpon = time when external pump is turned on (-1=never)
! ionoff = time when ions are frozen and their charge saved (-1=never)
      integer :: itpon = -1, ionoff = -1
! nsrand = (0,1) = (no,yes) randomize spatially positions locally
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
! nsrandi = (0,1) = (no,yes) randomize spatially ion positions locally
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: nsrand = 0, ndprof = 0, nsrandi = 0, ndprofi = 0
! ampdx/ampdx = amplitude of density compared to uniform in x/y
! scaledx/scaledx = scale length for spatial coordinate in x/y
! shiftdx/shiftdx = shift of spatial coordinate in x/y
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
! ampdxi/ampdxi = amplitude of ion density compared to uniform in x/y
! scaledxi/scaledxi = scale length for spatial ion coordinate in x/y
! shiftdxi/shiftdxi = shift of spatial ion coordinate in x/y
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
! modesxd/modesyd = number of modes in x/y to keep for ion density 
! diagnostic
      integer :: modesxd = 11, modesyd = 11
! modesxp/modesyp = number of modes in x/y to keep for potential 
! diagnostic
      integer :: modesxp = 11, modesyp = 11
! modesxa/modesya = number of modes in x/y to keep for vector potential
! diagnostic
      integer :: modesxa = 11, modesya = 11
! modesxe/modesye = number of modes in x/y to keep for electromagnetic
! diagnostic
      integer :: modesxe = 11, modesye = 11

! imbalance = load imbalance fraction repartition trigger
! (< 0.  to suppress repartion)
      real :: imbalance = .08
! mpimon = (0,1,2) = (suppress,display,display and log) mpi messages
      integer :: mpimon = 1
! nensemble = how many times to initialize the velocity distribution function
!             for ensemble averaging (default is 0)
      integer :: nensemble = 0

! b0 --> uniform magnetic field
      real b0

! define namelist
      namelist /pinput2/ idrun, idrun0, indx, indy, npx, npy, npxb, npyb&
     &, inorder, popt, dopt, djopt, nustrt, ntr, ntw, ntp, ntd, nta, ntv&
     &, nts, ntm, nte, ndw, ndp, ndd, nda, ndv, nds, ndm, nde, tend, dt,&
     & qme, vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, vtdy,    &
     &vtdz, psolve, relativity, omx, omy, omz, ci, ax, ay, ndc, movion, &
     &npxi, npyi, npxbi, npybi, qmi, rmass, rtempxi, rtempyi, rtempzi,  &
     &vxi0, vyi0, vzi0, vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, &
     &v0, w0, sortime, sortimi, nplot, idpal, ndstyle, sntasks, itpon,  &
     &ionoff, nsrand, ndprof, ampdx, scaledx, shiftdx, ampdy, scaledy,  &
     &shiftdy, nsrandi, ndprofi, ampdxi, scaledxi, shiftdxi, ampdyi,    &
     &scaledyi, shiftdyi, modesxd, modesyd, modesxp, modesyp, modesxa,  &
     &modesya, modesxe, modesye, imbalance, mpimon, nensemble,b0
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
! indian = (0,1) = architecture is (little-endian,big-endian)
! rlprec = (0,1) = default reals are (normal,double-precision)
! inprec = (0,1) = default integers are (normal,double-precision)
      integer :: indian = 1, rlprec = 1, inprec = 0
!
! Namelist output for ion density diagnostic
! ndrec = current record number for ion density writes
      integer :: ndrec = 0
! fdname = file name for potential diagnostic
      character(len=32) :: fdname
! define namelist
      namelist /pden2d/ idrun, indx, indy, ntd, modesxd, modesyd, psolve&
     &, ndrec, t0, tend, dt, ceng, indian, rlprec, inprec, fdname 
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname
! define namelist
      namelist /ppot2d/ idrun, indx, indy, ntp, modesxp, modesyp, psolve&
     &, omx, omy, omz, nprec, t0, tend, dt, ceng, indian, rlprec, inprec&
     &, fpname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
      integer :: narec = 0
! faname = file name for potential diagnostic
      character(len=32) :: faname
! define namelist
      namelist /pvpot2d/ idrun, indx, indy, nta, modesxa, modesya,      &
     &psolve, omx, omy, omz, ci, narec, t0, tend, dt, ceng, indian,     &
     &rlprec, inprec, faname
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
      integer :: nerec = 0
! fename = file name for potential diagnostic
      character(len=32) :: fename
! define namelist
      namelist /pem2d/ idrun, indx, indy, nte, modesxe, modesye, psolve,&
     &omx, omy, omz, ci, nerec, t0, tend, dt, ceng, indian, rlprec,     &
     &inprec, fename
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx&
     &,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,i&
     &dimp,npmax,nblok,ipbc,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, ipbc, kstrt
         integer :: nvp, ndv, nvrp, ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         end subroutine
      end interface
      interface
         subroutine PISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,n&
     &px,npy,nx,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, idps, ipbc
         integer :: ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,np&
     &y,nx,ny,idimp,npmax,nblok,ipbc,vranx,vrany,vranz,kstrt,nvp,ndv,nvr&
     &p,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, ipbc, kstrt
         integer :: nvp, ndv, nvrp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany, vranz
         end subroutine
      end interface
      interface
         subroutine PLDISTR2(part,nps,anlx,anly,npx,npy,nx,ny,idimp,npma&
     &x,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc, ierr
         real :: anlx, anly
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         end subroutine
      end interface
      interface
         subroutine PFDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,ar&
     &gy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PRDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,ar&
     &gy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ipbc
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PVRDISTR2(part,nps,fnx,argx1,argx2,argx3,fny,argy1,a&
     &rgy2,argy3,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt,nvp,n&
     &dv,nvrp,ipbc,ierr)
         implicit none
         integer :: npx, npy, nx, ny, idimp, npmax, nblok, kstrt, nvp
         integer :: ndv, nvrp, ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine PVDISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,&
     &npmax,nblok,kstrt,nvp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,&
     &npmax,nblok,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ndv, nvrp
         integer :: ierr
         real :: vtx, vty, vdx, vdy
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany
         end subroutine
      end interface
      interface
         subroutine PVDISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,idimp,npmax,nblok,kstrt,nvp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,idimp,npmax,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, idimp, npmax, nblok, kstrt, nvp, ndv, nvrp
         integer :: ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,nblok) :: part
         integer, dimension(nblok) :: npp, nps
         double precision, dimension(nvrp,nblok) :: vranx, vrany, vranz
         end subroutine
      end interface
      interface
         subroutine PBDISTR2L(part,bx,by,bz,npp,noff,qbm,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real, dimension(nxv,nypmx,nblok):: bx, by, bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBDISTR2L(part,bxy,npp,noff,qbm,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGBZDISTR2L(part,bz,npp,noff,qbm,nx,ny,idimp,npmax,n&
     &blok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm
         real, dimension(idimp,npmax,nblok) :: part
         real :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBDISTR2L(part,bxy,npp,noff,qbm,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real :: bxy
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine PGRBZDISTR2L(part,bz,npp,noff,qbm,ci,nx,ny,idimp,npm&
     &ax,nblok,nxv,nypmx,ipbc)
         implicit none
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc
         real :: qbm, ci
         real, dimension(idimp,npmax,nblok) :: part
         real :: bz
         integer, dimension(nblok) :: npp, noff
         end subroutine
      end interface
      interface
         subroutine FEDGES2(edges,noff,nyp,fny,arg1,arg2,arg3,ny,nypmin,&
     &nypmax,kstrt,nvp,nblok,idps,ipbc)
         implicit none
         integer :: ny, nypmin, nypmax, kstrt, nvp, nblok, idps, ipbc
         real :: arg1, arg2, arg3
         real, dimension(idps,nblok) :: edges
         integer, dimension(nblok) :: noff, nyp
         real, external :: fny
         end subroutine
      end interface
      
      interface
          real function waterbag(v)
          real v
          end function waterbag
          real function supergauss4(v)
          real v
          end function supergauss4
          real function supergauss5(v)
          real v
          end function supergauss5
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface distr
         module procedure ipistr2
         module procedure ipistrh2
         module procedure ipbpistr2
         module procedure iprbpistr2
      end interface
!
      interface pvdistr
         module procedure ipvistr2
         module procedure ipvistrh2
      end interface
!
      interface ldistr
         module procedure ipldistr2
      end interface
!
      interface fdistr
         module procedure ipfdistr2
      end interface
!
      interface vdistr
         module procedure ipvdistr2
         module procedure ipvdistr2_arb
         module procedure ipvdistrh2
      end interface
!
      interface vfdistr
         module procedure ipvfdistr2
      end interface
!
      interface vvdistr
         module procedure ipvvdistr2
         module procedure ipvvdistrh2
      end interface
!
      interface fedges
         module procedure ifedges2
      end interface
!
      interface sendnml
         module procedure sendnml2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine sendnml2()
! this subroutine packs 2d namelist variables into a double precision
! buffer and broadcasts them to other nodes
         integer, parameter :: lenml = 109
         double precision, dimension(lenml) :: ddata
! pack data
         ddata(1) = idrun; ddata(2) = idrun0
         ddata(3) = indx; ddata(4) = indy
         ddata(5) = npx; ddata(6) = npy
         ddata(7) = npxb; ddata(8) = npyb
         ddata(9) = inorder; ddata(10) = popt; ddata(11) = dopt
         ddata(12) = djopt;; ddata(13) = nustrt; ddata(14) = ntr
         ddata(15) = ntw; ddata(16) = ntp; ddata(17) = ntd
         ddata(18) = nta; ddata(19) = ntv; ddata(20) = nts
         ddata(21) = ntm; ddata(22) = nte
         ddata(23) = ndw; ddata(24) = ndp; ddata(25) = ndd
         ddata(26) = nda; ddata(27) = ndv; ddata(28) = nds
         ddata(29) = ndm; ddata(30) = nde
         ddata(31) = tend; ddata(32) = dt; ddata(33) = qme
         ddata(34) = vtx; ddata(35) = vty; ddata(36) = vtz
         ddata(37) = vx0; ddata(38) = vy0; ddata(39) = vz0
         ddata(40) = vdx; ddata(41) = vdy; ddata(42) = vdz
         ddata(43) = vtdx; ddata(44) = vtdy; ddata(45) = vtdz
         ddata(46) = psolve; ddata(47) = relativity
         ddata(48) = omx; ddata(49) = omy; ddata(50) = omz
         ddata(51) = ci; ddata(52) = ax; ddata(53) = ay
         ddata(54) = ndc; ddata(55) = movion
         ddata(56) = npxi; ddata(57) = npyi
         ddata(58) = npxbi; ddata(59) = npybi
         ddata(60) = qmi; ddata(61) = rmass
         ddata(62) = rtempxi; ddata(63) = rtempyi; ddata(64) = rtempzi
         ddata(65) = vxi0; ddata(66) = vyi0; ddata(67) = vzi0
         ddata(68) = vdxi; ddata(69) = vdyi; ddata(70) = vdzi
         ddata(71) = rtempdxi; ddata(72) = rtempdyi
         ddata(73) = rtempdzi; ddata(74) = v0; ddata(75) = w0
         ddata(76) = sortime; ddata(77) = sortimi
         ddata(78) = nplot; ddata(79) = idpal; ddata(80) = ndstyle
         ddata(81) = sntasks; ddata(82) = itpon; ddata(83) = ionoff
         ddata(84) = nsrand; ddata(85) = ndprof
         ddata(86) = ampdx; ddata(87) = scaledx; ddata(88) = shiftdx
         ddata(89) = ampdy; ddata(90) = scaledy; ddata(91) = shiftdy
         ddata(92) = nsrandi; ddata(93) = ndprofi
         ddata(94) = ampdxi; ddata(95) = scaledxi; ddata(96) = shiftdxi
         ddata(97) = ampdyi; ddata(98) = scaledyi; ddata(99) = shiftdyi
         ddata(100) = modesxd; ddata(101) = modesyd
         ddata(102) = modesxp; ddata(103) = modesyp
         ddata(104) = modesxa; ddata(105) = modesya
         ddata(106) = modesxe; ddata(107) = modesye
         ddata(108) = imbalance; ddata(109) = mpimon
! broadcast data
         call PBDCAST(ddata,lenml)
! unpack data
         idrun = ddata(1); idrun0 = ddata(2)
         indx = ddata(3); indy = ddata(4)
         npx = ddata(5); npy = ddata(6)
         npxb = ddata(7); npyb = ddata(8)
         inorder = ddata(9); popt = ddata(10); dopt = ddata(11)
         djopt = ddata(12); nustrt = ddata(13); ntr = ddata(14)
         ntw = ddata(15); ntp = ddata(16); ntd = ddata(17)
         nta = ddata(18); ntv = ddata(19); nts = ddata(20)
         ntm = ddata(21); nte = ddata(22)
         ndw = ddata(23); ndp = ddata(24); ndd = ddata(25)
         nda = ddata(26); ndv = ddata(27); nds = ddata(28)
         ndm = ddata(29); nde = ddata(30)
         tend = ddata(31); dt = ddata(32); qme = ddata(33)
         vtx = ddata(34); vty = ddata(35); vtz = ddata(36)
         vx0 = ddata(37); vy0 = ddata(38); vz0 = ddata(39)
         vdx = ddata(40); vdy = ddata(41); vdz = ddata(42)
         vtdx = ddata(43); vtdy = ddata(44); vtdz = ddata(45)
         psolve = ddata(46); relativity = ddata(47)
         omx = ddata(48); omy = ddata(49); omy = ddata(50)
         ci = ddata(51); ax = ddata(52); ay = ddata(53)
         ndc= ddata(54); movion = ddata(55)
         npxi = ddata(56); npyi = ddata(57)
         npxbi = ddata(58); npybi = ddata(59)
         qmi = ddata(60); rmass = ddata(61)
         rtempxi = ddata(62); rtempyi = ddata(63); rtempzi = ddata(64)
         vxi0 = ddata(65); vyi0 = ddata(66); vzi0 = ddata(67)
         vdxi = ddata(68); vdyi = ddata(69); vdzi = ddata(70)
         rtempdxi = ddata(71); rtempdyi = ddata(72)
         rtempdzi = ddata(73); v0 = ddata(74); w0 = ddata(75)
         sortime = ddata(76); sortimi = ddata(77)
         nplot = ddata(78); idpal = ddata(79); ndstyle = ddata(80)
         sntasks = ddata(81); itpon = ddata(82); ionoff = ddata(83)
         nsrand = ddata(84); ndprof = ddata(85)
         ampdx = ddata(86); scaledx = ddata(87); shiftdx = ddata(88)
         ampdy = ddata(89); scaledy = ddata(90); shiftdy = ddata(91)
         nsrandi = ddata(92); ndprofi = ddata(93)
         ampdxi = ddata(94); scaledxi = ddata(95); shiftdxi = ddata(96)
         ampdyi = ddata(97); scaledyi = ddata(98); shiftdyi = ddata(99)
         modesxd = ddata(100); modesyd = ddata(101)
         modesxp = ddata(102); modesyp = ddata(103)
         modesxa = ddata(104); modesya = ddata(105)
         modesxe = ddata(106); modesye = ddata(107)
         imbalance = ddata(108); mpimon = ddata(109)
         end subroutine sendnml2
!
         subroutine ipistr2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,n&
     &x,ny,ipbc)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: npx, npy, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, idps, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); idps = size(edges,1)
         call PISTR2(part,edges,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,id&
     &imp,npmax,nblok,idps,ipbc,ierr)
         end subroutine ipistr2
!
         subroutine ipistrh2(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,nx,ny,ipbc)
! calculates initial particle co-ordinates and velocities in 2-1/2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: npx, npy, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, idps, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3); idps = size(edges,1)
         call PISTR2H(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,nx,ny,idimp,npmax,nblok,idps,ipbc,ierr)
         end subroutine ipistrh2
!
         subroutine ipvistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,&
     &ipbc,kstrt,nvp)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, ipbc, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,nx,ny,idimp,n&
     &pmax,nblok,ipbc,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvistr2
!
         subroutine ipvistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,nx,ny,ipbc,kstrt,nvp)
! calculates initial particle co-ordinates and velocities in 2-1/2d
! with uniform density and maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, ipbc, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,nx,n&
     &y,idimp,npmax,nblok,ipbc,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr&
     &)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvistrh2
!
         subroutine ipldistr2(part,nps,anlx,anly,npx,npy,nx,ny,kstrt,nvp&
     &,ipbc)
! calculates initial particle co-ordinates in 2d
! with bi-linear density profile
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc
         real :: anlx, anly
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PLDISTR2(part,nps,anlx,anly,npx,npy,nx,ny,idimp,npmax,nblo&
     &k,kstrt,nvp,ipbc,ierr)
         end subroutine ipldistr2
!
         subroutine ipfdistr2(part,nps,ampx,scalex,shiftx,ampy,scaley,sh&
     &ifty,npx,npy,nx,ny,kstrt,nvp,ipbc,ndpro,nsran)
! calculates initial particle co-ordinates in 2d
! with various density profiles
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc, ndpro, nsran
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, nblok, ierr
         real :: sxi, syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zero&
     &,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,z&
     &ero,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if (nsran /= 0) then
               call PRDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,&
     &ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc)
            endif
         endif
         end subroutine ipfdistr2
!
         subroutine ipvdistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,kstrt&
     &,nvp)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PVDISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npmax,&
     &nblok,kstrt,nvp,ierr)
         end subroutine ipvdistr2
!
!23456789*123456789*123456789*123456789*123456789*123456789*123456789*123
         subroutine ipvdistr2_arb(part,dist_func,npp,nps,vxmin,vxmax,   &
     &vymin,vymax,vdx,vdy,npx,npy,kstrt,nvp)
     
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
  
         implicit none
         
         interface
	     real function dist_func(v)
	     real v
	     end function dist_func
	     end interface
	    
	     
         integer :: npx, npy, kstrt, nvp
	     real :: vxmin,vxmax,vymin,vymax
         real :: vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1) 
         npmax = size(part,2)
         nblok = size(part,3)
!23456789*123456789*123456789*123456789*123456789*123456789*123456789*123
         call PVDISTR2_ARB(part,dist_func,npp,nps,vxmin,vxmax,          &
     & vymin,vymax,vdx,vdy,npx,npy,idimp,npmax,nblok,kstrt,nvp,ierr)
     
         end subroutine ipvdistr2_arb
!
         subroutine ipvdistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,&
     &npy,kstrt,nvp)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         call PVDISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idi&
     &mp,npmax,nblok,kstrt,nvp,ierr)
         end subroutine ipvdistrh2
!
         subroutine ipvfdistr2(part,nps,ampx,scalex,shiftx,ampy,scaley,s&
     &hifty,npx,npy,nx,ny,kstrt,nvp,ipbc,ndpro,nsran)
! calculates initial particle co-ordinates in 2d
! with various density profiles
! using parallel random number generator
         implicit none
         integer :: npx, npy, nx, ny, kstrt, nvp, ipbc, ndpro
         integer :: nsran
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         real :: sxi, syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zero&
     &,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,&
     &zero,zero,zero,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt,n&
     &vp,ndv,nvrp,ipbc,ierr)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,amp&
     &y,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR2(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1&
     &,ampy,syi,shifty,npx,npy,nx,ny,idimp,npmax,nblok,vranx,vrany,kstrt&
     &,nvp,ndv,nvrp,ipbc,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvfdistr2
!
         subroutine ipvvdistr2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,kstr&
     &t,nvp)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVVISTR2(part,npp,nps,vtx,vty,vdx,vdy,npx,npy,idimp,npmax,&
     &nblok,vranx,vrany,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistr2
!
         subroutine ipvvdistrh2(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,npy,kstrt,nvp)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, kstrt, nvp
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nvrp = (ndv - 1)/nvp + 1
         call PVVISTR2H(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idi&
     &mp,npmax,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistrh2
!
         subroutine ipbpistr2(part,bxy,npp,noff,qbm,nx,ny,ipbc,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for 2d
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(bxy,2); nypmx = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PGBDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,nx,ny,idimp,n&
     &pmax,nblok,nxv,nypmx,ipbc)
         else
            call PGBDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,nx,ny,idimp,n&
     &pmax,nblok,nxv,nypmx,ipbc)
         endif
         end subroutine ipbpistr2
!
         subroutine iprbpistr2(part,bxy,npp,noff,qbm,ci,nx,ny,ipbc,inord&
     &er)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for relativistic 2d
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, ci
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(bxy,2); nypmx = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PGRBDISTR2L(part,bxy(1,1,1,1),npp,noff,qbm,ci,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
         else
            call PGRBDISTR2L(part,bxy(1,2,2,1),npp,noff,qbm,ci,nx,ny,idi&
     &mp,npmax,nblok,nxv,nypmx,ipbc)
         endif
         end subroutine iprbpistr2
!
         subroutine ifedges2(edges,noff,nyp,ampy,scaley,shifty,ny,kstrt,&
     &nvp,nypmx,ipbc,ndpro,nterg,ierr,inorder)
! finds new 1d partitions from initial analytic distribution function
         implicit none
         integer :: ny, kstrt, nvp, nypmx, ipbc, ndpro, nterg, ierr
         real :: ampy, scaley, shifty
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: noff, nyp
! local data
         integer :: idps, nblok, nypmin, nypmax, order
         real :: syi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idps = size(edges,1); nblok = size(edges,2)
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
! uniform density
         if (ndpro==0) then
            call FEDGES2(edges,noff,nyp,FLDISTR1,zero,zero,zero,ny,nypmi&
     &n,nypmax,kstrt,nvp,nblok,idps,ipbc)
! linear density
         else if (ndpro==1) then
            call FEDGES2(edges,noff,nyp,FLDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! sinusoidal density
         else if (ndpro==2) then
            call FEDGES2(edges,noff,nyp,FSDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! gaussian density
         else if (ndpro==3) then
            call FEDGES2(edges,noff,nyp,FGDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
! hyperbolic secant squared density
         else if (ndpro==4) then
            call FEDGES2(edges,noff,nyp,FHDISTR1,ampy,syi,shifty,ny,nypm&
     &in,nypmax,kstrt,nvp,nblok,idps,ipbc)
         endif
         if (order==LINEAR) then
            nypmax = nypmax + 1
         else
            nypmax = nypmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         nterg = nypmin - 1
         end subroutine ifedges2
!
      end module pinit2d
