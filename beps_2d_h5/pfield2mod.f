!-----------------------------------------------------------------------
!
      module pfield2d
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: july 25, 2008
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: PSCGUARD2, PSGUARD2, PSCGUARD2L, PSGUARD2L
      public :: PACGUARD2X, PAGUARD2X, PACGUARD2XL, PAGUARD2XL
      public :: cguard, bguard, sguard, aguard, zguard
      public :: pois_init, pois, pois3, cuperp, bpois, sbpois
      public :: ibpois, maxwel, emfield, emfieldr, apois, avpot
      public :: gtmodes, ptmodes, addqei, baddext
      public :: imoment, ipdivf2, ipgradf2, ipcurlf2, ipcurlf22
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PZGUARD2(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBGUARD2XL(bxy,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: bxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0, zj0
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
         implicit none
         real :: xj0, yj0
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
         implicit none
         real :: qi0
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(3,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(2,nxe,nypmx,nblok) :: cu
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PZGUARD2L(q,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer nx, nxe, nypmx, nblok
         real, dimension(nxe,nypmx,nblok) :: q
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,&
     &nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q, fx, fy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,ny&
     &v,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(2,nyv,kxp,jblok) :: fxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,ny&
     &v,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, we
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(nyv,kxp,jblok) :: q
         complex, dimension(2,nyv,kxp,jblok) :: fxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: f
         complex, dimension(nyv,kxp,jblok) :: df
         end subroutine
      end interface
      interface
         subroutine PGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, ndim, nyv, kxp, jblok
         complex, dimension(nyv,kxp,jblok) :: df
         complex, dimension(3,nyv,kxp,jblok) :: f
         end subroutine
      end interface
      interface
         subroutine PCURLF2(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: f, g
         end subroutine
      end interface
      interface
         subroutine PCURLF22(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: f
         complex, dimension(nyv,kxp,jblok) :: g
         end subroutine
      end interface
      interface
         subroutine PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: cu
         end subroutine
      end interface
      interface
         subroutine PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(2,nyv,kxp,jblok) :: cu
         end subroutine
      end interface
      interface
         subroutine PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,ks&
     &trt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: cu
         complex, dimension(3,nyv,kxp,jblok) :: bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PBPOISP22(cu,bxy,bz,isign,ffc,ax,ay,affp,ci,wm,nx,ny&
     &,kstrt,nyv,kxp,jblok,nyhd)
         implicit none
         real :: ax, ay, affp, ci, wm
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(2,nyv,kxp,jblok) :: cu
         complex, dimension(2,nyv,kxp,jblok) :: bxy
         complex, dimension(nyv,kxp,jblok) :: bz
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine IPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,jblo&
     &k,nyhd)
         implicit none
         real :: ci, wm
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: cu
         complex, dimension(3,nyv,kxp,jblok) :: bxy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         implicit none
         real :: affp, ci, dt, wf, wm
         integer :: nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: exy, bxy, cu
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,jblo&
     &k,nyhd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp, jblok, nyhd
         complex, dimension(3,nyv,kxp,jblok) :: fxy, exy
         complex, dimension(nyhd,kxp,jblok) :: ffc
         end subroutine
      end interface
      interface
         subroutine PEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,j2&
     &blok,nyd)
         implicit none
         integer :: isign, nx, ny, kstrt, nyv, kxp2, j2blok, nyd
         real, dimension(3,nyv,kxp2+1,j2blok) :: fxy, exy
         complex, dimension(nyd,kxp2,j2blok) :: ffd
         end subroutine
      end interface
      interface
         subroutine PAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp,jblok)
         implicit none
         integer :: nx, ny, kstrt, nyv, kxp, jblok
         complex, dimension(3,nyv,kxp,jblok) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine PGTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,&
     &kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
         integer :: nt, modesxpd, modesyd
         complex, dimension(nyv,kxp,jblok) :: pot
         complex, dimension(nt,modesyd,modesxpd,jblok) :: pott
         end subroutine
      end interface
      interface
         subroutine PPTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,&
     &kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, kstrt, nyv, kxp, jblok
         integer :: nt, modesxpd, modesyd
         complex, dimension(nyv,kxp,jblok) :: pot
         complex, dimension(nt,modesyd,modesxpd,jblok) :: pott
         end subroutine
      end interface
      interface
         subroutine PGTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,ks&
     &trt,nyv,kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp
         integer :: jblok, nt, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp,jblok) :: vpot
         complex, dimension(nt,ndim,modesyd,modesxpd,jblok) :: vpott
         end subroutine
      end interface
      interface
         subroutine PPTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,ks&
     &trt,nyv,kxp,jblok,nt,modesxpd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, kstrt, nyv, kxp
         integer :: jblok, nt, modesxpd, modesyd
         complex, dimension(ndim,nyv,kxp,jblok) :: vpot
         complex, dimension(nt,ndim,modesyd,modesxpd,jblok) :: vpott
         end subroutine
      end interface
      interface
         subroutine PADDQEI2(qe,qi,nyp,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
!        real, dimension(*) :: qe, qi
         real :: qe, qi
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PADDQEI2X(qe,qi,nyp,qbme,qbmi,wpmax,wpmin,nx,nxe,nyp&
     &mx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: qbme, qbmi ,wpmax, wpmin
!        real, dimension(*) :: qe, qi
         real :: qe, qi
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBADDEXT2(bxy,nyp,omx,omy,omz,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: omx, omy, omz
!        real, dimension(*) :: bxy
         real :: bxy
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PBADDEXT22(bz,nyp,omz,nx,nxe,nypmx,nblok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: omz
!        real, dimension(*) :: bz
         real :: bz
         integer,  dimension(nblok) :: nyp
         end subroutine
      end interface
      interface
         subroutine PIMOMENT2(qi,fxy,nyp,pxi,pyi,pzi,dt,nx,nxe,nypmx,nbl&
     &ok)
         implicit none
         integer :: nx, nxe, nypmx, nblok
         real :: pxi, pyi, pzi, dt
!        real, dimension(*) :: qi
!        real, dimension(*) :: fxy
         real :: qi, fxy
         integer, dimension(nblok) :: nyp
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface cguard
         module procedure ipcguard2x
         module procedure ipdguard2x
      end interface
!
      interface bguard
         module procedure ipbguard2x
      end interface
!
      interface sguard
         module procedure ipscguard2
         module procedure ipsguard2
      end interface
!      
      interface aguard
         module procedure ipacguard2x
         module procedure ipaguard2x
      end interface
!
      interface zguard
         module procedure ipzguard2
      end interface
!
       interface pois_init
         module procedure ippois22init
      end interface
! 
      interface pois
         module procedure ippois2
         module procedure ipspois2
         module procedure ippois22
      end interface
!
      interface pois3
         module procedure ippois23
      end interface
!
      interface cuperp
         module procedure ipcuperp2
      end interface
!
      interface bpois
         module procedure jpbpois23
      end interface
!
      interface sbpois
         module procedure ipsbpois23
      end interface
!
      interface apois
         module procedure ipapois23
      end interface
!
      interface ibpois
         module procedure iipbpois23
      end interface
!
      interface maxwel
         module procedure ipmaxwel2
      end interface
!
      interface emfield
         module procedure ipemfield2
      end interface
!
      interface emfieldr
         module procedure ipemfieldr2
      end interface
!
      interface avpot
         module procedure ipavpot23
      end interface
!
      interface gtmodes
         module procedure ipgtmodes2
      end interface
!
      interface ptmodes
         module procedure ipptmodes2
      end interface
!
      interface addqei
         module procedure ipaddqei2
         module procedure ipaddqei2x
      end interface
!
      interface baddext
         module procedure ipbaddext2
      end interface
!
      interface imoment
         module procedure ipimoment2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ipcguard2x(fxy,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(fxy,2); nypmx = size(fxy,3); nblok = size(fxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==LINEAR) then
               call PDGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PDGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         case (2)
            if (order==LINEAR) then
               call PCGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PCGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PBGUARD2XL(fxy,nyp,nx,nxe,nypmx,nblok)
            else
               call PBGUARD2X(fxy,nyp,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipcguard2x
!
         subroutine ipdguard2x(q,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PDGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PDGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipdguard2x
!
         subroutine ipbguard2x(bxy,nyp,nx,inorder)
! copy guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(bxy,2); nypmx = size(bxy,3); nblok = size(bxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PBGUARD2XL(bxy,nyp,nx,nxe,nypmx,nblok)
         else
            call PBGUARD2X(bxy,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipbguard2x
!
         subroutine ipscguard2(cu,nyp,xj0,yj0,zj0,nx,inorder)
! initialize periodic 2d vector field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call PSCGUARD22L(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
            else
               call PSCGUARD22(cu,nyp,xj0,yj0,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PSCGUARD2L(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
            else
               call PSCGUARD2(cu,nyp,xj0,yj0,zj0,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipscguard2
!
         subroutine ipsguard2(q,nyp,qi0,nx,inorder)
! initialize periodic 2d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PSGUARD2L(q,nyp,qi0,nx,nxe,nypmx,nblok)
         else
            call PSGUARD2(q,nyp,qi0,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipsguard2
!
         subroutine ipacguard2x(cu,nyp,nx,inorder)
! add guard cells in x for non-uniform, periodic 2d vector data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(cu,2); nypmx = size(cu,3); nblok = size(cu,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call PACGUARD22XL(cu,nyp,nx,nxe,nypmx,nblok)
            else
               call PACGUARD22X(cu,nyp,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PACGUARD2XL(cu,nyp,nx,nxe,nypmx,nblok)
            else
               call PACGUARD2X(cu,nyp,nx,nxe,nypmx,nblok)
            endif
         end select
         end subroutine ipacguard2x
!
         subroutine ipaguard2x(q,nyp,nx,inorder)
! add guard cells in x for non-uniform, periodic 2d scalar data
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PAGUARD2XL(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PAGUARD2X(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipaguard2x
!
         subroutine ipzguard2(q,nyp,nx,inorder)
! zero out guard cells in periodic 2d scalar field
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: q
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(q,1); nypmx = size(q,2); nblok = size(q,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PZGUARD2L(q,nyp,nx,nxe,nypmx,nblok)
         else
            call PZGUARD2(q,nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipzguard2
!
         subroutine ippois2init(ffc,ax,ay,affp,nx,ny,kstrt)
! initialize 2d periodic poisson solver
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 0, nyv, kxp, jblok, nyhd
         real :: we
         complex, dimension(1,1,1) :: q, fx, fy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ippois2init
!
         subroutine ippois2(q,fx,ffc,we,nx,ny,kstrt)
! poisson solver for periodic 2d potential
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, fx, ffc
! local data
         integer :: isign = 1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: fy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ippois2
!
         subroutine ipspois2(q,fy,ffc,nx,ny,kstrt)
! smoother for periodic 2d scalar field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: q, fy, ffc
! local data
         integer :: isign = 2, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, we
         complex, dimension(1,1,1) :: fx
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kx&
     &p,jblok,nyhd)
         end subroutine ipspois2
!
         subroutine ippois22init(ffc,ax,ay,affp,nx,ny,kstrt)
! initialize 2d periodic electric field solver
         implicit none
         integer :: nx, ny, kstrt
         real :: ax, ay, affp
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 0, nyv, kxp, jblok, nyhd
         real :: we
         complex, dimension(1,1,1) :: q
         complex, dimension(2,1,1,1) :: fxy
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
         end subroutine ippois22init
!
         subroutine ippois22(q,fxy,ffc,we,nx,ny,kstrt)
! poisson solver for periodic 2d electric field
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, ffc
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
         end subroutine ippois22
!
         subroutine ippois23(q,fxy,ffc,we,nx,ny,kstrt)
! poisson solver for periodic 2-1/2d electric field
         implicit none
         integer :: nx, ny, kstrt
         real :: we
         complex, dimension(:,:,:), pointer :: q, ffc
         complex, dimension(:,:,:,:), pointer :: fxy
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         nyv = size(q,1)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PPOIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nyv,kxp,&
     &jblok,nyhd)
         end subroutine ippois23
!
         subroutine ipdivf2(f,df,nx,ny,kstrt)
! calculates the divergence of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: df
! local data
         integer :: ndim, nyv, kxp, jblok
         ndim = size(f,1)
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PDIVF2(f,df,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         end subroutine ipdivf2
!
         subroutine ipgradf2(df,f,nx,ny,kstrt)
! calculates the gradient of periodic 2d scalar field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:), pointer :: df
         complex, dimension(:,:,:,:), pointer :: f
! local data
         integer :: ndim, nyv, kxp, jblok
         ndim = size(f,1)
         nyv = size(df,1); kxp = size(df,2); jblok = size(df,3)
         call PGRADF2(df,f,nx,ny,kstrt,ndim,nyv,kxp,jblok)
         end subroutine ipgradf2
!
         subroutine ipcurlf2(f,g,nx,ny,kstrt)
! calculates the curl of periodic 2-1/2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f, g
! local data
         integer :: nyv, kxp, jblok
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PCURLF2(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipcurlf2
!
         subroutine ipcurlf22(f,g,nx,ny,kstrt)
! calculates the curl of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: f
         complex, dimension(:,:,:), pointer :: g
! local data
         integer :: nyv, kxp, jblok
         nyv = size(f,2); kxp = size(f,3); jblok = size(f,4)
         call PCURLF22(f,g,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipcurlf22
!
         subroutine ipcuperp2(cu,nx,ny,kstrt)
! calculates the transverse part of periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu
! local data
         integer :: nyv, kxp, jblok
         nyv = size(cu,2); kxp = size(cu,3); jblok = size(cu,4)
         select case(size(cu,1))
         case (2)
            call PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
         case (3)
            call PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
         end select
         end subroutine ipcuperp2
!
         subroutine jpbpois23(cu,bxy,ffc,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = -1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         complex, dimension(2,1,1,1) :: bxy0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,bxy0,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,&
     &kstrt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
         end subroutine jpbpois23
!
         subroutine ipsbpois23(cu,bxy,ffc,nx,ny,kstrt)
! smoother for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: cu, bxy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 2, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp, ci, wm
         complex, dimension(1,1,1) :: bz0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,bxy,bz0,isign,ffc,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
         end subroutine ipsbpois23
!
         subroutine ipapois23(cu,axy,ffc,ci,wm,nx,ny,kstrt)
! calculates static vector potential for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:,:), pointer :: cu, axy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: isign = 1, nyv, kxp, jblok, nyhd
         real :: ax, ay, affp
         complex, dimension(1,1,1) :: bz0
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         select case(size(cu,1))
         case (2)
            call PBPOISP22(cu,axy,bz0,isign,ffc,ax,ay,affp,ci,wm,nx,ny,k&
     &strt,nyv,kxp,jblok,nyhd)
         case (3)
            call PBPOISP23(cu,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,kstrt&
     &,nyv,kxp,jblok,nyhd)
         end select
         end subroutine ipapois23
!
         subroutine iipbpois23(cu,bxy,ffc,ci,wm,nx,ny,kstrt)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: ci, wm
         complex, dimension(:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:), pointer :: cu, bxy
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call IPBPOISP23(cu,bxy,ffc,ci,wm,nx,ny,kstrt,nyv,kxp,jblok,nyhd&
     &)
         end subroutine iipbpois23
!
         subroutine ipmaxwel2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstr&
     &t)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, kstrt
         real :: affp, ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: ffc
         complex, dimension(:,:,:,:), pointer :: exy, bxy, cu
! local data
         integer :: nyv, kxp, jblok, nyhd
         nyv = size(cu,2)
         nyhd = size(ffc,1); kxp = size(ffc,2); jblok = size(ffc,3)
         call PMAXWEL2(exy,bxy,cu,ffc,affp,ci,dt,wf,wm,nx,ny,kstrt,nyv,k&
     &xp,jblok,nyhd)
         end subroutine ipmaxwel2
!
         subroutine ipemfield2(fxy,exy,ffc,isign,nx,ny,kstrt)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny, kstrt, nyhd
         complex, dimension(:,:,:,:), pointer :: fxy, exy
         complex, dimension(:,:,:), pointer :: ffc
! local data
         integer :: nyv, kxp, jblok
         nyv = size(fxy,2); kxp = size(fxy,3); jblok = size(fxy,4)
         nyhd = size(ffc,1) 
         call PEMFIELD2(fxy,exy,ffc,isign,nx,ny,kstrt,nyv,kxp,jblok,nyhd&
     &)
         end subroutine ipemfield2 
!
         subroutine ipemfieldr2(fxy,exy,ffd,isign,nx,ny,kstrt)
! combines and smooths real 2d vector fields
         implicit none
         integer :: isign, nx, ny, kstrt, nyd
         real, dimension(:,:,:,:), pointer :: fxy, exy
         complex, dimension(:,:,:), pointer :: ffd
! local data
         integer :: nyv, kxp2, j2blok
         nyv = size(fxy,2); j2blok = size(fxy,4)
         nyd = size(ffd,1); kxp2 = size(ffd,2); 
         call PEMFIELDR2(fxy,exy,ffd,isign,nx,ny,kstrt,nyv,kxp2,j2blok,n&
     &yd)
         end subroutine ipemfieldr2 
!
         subroutine ipavpot23(bxy,axy,nx,ny,kstrt)
! calculates periodic 2d vector potential from magnetic field
         implicit none
         integer :: nx, ny, kstrt
         complex, dimension(:,:,:,:), pointer :: bxy, axy
! local data
         integer :: nyv, kxp, jblok
         nyv = size(bxy,2); kxp = size(bxy,3); jblok = size(bxy,4)
         call PAVPOT23(bxy,axy,nx,ny,kstrt,nyv,kxp,jblok)
         end subroutine ipavpot23
!
         subroutine ipgtmodes2(pot,pott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes from periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:), pointer :: pot
         complex, dimension(:,:,:), pointer :: pott
! local data
         integer :: it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         nyv = size(pot,1); kxp = size(pot,2); jblok = size(pot,3)
         modesyd = size(pott,1); modesxpd = size(pott,2)
         call PGTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp,jb&
     &lok,nt,modesxpd,modesyd)
         end subroutine ipgtmodes2
!
         subroutine ipptmodes2(pot,pott,nx,ny,modesx,modesy,kstrt)
! extracts lowest order modes to periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy, kstrt
         complex, dimension(:,:,:), pointer :: pot
         complex, dimension(:,:,:), pointer :: pott
! local data
         integer :: it, nyv, kxp, jblok, nt, modesxpd, modesyd
         it = 1; nt = 1
         nyv = size(pot,1); kxp = size(pot,2); jblok = size(pot,3)
         modesyd = size(pott,1); modesxpd = size(pott,2)
         call PPTMODES2(pot,pott,nx,ny,it,modesx,modesy,kstrt,nyv,kxp,jb&
     &lok,nt,modesxpd,modesyd)
         end subroutine ipptmodes2
!
         subroutine ipaddqei2(qe,qi,nyp,nx,inorder)
! adds electron and ion densities
         implicit none
         integer :: nx
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: qe, qi
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(qe,1); nypmx = size(qe,2); nblok = size(qe,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PADDQEI2(qe(1,1,1),qi(1,1,1),nyp,nx,nxe,nypmx,nblok)
         else
            call PADDQEI2(qe(2,2,1),qi(2,2,1),nyp,nx,nxe,nypmx,nblok)
         endif
         end subroutine ipaddqei2
!
         subroutine ipaddqei2x(qe,qi,qbme,qbmi,wpmax,wpmin,nyp,nx,inorde&
     &r)
! adds electron and ion densities, and calculates maximum and minimum
! plasma frequency
         implicit none
         integer :: nx
         integer, optional :: inorder
         real :: qbme, qbmi, wpmax, wpmin
         real, dimension(:,:,:), pointer :: qe, qi
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         double precision, dimension(2) :: sum2, work2
         nxe = size(qe,1); nypmx = size(qe,2); nblok = size(qe,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PADDQEI2X(qe(1,1,1),qi(1,1,1),nyp,qbme,qbmi,wpmax,wpmin&
     &,nx,nxe,nypmx,nblok)
         else
            call PADDQEI2X(qe(2,2,1),qi(2,2,1),nyp,qbme,qbmi,wpmax,wpmin&
     &,nx,nxe,nypmx,nblok)
         endif
! maximum/minimum over the y direction
         sum2(1) = wpmax
         sum2(2) = -wpmin
         call PDMAX(sum2,work2,2,1)
         wpmax = sum2(1)
         wpmin = -sum2(2)
         end subroutine ipaddqei2x
!
         subroutine ipbaddext2(bxy,nyp,omx,omy,omz,nx,inorder)
! adds constant to magnetic field
         implicit none
         integer :: nx
         real :: omx, omy, omz
         integer, optional :: inorder
         real, dimension(:,:,:,:), pointer :: bxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         nxe = size(bxy,2); nypmx = size(bxy,3); nblok = size(bxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call PBADDEXT22(bxy(1,1,1,1),nyp,omz,nx,nxe,nypmx,nblok)
            else
               call PBADDEXT22(bxy(1,2,2,1),nyp,omz,nx,nxe,nypmx,nblok)
            endif
         case (3)
            if (order==LINEAR) then
               call PBADDEXT2(bxy(1,1,1,1),nyp,omx,omy,omz,nx,nxe,nypmx,&
     &nblok)
            else
               call PBADDEXT2(bxy(1,2,2,1),nyp,omx,omy,omz,nx,nxe,nypmx,&
     &nblok)
            endif
         end select
         end subroutine ipbaddext2
!
         subroutine ipimoment2(qi,fxy,nyp,id0,iunit,px,py,pz,dt,wx,wy,wz&
     &,nx,inorder)
! calculate ion momentum from integral of qi*fxy,
! and prints it, and adds it to total momentum, for 2 or 2-1/2d code
         implicit none
         integer :: nx, id0, iunit
         real :: px, py, pz, dt, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: qi
         real, dimension(:,:,:,:), pointer :: fxy
         integer, dimension(:), pointer :: nyp
! local data
         integer :: nxe, nypmx, nblok, order
         double precision, dimension(3) :: sum3, work3
  995    format (' ion momentum = ',3e14.7)
         nxe = size(fxy,2); nypmx = size(fxy,3); nblok = size(fxy,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! calculate and print ion momentum
         select case(size(fxy,1))
         case (2)
            if (order==LINEAR) then
               call PIMOMENT22(qi(1,1,1),fxy(1,1,1,1),nyp,px,py,dt,nx,nx&
     &e,nypmx,nblok)
            else
               call PIMOMENT22(qi(2,2,1),fxy(1,2,2,1),nyp,px,py,dt,nx,nx&
     &e,nypmx,nblok)
            endif
            pz = 0.0
         case (3)
            if (order==LINEAR) then
               call PIMOMENT2(qi(1,1,1),fxy(1,1,1,1),nyp,px,py,pz,dt,nx,&
     &nxe,nypmx,nblok)
            else
               call PIMOMENT2(qi(2,2,1),fxy(1,2,2,1),nyp,px,py,pz,dt,nx,&
     &nxe,nypmx,nblok)
            endif
         end select
! sum over the y and z directions
         sum3(1) = px
         sum3(2) = py
         sum3(3) = pz
         call PDSUM(sum3,work3,3,1)
         px = sum3(1)
         py = sum3(2)
         pz = sum3(3)
         if (id0==0) write (iunit,995) px, py, pz
! add to total momentum
         wx = wx + px
         wy = wy + py
         wz = wz + pz
         end subroutine ipimoment2
!
      end module pfield2d

