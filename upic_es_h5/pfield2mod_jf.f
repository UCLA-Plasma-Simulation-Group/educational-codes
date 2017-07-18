!-----------------------------------------------------------------------
!
      module pfield2d_jf
!
! Fortran90 interface to 2d parallel PIC Fortran77 library pfield2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: february 25, 2006
!
      implicit none
      private
      public :: cuperp_jf
      
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

      interface cuperp_jf
      	module procedure ipcuperp2_jf
      end interface

		contains
      
		subroutine ipcuperp2_jf(cu,nx,ny,kstrt)
! calculates the transverse part of periodic 2d vector field
		implicit none
		integer :: nx, ny, kstrt
		complex, dimension(:,:,:,:), pointer :: cu
! local data
		integer :: nyv, kxp, jblok
		nyv = size(cu,2); kxp = size(cu,3); jblok = size(cu,4)
		if (size(cu,1) == 3) then 
			call PCUPERP2(cu,nx,ny,kstrt,nyv,kxp,jblok)
		else if (size(cu,1) == 2) then
			call PCUPERP22(cu,nx,ny,kstrt,nyv,kxp,jblok)
		endif
		end subroutine ipcuperp2_jf
!

      end module pfield2d_jf