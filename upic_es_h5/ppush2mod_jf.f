! Modified push functions for doing v dot E stuff
module ppush2_jf
use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
use p0d, only: wtimer
use diag_jf
use pinit2d_jf
implicit none

private
public :: push_jf,djpost_jf,get_v_phi, gcjejpost_jf,gcjekejpost_jf,pgcjepost2_onlytrack_jf

interface
	subroutine PGSPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,&
		&nxv,nxyp,ipbc)
	 implicit none
	 integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc
	 real :: qbm, dt, ek
	 real, dimension(idimp,npmax,nblok) :: part
	 real, dimension(2,nxyp,nblok) :: fxy, vdotE
	 integer, dimension(nblok) :: npp, noff
	 end subroutine PGSPUSH2_jf
end interface

interface
	 subroutine PGPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
	 	&ipbc,nvexloc,nveyloc)
	 implicit none
	 ! for some stupid reason if I don't declare nvexloc and nveyloc with the n, then
	 ! they are passed as real and it doesn't work.  f77 is weird.
	 integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nvexloc,nveyloc
	 real :: qbm, dt, ek
	 real, dimension(idimp,npmax,nblok) :: part
	 real, dimension(2,nxv,nypmx,nblok) :: fxy,vdotE
	 integer, dimension(nblok) :: npp, noff
	 end subroutine
end interface

interface
	 subroutine GETVPHI(part,fxy,phi,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
	 	&ipbc,nvexloc,nveyloc)
	 implicit none
	 ! for some stupid reason if I don't declare nvexloc and nveyloc with the n, then
	 ! they are passed as real and it doesn't work.  f77 is weird.
	 integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nvexloc,nveyloc
	 real :: qbm, dt, ek
	 real, dimension(idimp,npmax,nblok) :: part
	 real, dimension(2,nxv,nypmx,nblok) :: fxy
	 real, dimension(nxv,nypmx,nblok) :: phi
	 integer, dimension(nblok) :: npp, noff
	 end subroutine
end interface


interface
	 subroutine PGPUSH2_jf_sum(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
	 	&ipbc,nvexloc,nveyloc)
	 implicit none
	 ! for some stupid reason if I don't declare nvexloc and nveyloc with the n, then
	 ! they are passed as real and it doesn't work.  f77 is weird.
	 integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, ipbc, nvexloc,nveyloc
	 real :: qbm, dt, ek
	 real, dimension(idimp,npmax,nblok) :: part
	 real, dimension(2,nxv,nypmx,nblok) :: fxy,vdotE
	 integer, dimension(nblok) :: npp, noff
	 end subroutine
end interface

interface
	 subroutine PGSJPOST22_jf(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,&
	 	&nvexloc,nveyloc)
	 implicit none
	 integer :: nx, ny, idimp, npmax, nblok, nxv, nxyp, ipbc,nvexloc,nveyloc
	 real :: qm, dt
	 real, dimension(idimp,npmax,nblok) :: part
	 real, dimension(3,nxyp,nblok) :: cu
	 integer, dimension(nblok) :: npp, noff
	 end subroutine
end interface

interface
	subroutine PGCJEJPOST2_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
	implicit none
	integer :: idimp, npmax, nblok, nxv, nypmx
	real :: qm, qbm, dt
	real, dimension(idimp,npmax,nblok) :: part
	real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, jE
	integer, dimension(nblok) :: npp, noff
	end subroutine
end interface

interface
	subroutine PGCJEJPOST2L_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
	implicit none
	integer :: idimp, npmax, nblok, nxv, nypmx
	real :: qm, qbm, dt
	real, dimension(idimp,npmax,nblok) :: part
	real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, jE
	integer, dimension(nblok) :: npp, noff
	end subroutine
end interface

interface
	subroutine PGCJEPOST2_onlytrackjf(part,fxy,npp,noff,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,kivx_loc)
	implicit none
	integer :: idimp, npmax, nblok, nxv, nypmx, kivx_loc
	real :: qm, qbm, dt
	real, dimension(idimp,npmax,nblok) :: part
	real, dimension(2,nxv,nypmx,nblok) :: fxy, jE
	integer, dimension(nblok) :: npp, noff
	end subroutine
end interface

interface
	subroutine PGCJEKEJPOST2_jf(part,fxy,npp,noff,cu,jE,kE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
	implicit none
	integer :: idimp, npmax, nblok, nxv, nypmx
	real :: qm, qbm, dt
	real, dimension(idimp,npmax,nblok) :: part
	real, dimension(2,nxv,nypmx,nblok) :: fxy, cu, jE
	real, dimension(nxv,nypmx,nblok) :: kE
	integer, dimension(nblok) :: npp, noff
	end subroutine
end interface


interface push_jf
	module procedure ipgpush2_jf
end interface
interface djpost_jf
	module procedure ipgjpost2_jf
end interface
interface gcjejpost_jf
	module procedure ipgcjejpost2_jf
end interface
interface gcjekejpost_jf
	module procedure ipgcjekejpost2_jf
end interface
interface pgcjepost2_onlytrack_jf
	module procedure ipgcjepost2_onlytrack_jf
end interface

 contains

      subroutine ipgpush2_jf(part,fxy,npp,noff,qbm,dt,ek,tpush,nx,ny,ipb&
     &c,inorder,popt)
! push particles with 2d electrostatic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, vdotE
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tpush,dtime,-1)
!For now, this push that does some fancy optimization doesn't work cause I can't figure it
!out and UCLA doesn't have the journal, I'll have to talk to Viktor about it.
!         call PGSPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npm&
!			     &ax,nblok,nxv,nxyp,ipbc)
!This push is straight forward but probably slower, but I can figure it out.
				 if (nvdotE_follow_part .eq. 0) then
					 call PGPUSH2_jf(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
							&ipbc,vEx_loc,vEy_loc)
				 else
					 call PGPUSH2_jf_sum(part,fxy,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
							&ipbc,vEx_loc,vEy_loc)
				 endif
! record time
         call wtimer(tpush,dtime)
       end subroutine ipgpush2_jf
       
! Current deposit function, but modified here to deposit not just the current but
! and quantity in the particle array, the point is to deposit the vdotE in the same way
! that the current is deposited.  
! Can use the lookahead functions here even though I don't understand them because only need 
! to change part(3,i,l) and part(4,i,l) to part(vEx_loc,i,l) and part(vEy_loc,i,l)
! This function is usually in pbpush2mod.f but put it here to reduce keep the number of object files
! low.

         subroutine ipgjpost2_jf(part,cu,npp,noff,qm,dt,tdjpost,nx,ny,ipbc,&
     &inorder,djopt)
! deposit current, 1d partition
         implicit none
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: cu
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: nx, ny, idimp, npmax, nblok, nxv, nypmx, nxyp, ipbc
         integer :: order, opt, nvdim
! npd = size of scratch buffers for vectorized charge deposition
         integer, parameter :: npd = 128, n12 = 12, n27 = 27
         integer, dimension(n27,npd,size(part,3)) :: nn
         real, dimension(n27,npd,size(part,3)) :: amxy
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(cu,2); nypmx = size(cu,3); nxyp = nxv*nypmx
         nvdim = size(cu,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdjpost,dtime,-1)
         if (nvdim == 3) then
         	print*,"Can't do 3d in velocity with vdotE_part diagnostic. Exiting."
         	stop
				if (order==LINEAR) then
					if (opt==LOOKAHEAD) then
						call PGSJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,&
		  &nblok,nxv,nxyp,ipbc)
					else if (opt==VECTOR) then
						call PGSJOST2XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc)
					else
						call PGJPOST2L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
					endif
				else
					if (opt==LOOKAHEAD) then
						call PGSJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc)
					else if (opt==VECTOR) then
						call PGSJOST2X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp&
		  &,npmax,nblok,nxv,nxyp,npd,n27,ipbc)
					else
						call PGJPOST2(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
					endif
				endif			
			else 
				if (order==LINEAR) then
         	print*,"Can't use linear interpolation with vdotE_part diagnostic. Exiting."
         	stop
					if (opt==LOOKAHEAD) then
						call PGSJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,&
		  &nblok,nxv,nxyp,ipbc)
					else if (opt==VECTOR) then
						call PGSJOST22XL(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,npd,n12,ipbc)
					else
						call PGJPOST22L(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
					endif
				else
					if (opt==LOOKAHEAD) then
					!This is the only one I modified for now.  Maybe someday I won't be lazy.
						call PGSJPOST22_jf(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nxyp,ipbc,&
							&vEx_loc,vEy_loc)
					else if (opt==VECTOR) then
         	print*,"Can't use vector with vdotE_part diagnostic. Exiting."
         	stop
						call PGSJOST22X(part,cu,npp,noff,nn,amxy,qm,dt,nx,ny,idimp&
		  &,npmax,nblok,nxv,nxyp,npd,n27,ipbc)
					else
         	print*,"Can't use this deposit with vdotE_part diagnostic. Exiting."
         	stop
						call PGJPOST22(part,cu,npp,noff,qm,dt,nx,ny,idimp,npmax,nblok,nxv,nypmx,ipbc)
					endif
				endif				
         endif
! record time
         call wtimer(tdjpost,dtime)
         end subroutine ipgjpost2_jf
         
      subroutine get_v_phi(part,fxy,phi,npp,noff,qbm,dt,ek,tpush,nx,ny,ipb&
     &c,inorder,popt)
! push particles with 2d electrostatic fields, 1d partition
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy
         real, dimension(:,:,:), pointer :: phi
         integer, dimension(:), pointer :: npp, noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nxyp, order, opt
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(phi,1); nypmx = size(phi,2); nxyp = nxv*nypmx
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tpush,dtime,-1)

				 call GETVPHI(part,fxy,phi,npp,noff,qbm,dt,ek,nx,ny,idimp,npmax,nblok,nxv,nypmx,&
						&ipbc,vEx_loc,vEy_loc)
! record time
         call wtimer(tpush,dtime)
       end subroutine get_v_phi

         subroutine ipgcjejpost2_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,tdcjpost,inorder)
! deposit current density with 2d electrostatic fields, 1d partition
! Modified to also deposit j.E into jE
         implicit none
         real :: qm, qbm, dt, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, cu, jE
         integer, dimension(:), pointer :: npp, noff
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tdc,dtime,-1)
         if (order==LINEAR) then
            call PGCJEJPOST2L_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
         else
            call PGCJEJPOST2_jf(part,fxy,npp,noff,cu,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgcjejpost2_jf

         subroutine ipgcjekejpost2_jf(part,fxy,npp,noff,cu,jE,kE,qm,qbm,dt,tdcjpost,inorder)
! deposit current density with 2d electrostatic fields, 1d partition
! Modified to also deposit j.E into jE
         implicit none
         real :: qm, qbm, dt, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, cu, jE
         real, dimension(:,:,:), pointer :: kE
         integer, dimension(:), pointer :: npp, noff
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tdc,dtime,-1)
         call PGCJEKEJPOST2_jf(part,fxy,npp,noff,cu,jE,kE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx)
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgcjekejpost2_jf

         subroutine ipgcjepost2_onlytrack_jf(part,fxy,npp,noff,jE,qm,qbm,dt,tdcjpost,inorder)
! deposit current density with 2d electrostatic fields, 1d partition
! Modified to also deposit j.E into jE
         implicit none
         real :: qm, qbm, dt, tdcjpost
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:), pointer :: fxy, jE
         integer, dimension(:), pointer :: npp, noff
         integer, optional :: inorder
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, order
         real :: tdc
         double precision :: dtime
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)
         nxv = size(fxy,2); nypmx = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tdc,dtime,-1)
         if (order==LINEAR) then
            print*,"no linear PGCJEPOST2_jf"
         else
            call PGCJEPOST2_onlytrackjf(part,fxy,npp,noff,jE,qm,qbm,dt,idimp,npmax,nblok,nxv,nypmx,kivx_loc)
         endif
! record time
         call wtimer(tdc,dtime)
         tdcjpost = tdcjpost + tdc
         end subroutine ipgcjepost2_onlytrack_jf


end module ppush2_jf
