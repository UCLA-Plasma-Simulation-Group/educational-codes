!-----------------------------------------------------------------------
!
      module diag_jf
      
! Has diagnostic subroutines written by Jay Fahlen for use with 
! Viktor Decyk's 2d BEPS code
! written Jan. 2008
			
      use pinit2d, only: idrun, indx, indy, ntp, psolve, tend, dt,npx&
     &,npy
      use p2d_jf
      use p0d
      use pinit2d_jf
      use hdf_write_jf
      use par_track_jf
      implicit none
			include 'mpif.h'
      private
      public :: idrun, indx, indy, ntp!, modesx, modesy
      public :: tend, dt
      public :: ntfield
      public :: write_jf_diag_file,write_Eoft_aty,mkdir_structure
      public :: phase_slices_x,phase_slices_y,phase_space_vxy_vs_x_and_y
		! name convention is xx mean vx vs x, xy means vx vs y
     	public :: nlines
      public :: j_diag_deposit,vdotE
      public :: addtag
      public :: store_init_vel,calc_ene_dist,calc_vx_dist,calc_vy_dist
      public :: kivx_loc,kivy_loc, vEx_loc,vEy_loc
      public :: set_part_array_vdotE,n_jE_u_poynt
      public :: phase_vx_vy
      public :: phase_space_vxy_vs_x,phase_space_vxy_vs_y
      public :: phase_space_vxy_vs_x_ion, phase_space_vxy_vs_y_ion
      public :: get_ene_in_driver_region,calc_ene_three_regions,get_Py_sumover_x,get_1D_sumover_x
      public :: get_2D_sumover_x,get_2D_sumover_x_xcomp
      public :: bin_Ex_and_Ey
      public :: phase_vx_vy_x_y_movie_at_each_dump
      public :: linesfunit,nlinerec
      public :: nt_write_jE_sumover_x_funit,nt_write_U_sumover_x_funit,nt_write_div_P_vs_y_funit
      public :: nt_write_Py_sumover_x_funit,nt_write_Ux_sumover_x_funit,nt_write_U_sumover_x_fromE_funit
      public :: nt_write_S_sumover_x_funit,nt_write_kE_sumover_x_funit
      public :: nt_write_jE_onlytracked_sumover_x_funit
      public :: init_through_wave,calc_through_wave_dist,div_finite_difference
      
      save
      integer :: addtag = 0
      integer :: kivx_loc=0,kivy_loc=0, vEx_loc=0,vEy_loc=0
			integer :: set_part_array_vdotE = 0, n_jE_u_poynt = 0
      integer,dimension(:),allocatable :: nlinerec
      integer,dimension(5) :: linesfunit
      integer :: nlines
      integer :: nt_write_jE_sumover_x_funit,nt_write_U_sumover_x_funit,nt_write_div_P_vs_y_funit
      integer :: nt_write_Py_sumover_x_funit,nt_write_Ux_sumover_x_funit,nt_write_U_sumover_x_fromE_funit
      integer :: nt_write_S_sumover_x_funit,nt_write_kE_sumover_x_funit
      integer :: nt_write_jE_onlytracked_sumover_x_funit

!     HDF Function declaration.
			integer,external :: sfstart, sfcreate, sfwdata, sfsdtstr, sfsdmstr
			integer,external :: sfdimid, sfsdmname, sfsdscale, sfsblsz 
			integer,external :: sfendacc, sfend, sfsnatt, sfscompress

!     HDF Constant declaration
			integer, parameter :: DFACC_CREATE = 4
			integer, parameter :: DFACC_WRITE = 2
			integer, parameter :: DFNT_CHAR = 4
			integer, parameter :: DFNT_FLOAT32 = 5
			integer, parameter :: DFNT_FLOAT64 = 6
			integer, parameter :: DFNT_INT32 = 24
			integer, parameter :: COMP_CODE_DEFLATE = 4
			integer, parameter :: SD_UNLIMITED = 0
			integer, parameter :: SD_FILL = 0

            
!      character(len=*)::idrun
!      integer::modesx=11,modesy=11
!      real::t0,tend,dt,ceng
      
      contains
      	subroutine mkdir_structure(idproc)
      		implicit none
      		integer :: idproc
      		character(len=4) :: temp
      		
				integer :: ierr, i
      		
      		if (idproc == 0) then
					
						if (ntden .ne. 0) then
							call mkdir_f('./DIAG/Den'//char(0),ierr)
							call mkdir_f('./DIAG/IDen'//char(0),ierr)
						endif
						if (ntfield .ne. 0) then
							call mkdir_f('./DIAG/Ex'//char(0),ierr)
							call mkdir_f('./DIAG/Ey'//char(0),ierr)
							call mkdir_f('./DIAG/Ez'//char(0),ierr)
						endif
						if (nt_b_field .ne. 0) then
							call mkdir_f('./DIAG/Bx'//char(0),ierr)
							call mkdir_f('./DIAG/By'//char(0),ierr)
							call mkdir_f('./DIAG/Bz'//char(0),ierr)
						endif
						if (nt_write_grad_phi .ne. 0) then
							call mkdir_f('./DIAG/grad_phi_x'//char(0),ierr)
							call mkdir_f('./DIAG/grad_phi_y'//char(0),ierr)
						endif
						if (ntESPoynt .ne. 0) then
							call mkdir_f('./DIAG/ESPoynt_x/'//char(0),ierr)
							call mkdir_f('./DIAG/ESPoynt_y/'//char(0),ierr)
						endif
						if (ntESPoynt_int .ne. 0) then
							call mkdir_f('./DIAG/ESPoynt_int_x/'//char(0),ierr)
							call mkdir_f('./DIAG/ESPoynt_int_y/'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_ESPoynt_int/'//char(0),ierr)							
						endif
						if (ntj .ne. 0) then
							call mkdir_f('./DIAG/jx'//char(0),ierr)
							call mkdir_f('./DIAG/jy'//char(0),ierr)
						endif
						if (nttrack .ne. 0) then
							call mkdir_f('./DIAG/TRACKS'//char(0),ierr)
						endif
						if (keep_init_vel .ne. 0) then
							call mkdir_f('./DIAG/EneBin'//char(0),ierr)
							call mkdir_f('./DIAG/VxBin'//char(0),ierr)
							call mkdir_f('./DIAG/VyBin'//char(0),ierr)
						endif
						if (nvdotE_int .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_int'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_int'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_vdotE_int'//char(0),ierr)
						endif
						if (nvdotE_part .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_part'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_part'//char(0),ierr)
						endif
						if (ntvdotE .ne. 0) then
							call mkdir_f('./DIAG/vdotEx'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy'//char(0),ierr)
						endif
						if (ntphxy .ne. 0) then
							call mkdir_f('./DIAG/fvx_xy'//char(0),ierr)
							call mkdir_f('./DIAG/fvy_xy'//char(0),ierr)
						endif
						if (nvdotE_follow_part .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_follow_part'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_follow_part'//char(0),ierr)
						endif
						if (ntfield_ene_int .ne. 0) then
							call mkdir_f('./DIAG/u_x_int'//char(0),ierr)
							call mkdir_f('./DIAG/u_y_int'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_u_int'//char(0),ierr)
						endif
						if (ntraw .ne. 0) then
							call mkdir_f('./DIAG/RAW'//char(0),ierr)
						endif
						if (nt_div_ESPoynt_int .ne. 0) then
							call mkdir_f('./DIAG/div_ESPoynt_int'//char(0),ierr)
						endif
						if (nt_div_ESPoynt .ne. 0) then
							call mkdir_f('./DIAG/div_ESPoynt'//char(0),ierr)
						endif
						if (nt_phase_vx_vy .ne. 0) then
							call mkdir_f('./DIAG/VxVy'//char(0),ierr)
						endif
						if (nt_bin_E .ne. 0) then
							call mkdir_f('./DIAG/bin_Ex'//char(0),ierr)
							call mkdir_f('./DIAG/bin_Ey'//char(0),ierr)
						endif
						if (nt_vx_vy_speed .ne. 0) then
							call mkdir_f('./DIAG/vx_vy_speed'//char(0),ierr)
						endif
						if (nt_dEdt .ne. 0) then
							call mkdir_f('./DIAG/dExdt'//char(0),ierr)
							call mkdir_f('./DIAG/dEydt'//char(0),ierr)
						endif
						if (nt_dEdt .ne. 0) then
							call mkdir_f('./DIAG/grad_phi_x'//char(0),ierr)
							call mkdir_f('./DIAG/grad_phi_y'//char(0),ierr)
						endif
						if (nt_kE .ne. 0) then
							call mkdir_f('./DIAG/kE'//char(0),ierr)
						endif

						if (nphxx .ne. 0) then
							call mkdir_f('./DIAG/Vx_x'//char(0),ierr)
							call mkdir_f('./DIAG/IVx_x'//char(0),ierr)
						endif
						if (nphyx .ne. 0) then
							call mkdir_f('./DIAG/Vy_x'//char(0),ierr)
                                                        call mkdir_f('./DIAG/IVy_x'//char(0),ierr)
						endif
						if (nphxy .ne. 0) then
							call mkdir_f('./DIAG/Vx_y'//char(0),ierr)
						endif
						if (nphyy .ne. 0) then
							call mkdir_f('./DIAG/Vy_y'//char(0),ierr)
						endif
						if (nt_write_Py_vs_y .ne. 0) then
							call mkdir_f('./DIAG/Py_sumover_x'//char(0),ierr)
						endif
						if (nt_write_div_P_vs_y .ne. 0) then
							call mkdir_f('./DIAG/div_P_sumover_x'//char(0),ierr)
						endif
						if (ntp .ne. 0) then
							call mkdir_f('./DIAG/pot'//char(0),ierr)
						endif
						if (nt_through_wave .ne. 0) then
							call mkdir_f('./DIAG/ThroughWave/Vx_Vx'//char(0),ierr)
							call mkdir_f('./DIAG/ThroughWave/Vx_Vy'//char(0),ierr)
							call mkdir_f('./DIAG/ThroughWave/Vy_Vx'//char(0),ierr)
							call mkdir_f('./DIAG/ThroughWave/Vy_Vy'//char(0),ierr)
							call mkdir_f('./DIAG/ThroughWave/dEne'//char(0),ierr)
						endif
						
						do i = 1, num_phsl_x
							write(temp,'(i4)') phsl_x_pos(i)
							temp = adjustl(temp)
							call mkdir_f('./DIAG/vx_x/'//trim(temp)//char(0),ierr)
							call mkdir_f('./DIAG/vy_x/'//trim(temp)//char(0),ierr)
						enddo
						do i = 1, num_phsl_y
							write(temp,'(i4)') phsl_y_pos(i)
							temp = adjustl(temp)
							call mkdir_f('./DIAG/vx_y/'//trim(temp)//char(0),ierr)
							call mkdir_f('./DIAG/vy_y/'//trim(temp)//char(0),ierr)
						enddo
				endif
      		
      	end subroutine mkdir_structure
      	
      	subroutine write_jf_diag_file(id0,cdrun)
      		implicit none
      		integer::id0
      		character(len=*) :: cdrun
      		character(len=32)::filename
      		integer :: funitnum
      		!integer::modesx=11,modesy=11

      		if (id0 == 0) then
      			filename='diagparams_jf'
      			funitnum = get_funit(20)
      			open(unit=funitnum,file=trim(filename),form='unformatted',status=&
      &                 'replace')
!234567
					write (funitnum) idrun,indx,indy,ntfield,psolve,tend,dt,npx,npy,&
						&fvxmax,fvymax,nphbx,nphby,nphxx,nphyx,ntlines,nlines,linepos,phsl_x_pos,&
						&phsl_x_thick,ntphsl_x,num_phsl_x,ntden,ntpene,ntj,ntvdotE, &
                                                &fvxmax_ion,fvymax_ion

				endif
			end subroutine write_jf_diag_file
			
			subroutine get_1D_sumover_x(div_ESPoynt,div_P_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:) :: div_ESPoynt
				real,dimension(:) :: div_P_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					div_P_sumover_x(cury-inorder+1) = sum(div_ESPoynt(inorder:inorder+nx-1,j,1))
				enddo
				call plsum(div_P_sumover_x)
			end subroutine get_1D_sumover_x

			subroutine get_2D_sumover_x(jE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: jE
				real,dimension(:) :: jE_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo

      	if (size(jE,1) .eq. 2) then
					do j = inorder, nyproc+inorder-1
						cury = proc_pos(idproc+1) + j
						jE_sumover_x(cury-inorder+1) = sum(jE(1,inorder:inorder+nx-1,j,1)) + &
																				&	 sum(jE(2,inorder:inorder+nx-1,j,1))
					enddo
				endif
      	if (size(jE,1) .eq. 3) then
					do j = inorder, nyproc+inorder-1
						cury = proc_pos(idproc+1) + j
						jE_sumover_x(cury-inorder+1) = sum(jE(1,inorder:inorder+nx-1,j,1)) + &
																				&	 sum(jE(2,inorder:inorder+nx-1,j,1)) + &
																				&	 sum(jE(3,inorder:inorder+nx-1,j,1))
					enddo
				endif
				
				call plsum(jE_sumover_x)
			end subroutine get_2D_sumover_x

			!same as above function but just sums of the x-component, not x and y
			subroutine get_2D_sumover_x_xcomp(jE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: jE
				real,dimension(:) :: jE_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					jE_sumover_x(cury-inorder+1) = sum(jE(1,inorder:inorder+nx,j,1))
				enddo
				call plsum(jE_sumover_x)
			end subroutine get_2D_sumover_x_xcomp
			
			subroutine get_Py_sumover_x(ESPoynt,Py_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: ESPoynt
				real,dimension(:) :: Py_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					Py_sumover_x(cury-inorder+1) = sum(ESPoynt(2,inorder:inorder+nx,j,1))
				enddo
				call plsum(Py_sumover_x)
			end subroutine get_Py_sumover_x
			
			subroutine calc_ene_three_regions(grad_phi_int,ESPoynt_int,jE_sum,div_ESPoynt_sum,totU,Pflow,&
				&totjE,div_tot,num_par_cell,dt,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: grad_phi_int,ESPoynt_int,jE_sum
				real,dimension(:,:,:) :: div_ESPoynt_sum
				real,dimension(:,:) :: Pflow
				real,dimension(:) :: totjE, totU, div_tot
				real :: num_par_cell,dt
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder
				
				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
!				print*,nyproc,proc_pos(idproc+1),idproc,inorder
				do j = inorder, nyproc+inorder-1
!					print*,idproc, proc_pos(idproc+1)+j
					cury = proc_pos(idproc+1) + j
					!Do bottom box
					if (cury .ge. int_pos(1) .AND. cury .lt. int_pos(2)) then
						do i = inorder, inorder + nx
							totU(1) = totU(1) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
							totjE(1) = totjE(1) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
							div_tot(1) = div_tot(1) + div_ESPoynt_sum(i,j,1)
						enddo
						if (cury .eq. int_pos(1)) then
							do i = inorder, inorder + nx
								Pflow(1,1) = Pflow(1,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(1,1) = Pflow(1,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
						if (cury .eq. int_pos(2)-1) then
							do i = inorder, inorder + nx
								Pflow(1,3) = Pflow(1,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(1,3) = Pflow(1,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
					endif
					!Do middle box
					if (cury .ge. int_pos(2) .AND. cury .lt. int_pos(3)) then
						do i = inorder, inorder + nx
							totU(2) = totU(2) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
							totjE(2) = totjE(2) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
							div_tot(2) = div_tot(2) + div_ESPoynt_sum(i,j,1)
						enddo
						if (cury .eq. int_pos(2)) then
							do i = inorder, inorder + nx
								Pflow(2,1) = Pflow(2,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(2,1) = Pflow(2,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
						if (cury .eq. int_pos(3)-1) then
							do i = inorder, inorder + nx
								Pflow(2,3) = Pflow(2,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(2,3) = Pflow(2,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
					endif
					!Do top box
					if (cury .ge. int_pos(3) .AND. cury .lt. int_pos(4)) then
						do i = inorder, inorder + nx
							totU(3) = totU(3) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
							totjE(3) = totjE(3) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
							div_tot(3) = div_tot(3) + div_ESPoynt_sum(i,j,1)
						enddo
						if (cury .eq. int_pos(3)) then
							do i = inorder, inorder + nx
								Pflow(3,1) = Pflow(3,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(3,1) = Pflow(3,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
						if (cury .eq. int_pos(4)-1) then
							do i = inorder, inorder + nx
								Pflow(3,3) = Pflow(3,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!								Pflow(3,3) = Pflow(3,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
							enddo
						endif
					endif
				enddo
				
				call plsum(totU)
				call plsum(totjE)
				call plsum(Pflow)
				call plsum(div_tot)
				totU = totU*num_par_cell
				totjE = totje
				Pflow = Pflow
				div_tot = div_tot
			end subroutine calc_ene_three_regions
			
			!This returns the field energy in the rectagle -yrise_fall < y - ny/2 < yrise_fall and all x
			subroutine get_ene_in_driver_region(efield,field_ene,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				integer :: nx,nxe,ny,nypmx,nvp,idproc,inorder
				real :: field_ene
				real,dimension(:,:,:,:) :: efield
				
				integer :: i, j, nyproc,nylow,nyhigh
      	integer, dimension(nvp) :: proc_pos
      	real,dimension(1) :: temp
				
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
      	field_ene = 0.
      	
!      	print*,idproc,"a"
      	do j = inorder, nyproc+inorder
      		if ( (proc_pos(idproc+1) + j-1 >= ny/2-yrise_fall) .AND. &
      			&  (proc_pos(idproc+1) + j-1 <= ny/2+yrise_fall) ) then
!						if (idproc ==0) then
!						print*,j
!						endif
						do i = inorder, inorder+nx
							field_ene = field_ene + 0.5*(efield(1,i,j,1)**2 + efield(2,i,j,1)**2)
						enddo
					endif      		
      	enddo
!      	print*,idproc,"b"
      	
      	temp(1) = field_ene
      	call plsum(temp)
      	field_ene = temp(1)

			end subroutine get_ene_in_driver_region
			
			subroutine store_init_vel(part,npp)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				
				integer :: i
				
				do i = 1, npp(1)
					part(kivx_loc,i,1) = part(3,i,1)
					part(kivy_loc,i,1) = part(4,i,1)
				enddo
			end subroutine store_init_vel
			
			!This function puts all particles in a grid of their initial energy
			! by the final energy and writes it to a file
			subroutine calc_ene_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nEneBin_init,nEneBin_final) :: ene_bins
				real,dimension(nEneBin_init,nEneBin_final,1) :: ene_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_ene,f_ene
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = initEneRange / real(nEneBin_init)
				dfin = finalEneRange / real(nEneBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				ene_bins = 0.
				
				
				do i = 1, npp(1)
					i_ene = 0.5*(part(kivx_loc,i,1)*part(kivx_loc,i,1) + part(kivy_loc,i,1)*part(kivy_loc,i,1))
					f_ene = 0.5*(part(3,i,1)*part(3,i,1) + part(4,i,1)*part(4,i,1))
					
					i_arr = floor( i_ene * inv_dini ) + 1
					f_arr = floor( f_ene * inv_dfin ) + 1
					
					if (i_arr > nEneBin_init) then 
						i_arr = nEneBin_init
					endif
					if (f_arr > nEneBin_final) then 
						f_arr = nEneBin_final
					endif
					
					ene_bins(i_arr,f_arr) = ene_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(ene_bins)
				
!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				
				if (idproc == 0) then
					write (nlabel,'(i5)') it + 10000
					name = './DIAG/EneBin/eneBin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					ene_bins3d(:,:,1) = ene_bins
					call write_jf(ene_bins3d,nEneBin_init,nEneBin_final,name,FIN_VS_INIT_ENE,&
						&time,it,dini,dfin,0.,0.)
				endif

			end subroutine calc_ene_dist

			!This function puts all particles in a grid of their initial vx
			! by the final vx and writes it to a file
			subroutine calc_vx_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nVxBin_init,nVxBin_final) :: vx_bins
				real,dimension(nVxBin_init,nVxBin_final,1) :: vx_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_vx,f_vx
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = (initVxHigh-initVxLow) / real(nVxBin_init)
				dfin = (finalVxHigh-finalVxLow) / real(nVxBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				vx_bins = 0.
				
				
				do i = 1, npp(1)
					i_vx = part(kivx_loc,i,1)
					f_vx = part(3,i,1)
					
					i_arr = floor( (i_vx - initVxLow) * inv_dini ) + 1
					f_arr = floor( (f_vx - finalVxLow) * inv_dfin ) + 1
					
					if (i_arr > nVxBin_init) then 
						i_arr = nVxBin_init
					endif
					if (i_arr < 1) then 
						i_arr = 1
					endif
					if (f_arr > nVxBin_final) then 
						f_arr = nVxBin_final
					endif
					if (f_arr < 1) then 
						f_arr = 1
					endif
					
					vx_bins(i_arr,f_arr) = vx_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(vx_bins)
				
!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				if (idproc == 0) then
					write (nlabel,'(i5)') it + 10000
					name = './DIAG/VxBin/vxbin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					vx_bins3d(:,:,1) = vx_bins
					call write_jf(vx_bins3d,nVxBin_init,nVxBin_final,name,FIN_VS_INIT_VX,time,it,&
						&dini,dfin,initVxLow,finalVxLow)
				endif

			end subroutine calc_vx_dist

			!This function puts all particles in a grid of their initial vy
			! by the final vy and writes it to a file
			subroutine calc_vy_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nVyBin_init,nVyBin_final) :: vy_bins
				real,dimension(nVyBin_init,nVyBin_final,1) :: vy_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_vy,f_vy
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = (initVyHigh-initVyLow) / real(nVyBin_init)
				dfin = (finalVyHigh-finalVyLow) / real(nVyBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				vy_bins = 0.
								
				do i = 1, npp(1)
					i_vy = part(kivy_loc,i,1)
					f_vy = part(4,i,1)
					
					i_arr = floor( (i_vy - initVyLow) * inv_dini ) + 1
					f_arr = floor( (f_vy - finalVyLow) * inv_dfin ) + 1
					
					if (i_arr > nVyBin_init) then 
						i_arr = nVyBin_init
					endif
					if (i_arr < 1) then 
						i_arr = 1
					endif
					if (f_arr > nVyBin_final) then 
						f_arr = nVyBin_final
					endif
					if (f_arr < 1) then 
						f_arr = 1
					endif
					
					vy_bins(i_arr,f_arr) = vy_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(vy_bins)

!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				if (idproc == 0) then
					write (nlabel,'(i5)') it + 10000
					name = './DIAG/VyBin/vybin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					vy_bins3d(:,:,1) = vy_bins
					call write_jf(vy_bins3d,nVyBin_init,nVyBin_final,&
						&name,FIN_VS_INIT_VY,time,it,dini,dfin,initVyLow,finalVyLow)
				endif
			end subroutine calc_vy_dist

			! This function sets part(kivx_loc) = 10001. for all particles not in init_range
			! all other particles get their init coords
			subroutine init_through_wave(part,npp)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				real :: x,y,vx,vy
				
				integer :: i

				do i = 1, npp(1)
					x = part(1,i,1)
					y = part(2,i,1)
					vx = part(3,i,1)
					vy = part(4,i,1)
					
					if ( (x > init_range(1,1)) .AND. (x < init_range(1,2)) .AND. &
						&  (y > init_range(2,1)) .AND. (y < init_range(2,2)) .AND. &
						&  (vx > init_range(3,1)) .AND. (vx < init_range(3,2)) .AND. &
						&  (vy > init_range(4,1)) .AND. (vy < init_range(4,2)) ) then
						part(kivx_loc,i,1) = part(3,i,1)
						part(kivy_loc,i,1) = part(4,i,1)
					else
						part(kivx_loc,i,1) = 10001.
					endif
				enddo
			end subroutine init_through_wave


			! This function puts all particles in grids of their initial velocities
			! by the final velocities only for those particles that have crossed final_y
			! This function must be called at then end of simulation with write_files = 1
			! all other calls should have write_files = 0
			!****************************************************************************
			! *********ONLY WORKS WHEN INIT_RANGE(2,1) AND INIT_RANGE(2,2) ARE BOTH > 0
			! *********CAN ONLY DO PARTICLES MOVING UP INITIALLY, NOT DOWN
			!****************************************************************************
			subroutine calc_through_wave_dist(part,nx,npp,time,it,idproc,tot_dEne,tot_part_count,tot_mag_dEne)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc,nx
				real :: time
				real,dimension(1) :: tot_dEne, tot_mag_dEne
				integer,dimension(1) :: tot_part_count

				real,dimension(:,:,:),allocatable,save :: vx_vy_bins,vy_vy_bins,vx_vx_bins,vy_vx_bins,temp
				real,dimension(:,:,:),allocatable,save :: dEne
				real,dimension(:,:,:),allocatable,save :: tempxx,tempxy,tempyx,tempyy,tempdEne
				real :: dinix,dfinx,inv_dinix,inv_dfinx,diniEne,dfinEne
				real :: diniy,dfiny,inv_diniy,inv_dfiny,inv_diniEne,inv_dfinEne
				real :: i_vx,f_vx, i_vy,f_vy,i_ene,f_ene
				real :: up, vph, temp_tot_dEne, temp_tot_mag_dEne
				integer :: i, ix_arr,fx_arr, iy_arr,fy_arr, iene_arr,fene_arr
				integer :: funitnum
				integer :: particle_count
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				if (.not. allocated(vx_vy_bins)) then
					allocate(vx_vy_bins(nVxBin_init,nVyBin_final,1))
					allocate(vy_vy_bins(nVyBin_init,nVyBin_final,1))
					allocate(vx_vx_bins(nVxBin_init,nVxBin_final,1))
					allocate(vy_vx_bins(nVyBin_init,nVxBin_final,1))
					allocate(dEne(nEneBin_init,nEneBin_final,1))
					allocate(tempxy(nVxBin_init,nVyBin_final,1))
					allocate(tempyy(nVyBin_init,nVyBin_final,1))
					allocate(tempxx(nVxBin_init,nVxBin_final,1))
					allocate(tempyx(nVyBin_init,nVxBin_final,1))
					allocate(tempdEne(nEneBin_init,nEneBin_final,1))
					vx_vy_bins = 0.
					vx_vx_bins = 0.
					vy_vy_bins = 0.
					vy_vx_bins = 0.
					dEne = 0.
				endif

				tempxx = 0.
				tempxy = 0.
				tempyx = 0.
				tempyy = 0.

				dinix = (initVxHigh-initVxLow) / real(nVxBin_init)
				dfinx = (finalVxHigh-finalVxLow) / real(nVxBin_final)
				inv_dinix = 1./ dinix
				inv_dfinx = 1./ dfinx

				diniy = (initVyHigh-initVyLow) / real(nVyBin_init)
				dfiny = (finalVyHigh-finalVyLow) / real(nVyBin_final)
				inv_diniy = 1./ diniy
				inv_dfiny = 1./ dfiny
				
				diniEne = real(initEneRange) / real(nEneBin_init)
				dfinEne = real(finalEneRange) / real(nEneBin_final)
				inv_diniEne = 1. / diniEne
				inv_dfinEne = 1. / dfinEne
				
				vph = real(nx) / real(wavemode)
				vph = 6.283185307/vph !vph = k after this line
				vph = wavew / vph

				particle_count = 0
				temp_tot_dEne = 0.
				temp_tot_mag_dEne = 0.
				
				do i = 1, npp(1)
					! Do x component
					i_vx = part(kivx_loc,i,1)
					if (i_vx < 10000.) then
						if (part(2,i,1) > final_y) then
							particle_count = particle_count + 1
							part(kivx_loc,i,1) = 10001.
	
							f_vx = part(3,i,1)
							
							ix_arr = floor( (i_vx - initVxLow) * inv_dinix ) + 1
							fx_arr = floor( (f_vx - finalVxLow) * inv_dfinx ) + 1
							
							if (ix_arr > nVxBin_init) then 
								ix_arr = nVxBin_init
							endif
							if (ix_arr < 1) then 
								ix_arr = 1
							endif
							if (fx_arr > nVxBin_final) then 
								fx_arr = nVxBin_final
							endif
							if (fx_arr < 1) then 
								fx_arr = 1
							endif
							
							! Do y component
							i_vy = part(kivy_loc,i,1)
							f_vy = part(4,i,1)
							
							iy_arr = floor( (i_vy - initVyLow) * inv_diniy ) + 1
							fy_arr = floor( (f_vy - finalVyLow) * inv_dfiny ) + 1
							
							if (iy_arr > nVyBin_init) then 
								iy_arr = nVyBin_init
							endif
							if (iy_arr < 1) then 
								iy_arr = 1
							endif
							if (fy_arr > nVyBin_final) then 
								fy_arr = nVyBin_final
							endif
							if (fy_arr < 1) then 
								fy_arr = 1
							endif
							i_ene = 0.5 * ( (i_vx-vph)**2 + i_vy**2)
							f_ene = 0.5 * ( (f_vx-vph)**2 + f_vy**2)
							
							iene_arr = floor( (i_ene - 0.) * inv_diniEne ) + 1  !minimum energy is 0
							fene_arr = floor( (f_ene - 0.) * inv_dfinEne ) + 1
							
							if (iene_arr > nEneBin_init) then 
								iene_arr = nEneBin_init
							endif
							if (iene_arr < 1) then 
								iene_arr = 1
							endif
							if (fene_arr > nEneBin_final) then 
								fene_arr = nEneBin_final
							endif
							if (fene_arr < 1) then 
								fene_arr = 1
							endif
		
							vx_vx_bins(ix_arr,fx_arr,1) = vx_vx_bins(ix_arr,fx_arr,1) + 1.
							vx_vy_bins(ix_arr,fy_arr,1) = vx_vy_bins(ix_arr,fy_arr,1) + 1.
							vy_vx_bins(iy_arr,fx_arr,1) = vy_vx_bins(iy_arr,fx_arr,1) + 1.
							vy_vy_bins(iy_arr,fy_arr,1) = vy_vy_bins(iy_arr,fy_arr,1) + 1.
							dEne(iene_arr,fene_arr,1) = dEne(iene_arr,fene_arr,1) + 1!f_ene - i_ene
							
							temp_tot_dEne = temp_tot_dEne + f_ene - i_ene
							temp_tot_mag_dEne = temp_tot_mag_dEne + abs(f_ene - i_ene)
						endif
					endif
				enddo
				
				if (particle_count .ne. 0) then
					tot_dEne(1) = tot_dEne(1) + temp_tot_dEne
					tot_mag_dEne(1) = tot_mag_dEne(1) + temp_tot_mag_dEne
					tot_part_count(1) = tot_part_count(1) + particle_count
				endif

				tempxx = vx_vx_bins
				tempxy = vx_vy_bins
				tempyx = vy_vx_bins
				tempyy = vy_vy_bins
				tempdEne = dEne
				call plsum(tempxx)
				call plsum(tempxy)
				call plsum(tempyx)
				call plsum(tempyy)
				call plsum(tempdEne)
		
!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				if (idproc == 0) then
					
					write (nlabel,'(i5)') it + 10000
					name = './DIAG/ThroughWave/Vx_Vx/vx_vx_bin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					call write_jf(tempxx,nVxBin_init,nVxBin_final,name,THROUGH_VX_VX,time,it,&
						&dinix,dfinx,initVxLow,finalVxLow)

					name = './DIAG/ThroughWave/Vx_Vy/vx_vy_bin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					call write_jf(tempxy,nVxBin_init,nVyBin_final,name,THROUGH_VX_VY,time,it,&
						&dinix,dfiny,initVxLow,finalVyLow)

					name = './DIAG/ThroughWave/Vy_Vx/vy_vx_bin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					call write_jf(tempyx,nVyBin_init,nVxBin_final,name,THROUGH_VY_VX,time,it,&
						&diniy,dfinx,initVyLow,finalVxLow)

					name = './DIAG/ThroughWave/Vy_Vy/vy_vy_bin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					call write_jf(tempyy,nVyBin_init,nVyBin_final,name,THROUGH_VY_VY,time,it,&
						&diniy,dfiny,initVyLow,finalVyLow)

					name = './DIAG/ThroughWave/dEne/dEne_bin'
					fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
					call write_jf(tempdEne,nEneBin_init,nEneBin_final,name,THROUGH_DENE,time,it,&
						&diniEne,dfinEne,0.,0.)
				endif

			end subroutine calc_through_wave_dist


			! This will calculate the current from the particle data 
			subroutine j_diag_deposit(part,arr,npp,xory,ny,nvp,idproc)
				implicit none
				real,dimension(:,:,:), pointer :: part, arr
				integer,dimension(:) :: npp
				integer :: xory,ny,nvp,idproc
				
				integer :: i,xpos,ypos,dir,procpos
				
				procpos = ny /nvp * idproc
				procpos = procpos - 2 !to get correct array pos with guard cells

				dir = 2+xory
				do i = 1, npp(1)
					xpos = floor(part(1,i,1)) + 2
					ypos = floor(part(2,i,1)) - procpos
					arr(xpos,ypos,1) = arr(xpos,ypos,1) + part(dir,i,1)
				enddo
			end subroutine j_diag_deposit
			
			!This will calculate v_i E_i, where i is specified by xory
			subroutine vdotE(part,fxye,arr,npp,xory,ny,nvp,idproc)
				implicit none
				real,dimension(:,:,:), pointer :: part,arr
				real,dimension(:,:,:,:) :: fxye
				integer,dimension(:) :: npp
				integer :: xory,ny,nvp,idproc
				
				integer :: i,xpos,ypos,dir,procpos

				procpos = ny /nvp * idproc
				procpos = procpos - 2 !to get correct array pos with guard cells

				dir = 2+xory
				do i = 1, npp(1)
					xpos = floor(part(1,i,1)) + 2
					ypos = floor(part(2,i,1)) - procpos
					arr(xpos,ypos,1) = arr(xpos,ypos,1) + part(dir,i,1) * fxye(xory,xpos,ypos,1)
				enddo
			end subroutine vdotE

! PHASE SPACE  -- > electrons
			! This calculates the phase space at a given time for 1 direction
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_x(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname
				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
				
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
						
							pos = floor(part(1,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbx) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(1,i,1))+1
								v = part(3,i,1)
								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphbx) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then					
						fname = './DIAG/Vx_x/vx_x'
						call write_jf(fvxy,nxy,nphbx,trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(1,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(1,i,1))+1
								v = part(4,i,1)
								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphby) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/Vy_x/vy_x'
						call write_jf(fvxy,nxy,nphby,trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
					endif
				endif

			end subroutine phase_space_vxy_vs_x

			! This calculates the phase space at a given time for 1 direction
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_y(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname

				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
								
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(2,i,1))+1
								v = part(3,i,1)
								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
								
								if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/Vx_y/vx_y'
						call write_jf(fvxy,nxy,nphbx,trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(2,i,1))+1
								v = part(4,i,1)
								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphby) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/Vy_y/vy_y'
						call write_jf(fvxy,nxy,nphby,trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
					endif
				endif

			end subroutine phase_space_vxy_vs_y
			
! PHASE SPACE --> ions
			! This calculates the phase space at a given time for 1 direction
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_x_ion(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname
				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
				
				if (xory == 1) then
					vrange = 2.*fvxmax_ion
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -fvxmax_ion
					high = fvxmax_ion
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
						
							pos = floor(part(1,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v - low)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbx) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(1,i,1))+1
								v = part(3,i,1)
								varrpos = floor( (v - low)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphbx) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then					
						fname = './DIAG/IVx_x/vx_x'
						call write_jf(fvxy,nxy,nphbx,trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
					endif
				else
					vrange = 2.*fvymax_ion
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax_ion
					high = fvymax_ion
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(1,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v - low)*dvinv + 0.5) + 1
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(1,i,1))+1
								v = part(4,i,1)
								varrpos = floor( (v - low)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphby) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/IVy_x/vy_x'
						call write_jf(fvxy,nxy,nphby,trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
					endif
				endif

			end subroutine phase_space_vxy_vs_x_ion

			! This calculates the phase space at a given time for 1 direction
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_y_ion(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname

				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
								
				if (xory == 1) then
					vrange = 2.*fvxmax_ion
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(2,i,1))+1
								v = part(3,i,1)
								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
								
								if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/IVx_y/vx_y'
						call write_jf(fvxy,nxy,nphbx,trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
					endif
				else
					vrange = 2.*fvymax_ion
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax_ion
					high = fvymax_ion
	
					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
							endif
						enddo
					else
						do i = 1, npp(1)
							real_tag = part(addtag_loc,i,1)
							if (int_tag(1) == -1) then
								pos = floor(part(2,i,1))+1
								v = part(4,i,1)
								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
								
								if (varrpos .le. nphby) then 
									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
								endif
							endif
						enddo
					endif
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						fname = './DIAG/IVy_y/vy_y'
						call write_jf(fvxy,nxy,nphby,trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
					endif
				endif

			end subroutine phase_space_vxy_vs_y_ion
			
			subroutine bin_Ex_and_Ey(fxye,Ex_bin,Ey_bin,nx,ny,kyp,it,time,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: fxye
				real,dimension(:,:,:) :: Ex_bin, Ey_bin
				integer :: nx, ny, it, nvp, idproc, inorder,kyp
				real :: time
				
				integer :: i, j, nyproc, arr, ypos,xrange, yrange
				real :: Exrange, Eyrange, dEx, dEy,invdEx,invdEy
				integer, dimension(nvp) :: proc_pos
				character(len=32) :: fname
				
				xrange = bin_E_xrange(2) - bin_E_xrange(1)
				yrange = bin_E_yrange(2) - bin_E_yrange(1)
				Exrange = bin_E_Exrange(2) - bin_E_Exrange(1)
				Eyrange = bin_E_Eyrange(2) - bin_E_Eyrange(1)
				dEx = Exrange / real(bin_E)
				dEy = Exrange / real(bin_E)
				invdEx = 1./dEx
				invdEy = 1./dEy
				
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo
				do i = 0, kyp-1
					do j = 0, nx-1
	
						ypos = proc_pos(idproc+1) + i
						
						if (ypos >= bin_E_yrange(1) .AND. ypos < bin_E_yrange(2)) then
							if (j >= bin_E_xrange(1) .AND. j < bin_E_xrange(2)) then
								arr = floor( (fxye(1,j+inorder+1,i+inorder+1,1) - bin_E_Exrange(1))*invdEx + 0.5) + 1
								if (arr > 0 .AND. arr < bin_E) then
									Ex_bin(j-bin_E_xrange(1)+1,ypos-bin_E_yrange(1)+1,arr) = 1.
								endif

								arr = floor( (fxye(2,j+inorder+1,i+inorder+1,1) - bin_E_Eyrange(1))*invdEy + 0.5) + 1
								if (arr > 0 .AND. arr < bin_E) then
									Ey_bin(j-bin_E_xrange(1)+1,ypos-bin_E_yrange(1)+1,arr) = 1.
								endif
							endif
						endif
					enddo
				enddo
				
				call plsum(Ex_bin)
				call plsum(Ey_bin)

				if (idproc == 0) then
					fname = './DIAG/bin_Ex/bin_Ex'
					!Use LINEAR because no gaurd cells in fvy
					call write_jf(Ex_bin,xrange,yrange,bin_E,fname,BIN_EX_LABEL,&
						&time,it,1.,1.,dEx,real(bin_E_xrange(1)),real(bin_E_yrange(1)),bin_E_Exrange(1))

					fname = './DIAG/bin_Ey/bin_Ey'
					!Use LINEAR because no gaurd cells in fvy
					call write_jf(Ey_bin,xrange,yrange,bin_E,&
						&fname,BIN_EY_LABEL,time,it,1.,1.,dEy,real(bin_E_xrange(1)),&
						&real(bin_E_yrange(1)),bin_E_Eyrange(1))
				endif

			end subroutine bin_Ex_and_Ey

			! This calculates the phase space for vx or vy vs. x and y
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_x_and_y(part,fvxy,xory,nx,ny,np,npp,it,time,idproc)
				implicit none
				integer :: np,xory,nx,ny,idproc,it
				real :: time
      	real,dimension(:,:,:) :: part
				real,dimension(:,:,:),pointer :: fvxy
				integer,dimension(:),pointer :: npp
				integer :: i,varrpos,xarrpos
!				real :: dv,low,high,v,vrange,dvinv
				real :: v
				character(len=32) :: fname
				real,dimension(2) :: fmax,vrange,low,dv,dvinv,r
				real :: dx,dy,invdx,invdy
				integer,dimension(2) :: pos
				
				fmax(1) = fvx_xy_max
				fmax(2) = fvy_xy_max
				vrange(1) = 2.*fvx_xy_max
				vrange(2) = 2.*fvy_xy_max
				low(1) = -1.*fmax(1)
				low(2) = -1.*fmax(2)
				dv(1) = vrange(1) / real(nphb_xy)
				dv(2) = vrange(2) / real(nphb_xy)
				dvinv = 1./dv
				
				dx = (fvxy_xy_xrange(2) - fvxy_xy_xrange(1)) / real(nfvxy_xy_xrange)
				dy = (fvxy_xy_yrange(2) - fvxy_xy_yrange(1)) / real(nfvxy_xy_yrange)
				invdx = 1./dx
				invdy = 1./dy
				
				do i = 1, npp(1)
					r = part(1:2,i,1)
					if (r(1) > fvxy_xy_xrange(1) .AND. r(1) < fvxy_xy_xrange(2)) then
						if (r(2) > fvxy_xy_yrange(1) .AND. r(2) < fvxy_xy_yrange(2)) then
							
							pos(1) = floor( (r(1) - fvxy_xy_xrange(1))*invdx + 0.5 ) + 1
							pos(2) = floor( (r(2) - fvxy_xy_yrange(1))*invdy + 0.5 ) + 1
							
							if (pos(1) > 0 .AND. pos(1) < nfvxy_xy_xrange+1) then
								if (pos(2) > 0 .AND. pos(2) < nfvxy_xy_yrange+1) then
									
									v = part(2+xory,i,1)
									varrpos = floor( (v + fmax(xory))*dvinv(xory) + 0.5) + 1
									if (varrpos .le. nphb_xy) then 
										fvxy(pos(1),pos(2),varrpos) = fvxy(pos(1),pos(2),varrpos) + 1.
									endif
								endif
							endif
						endif
					endif
				enddo
				
				call plsum(fvxy(:,:,:))				
				
				if (xory == 1) then
					if (idproc == 0) then
						fname = './DIAG/fvx_xy/fvx_xy'
						!Use LINEAR because no gaurd cells in fvy
						call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,&
							&fname,FVX_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
							&fvxy_xy_yrange(1),low(xory))
					endif
				else
					if (idproc == 0) then
						fname = './DIAG/fvy_xy/fvy_xy'
						!Use LINEAR because no gaurd cells in fvxy
						call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,&
							&fname,FVY_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
							&fvxy_xy_yrange(1),low(xory))
					endif
				endif
			end subroutine phase_space_vxy_vs_x_and_y

			!make vx vs. vy phase space plot
			subroutine phase_vx_vy(part,fvxvy,np,npp,time,it,idimp,nblok,idproc)
				implicit none
				integer :: np, idimp, nblok,idproc,it
      	real,dimension(:,:,:) :: part
				real,dimension(:,:,:),pointer :: fvxvy
				integer,dimension(nblok) :: npp
				real :: time
				integer :: i,vxarrpos,vyarrpos
				real :: dvx,lowx,highx,vx,vxrange,dvxinv
				real :: dvy,lowy,highy,vy,vyrange,dvyinv
				character(len=32) :: fname
				
				vxrange = 2.*fvxmax
				dvx = vxrange / real(nphbx)
				dvxinv = 1. / dvx
				lowx = -1.*fvxmax
				highx = fvxmax
				vyrange = 2.*fvymax
				dvy = vyrange / real(nphby)
				dvyinv = 1. / dvy
				lowy = -1.*fvymax
				highy = fvymax

				do i = 1, npp(1)
					vx = part(3,i,1)
					vxarrpos = floor( (vx + fvxmax)*dvxinv + 0.5) + 1
					vy = part(4,i,1)
					vyarrpos = floor( (vy + fvymax)*dvyinv + 0.5) + 1
					
					if ((vxarrpos .le. nphbx) .AND. (vyarrpos .le. nphby) .AND. &
						&	(vxarrpos .gt. 0) .AND. (vyarrpos .gt. 0)) then 
						fvxvy(vxarrpos,vyarrpos,1) = fvxvy(vxarrpos,vyarrpos,1) + 1.
					endif
				enddo
				call plsum(fvxvy(:,:,1))
				fname = './DIAG/VxVy/vx_vy'
				call write_jf(fvxvy,nphby,nphbx,trim(fname),VXVY,time,it,dvx,dvy,lowx,lowy)					
				
			end subroutine phase_vx_vy
			
			!This function makes, at each time that it is to be dumped, a number of 
			!vx vs. vy phase space plots that are taken over a range of x and y, with
			!the y moving at y+vx_vy_speed*n, where n ranges from 0 to nvx_vy_speed-1.
			!The purpose is to see how the dist fxn varies over a range of space at each time
			subroutine phase_vx_vy_x_y_movie_at_each_dump(part,fvxvy,np,npp,time,it,idimp,nblok,idproc)
				implicit none
				integer :: np, idimp, nblok,idproc,it
      	real,dimension(:,:,:) :: part
				real,dimension(:,:,:),pointer :: fvxvy
				integer,dimension(nblok) :: npp
				real :: time
				integer :: i,vxarrpos,vyarrpos, j, ierr, rec
				real :: dvx,lowx,highx,vx,vxrange,dvxinv
				real :: dvy,lowy,highy,vy,vyrange,dvyinv
				real :: yshift
				character(len=100) :: fname,name,dname
	      character(len=5) :: nlabel
				
				vxrange = 2.*fvxmax
				dvx = vxrange / real(nphbx)
				dvxinv = 1. / dvx
				lowx = -1.*fvxmax
				highx = fvxmax
				vyrange = 2.*fvymax
				dvy = vyrange / real(nphby)
				dvyinv = 1. / dvy
				lowy = -1.*fvymax
				highy = fvymax
				
				write (nlabel,'(i5)') it + 10000
				dname = './DIAG/vx_vy_speed/vx_vy_'//trim(adjustl(nlabel))
				call mkdir_f(trim(adjustl(dname))//char(0),ierr)

				do j = 1, nvx_vy_speed
					fvxvy = 0.

					yshift = real(j-1)*vx_vy_speed
					do i = 1, npp(1)
						if ((part(2,i,1) >= vx_vy_yrange(1)+yshift) .AND. (part(2,i,1) < vx_vy_yrange(2)+yshift)) then
						if (part(1,i,1) >= vx_vy_xrange(1) .AND. part(1,i,1) < vx_vy_xrange(2)) then
							vx = part(3,i,1)
							vxarrpos = floor( (vx + fvxmax)*dvxinv + 0.5) + 1
							vy = part(4,i,1)
							vyarrpos = floor( (vy + fvymax)*dvyinv + 0.5) + 1
							
							if ((vxarrpos .le. nphbx) .AND. (vyarrpos .le. nphby) .AND. &
								&	(vxarrpos .gt. 0) .AND. (vyarrpos .gt. 0)) then 
								fvxvy(vxarrpos,vyarrpos,1) = fvxvy(vxarrpos,vyarrpos,1) + 1.
							endif
						endif
						endif
					enddo
					call plsum(fvxvy(:,:,1))
					fname = ''
					fname = trim(adjustl(dname))//'/vx_vy'
					call write_jf(fvxvy,nphbx,nphby,trim(adjustl(fname)),VXVY,vx_vy_yrange(1)+yshift,j,dvx,dvy,lowx,lowy)					

				enddo					
			end subroutine phase_vx_vy_x_y_movie_at_each_dump
					
			subroutine write_Eoft_aty(lines,fxye,nx,ny,nvp,idproc)
				implicit none
				integer :: nx,ny,nvp,idproc
				real,dimension(:,:,:,:) :: fxye
				real,dimension(:,:) :: lines
				integer :: nyproc,i
				integer, dimension(nvp) :: proc_pos
				real,dimension(:),allocatable,save :: sfieldtemp1
				character(len=10) :: num
				character(len=30) :: fname
				
				if (.not. allocated(sfieldtemp1)) then
					allocate(sfieldtemp1(nx))
				endif
				
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo
				lines = 0.
				do i = 1, nlines
!					linepos = ny / nlines * (i-1)
					if (linepos(i) < proc_pos(idproc+1) + nyproc .AND. &
						&	linepos(i) >= proc_pos(idproc+1) ) then
						lines(i,:) = fxye(1,2:nx+1,linepos(i)-proc_pos(idproc+1)+2,1)
					endif
				enddo
				call plsum(lines)
				if (idproc == 0) then
					do i = 1, nlines
!						linepos = ny / nlines * (i-1)
						if (nlinerec(i) == 0) then
							nlinerec(i) = -1
							write (num,'(i10)') linepos(i)
							num = adjustl(num)
							fname = 'line.'//num
							linesfunit(i) = get_funit(20)
							sfieldtemp1 = lines(i,:)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtemp1,nx,linesfunit(i),nlinerec(i),fname,LINEAR)
						else 
							write (num,'(i10)') linepos(i)
							num = adjustl(num)
							fname = 'line.'//num
							sfieldtemp1 = lines(i,:)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtemp1,nx,linesfunit(i),nlinerec(i),fname,LINEAR)
						endif
					enddo
				endif
			end subroutine	write_Eoft_aty

			subroutine phase_slices_y(part,fv,xory,nx,ny,np,npp,time,it,idimp,nblok,nvp,idproc)
				implicit none
				integer :: np, idimp, nblok,xory,nx,ny,nvp,idproc,it
				real :: time
	      	real,dimension(:,:,:) :: part
!	      	real,dimension(idimp,np,nblok) :: part
				real,dimension(:,:,:),pointer :: fv
				real,dimension(:,:,:),allocatable,save :: sfieldtempx,sfieldtempy
				integer,dimension(nblok) :: npp
				integer, dimension(nvp) :: proc_pos
				integer :: i,j,varrpos,xarrpos,pos,nyproc
				real :: dv,low,high,posx,v,vrange,dvinv
				character(len=100) :: fname,num
				real :: pposhigh,pposlow
				
				if (.not. allocated(sfieldtempx)) then
					allocate(sfieldtempx(ny,nphbx,1))
				endif
				if (.not. allocated(sfieldtempy)) then
					allocate(sfieldtempy(ny,nphby,1))
				endif
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo

				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					do i = 1, npp(1)
						pos = floor(part(2,i,1))+1
						v = part(3,i,1)
						do j = 1, num_phsl_y
							pposlow = real(phsl_y_pos(j))
							pposhigh = real(phsl_y_pos(j)+phsl_y_thick)
							if ((part(1,i,1) .ge. pposlow) .AND. &
								& (part(1,i,1) .lt. pposhigh) ) then
								
								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
								if (varrpos .le. nphbx) then 
									fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
								endif
							endif
						enddo
					enddo
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_y
						call plsum(fv(:,:,i))
					enddo
					if (idproc == 0) then
						do i = 1, num_phsl_y
							write (num,'(i10)') phsl_y_pos(i)
							num = adjustl(num)
							fname = './DIAG/vx_y/'//trim(num)//'/ph_xy_slice_'//trim(num)
							sfieldtempx(:,:,1) = fv(:,:,i)
							!Use LINEAR because no gaurd cells in fvxy
! HDF write
							call write_jf(sfieldtempx,ny,nphbx,trim(fname),PHSL_XY,time,it,1.,1.,0.,0.)
						enddo
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
					
					do i = 1, npp(1)
						pos = floor(part(2,i,1))+1
						v = part(4,i,1)
						do j = 1, num_phsl_y
							pposlow = real(phsl_y_pos(j))
							pposhigh = real(phsl_y_pos(j)+phsl_y_thick)
							if ((part(1,i,1) .ge. pposlow) .AND. &
								& (part(1,i,1) .lt. pposhigh) ) then
								
								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
								if (varrpos .le. nphby) then 
									fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
								endif
							endif
						enddo
					enddo
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_y
						call plsum(fv(:,:,i))
					enddo
					if (idproc == 0) then
						do i = 1, num_phsl_y
							write (num,'(i10)') phsl_y_pos(i)
							num = adjustl(num)
							fname = './DIAG/vy_y/'//trim(num)//'/ph_yy_slice_'//trim(num)
							sfieldtempy(:,:,1) = fv(:,:,i)
							call write_jf(sfieldtempy,ny,nphby,trim(fname),PHSL_YY,time,it,1.,1.,0.,0.)
						enddo
					endif
				endif
			end subroutine phase_slices_y
			
			subroutine phase_slices_x(part,fv,xory,nx,ny,np,npp,time,it,idimp,nblok,nvp,idproc)
				implicit none
				integer :: np, idimp, nblok,xory,nx,ny,nvp,idproc,it
				real :: time
	      	real,dimension(:,:,:) :: part
!	      	real,dimension(idimp,np,nblok) :: part
				real,dimension(:,:,:),pointer :: fv
				real,dimension(:,:,:),allocatable,save :: sfieldtempx,sfieldtempy
				integer,dimension(nblok) :: npp
				integer, dimension(nvp) :: proc_pos
				integer :: i,j,varrpos,xarrpos,pos,nyproc
				real :: dv,low,high,posx,v,vrange,dvinv
				character(len=100) :: fname,num
				real :: pposhigh,pposlow
				
				if (.not. allocated(sfieldtempx)) then
					allocate(sfieldtempx(nx,nphbx,1))
				endif
				if (.not. allocated(sfieldtempy)) then
					allocate(sfieldtempy(nx,nphby,1))
				endif
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
	
					do j = 1, num_phsl_x
						! If the jth slice is on this processor then do it
						if (phsl_x_pos(j) >= proc_pos(idproc+1) .AND. &
							&	phsl_x_pos(j) < proc_pos(idproc+1) + nyproc) then
							
							pposlow = real(phsl_x_pos(j))
							pposhigh = real(phsl_x_pos(j)+phsl_x_thick)
							do i = 1, npp(1)
								! If this particle falls in slice thickness
								if (part(2,i,1) < pposhigh .AND. part(2,i,1) >= pposlow) then
									pos = floor(part(1,i,1))+1
									v = part(3,i,1)
									varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
									if (varrpos .le. nphbx) then 
										fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
									endif
								endif
							enddo
						endif
					enddo
					
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_x
						call plsum(fv(:,:,i))
					enddo

					if (idproc == 0) then
						do i = 1, num_phsl_x
							write (num,'(i10)') phsl_x_pos(i)
							num = adjustl(num)
!								fname = './DIAG/vx_x/ph_xx_slice.'//num
							fname = './DIAG/vx_x/'//trim(num)//'/ph_xx_slice_'//trim(num)
							sfieldtempx(:,:,1) = fv(:,:,i)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtempx,nx,nphbx,trim(fname),PHSL_XX,time,it,1.,1.,0.,0.)
						enddo
					endif

				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
					do j = 1, num_phsl_x
						! If the jth slice is on this processor then do it
						if (phsl_x_pos(j) >= proc_pos(idproc+1) .AND. &
							&	phsl_x_pos(j) < proc_pos(idproc+1) + nyproc) then							
							pposlow = real(phsl_x_pos(j))
							pposhigh = real(phsl_x_pos(j)+phsl_x_thick)

							do i = 1, npp(1)
								! If this particle falls in slice thickness
								if (part(2,i,1) < pposhigh .AND. part(2,i,1) >= pposlow) then
									pos = floor(part(1,i,1))+1
									v = part(3,i,1)
									varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
									if (varrpos .le. nphby) then 
										fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
									endif
								endif
							enddo
						endif
					enddo

					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_x
						call plsum(fv(:,:,i))
					enddo
	
					if (idproc == 0) then
						do i = 1, num_phsl_x
							write (num,'(i10)') phsl_x_pos(i)
							num = adjustl(num)
							fname = './DIAG/vy_x/'//trim(num)//'/ph_yx_slice_'//trim(num)
							sfieldtempy(:,:,1) = fv(:,:,i)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtempy,nx,nphby,trim(fname),PHSL_YX,time,it,1.,1.,0.,0.)
						enddo
					endif
				endif

			end subroutine phase_slices_x

!-------------------------------------------------------------------------------

			subroutine div_finite_difference(f,df,nx,ny)
				real,dimension(:,:,:,:) :: f
				real,dimension(:,:,:) :: df
				integer :: nx, ny
				integer :: i, j

				df = 0.
				do i = 2, nx + 2
					df(i,:,1) = f(1,i+1,:,1)-f(1,i,:,1)
				enddo
				do i = 2, ny + 2
					df(:,i,1) = df(:,i,1) + f(2,:,i+1,1)-f(2,:,i,1)
				enddo
				
			end subroutine div_finite_difference
				
				
				
				
				





			end module diag_jf
