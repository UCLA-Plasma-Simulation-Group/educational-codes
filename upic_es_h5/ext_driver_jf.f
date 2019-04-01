!-----------------------------------------------------------------------
!
module ext_driver_jf
	use pinit2d_jf
	implicit none
	private
	public :: plane_wave,plane_wave_in_y,plane_finite_wavelen,gauss_tran_finite_wavelen
	public :: angle_wave,cross_wave,rect_tran_finite_wavelen
	public :: gauss_tran_per_wavelen,doub_gauss_tran_finite_wavelen
	public :: gauss_tran_per_wavelen_circle_wavefronts
	public :: gauss_tran_per_wavelen_pois,rect_tran_per_wavelen,supergauss_tran_per_wavelen
	public :: plane_finite_wavelen_supergauss
	
	contains
	
	subroutine cross_wave(fxye,t,nx,nxe,ny,nypmx,&
		&nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos, yhalf
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f, phase
		real :: tan_value
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo

		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
				
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			yhalf = real(ny/2)
			tan_value = tan(angle)
			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					if (ypos > ny/2) then
						ytemp = 1.
						! phase should be tan(angle) * (y - y_0)
						phase = -1. * tan_value * (ypos - yhalf)
					else
						ytemp = 1.
						phase = tan_value * (ypos - yhalf)
					endif
				else
					ytemp = 0.
				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
					fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp2*sin(wavek*xpos-wavew*t + phase)
				enddo
			enddo
		endif
	end subroutine cross_wave
	
	subroutine angle_wave(fxye,t,nx,nxe,ny,nypmx,&
		&nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f, phase,tan_value

		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
		!				tan_value = tan(1.5707963 - angle)
			tan_value = tan(angle)
			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					ytemp = 1.
					! phase should be tan(angle) * (y - y_0)
					phase = tan_value * (ypos - ylow)
				else
					ytemp = 0.
				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
					fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp2*sin(wavek*xpos-wavew*t - phase)
				enddo
			enddo
		endif
	end subroutine angle_wave
	
	subroutine plane_wave(fxye,t,nx,nxe,nvp,idproc)
		integer :: nx,nxe,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i
		real :: xstart,xpos,wavek,tempamp,fac, tfac

		real :: wavek2, t_local
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		fac = 1.
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		

		tfac = tfac * amp
		
		if (t < timerise + timeflat + timefall) then
		
			! Array position 2 corresponds to x=0 cause of gaurdcells
			do i = 1, nxe
				xpos = real(i)-2.
				fxye(1,i,:,1) = fxye(1,i,:,1) + tfac * sin(wavek*xpos-wavew*t)
			enddo
		endif
		
! second driver

        if(amp2 .ne. 0) then
	wavek2 = real(nx)/wavemode2
        wavek2=6.283185307/wavek2
        fac=1.
        t_local=(t-time_delay)
        if (t_local < 0) then
                tfac = 0
        else if(t_local < timerise2) then
        	tfac = t_local/timerise2
        else if (t_local < (timerise2 + timeflat2)) then
        	tfac = 1.0
        else if (t_local < (timerise2 + timeflat2 + timefall2)) then
        	tfac = 1 - (t_local - (timerise2+timeflat2))/timefall2
        endif

        tfac = tfac * amp2
        if (t < timerise2 + timeflat2 + timefall2) then
		
			! Array position 2 corresponds to x=0 cause of gaurdcells
			do i = 1, nxe
				xpos = real(i)-2.
				fxye(1,i,:,1) = fxye(1,i,:,1) + tfac * sin(wavek2*xpos-wavew2*t)
			enddo
	endif
        endif
	end subroutine plane_wave

	subroutine plane_wave_in_y(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i, j, nyproc
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,ypos,wavek,tempamp,fac, tfac
		
		wavek = real(ny)/wavemode
		wavek = 6.283185307/wavek
		fac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		tfac = tfac * amp
		
		if (t < timerise + timeflat + timefall) then
		
			! Array position 2 corresponds to x=0 cause of gaurdcells
!			do i = 1, nxe
!				xpos = real(i)-2.
!				fxye(1,i,:,1) = fxye(1,i,:,1) + tfac * sin(wavek*xpos-wavew*t)
!			enddo
!			print *, 'at proc_pos ', proc_pos, ' at time t ', t, ' with nypmx ',nypmx
			do j = 1, nypmx
				ypos = real(proc_pos(idproc+1) + j - 2)
!				print *, ypos
				do i = 2, nx+1
!					xpos = real(i)-2.
				!do i = 2, nx
					fxye(2,i,j,1) = fxye(2,i,j,1) + tfac * sin(wavek*ypos-wavew*t)
				enddo
!				print *, fxye(2,1,j,1)
			enddo
			!print *, 'at proc_pos ', proc_pos, ' at time t ', t, ' with nypmx ',nypmx
			!print *, fxye(2,
		endif
		
	end subroutine plane_wave_in_y
	
	subroutine plane_finite_wavelen(fxye,t,nx,nxe,nvp,idproc)
		integer :: nx,nxe,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,stoppos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: tempamp, tempamp2
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
			tempamp = amp * tfac
		
			! Array position 2 corresponds to x=0 cause of gaurdcells      	
			stoppos = floor(r_f_f)+3
			do i = 2, stoppos-2
				xpos = real(i)-2.
				if (xpos < rise) then
					xfac = xpos * riseinv
				else if (xpos < r_f) then
					xfac = 1.
				else if (xpos < r_f_f) then
					xfac = mfallinv*xpos + falloffset
				endif
				
				tempamp2 = tempamp*xfac      		
				fxye(1,i,:,1) = fxye(1,i,:,1) + &
					& tempamp2*sin(wavek*xpos-wavew*t)
			enddo
		endif
	end subroutine plane_finite_wavelen

	subroutine plane_finite_wavelen_supergauss(fxye,t,nx,nxe,nvp,idproc)
		integer :: nx,nxe,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,stoppos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: tempamp, tempamp2
		real :: mfallinv,falloffset,r_f,r_f_f,invW_n
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.

		!treat rise as W, and flat at x offset for center of wave
		invW_n = -0.5 * (1. / rise)**superGauss

		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
			tempamp = amp * tfac
		
			! Array position 2 corresponds to x=0 cause of gaurdcells      	
			stoppos = floor(r_f_f)+3
			do i = 2, 2+nx
				xpos = real(i)-2. - flat
				
				xfac = exp( invW_n * xpos**superGauss )
				
				tempamp2 = tempamp*xfac      		
				fxye(1,i,:,1) = fxye(1,i,:,1) + &
					& tempamp2*sin(wavek*xpos-wavew*t)
			enddo
		endif
	end subroutine plane_finite_wavelen_supergauss

	
	subroutine gauss_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f,yrise_fall_inv,a1
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		yrise_fall_inv = 1. / real(yrise_fall)
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
			a1 = -.5 * superGauss * (yrise_fall_inv**superGauss) / wavek
			ylow = -4.*real(yrise_fall)
			yhigh = 4.*real(yrise_fall)

			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2 - ny/2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					ytemp = exp( -0.5 * (ypos * yrise_fall_inv)**2 )
				else
					ytemp = 0.
				endif


!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					if (ypos <= ny / 2) then
!						ytemp = ( ypos - ylow ) / yrise_fall
!						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!					else 
!						ytemp = ( ypos - real(ny/2) ) / yrise_fall
!						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
!					endif
!				else
!					ytemp = 0.
!				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
		
					fxye(1,i,j,1) = fxye(1,i,j,1) + &
						& tempamp2*sin(wavek*xpos-wavew*t)
				enddo
			enddo
		endif
	end subroutine gauss_tran_finite_wavelen
	
	subroutine doub_gauss_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
			!Do driver centered around center1 first
			!both driver one and two have the same yrise_fall
			ylow = real(center1-yrise_fall)
			yhigh = real(center1+yrise_fall)
			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					if (ypos <= center1) then
						ytemp = ( ypos - ylow ) / yrise_fall
						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
					else 
						ytemp = ( ypos - real(center1) ) / yrise_fall
						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
					endif
				else
					ytemp = 0.
				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
		
					fxye(1,i,j,1) = fxye(1,i,j,1) + &
						& tempamp2*sin(wavek*xpos-wavew*t)
				enddo
			enddo
		
			!Now do driver centered around center2
			!both driver one and two have the same yrise_fall
			ylow = real(center2-yrise_fall)
			yhigh = real(center2+yrise_fall)
			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					if (ypos <= center2) then
						ytemp = ( ypos - ylow ) / yrise_fall
						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
					else 
						ytemp = ( ypos - real(center2) ) / yrise_fall
						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
					endif
				else
					ytemp = 0.
				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
		
					fxye(1,i,j,1) = fxye(1,i,j,1) + &
						& tempamp2*sin(wavek*xpos-wavew*t)
				enddo
			enddo
		
		endif
	end subroutine doub_gauss_tran_finite_wavelen
	
	
	subroutine gauss_tran_per_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
		ylow = real(ny/2-yrise_fall)
		yhigh = real(ny/2+yrise_fall)
		do j = 1,nypmx
	
			tempamp = amp * tfac
	
			ypos = real(proc_pos(idproc+1) + j - 2)
			if (ypos >= ylow .AND. ypos <= yhigh) then
				if (ypos <= ny / 2) then
					ytemp = ( ypos - ylow ) / yrise_fall
					ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
				else 
					ytemp = ( ypos - real(ny/2) ) / yrise_fall
					ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
				endif
			else
				ytemp = 0.
			endif
	
			tempamp = tempamp * ytemp
	
			! Array position 2 corresponds to x=0 cause of gaurdcells
			do i = 2, nx
				xpos = real(i)-2.
				
				fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t)
			enddo
		enddo
	endif
	end subroutine gauss_tran_per_wavelen
	
	subroutine gauss_tran_per_wavelen_circle_wavefronts(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend,yarr
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		real,allocatable,save,dimension(:) :: profile,phase_off
		real :: r, phi
						
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		!Initialize phase_offset array
		if (.not. allocated(phase_off)) then
			allocate(phase_off(0:int(yrise_fall)))
			r = ( phase_offset**2 + yrise_fall**2 ) / 2. / phase_offset
		!				print*,"r",r
			phase_off(0) = -1.*phase_offset
			do j = 1, int(yrise_fall)
				phi = asin( (real(j)) / r)
				phase_off(j) = -((real(j)) / tan(phi)-(r-phase_offset))
		!					print*,"phase_off",j,phase_off(j)-(r-phase_offset)
			enddo
		endif
		
		!Initialize transverse profile array
		if (.not. allocated(profile)) then
			allocate(profile(0:int(yrise_fall)))
			do j = 0,int(yrise_fall)
				ypos = real(j)
				ytemp = (yrise_fall-ypos) / yrise_fall
				profile(j) = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
			enddo
		endif
					
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			do j = 1,nypmx
				
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					yarr = abs( ny/2 - int(ypos))
		
					tempamp = amp * tfac * profile(yarr)
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
					do i = 2, nx
						xpos = real(i)-2.
						
						fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t + phase_off( yarr ) )
					enddo
				endif
			enddo
		endif
	end subroutine gauss_tran_per_wavelen_circle_wavefronts
	
	!This driver is like gauss_tran_per_wavelen but also adds a component in Ey to make wave like one
	! that would come from a potential with gaussian profile
	subroutine gauss_tran_per_wavelen_pois(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend,yarr
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp, kW2inv, yfac
		real :: mfallinv,falloffset,r_f,r_f_f
		real,allocatable,save,dimension(:) :: profile
		real :: r, phi
						
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		kW2inv = 1./wavek/yrise_fall/yrise_fall
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		!Initialize transverse profile array
		if (.not. allocated(profile)) then
			allocate(profile(0:int(yrise_fall)))
			do j = 0,int(yrise_fall)
				ypos = real(j)
				ytemp = (yrise_fall-ypos) / yrise_fall
				profile(j) = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
			enddo
		endif
					
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			do j = 1,nypmx
				
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					yarr = abs( ny/2 - int(ypos))
		
					tempamp = amp * tfac * profile(yarr)
					yfac = tempamp*(ypos-real(ny/2))*kW2inv
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
					do i = 2, nx
						xpos = real(i)-2.
						
						fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t)
						fxye(2,i,j,1) = fxye(2,i,j,1) + yfac*cos(wavek*xpos-wavew*t)
					enddo
				endif
			enddo
		endif
	end subroutine gauss_tran_per_wavelen_pois
	
	subroutine rect_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		tfac = 1.
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
			
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			do j = 1,nypmx
		
				tempamp = amp * tfac
		
				ypos = real(proc_pos(idproc+1) + j - 2)
				if (ypos >= ylow .AND. ypos <= yhigh) then
					ytemp = 1.
				else
					ytemp = 0.
				endif
		
				tempamp = tempamp * ytemp
		
				! Array position 2 corresponds to x=0 cause of gaurdcells
				stoppos = floor(r_f_f)+3
				do i = 2, stoppos-2
					xpos = real(i)-2.
					if (xpos < rise) then
						xfac = xpos * riseinv
					else if (xpos < r_f) then
						xfac = 1.
					else if (xpos < r_f_f) then
						xfac = mfallinv*xpos + falloffset
					endif
					
					tempamp2 = tempamp*xfac  
		
					fxye(1,i,j,1) = fxye(1,i,j,1) + &
						& tempamp2*sin(wavek*xpos-wavew*t)
				enddo
			enddo
		endif
	end subroutine rect_tran_finite_wavelen
	
	subroutine rect_tran_per_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
	integer :: nx,nxe,ny,nypmx,nvp,idproc
	real :: t
	real,dimension(:,:,:,:) :: fxye
	integer :: i,j,stoppos, nyproc,ystart,yend
	integer, dimension(nvp) :: proc_pos
	real :: xstart,xpos,wavek,riseinv,tfac,xfac
	real :: ylow, yhigh, ypos
	real :: tempamp, tempamp2,ytemp
	real :: mfallinv,falloffset,r_f,r_f_f
	
	wavek = real(nx)/wavemode
	wavek = 6.283185307/wavek
	r_f = rise + flat
	r_f_f = rise + flat + fall
	riseinv = 1. / rise
	mfallinv = -1. / fall
	falloffset = rise/fall + flat/fall + 1.
	tfac = 1.
	
	nyproc = ny / nvp
	do i = 0, nvp-1
		proc_pos(i+1) = ny / nvp * (i)
	enddo
	
	if (t < timerise) then
		tfac = t / timerise
	else if (t < timerise + timeflat) then
		tfac = 1.
	else if (t < timerise + timeflat + timefall) then
		tfac = 1. - (t - (timerise+timeflat))/timefall
	endif
	
	if (t < timerise + timeflat + timefall) then
		
		ylow = real(ny/2-yrise_fall)
		yhigh = real(ny/2+yrise_fall)
		do j = 1,nypmx
	
			tempamp = amp * tfac
	
			ypos = real(proc_pos(idproc+1) + j - 2)
			if (ypos >= ylow .AND. ypos <= yhigh) then
				ytemp = 1.
			else
				ytemp = 0.
			endif
	
			tempamp = tempamp * ytemp
	
			! Array position 2 corresponds to x=0 cause of gaurdcells
			stoppos = floor(r_f_f)+3
			do i = 1, nxe
				xpos = real(i)-2.						
				fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t)
			enddo
		enddo
	endif
	end subroutine rect_tran_per_wavelen

	subroutine supergauss_tran_per_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac, pshift, Ey, a1,p, a2
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp, yrise_fall_inv
		real :: mfallinv,falloffset,r_f,r_f_f
		
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		yrise_fall_inv = 1. / real(yrise_fall)
		
		nyproc = ny / nvp
		do i = 0, nvp-1
			proc_pos(i+1) = ny / nvp * (i)
		enddo
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
		
		a1 = -.5 * superGauss * (yrise_fall_inv**superGauss) / wavek
		ylow = -4.*real(yrise_fall)
		yhigh = 4.*real(yrise_fall)
		do j = 1,nypmx
	
			tempamp = amp * tfac
	
			ypos = real(proc_pos(idproc+1) + j - 2 - ny/2)
			if (ypos >= ylow .AND. ypos <= yhigh) then
				ytemp = exp( -0.5 * (ypos * yrise_fall_inv)**superGauss )
			else
				ytemp = 0.
			endif
	
			tempamp = tempamp * ytemp
			Ey = a1 * (ypos**(superGauss-1)) * tempamp
			a2 = bow*bow_power * ytemp**bow_power
			
			pshift = -1.* bow*(ytemp**bow_power)
	
			! Array position 2 corresponds to x=0 cause of gaurdcells
			do i = 2, nx+1
				xpos = real(i)-2.
				p = wavek*xpos-wavew*t + pshift
				fxye(1,i,j,1) = fxye(1,i,j,1) - tempamp*sin(p)
				fxye(2,i,j,1) = fxye(2,i,j,1) + Ey * (a2*sin(p) - cos(p))
			enddo
		enddo
	endif
	end subroutine supergauss_tran_per_wavelen


end module ext_driver_jf
