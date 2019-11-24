!-----------------------------------------------------------------------
!
      module pinit2d_jf
!

      implicit none
      public :: ntfield
      public :: sendnml_jf
      public :: amp,wavemode,wavew,angle
!
! FST -> BVP
! August, 2019
! additional input parameters for boundary value problems (BVP)
! antenna is located @ the middle (n/2) and slight right (n/2+1) 
! of the box.  The value @ mid-point is exactly as specified by 
! the input deck and the n/2+1 cell is the exact opposite.  The 
! idea is to mimic a shielded potential/field source at some 
! location x0. 
!

      public :: ant_amp, ant_omega, ant_trise, ant_tflat, ant_tfall

!
! FST -> Jan 2019
! ECHO
! second antenna for the echo problem
      public :: time_delay, timerise2,timeflat2,timefall2
      public :: amp2, wavemode2, wavew2
! ECHO
      public :: fvxmax,fvymax,nphbx,nphby,nphxx,nphyx,nphxy,nphyy
      public :: fvxmax_ion,fvymax_ion
      public :: driver_select,rise,flat,fall,timerise,timeflat,timefall,phase_offset
      public :: yrise_fall,linepos,ntlines,center1,center2
      public :: phsl_x_pos,phsl_x_thick,ntphsl_x,num_phsl_x
      public :: phsl_y_pos,phsl_y_thick,ntphsl_y,num_phsl_y
      public :: ntden,ntpene,ntj,ntvdotE
      public :: ntraw,nttrack,nt_dump_track,raw_fraction,raw_range,raw_dump_tagged,track_no_h5
      public :: ntESPoynt,ntESPoynt_int
      public :: keep_init_vel,initEneRange,finalEneRange,nEneBin_init,nEneBin_final
      public :: initVxLow,initVxHigh,finalVxLow,finalVxHigh,nVxBin_init,nVxBin_final
      public :: initVyLow,initVyHigh,finalVyLow,finalVyHigh,nVyBin_init,nVyBin_final
      public :: nvdotE_part,nvdotE_int,fvx_xy_max,nphb_xy,ntphxy,fvy_xy_max
      public :: fvxy_xy_xrange,fvxy_xy_yrange,nfvxy_xy_xrange,nfvxy_xy_yrange
      public :: nvdotE_follow_part,ntfield_ene_int
      public :: dump_start,dump_end,nt_div_ESPoynt_int,nt_div_ESPoynt,nt_phase_vx_vy
      public :: phase_keep_only_tracked,nt_write_field_ene_in_driver_region,int_pos,nt_write_ene_int
      public :: nt_write_Py_vs_y,nt_write_div_P_vs_y,nt_write_jE_sumover_x,nt_write_Py_sumover_x
      public :: nt_write_U_sumover_x, npfac, nt_write_Ux_sumover_x
      public :: nt_bin_E,bin_E,bin_E_xrange,bin_E_yrange,bin_E_Exrange,bin_E_Eyrange,raw_cycle_x_move
      public :: superGauss,bow,bow_power
      public :: nt_vx_vy_speed, nvx_vy_speed, vx_vy_xrange, vx_vy_yrange, vx_vy_speed
      public :: nt_write_U_sumover_x_fromE,nt_write_grad_phi,ampere_k0,nt_write_S_sumover_x,nt_b_field
      public :: turn_off_self_con
      public :: nt_through_wave, init_range, final_y,nt_dEdt
      public :: nt_write_kE_sumover_x, nt_kE,nt_write_jE_onlytracked_sumover_x
      
      save
      
      integer :: ntfield = 0, nt_b_field = 0
			! ntfield writes the fxye array
			! nt_write_grad_phi writes the grad_phi array
      integer :: nt_write_grad_phi = 0
      !wavemode is the mode number for the external driver, the k of the wave is 2*pi/(nx/wavemode)
      ! wavew is the frequency of the driver
      integer :: wavemode=0
      real :: amp = 0., wavew = 0.,angle = 0.
! BVP
      real :: ant_amp, ant_omega, ant_trise, ant_tflat, ant_tfall
      real :: ant_amp2, ant_omega2, ant_trise2, ant_tflat2, ant_tfall2
! BVP 
! ECHO
      integer :: wavemode2 = 0
      real :: amp2 = 0., wavew2 = 1.0 
! ECHO
      ! these two are the ranges for the phase space velocities
      ! the phase space will be from v = -fvxmax to fvxmax
      real :: fvxmax = 10., fvymax = 10.
      real :: fvxmax_ion = 1., fvymax_ion = 1.
      ! nphbx,nphby = # phase space bins in vx,vy
      ! nphxx,nphyx,nphxy,nphyy = dump factor for phase space in vx vs. x, vx vs. y, and so on
      integer :: nphbx = 10, nphby = 10, nphxx = 0, nphyx = 0,nphxy=0,nphyy=0
      ! see main loop for list of drivers and their code
      integer :: driver_select = 0
	  	! rise, flat, and fall times for finite wavelength driver
			! these are given in # of wavelengths!
			! so, for example, rise = 3 means the rise distance is 3*lambda
      real :: rise=1.,flat=1.,fall=1.
			! rise time flat time and fall time
! ECHO
      real :: timerise = 0.,timeflat = 1000000., timefall = 0.
      real :: timerise2 = 0., timeflat2 = 1.0, timefall2 = 0. 
      real :: time_delay = 1.0
! ECHO
      ! yrise_fall =  spatial rise and fall for the transverse profile
      ! phase_offset is the distance that the phase of the center of the transverse profile is 
      !  offset from the top and bottom edges, for use with the circular wavefront driver
      real :: yrise_fall, phase_offset
      !these two specify the center of the two drivers for doub_gauss_tran_finite_wavelen
      integer :: center1,center2
      ! ntlines specifies how often to dump the lines data
      integer,dimension(5) :: linepos
      integer :: ntlines = 0
      ! phsl_x_pos is an array that lists the positions for each thin phasespace dump
      ! phsl_x_thick is the thickness of the sum over the y direction for the
      ! narrow phasespace.  If the range over y to sum falls across processors then 
      ! it will only get those particles in the lowest processor.  Only five places
      ! work for now, if values is set to -1 then don't do it.
      ! those with x in the name are a lineout in y, with y in the name a lineout in x
      integer,dimension(5) :: phsl_x_pos,phsl_y_pos
      integer :: phsl_x_thick=2,ntphsl_x = 0,num_phsl_x=0
      integer :: phsl_y_thick=2,ntphsl_y = 0,num_phsl_y=0
      ! ntden is how often to dump density, ntpene is particle energy dump,
      ! ntj is current dump
      ! ntvdotE is how often to dumb v dot E
      integer :: ntden=0,ntpene=0,ntj=0,ntvdotE=0
      ! ntraw is how often to dump raw particle data, nt_dump_track is how often to
      ! flush the track data from the buffer to disk
      ! raw_fraction is fraction of
      ! particles to dump data for, nttrack adds the tag data to track the particles
      ! raw_dump_tagged means dump to raw file only tagged particles, 0 off, 1 on
      ! it cannot be used with nttrack
      ! track_no_h5 writes track data in plain Fortran unformatted to be converted to H5 in postprocess
      !   using the program convert_h5.e contained in the convert_h5 folder.  This is because
      !   I had trouble with the h5 file writes taking forever, so I just write unformatted fortran
      !   and convert it to h5 later.
      integer :: ntraw = 0, nttrack = 0, nt_dump_track = 0,track_no_h5=1
      ! raw_dump_tagged means dump to raw file only tagged particles, 0 off, 1 on
      ! it cannot be used with nttrack, file input_tags must be present to read list of
      ! particles to dump.
      ! raw_cycle_x_move moves the window specified by raw_range at the specified speed 
      ! in the x direction
      real :: raw_fraction=0., raw_dump_tagged = 0, raw_cycle_x_move = -1.
      ! raw_range is the 4D cube inside of which the raw data is stored for the particles
      ! raw_range(1,*) is the lower and upper boundary in x
      ! raw_range(2,*) is the lower and upper boundary in y
      ! raw_range(3,*) is the lower and upper boundary in vx
      ! raw_range(4,*) is the lower and upper boundary in vy
      real,dimension(4,2) :: raw_range = -1.
      ! ntESPoynt is the electrostatic component of the Poynting vector 
      ! See VK Decyk Phys. Fluids 25, 1205 (1982)
      integer :: ntESPoynt = 0, ntESPoynt_int = 0
      ! keep_init_vel = 1 means to add two more arrays to particle array that 
      ! store the initial velocity of the particles so comparisons and statistics
      ! can be calculated at the end of the run
      integer :: keep_init_vel = 0
      ! the EneRanges are the maximum energies to bin for the initial and final energy diagnostic
      real :: initEneRange,finalEneRange
      ! nEneBin_init and nEneBin_final are the number of bins for the diagnostic above
      integer :: nEneBin_init,nEneBin_final
      ! same as Ene stuff above, but for vx and must give high and low values of vx to save
      real :: initVxLow,initVxHigh,finalVxLow,finalVxHigh
      real :: initVyLow,initVyHigh,finalVyLow,finalVyHigh
      integer :: nVxBin_init,nVxBin_final
      integer :: nVyBin_init,nVyBin_final
      ! These two calculate the v dot E for each particle and store it in the particle array.
      ! For diagnostics, this is then deposited just like the current would be deposited.  
      ! These are better diagnostics than the other vdotE (ntvdotE) and are meant to replace it.
      ! nvdotE_part dumps vdot E at specified interval, while nvdotE_int sums the vdotE at every
      ! timestep  and then dumps the running sum at spedified intervals.
      ! nvdotE_part MUST be specified to use nvdotE_int, but if one is nonzero, then BOTH must be.
      integer :: nvdotE_part=0,nvdotE_int=0
      real :: fvx_xy_max=10.,fvy_xy_max=10.
      ! These specifiy the spatial extent of the phase space
      real,dimension(2) :: fvxy_xy_xrange=(/-1.,-1./),fvxy_xy_yrange=(/-1.,-1./)
      ! These are the number of velocity bins, timestep skip, number of bins in x and y
      integer :: nphb_xy,ntphxy=0,nfvxy_xy_xrange=0,nfvxy_xy_yrange=0
			! if nvdotE_follow_part is set, then the particles will carry their total vdotE rather than
			! their instantaneous vdotE.  This cannot be used with either nvdotE_int or nvdotE_part
      integer :: nvdotE_follow_part=0
      integer :: ntfield_ene_int = 0		!this will write out the integrated field energy
      ! These two are to specify at what time to start and stop dumping diagnostic data
      ! So far only works for e-field and raw, negative means dump from beginning to end
      real :: dump_start=-1.,dump_end=-1.
      ! This dumps the divergence and integrated divergence of the ES Poynting vector
      integer :: nt_div_ESPoynt = 0, nt_div_ESPoynt_int = 0
      ! nt_phase_vx_vy plots vx vs. vy phase space
      integer :: nt_phase_vx_vy
      ! phase_keep_only_tracked, this deposits only the tracked particles to the phase space 
      ! diagnostics.  Therefore, must be used with the tracking stuff on, but just set them
      ! so that they are too big to ever actually be written.
      integer :: phase_keep_only_tracked
      ! nt_write_field_ene_in_driver_region writes electric field energy in the region:
      ! -yrise_fall < y - ny/2 < yrise_fall and all x to the formatted file on unit 77
      integer :: nt_write_field_ene_in_driver_region = 0
      ! int_pos are the y positions of the three boxes used in integrating the Poynting flow energy diagnostic
      ! int_pos(1) is y pos of bottom of bottom box, int_pos(2) is top of bottom box and bottom of middle box
      ! and so on
      integer,dimension(4) :: int_pos=(/-1.,-1.,-1.,-1./)
      ! nt_write_ene_int writes to file the integrated energies for the boxes indicated in int_pos
      integer :: nt_write_ene_int = 0
      ! nt_write_Py_vs_y sums Py over all x and writes it for each y position
      integer :: nt_write_Py_vs_y = 0, nt_write_div_P_vs_y = 0
      ! nt_write_jE_sumover_x sums jE over all x and writes it for each y position
      ! nt_write_Py_sumover_x sums Py over all x and writes it for each y position
      ! nt_write_Ux_sumover_x sums U over just the x component and ignores the part from Ey
      ! nt_write_U_sumover_x_fromE calculates U from the fxye array, not from grad_phi
      integer :: nt_write_jE_sumover_x = 0,nt_write_Py_sumover_x = 0, nt_write_U_sumover_x = 0
      integer :: nt_write_Ux_sumover_x = 0, nt_write_U_sumover_x_fromE
      ! npfac is the percentage increase to np for allocating the part array to ensure that when
      ! new particles come onto the current processor the part array isn't over filled.
      real :: npfac = 0.05
      ! nt_bin_E writes Ex and Ey as 3d data, bin_E is the number of bins in Ex and Ey
      ! bin_E_xrange and bin_E_yrange are the space ranges over which to do the binning, must be
      ! 	cells, cannot do finer or coarser grid than the simulation grid
      ! bin_E_Exrange,bin_E_Eyrange are the ranges in Ex and Ey to do it
      integer :: nt_bin_E = 0, bin_E = 100
      integer, dimension(2) :: bin_E_xrange = (/0,1/), bin_E_yrange = (/0,1/)
      real, dimension(2) :: bin_E_Exrange = (/-1.,1./), bin_E_Eyrange = (/-1.,1./)
      ! superGauss is the integer exponent of the super Gaussian used for the supergauss external driver
      integer :: superGauss = 4
      ! bow is the coefficient to the amplitude dependent phase shift in the driver, bow_power is the
      ! exponent of the amplitude in the phase shift
      real :: bow = 0., bow_power = 0.
      ! nt_vx_vy_speed is the number of timesteps between dumps for making nvx_vy_speed of vx_vy plots
      ! over the range of space given by vx_vy_xrange and vx_vy_yrange, with y shifting by y-vx_vy_speed*n,
      ! where n goes from 0 to nvx_vy_speed-1.
      integer :: nt_vx_vy_speed = 0, nvx_vy_speed = 5
      real, dimension(2) :: vx_vy_xrange = (/0.,1./), vx_vy_yrange = (/0.,1./)
      real :: vx_vy_speed = 1.
      ! ampere_k0 if set to one then ampere's law is used to calculate field for the k=0 modes
      integer :: ampere_k0 = 0
      ! nt_write_S_sumover_x like nt_write_jE_sumover_x except it's the EM poynting vector
      integer :: nt_write_S_sumover_x = 0
      ! turn_off_self_con turns off the self-consistent field.  The particles will only feel the driver
      ! force, there is no interaction between particles.  If there is no driver, then the particles will
      ! simply free stream
      integer :: turn_off_self_con = 0
      ! nt_through_wave (CANNOT be used with keep_init_vel) finds all particles that start within init_range
      ! that cross final_y at some time.  Those that never cross are collected with function collect_not_crossed
      ! init_range(4,2), first dim is 1=x,2=y,3=vx,4=vy (like part array), second coord is lower and upper boundary
      ! diag works by setting part(kivx_loc,i,1) = -1000000. for all particles not initially located within init_range
      ! once a particle crosses final_y, then it's coords are added to the diag arrays and its 
      ! part(kivx_loc,i,1) is set to -10000.  The diag arrays (delta v vs. vinit for each velocity) are written at
      ! then end of the simulation.
      ! nt_write_jE_onlytracked_sumover_x writes the jE of only those particles that would be dumped above
      integer :: nt_through_wave = 0, nt_write_jE_onlytracked_sumover_x = 0
      real :: final_y = 0.
      real,dimension(4,2) :: init_range
      ! nt_dEdt dumps the time derivative of E one timestep back
      integer :: nt_dEdt = 0
      ! nt_write_kE_sumover_x, nt_kE show the deposited kinetic energy
      integer :: nt_write_kE_sumover_x = 0, nt_kE = 0
      
      namelist /pinput2_jf/ ntfield, amp,wavemode,wavew,fvxmax,fvymax,&
        &amp2,wavemode2,wavew2,timerise2,timeflat2,timefall2,time_delay,&
      	&nphbx,nphby,nphxx,nphyx,rise,flat,fall,driver_select,timerise,&
      	&timeflat,timefall,yrise_fall,linepos,ntlines,phsl_x_pos,phsl_x_thick,ntphsl_x,&
      	&phsl_y_pos,phsl_y_thick,ntphsl_y,ntden,ntpene,ntj,ntvdotE, ntraw,nttrack, &
      	&raw_fraction, nt_dump_track,ntESPoynt,angle, keep_init_vel,initEneRange,finalEneRange,&
      	&nEneBin_init,nEneBin_final,initVxLow,initVxHigh,finalVxLow,finalVxHigh,nVxBin_init,&
      	&nVxBin_final,initVyLow,initVyHigh,finalVyLow,finalVyHigh,nVyBin_init,nVyBin_final,&
      	&nvdotE_part,nvdotE_int,fvx_xy_max,nphb_xy,ntphxy,fvy_xy_max,fvxy_xy_xrange,fvxy_xy_yrange,&
      	&nfvxy_xy_xrange,nfvxy_xy_yrange,nvdotE_follow_part,raw_range,ntESPoynt_int,ntfield_ene_int,&
      	&center1,center2,raw_dump_tagged,dump_start,dump_end,nt_div_ESPoynt_int,track_no_h5,phase_offset,&
      	&nt_phase_vx_vy,nphxy,nphyy,phase_keep_only_tracked,nt_write_field_ene_in_driver_region,int_pos,&
      	&nt_write_ene_int,nt_write_Py_vs_y,nt_div_ESPoynt,nt_write_div_P_vs_y,nt_write_jE_sumover_x,&
      	&nt_write_Py_sumover_x, nt_write_U_sumover_x, npfac, nt_bin_E,bin_E,bin_E_xrange,bin_E_yrange,&
      	&bin_E_Exrange,bin_E_Eyrange,raw_cycle_x_move,superGauss,bow,bow_power,&
      	&nt_vx_vy_speed, nvx_vy_speed, vx_vy_xrange, vx_vy_yrange, vx_vy_speed,nt_write_Ux_sumover_x,&
      	&nt_write_U_sumover_x_fromE,nt_write_grad_phi,ampere_k0,nt_write_S_sumover_x, nt_b_field,&
      	&turn_off_self_con,nt_through_wave,final_y,init_range,nt_dEdt,nt_write_kE_sumover_x, nt_kE,&
      	&nt_write_jE_onlytracked_sumover_x, fvxmax_ion,fvymax_ion,&
        &ant_amp,ant_omega,ant_trise,ant_tflat,ant_tfall,&
        &ant_amp2,ant_omega2,ant_trise2,ant_tflat2,ant_tfall2
      
      contains
      	subroutine sendnml_jf()
      		implicit none
                  ! LENML -> Namelist Length
				integer,parameter :: lenml = 174
                  ! LENML
				double precision, dimension(lenml) :: ddata
				ddata(1) = ntfield
				ddata(2) = amp
				ddata(3) = wavemode
				ddata(4) = wavew
				ddata(5) = fvxmax
				ddata(6) = fvymax
				ddata(7) = nphbx
				ddata(8) = nphby
				ddata(9) = nphxx
				ddata(10) = nphyx
				ddata(11) = driver_select
				ddata(12) = rise
				ddata(13) = flat
				ddata(14) = fall
				ddata(15) = timerise
				ddata(16) = timeflat
				ddata(17) = yrise_fall
				ddata(18:22) = linepos
				ddata(23) = ntlines
				ddata(24:28) = phsl_x_pos
				ddata(29) = phsl_x_thick
				ddata(30) = ntphsl_x
				ddata(31) = ntden
				ddata(32) = ntpene
				ddata(33) = ntj
				ddata(34) = ntvdotE
				ddata(35:39) = phsl_y_pos
				ddata(40) = phsl_y_thick
				ddata(41) = ntphsl_y
				ddata(42) = ntraw
				ddata(43) = nttrack
				ddata(44) = raw_fraction
				ddata(45) = nt_dump_track
				ddata(46) = ntESPoynt
				ddata(47) = angle
				ddata(48) = keep_init_vel
				ddata(49) = initEneRange
				ddata(50) = finalEneRange
				ddata(51) = nEneBin_init
				ddata(52) = nEneBin_final
				ddata(53) = initVxLow
				ddata(54) = initVxHigh
				ddata(55) = finalVxLow
				ddata(56) = finalVxHigh
				ddata(57) = nVxBin_init
				ddata(58) = nVxBin_final
				ddata(59) = initVyLow
				ddata(60) = initVyHigh
				ddata(61) = finalVyLow
				ddata(62) = finalVyHigh
				ddata(63) = nVyBin_init
				ddata(64) = nVyBin_final
				ddata(65) = nvdotE_part
				ddata(66) = nvdotE_int
				ddata(67) = fvx_xy_max
				ddata(68) = nphb_xy
				ddata(69) = ntphxy
				ddata(70) = fvy_xy_max
				ddata(71:72) = fvxy_xy_xrange
				ddata(73:74) = fvxy_xy_yrange
				ddata(75) = nfvxy_xy_xrange
				ddata(76) = nfvxy_xy_yrange
				ddata(77) = nvdotE_follow_part
				ddata(78:79) = raw_range(1,:)
				ddata(80:81) = raw_range(2,:)
				ddata(82:83) = raw_range(3,:)
				ddata(84:85) = raw_range(4,:)
				ddata(86) = ntESPoynt_int
				ddata(87) = ntfield_ene_int
				ddata(88) = center1
				ddata(89) = center2
				ddata(90) = raw_dump_tagged
				ddata(91) = dump_start
				ddata(92) = dump_end
				ddata(93) = nt_div_ESPoynt_int
				ddata(94) = track_no_h5
				ddata(95) = phase_offset
				ddata(96) = nt_phase_vx_vy
				ddata(97) = nphxy
				ddata(98) = nphyy
				ddata(99) = phase_keep_only_tracked
				ddata(100) = nt_write_field_ene_in_driver_region
				ddata(101:104) = int_pos
				ddata(105) = nt_write_ene_int
				ddata(106) = nt_write_Py_vs_y
				ddata(107) = nt_div_ESPoynt
				ddata(108) = nt_write_div_P_vs_y
				ddata(109) = nt_write_jE_sumover_x
				ddata(110) = nt_write_Py_sumover_x
				ddata(111) = nt_write_U_sumover_x
				ddata(112) = timefall
				ddata(113) = npfac
				ddata(114) = nt_bin_E
				ddata(115) = bin_E
				ddata(116:117) = bin_E_xrange
				ddata(118:119) = bin_E_yrange
				ddata(120:121) = bin_E_Exrange
				ddata(122:123) = bin_E_Eyrange
				ddata(124) = raw_cycle_x_move
				ddata(125) = superGauss
				ddata(126) = bow
				ddata(127) = bow_power
				ddata(128) = nt_vx_vy_speed
				ddata(129) = nvx_vy_speed
				ddata(130:131) = vx_vy_xrange
				ddata(132:133) = vx_vy_yrange
				ddata(134) = vx_vy_speed
				ddata(135) = nt_write_Ux_sumover_x
				ddata(136) = nt_write_U_sumover_x_fromE
				ddata(137) = nt_write_grad_phi
				ddata(138) = ampere_k0
				ddata(139) = nt_write_S_sumover_x
				ddata(140) = nt_b_field
				ddata(141) = turn_off_self_con
				ddata(142) = nt_through_wave
				ddata(143) = final_y
				ddata(144:147) = init_range(:,1)
				ddata(148:151) = init_range(:,2)
				ddata(152) = nt_dEdt
				ddata(153) = nt_write_kE_sumover_x
				ddata(154) = nt_kE
				ddata(155) = nt_write_jE_onlytracked_sumover_x				
                        ddata(156) = fvxmax_ion
                        ddata(157) = fvymax_ion
                        ! ECHO
                        ddata(158) = amp2
                        ddata(159) = wavemode2
                        ddata(160) = wavew2
                        ddata(161) = time_delay
                        ddata(162) = timerise2
                        ddata(163) = timeflat2
                        ddata(164) = timefall2
                        ! ECHO
                        !BVP
                        ddata(165) = ant_amp
                        ddata(166) = ant_omega
                        ddata(167) = ant_trise
                        ddata(168) = ant_tflat
                        ddata(169) = ant_tfall

                        ddata(170) = ant_amp
                        ddata(171) = ant_omega
                        ddata(172) = ant_trise
                        ddata(173) = ant_tflat
                        ddata(174) = ant_tfall
                        !BVP

				call PBCAST(ddata,lenml)
				ntfield = ddata(1)
				amp = ddata(2)
				wavemode = ddata(3)
				wavew = ddata(4)
				fvxmax = ddata(5)
				fvymax = ddata(6)
				nphbx = ddata(7)
				nphby = ddata(8)
				nphxx = ddata(9)
				nphyx = ddata(10)
				driver_select = ddata(11)
				rise = ddata(12)
				flat = ddata(13)
				fall = ddata(14)
				timerise = ddata(15)
				timeflat = ddata(16)
				yrise_fall = ddata(17)
				linepos = ddata(18:22)
				ntlines = ddata(23)
				phsl_x_pos = ddata(24:28)
				phsl_x_thick = ddata(29)
				ntphsl_x = ddata(30)
				ntden = ddata(31)
				ntpene = ddata(32)
				ntj = ddata(33)
				ntvdotE = ddata(34)
				phsl_y_pos = ddata(35:39)
				phsl_y_thick = ddata(40)
				ntphsl_y = ddata(41)
				ntraw = ddata(42)
				nttrack = ddata(43)
				raw_fraction = ddata(44)
				nt_dump_track = ddata(45)
				ntESPoynt = ddata(46)
				angle = ddata(47)
				keep_init_vel = ddata(48)
				initEneRange = ddata(49)
				finalEneRange = ddata(50)
				nEneBin_init = ddata(51)
				nEneBin_final = ddata(52)
				initVxLow = ddata(53)
				initVxHigh = ddata(54)
				finalVxLow = ddata(55)
				finalVxHigh = ddata(56)
				nVxBin_init = ddata(57)
				nVxBin_final = ddata(58)
				initVyLow = ddata(59)
				initVyHigh = ddata(60)
				finalVyLow = ddata(61)
				finalVyHigh = ddata(62)
				nVyBin_init = ddata(63)
				nVyBin_final = ddata(64)
				nvdotE_part = ddata(65)
				nvdotE_int = ddata(66)
				fvx_xy_max = ddata(67)
				nphb_xy = ddata(68)
				ntphxy = ddata(69)
				fvy_xy_max = ddata(70)
				fvxy_xy_xrange = ddata(71:72)
				fvxy_xy_yrange = ddata(73:74)
				nfvxy_xy_xrange = ddata(75)
				nfvxy_xy_yrange = ddata(76)
				nvdotE_follow_part = ddata(77)
				raw_range(1,:) = ddata(78:79)
				raw_range(2,:) = ddata(80:81)
				raw_range(3,:) = ddata(82:83)
				raw_range(4,:) = ddata(84:85)
				ntESPoynt_int = ddata(86)
				ntfield_ene_int = ddata(87)
				center1 = ddata(88)
				center2 = ddata(89)
				raw_dump_tagged = ddata(90)
				dump_start = ddata(91)
				dump_end = ddata(92)
				nt_div_ESPoynt_int = ddata(93)
				track_no_h5 = ddata(94)
				phase_offset = ddata(95)
				nt_phase_vx_vy = ddata(96)
				nphxy = ddata(97)
				nphyy = ddata(98)
				phase_keep_only_tracked = ddata(99)
				nt_write_field_ene_in_driver_region = ddata(100)
				int_pos = ddata(101:104)
				nt_write_ene_int = ddata(105)
				nt_write_Py_vs_y = ddata(106)
				nt_div_ESPoynt = ddata(107)
				nt_write_div_P_vs_y = ddata(108)
				nt_write_jE_sumover_x = ddata(109)
				nt_write_Py_sumover_x = ddata(110)
				nt_write_U_sumover_x = ddata(111)
				timefall = ddata(112)
				npfac = ddata(113)
				nt_bin_E = ddata(114)
				bin_E = ddata(115)
				bin_E_xrange = ddata(116:117)
				bin_E_yrange = ddata(118:119)
				bin_E_Exrange = ddata(120:121)
				bin_E_Eyrange = ddata(122:123)
				raw_cycle_x_move = ddata(124)
				superGauss = ddata(125)
				bow = ddata(126)
				bow_power = ddata(127)
				nt_vx_vy_speed = ddata(128)
				nvx_vy_speed = ddata(129)
				vx_vy_xrange = ddata(130:131)
				vx_vy_yrange = ddata(132:133)
				vx_vy_speed = ddata(134)
				nt_write_Ux_sumover_x = ddata(135)
				nt_write_U_sumover_x_fromE = ddata(136)
				nt_write_grad_phi = ddata(137)
				ampere_k0 = ddata(138)
				nt_write_S_sumover_x = ddata(139)
				nt_b_field = ddata(140)
				turn_off_self_con = ddata(141)
				nt_through_wave = ddata(142)
				final_y = ddata(143)
				init_range(:,1) = ddata(144:147)
				init_range(:,2) = ddata(148:151)
				nt_dEdt = ddata(152)
				nt_write_kE_sumover_x = ddata(153)
				nt_kE = ddata(154)
				nt_write_jE_onlytracked_sumover_x = ddata(155)
                                fvymax_ion = ddata(156)
                                fvymax_ion = ddata(157)
                                ! ECHO
                                amp2 = ddata(158)
                                wavemode2 = ddata(159)
                                wavew2 = ddata(160)
                                time_delay = ddata(161)
                                timerise2 = ddata(162)
                                timeflat2 = ddata(163)
                                timefall2 = ddata(164)
                                ! ECHO
                        ! BVP
                        ant_amp = ddata(165)
                        ant_omega = ddata(166)
                        ant_trise = ddata(167)
                        ant_tflat = ddata(168)
                        ant_tfall = ddata(169)

                        ant_amp2 = ddata(170)
                        ant_omega2 = ddata(171)
                        ant_trise2 = ddata(172)
                        ant_tflat2 = ddata(173)
                        ant_tfall2 = ddata(174)
                        ! BVP
 

      	end subroutine sendnml_jf
      	
      end module pinit2d_jf
