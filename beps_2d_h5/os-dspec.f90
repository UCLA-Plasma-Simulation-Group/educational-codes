!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     species diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module m_diag_species 

        use m_system
!        use m_debug
        use m_parameters
        use m_file_system

        use m_node_conf
        use m_time_avg

        use m_diagnostic_utilities
        
        use m_utilities
        use m_units
        
        use m_space
        use m_species_define
        use m_species_tracks
        
        use m_error_handling

        implicit none

!       restrict access to things explicitly declared public
        private



!       HDF Function declaration.

        integer,external :: sfstart, sfcreate, sfwdata, sfsdtstr, sfsdmstr
        integer,external :: sfdimid, sfsdmname, sfsdscale, sfsblsz 
        integer,external :: sfendacc, sfend, sfsnatt, sfscompress

!       HDF Constant declaration

        integer, parameter :: DFACC_CREATE = 4
        integer, parameter :: DFACC_WRITE = 2
        integer, parameter :: DFNT_CHAR = 4
        integer, parameter :: DFNT_FLOAT32 = 5
        integer, parameter :: DFNT_FLOAT64 = 6
        integer, parameter :: DFNT_INT32 = 24
        integer, parameter :: COMP_CODE_DEFLATE = 4
        integer, parameter :: SD_UNLIMITED = 0
        integer, parameter :: SD_FILL = 0
        
!       End HDF Declarations 

        ! maximum number of bins for ene_bin diagnostics
        integer(p_k_is), parameter :: max_n_ene_bins = 32
 
        ! size of particle cache for phasespaces
        integer, parameter :: p_par_buf_size = 4096

!        type :: t_diag_species
!
!!         allow access to type components only to module procedures
!          private
!
!!         frequency of species dianostic data dumps
!          integer(p_k_is) :: ndump_fac_pha, ndump_fac_ene, ndump_fac_raw
!          integer(p_k_is) :: ndump_fac_charge
!
!!         file_name prefix for data files
!          character(80) :: file_name_pha, file_name_ene
!
!!         physical range for phasespace data dumps
!          real(p_k_rdoub), dimension(p_x_dim) :: ps_xmin, ps_xmax
!          real(p_k_rdoub), dimension(p_p_dim) :: ps_pmin, ps_pmax
!          real(p_k_rdoub), dimension(p_p_dim) :: ps_pmin_r, ps_pmax_r
!
!          
!!         switch for autorange of ps_p
!          logical, dimension(p_p_dim) :: if_ps_p_auto
!
!!         physical range for 1D momenta phasespace data dumps
!          real(p_k_rdoub) :: ps_gammamin, ps_gammamax
!          real(p_k_rdoub) :: ps_gammamin_r, ps_gammamax_r
!
!!         switch for autorange of ps_gamma
!          logical :: if_ps_gamma_auto
!
!!         switch for log deposition of gamma
!          logical :: if_ps_gamma_log
!          
!
!!         resolutions for phasespace data dumps
!          integer(p_k_is), dimension(p_x_dim) :: ps_nx, ps_nx_3D
!          integer(p_k_is), dimension(p_p_dim) :: ps_np, ps_np_3D
!
!!         resolution for 1D momenta phasespace data dump
!          integer(p_k_is) :: ps_ngamma
!
!!         deposition scheme for phase space calculations
!          integer(p_k_is) :: dep_sch
!
!!         switch to turn on and off species 1D momenta phase space diags
!          logical :: if_gamma
!
!!         switches to turn on and off species phase space diagnostics
!          logical :: if_x2x1, if_x3x1, if_p1x1, if_p2x1, if_p3x1
!          logical ::          if_x3x2, if_p1x2, if_p2x2, if_p3x2 
!          logical ::                   if_p1x3, if_p2x3, if_p3x3
!          logical ::                            if_p2p1, if_p3p1
!          logical ::                                     if_p3p2
!
!!         switches to turn on and off spatially resolved energy diagnoestics
!          logical :: if_x2x1_ene
!          logical :: if_x3x2x1_ene, if_x3x1_ene, if_x3x2_ene 
!          
!!         switches and parameters to turn on and off energy binned species diagnostics
!          logical :: if_x2x1_bin_ene
!          logical :: if_x3x2x1_bin_ene, if_x3x1_bin_ene, if_x3x2_bin_ene
!          integer(p_k_is) :: n_ene_bins
!
!          real(p_k_rsing), dimension(max_n_ene_bins) :: ene_bins  
!          
!!         switchs for spectrum diagnostics
!          logical :: if_spectrum
!          integer :: iter_spectrum
!          
!          
!!!         switches to turn on and off energy binned species diagnostics
!!          logical ::  if_nv_ene, if_nv_gamma_ene
!!          real(p_k_rdoub), dimension(2) :: ene_bin
!
!!         switches to turn on and off 3D species phase space diagnostics
!          logical :: if_x3x2x1
!          logical :: if_p1x2x1, if_p2x2x1, if_p3x2x1
!
!!         parameters for raw data dump
!          real(p_k_rdoub) :: gamma_limit
!          real(p_k_rdoub) :: particle_fraction
!
!!         switches to turn on and off other species diagnostics 
!          logical :: iforigin, iftestpar ! currently not used
! 
!!         for diagnostic purposes : - currently not used
!!         original position, momentum, and charge of particles
!          real(p_k_rdoub), dimension(:,:), pointer :: xor => NULL()
!          real(p_k_rdoub), dimension(:,:), pointer :: por => NULL()
!          real(p_k_rdoub), dimension(:),   pointer :: qor => NULL()
!!         particle identifier
!          integer(p_k_is), dimension(:),   pointer :: pid => NULL()
!
!        end type t_diag_species


        interface new
          module procedure new_diag_species
        end interface

        interface delete
          module procedure delete_diag_species
        end interface

        interface clear
          module procedure clear_diag_species
        end interface

        interface read_nml
          module procedure read_nml_diag_species
        end interface

        interface write_dbg
          module procedure write_dbg_diag_species
        end interface
  
        interface setup
          module procedure setup_diag_species
        end interface

        interface restart_write
          module procedure restart_write_diag_species
        end interface

        interface restart_read
          module procedure restart_read_diag_species
        end interface

        interface report
          module procedure report_diag_species_new
!          module procedure report_diag_species
!          module procedure report_diag_species_pE
!          module procedure report_diag_species_pE_xid
        end interface

        interface write_phasespace
          module procedure write_phasespace_to_hdf
        end interface

        interface write_spatial_ene
          module procedure write_spatial_ene_to_hdf         
        end interface

        interface write_bin_ene
          module procedure write_bin_ene_to_hdf         
        end interface


        interface write_3D_phasespace
          module procedure write_3D_phasespace_to_dx
        end interface

        interface write_1D_gamma_phasespace
          module procedure write_1D_gamma_phasespace_to_dx
        end interface

        interface write_particle_energies
          module procedure write_particle_energies_to_file
        end interface

        interface write_raw
          module procedure write_raw_new
          module procedure write_raw_
          module procedure write_raw_pE
          module procedure write_raw_pE_par_xid
        end interface
        
! OS 2.0
        interface if_report_charge
           module procedure if_report_charge_dspec
        end interface        

!       declare things that should be public
        public :: t_diag_species, t_node_conf, new, delete, clear
        public :: read_nml, write_dbg, setup
        public :: restart_write, restart_read
        public :: report
        public :: if_report_charge


       contains 


!---------------------------------------------------
        subroutine new_diag_species( this )
!---------------------------------------------------
!       allocate space to this pointer
!---------------------------------------------------

          implicit none

!       dummy variables

          type( t_diag_species ), pointer :: this
          integer :: ierr

!       local variables - none

!       executable statements

          if ( associated( this ) ) then
             call delete( this )
          end if

          allocate( this , stat=ierr)
     
        if (ierr .ne. 0) then
          write(file_id_msg,*) 'allocation problem in'
          write(file_id_msg,*) 'new_diag_species'
          write(file_id_msg,*) 'STOP'
      call file_flush(file_id_msg)
          call abort_program( -1_p_k_is )
        endif

        end subroutine new_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine delete_diag_species( this )
!---------------------------------------------------
!       release memory ocupied by this variable       
!---------------------------------------------------

          implicit none

!       dummy variables

          type( t_diag_species ), pointer :: this

!       local variables

          integer(p_k_is) :: result

!       executable statements

!          if ( associated( this ) ) then
             deallocate( this, stat=result )
!          endif
          nullify(this)

        end subroutine delete_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine clear_diag_species( this )
!---------------------------------------------------
!       release memory ocupied by this variable       
!---------------------------------------------------

          implicit none

!       dummy variables

          type( t_diag_species ), intent( inout ) :: this

!       local variables

          integer(p_k_is) :: result

!       executable statements

!          if ( associated( this%xor ) ) then
!!             dbg_tempo = size( this%xor )
!             deallocate(   this%xor, stat=result )
!!             if ( result .eq. 0 ) then
!!                dbg_dec_rdoub = dbg_dec_rdoub + dbg_tempo
!!             endif
!          endif
!          nullify(this%xor)
!
!          if ( associated( this%por ) ) then
!!             dbg_tempo = size( this%por )
!             deallocate(   this%por, stat=result )
!!             if ( result .eq. 0 ) then
!!                dbg_dec_rdoub = dbg_dec_rdoub + dbg_tempo
!!             endif
!          endif
!          nullify(this%por)
!
!          if ( associated( this%qor ) ) then
!!             dbg_tempo = size( this%qor )
!             deallocate(   this%qor, stat=result )
!!             if ( result .eq. 0 ) then
!!                dbg_dec_rdoub = dbg_dec_rdoub + dbg_tempo
!!             endif
!          endif
!          nullify(this%qor)
!
!          if ( associated( this%pid ) ) then
!!             dbg_tempo = size( this%pid )
!             deallocate(   this%pid, stat=result )
!!             if ( result .eq. 0 ) then
!!                dbg_dec_is = dbg_dec_is + dbg_tempo
!!             endif
!          endif
!          nullify(this%pid)

        end subroutine clear_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine read_nml_diag_species( this )
!---------------------------------------------------
!       read necessary information from inputdec
!---------------------------------------------------

          implicit none

!       dummy variables

          type( t_diag_species ), intent(out) :: this

!       local variables

          integer(p_k_is) :: ndump_fac_pha, ndump_fac_ene, ndump_fac_raw
! OS 2.0          
          integer(p_k_is) :: ndump_fac_charge
! OS 2.0

!         file_name prefix for data files
          character(80) :: file_name_pha, file_name_ene

!         physical range for phasespace data dumps
          real(p_k_rdoub), dimension(p_x_dim) :: ps_xmin, ps_xmax
          real(p_k_rdoub), dimension(p_p_dim) :: ps_pmin, ps_pmax

!         switch for autorange of ps_p
          logical, dimension(p_p_dim) :: if_ps_p_auto

!         physical range for 1D momenta phasespace data dumps
          real(p_k_rdoub) :: ps_gammamin, ps_gammamax

!         switch for log deposition of gamma
          logical :: if_ps_gamma_log

!         switch for autorange of ps_p
          logical :: if_ps_gamma_auto

!         resolutions for phasespace data dumps
          integer(p_k_is), dimension(p_x_dim) :: ps_nx, ps_nx_3D
          integer(p_k_is), dimension(p_p_dim) :: ps_np, ps_np_3D

!         resolution for 1D momenta phasespace data dump
          integer(p_k_is) :: ps_ngamma

!         deposition scheme for phase space calculations
          integer(p_k_is) :: dep_sch

!         switch to turn on and off species 1D momenta phase space diags
          logical :: if_gamma

!         switches to turn on and off species phase space diagnostics
          logical :: if_x2x1, if_x3x1, if_p1x1, if_p2x1, if_p3x1
          logical ::          if_x3x2, if_p1x2, if_p2x2, if_p3x2
          logical ::                   if_p1x3, if_p2x3, if_p3x3
          logical ::                            if_p2p1, if_p3p1
          logical ::                                     if_p3p2

!         parameters for raw data dump
          real(p_k_rdoub)  :: gamma_limit
          real(p_k_rdoub)  :: particle_fraction

!         switches for energy density diagnostics
          logical :: if_x2x1_ene
          logical :: if_x3x1_ene, if_x3x2_ene, if_x3x2x1_ene

!         switches and parameters for energy binned diagnostics
          logical :: if_x2x1_bin_ene
          logical :: if_x3x1_bin_ene, if_x3x2_bin_ene, if_x3x2x1_bin_ene
          integer(p_k_is) :: n_ene_bins
          real(p_k_rsing), dimension(max_n_ene_bins) :: ene_bins  
          
!         tags and tracks
          integer :: ndump_fac_tracks, niter_tracks
          character(len = p_max_filename_len) :: file_tags
  
          logical :: if_spectrum
          integer :: iter_spectrum
          

!          logical :: if_nv_ene, if_nv_gamma_ene
!          real(p_k_rdoub), dimension(2) :: ene_bin

          logical :: if_x3x2x1
          logical :: if_p1x2x1, if_p2x2x1, if_p3x2x1

          logical  ::  iforigin, iftestpar

          namelist /nl_diag_species/ &
     &             ndump_fac_pha, file_name_pha, &
     &             ndump_fac_ene, file_name_ene, &
     &             ndump_fac_raw, &
     &             ndump_fac_charge, &
     &             ps_xmin, ps_xmax, ps_pmin, ps_pmax, &
     &             if_ps_p_auto, &
     &             ps_gammamin, ps_gammamax, if_ps_gamma_log, if_ps_gamma_auto, &
     &             ps_nx, ps_nx_3D, ps_np, ps_np_3D, ps_ngamma, dep_sch, &
     &             if_gamma, &
     &             if_spectrum, &
     &             iter_spectrum, &
     &             if_x2x1, if_x3x1, if_p1x1, if_p2x1, if_p3x1, &
     &                      if_x3x2, if_p1x2, if_p2x2, if_p3x2, &
     &                               if_p1x3, if_p2x3, if_p3x3, &
     &                                        if_p2p1, if_p3p1, &
     &                                                 if_p3p2, &
     &             if_x3x2x1, &
     &             if_p1x2x1, if_p2x2x1, if_p3x2x1, &
     &             gamma_limit, particle_fraction, &
     &             iforigin, iftestpar, &
     &             if_x2x1_ene, &
     &             if_x3x1_ene, if_x3x2_ene, if_x3x2x1_ene, & 
     &             if_x2x1_bin_ene, &
     &             if_x3x1_bin_ene, if_x3x2_bin_ene, if_x3x2x1_bin_ene, &
     &             n_ene_bins, ene_bins, &
     &             ndump_fac_tracks, niter_tracks, file_tags
          
          integer :: ierr, nml_present,i

!       executable statements

!         write (*,*) 'Before read_nml_diag_species'

          ndump_fac_pha = 0
          file_name_pha = ' '

          ndump_fac_ene = 0
          file_name_ene = ' '

          ndump_fac_raw = 0
          
          ndump_fac_charge = 0

          ps_xmin       = 0.0_p_k_rdoub
          ps_xmax       = 0.0_p_k_rdoub
          ps_pmin       = 0.0_p_k_rdoub
          ps_pmax       = 0.0_p_k_rdoub
          if_ps_p_auto  = .false.
          ps_gammamin   = 0.0_p_k_rdoub
          ps_gammamax   = 0.0_p_k_rdoub
          if_ps_gamma_log  = .false.
          if_ps_gamma_auto  = .false.
          ps_nx         = 64
          ps_nx_3D      = -1
          ps_np         = 64
          ps_np_3D      = -1
          ps_ngamma     = 64
          dep_sch       = 1
          if_gamma      = .false.
          if_x2x1       = .false.
          if_x3x1       = .false.
          if_p1x1       = .false.
          if_p2x1       = .false.
          if_p3x1       = .false.
          if_x3x2       = .false.
          if_p1x2       = .false.
          if_p2x2       = .false.
          if_p3x2       = .false.
          if_p1x3       = .false.
          if_p2x3       = .false.
          if_p3x3       = .false.
          if_p2p1       = .false.
          if_p3p1       = .false.
          if_p3p2       = .false.

          if_x3x2x1     = .false.
          if_p1x2x1     = .false.
          if_p2x2x1     = .false.
          if_p3x2x1     = .false.

          iforigin      = .false.
          iftestpar     = .false.

          gamma_limit       = 0.0_p_k_rdoub
          particle_fraction = 1.0_p_k_rdoub
          
          if_spectrum=.false.
          iter_spectrum=0

          if_x2x1_ene   = .false.
          if_x3x1_ene   = .false.
          if_x3x2_ene   = .false.
          if_x3x2x1_ene   = .false.
 
          if_x2x1_bin_ene   = .false.
          if_x3x1_bin_ene   = .false.
          if_x3x2_bin_ene   = .false.
          if_x3x2x1_bin_ene   = .false.
 
          n_ene_bins = 0
          ene_bins = 0.0_p_k_rsing

          ndump_fac_tracks = 0
          niter_tracks = 1
          file_tags = ''
  
          
!          if_nv_ene     = .false. 
!          if_nv_gamma_ene = .false.
!          ene_bin       = 0.0_p_k_rdoub

          nml_present = nml_is_present(file_id_nml,"nl_diag_species")
          if (nml_present == 1) then
            read (file_id_nml,nl_diag_species, iostat = ierr)
            if (ierr /= 0) then
              write(*,*) ""
              write(*,*) "   Error reading diag_species parameters "
              write(*,*) "   aborting..."
              stop
            endif
          else 
            write(*,*) "   - no diagnostics specified"
          endif



          this%ndump_fac_pha = ndump_fac_pha
          this%file_name_pha = file_name_pha
          this%ndump_fac_ene = ndump_fac_ene
          this%file_name_ene = file_name_ene
          
          this%ndump_fac_raw = ndump_fac_raw
          this%ndump_fac_charge = ndump_fac_charge
          
          this%ps_xmin       = ps_xmin
          this%ps_xmax       = ps_xmax
          this%ps_pmin       = ps_pmin
          this%ps_pmax       = ps_pmax
          this%ps_pmin_r     = this%ps_pmin
          this%ps_pmax_r     = this%ps_pmax
          this%if_ps_p_auto  = if_ps_p_auto
          this%if_ps_gamma_auto  = if_ps_gamma_auto
          this%ps_gammamin   = ps_gammamin
          this%ps_gammamax   = ps_gammamax
          this%ps_gammamin_r = this%ps_gammamin
          this%ps_gammamax_r = this%ps_gammamax
          this%if_ps_gamma_log  = if_ps_gamma_log
          this%ps_nx         = ps_nx
          this%ps_np         = ps_np
          this%ps_ngamma     = ps_ngamma
          this%dep_sch       = dep_sch
          this%if_gamma      = if_gamma
          this%if_x2x1       = if_x2x1 
          this%if_x3x1       = if_x3x1    
          this%if_p1x1       = if_p1x1    
          this%if_p2x1       = if_p2x1    
          this%if_p3x1       = if_p3x1    
          this%if_x3x2       = if_x3x2    
          this%if_p1x2       = if_p1x2    
          this%if_p2x2       = if_p2x2    
          this%if_p3x2       = if_p3x2    
          this%if_p1x3       = if_p1x3   
          this%if_p2x3       = if_p2x3    
          this%if_p3x3       = if_p3x3    
          this%if_p2p1       = if_p2p1    
          this%if_p3p1       = if_p3p1   
          this%if_p3p2       = if_p3p2  
      
          this%if_x3x2x1     = if_x3x2x1
          this%if_p1x2x1     = if_p1x2x1
          this%if_p2x2x1     = if_p2x2x1
          this%if_p3x2x1     = if_p3x2x1

          this%raw_gamma_limit       = gamma_limit
          this%raw_fraction = particle_fraction
          
          this%if_spectrum=if_spectrum
          this%iter_spectrum=iter_spectrum

          this%iforigin      = iforigin 
          this%iftestpar     = iftestpar 

          this%if_x2x1_ene   = if_x2x1_ene
          this%if_x3x1_ene   = if_x3x1_ene
          this%if_x3x2_ene   = if_x3x2_ene
          this%if_x3x2x1_ene   = if_x3x2x1_ene

          this%if_x2x1_bin_ene   = if_x2x1_bin_ene
          this%if_x3x1_bin_ene   = if_x3x1_bin_ene
          this%if_x3x2_bin_ene   = if_x3x2_bin_ene
          this%if_x3x2x1_bin_ene   = if_x3x2x1_bin_ene
          
          this%n_ene_bins = n_ene_bins
          if (n_ene_bins > max_n_ene_bins) then
            write(*,*) "** n_ene_bins must be <= ", max_n_ene_bins
            stop
          endif
          
          do i=1, n_ene_bins
              this%ene_bins(i) = ene_bins(i)
              !write(*,*) "ene_bins(",i,") = ",ene_bins(i) 
          enddo

  ! store particle tracking data
  this%ndump_fac_tracks = ndump_fac_tracks
  if ( this%ndump_fac_tracks > 0 ) then
    this%tracks%niter = niter_tracks
    this%tracks%file_tags = file_tags
    if ( .not. file_exists(file_tags)) then
	   write(*,*) "(*error*) ndump_fac_tracks is set but tags file '"//trim(file_tags)// &
				  "' cannot be found."
	   stop
    endif
  endif

!          this%if_nv_ene     = if_nv_ene
!          this%if_nv_gamma_ene = if_nv_gamma_ene
!          this%ene_bin       = ene_bin
          
          if ( ps_nx_3D(1) .eq. -1 ) then
            this%ps_nx_3D    = ps_nx
          else
            this%ps_nx_3D    = ps_nx_3D
          endif
          if ( ps_np_3D(1) .eq. -1 ) then
            this%ps_np_3D    = ps_np
          else
            this%ps_np_3D    = ps_np_3D
          endif


!         write (*,*) 'After read_nml_diag_species'

        end subroutine read_nml_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_dbg_diag_species( this, file_number )
!---------------------------------------------------
!       write for debugging purposes
!---------------------------------------------------

          implicit none

!       dummy variables

          type( t_diag_species ),    intent( in )  :: this
          integer(p_k_is), optional, intent( in )  :: file_number 

!       local variables - none

!       executable statements

          if ( .not. ( present( file_number ) ) ) then

             write(*,*) 'diag_species :' 

             write(*,*) this%ndump_fac_pha
             write(*,*) this%file_name_pha
             write(*,*) this%ndump_fac_ene
             write(*,*) this%file_name_ene
             write(*,*) this%ps_xmin
             write(*,*) this%ps_xmax
             write(*,*) this%ps_pmin
             write(*,*) this%ps_pmax
             write(*,*) "ps_pmin_r: ", this%ps_pmin_r
             write(*,*) "ps_pmax_r: ", this%ps_pmax_r
             write(*,*) "if_ps_p_auto: ", this%if_ps_p_auto
             write(*,*) this%ps_gammamin
             write(*,*) this%ps_gammamax
             write(*,*) "if_ps_gamma_log: ", this%if_ps_gamma_log
             write(*,*) "if_ps_gamma_auto: ", this%if_ps_gamma_auto
             write(*,*) this%ps_nx
             write(*,*) "ps_nx_3D: ", this%ps_nx_3D
             write(*,*) this%ps_np
             write(*,*) "ps_np_3D: ", this%ps_np_3D
             write(*,*) this%ps_ngamma
             write(*,*) this%dep_sch
             write(*,*) this%if_gamma
             write(*,*) this%if_x2x1
             write(*,*) this%if_x3x1
             write(*,*) this%if_p1x1
             write(*,*) this%if_p2x1
             write(*,*) this%if_p3x1
             write(*,*) this%if_x3x2
             write(*,*) this%if_p1x2
             write(*,*) this%if_p2x2
             write(*,*) this%if_p3x2
             write(*,*) this%if_p1x3
             write(*,*) this%if_p2x3
             write(*,*) this%if_p3x3
             write(*,*) this%if_p2p1
             write(*,*) this%if_p3p1
             write(*,*) this%if_p3p2

             write(*,*) "if_x3x2x1: ", this%if_x3x2x1
             write(*,*) "if_p1x2x1: ", this%if_p1x2x1
             write(*,*) "if_p2x2x1: ", this%if_p2x2x1
             write(*,*) "if_p3x2x1: ", this%if_p3x2x1

             write(*,*) this%iforigin
             write(*,*) this%iftestpar


             write(*,*) 'xor, por, qor, pid - ', &
     &                  'write is not impemented'

             write(*,*) 'end of diag_species'
             write(*,*) ' '

          else

             write( file_number, *) 'diag_species :' 

             write( file_number, *) this%ndump_fac_pha
             write( file_number, *) this%file_name_pha
             write( file_number, *) this%ndump_fac_ene
             write( file_number, *) this%file_name_ene
             write( file_number, *) this%ps_xmin
             write( file_number, *) this%ps_xmax
             write( file_number, *) this%ps_pmin
             write( file_number, *) this%ps_pmax
             write( file_number, *) this%ps_pmin_r
             write( file_number, *) this%ps_pmax_r
             write( file_number, *) this%if_ps_p_auto
             write( file_number, *) this%ps_gammamin
             write( file_number, *) this%ps_gammamax
             write( file_number, *) this%if_ps_gamma_log
             write( file_number, *) this%if_ps_gamma_auto
             write( file_number, *) this%ps_nx
             write( file_number, *) this%ps_nx_3D
             write( file_number, *) this%ps_np
             write( file_number, *) this%ps_np_3D
             write( file_number, *) this%ps_ngamma
             write( file_number, *) this%dep_sch
             write( file_number, *) this%if_gamma
             write( file_number, *) this%if_x2x1
             write( file_number, *) this%if_x3x1
             write( file_number, *) this%if_p1x1
             write( file_number, *) this%if_p2x1
             write( file_number, *) this%if_p3x1
             write( file_number, *) this%if_x3x2
             write( file_number, *) this%if_p1x2
             write( file_number, *) this%if_p2x2
             write( file_number, *) this%if_p3x2
             write( file_number, *) this%if_p1x3
             write( file_number, *) this%if_p2x3
             write( file_number, *) this%if_p3x3
             write( file_number, *) this%if_p2p1
             write( file_number, *) this%if_p3p1
             write( file_number, *) this%if_p3p2

             write( file_number, *) this%if_x3x2x1
             write( file_number, *) this%if_p1x2x1
             write( file_number, *) this%if_p2x2x1
             write( file_number, *) this%if_p3x2x1

             write( file_number, *) this%iforigin
             write( file_number, *) this%iftestpar

             write( file_number, *) 'xor, por, qor, pid - ', &
     &                              'write not implemented'

             write( file_number, *) 'end of diag_species'
             write( file_number, *) ' '

          endif


        end subroutine write_dbg_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine setup_diag_species(this,ndump)
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type( t_diag_species ), intent( inout )  ::  this
          integer::ndump

!       local variables
          integer::ierr

!       executable statements

          call clear(this)
          
          ! correct gamma phasespace limits if necessary
          
          if (this%if_gamma) then
            
            if (this%ps_gammamin .lt. 1.0) this%ps_gammamin = 1.0
            
            if (this%ps_gammamax .le. this%ps_gammamin) then
              if (this%if_ps_gamma_log) then
                 this%ps_gammamax = 10.0 * this%ps_gammamin
              else 
                 this%ps_gammamax = 10.0 + this%ps_gammamin
              endif
            endif
            
          endif

	! set variable list for raw_math_expr
	if (this%raw_math_expr /= '') then
	   
	   select case (p_x_dim)
		 case (1)  
		   call setup(this%raw_func, trim(this%raw_math_expr), &
						 (/'x1', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
						 
		 case (2)
		   call setup(this%raw_func, trim(this%raw_math_expr), &
						 (/'x1', 'x2', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
						 
		 case (3)
		   call setup(this%raw_func, trim(this%raw_math_expr), &
						 (/'x1', 'x2', 'x3', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
					  
	   end select
	   
	   ! check if function compiled ok
	   if (ierr /= 0) then
             write(*,*)"Error compiling supplied function :"
             write(*,*) trim(this%raw_math_expr)
             call abort_program(-1)
	   endif
        endif	 
	
	! setup tracking diagnostic
!    call setup( this%tracks, ndump_fac*this%ndump_fac_tracks, restart )
    call setup( this%tracks, ndump*this%ndump_fac_tracks, .false. )



        end subroutine setup_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine restart_write_diag_species( this, file_id )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type( t_diag_species ), intent(in) :: this

          integer(p_k_is), intent(in) :: file_id

!       local variables

          logical :: if_xor, if_por, if_qor, if_pid

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'Now, before restart_write_diag_species'
!          call file_flush(file_id_msg)

!          if_xor = associated( this%xor )
!          if_por = associated( this%por )
!          if_qor = associated( this%qor )
!          if_pid = associated( this%pid )
          if_xor = .false.
          if_por = .false.
          if_qor = .false.
          if_pid = .false.

!      if_ps_gamma_auto was deliberately left out of the following
!      if_ps_gamma_log was deliberately left out of the following

          write(file_id) this%ndump_fac_pha, this%file_name_pha, &
     &                   this%ndump_fac_ene, this%file_name_ene, &
     &                   this%ps_xmin, this%ps_xmax, &
     &                   this%ps_pmin, this%ps_pmax, &
     &                   this%ps_pmin_r, this%ps_pmax_r, &
     &                   this%if_ps_p_auto, &
     &                   this%ps_gammamin, this%ps_gammamax, &
     &                   this%ps_nx, this%ps_nx_3D, &
     &                   this%ps_np, this%ps_np_3D, &
     &                   this%ps_ngamma, &
     &                   this%dep_sch, &
     &                   this%if_gamma, &
     &                   this%if_x2x1, this%if_x3x1, &
     &                   this%if_p1x1, this%if_p2x1, this%if_p3x1, &
     &                   this%if_x3x2, this%if_p1x2, &
     &                   this%if_p2x2, this%if_p3x2, &
     &                   this%if_p1x3, this%if_p2x3, this%if_p3x3, &
     &                   this%if_p2p1, this%if_p3p1, this%if_p3p2, &
     &                   this%if_x3x2x1, &
     &                   this%if_p1x2x1, &
     &                   this%if_p2x2x1, &
     &                   this%if_p3x2x1, &
     &                   this%iforigin, this%iftestpar, &
     &                   if_xor, if_por, if_qor, if_pid

!          if ( if_xor ) then
!             write(file_id) size(this%xor,1), size(this%xor,2)
!             write(file_id) this%xor
!          endif
!          if ( if_por ) then
!             write(file_id) size(this%por,1), size(this%por,2)
!             write(file_id) this%por           
!          endif
!          if ( if_qor ) then
!             write(file_id) size(this%qor,1)     
!             write(file_id) this%qor           
!          endif
!          if ( if_pid ) then
!             write(file_id) size(this%pid,1)      
!             write(file_id) this%pid           
!          endif
  ! only tracks have restart data
!  call restart_write( this%tracks, file_id )

        end subroutine restart_write_diag_species
!---------------------------------------------------


!---------------------------------------------------
        subroutine restart_read_diag_species( this, file_id )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type( t_diag_species ), intent(inout) :: this

          integer(p_k_is), intent(in) :: file_id

!       local variables

          integer :: s1, s2, ierr
          logical :: if_xor, if_por, if_qor, if_pid

!       executable statements

!          write(file_id_msg, *) ' '
          write(file_id_msg, *) 'Now, before restart_read_diag_species'
!          call file_flush(file_id_msg)

          call clear(this)

!      if_ps_gamma_auto was deliberately left out of the following
!      if_ps_gamma_log was deliberately left out of the following

          read(file_id) this%ndump_fac_pha, this%file_name_pha, &
     &                  this%ndump_fac_ene, this%file_name_ene, &
     &                  this%ps_xmin, this%ps_xmax, &
     &                  this%ps_pmin, this%ps_pmax, &
     &                  this%ps_pmin_r, this%ps_pmax_r, &
     &                  this%if_ps_p_auto, &
     &                  this%ps_gammamin, this%ps_gammamax, &
     &                  this%ps_nx, this%ps_nx_3D, &
     &                  this%ps_np, this%ps_np_3D, &
     &                  this%ps_ngamma, &
     &                  this%dep_sch, &
     &                  this%if_gamma, &
     &                  this%if_x2x1, this%if_x3x1, &
     &                  this%if_p1x1, this%if_p2x1, this%if_p3x1, &
     &                  this%if_x3x2, this%if_p1x2, &
     &                  this%if_p2x2, this%if_p3x2, &
     &                  this%if_p1x3, this%if_p2x3, this%if_p3x3, &
     &                  this%if_p2p1, this%if_p3p1, this%if_p3p2, &
     &                  this%if_x3x2x1, &
     &                  this%if_p1x2x1, &
     &                  this%if_p2x2x1, &
     &                  this%if_p3x2x1, &
     &                  this%iforigin, this%iftestpar, &
     &                  if_xor, if_por, if_qor, if_pid

!          if ( if_xor ) then
!             read(file_id) s1, s2
!             allocate( this%xor(s1,s2) , stat=ierr)
!     
!        if (ierr .ne. 0) then
!          write(file_id_msg,*) 'allocation problem in'
!          write(file_id_msg,*) 'restart_read_diag_species'
!          write(file_id_msg,*) 'STOP'
!          call file_flush(file_id_msg)
!          call abort_program( -1_p_k_is )
!        endif
!             read(file_id) this%xor
!          endif
!          if ( if_por ) then  
!             read(file_id) s1, s2
!             allocate( this%por(s1,s2) , stat=ierr)
!     
!        if (ierr .ne. 0) then
!          write(file_id_msg,*) 'allocation problem in'
!          write(file_id_msg,*) 'restart_read_diag_species'
!          write(file_id_msg,*) 'STOP'
!          call file_flush(file_id_msg)
!          call abort_program( -1_p_k_is )
!        endif       
!             read(file_id) this%por
!          endif
!          if ( if_qor ) then  
!             read(file_id) s1
!             allocate( this%qor(s1) , stat=ierr)
!     
!        if (ierr .ne. 0) then
!          write(file_id_msg,*) 'allocation problem in'
!          write(file_id_msg,*) 'restart_read_diag_species'
!          write(file_id_msg,*) 'STOP'
!          call file_flush(file_id_msg)
!          call abort_program( -1_p_k_is )
!        endif     
!             read(file_id) this%qor
!          endif
!          if ( if_pid ) then  
!             read(file_id) s1
!             allocate( this%pid(s1) , stat=ierr)
!     
!        if (ierr .ne. 0) then
!          write(file_id_msg,*) 'allocation problem in'
!          write(file_id_msg,*) 'restart_read_diag_species'
!          write(file_id_msg,*) 'STOP'
!          call file_flush(file_id_msg)
!          call abort_program( -1_p_k_is )
!        endif       
!             read(file_id) this%pid
!          endif

        end subroutine restart_read_diag_species
!---------------------------------------------------
!---------------------------------------------------
        subroutine report_diag_species_new( this,spec, g_space,no_co, &
     &                                  n, ndump, t,dt, dx, coordinates )
!---------------------------------------------------
!       report on physical field - diagnostic
!---------------------------------------------------

!         <use-statements for this subroutine>
          implicit none

!       passed values
!         the diagnostics class of the current species
          type( t_diag_species ), intent(inout) :: this
          
          type( t_species ),             intent(inout) :: spec
          
          type( t_space ),               intent(in) :: g_space

!         the system node configuration
          type( t_node_conf ),             intent(in) :: no_co
!         the iteraction number and ndump number
          integer(p_k_is),                 intent(in) :: n, ndump
!         the simulation time
          real(p_k_rdoub),                 intent(in) :: t
          real(p_k_rdoub),                 intent(in) :: dt
!         the dimensions of the simulation cell
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
!         gives information on which kind of coordinates system is used
          integer(p_k_is),                 intent(in) :: coordinates

!       local variables
!         position and momenta of the particles
          real(p_k_rdoub), dimension(:,:),pointer :: x, p
!         rqm of the individual particles
          real(p_k_rdoub), dimension(:),pointer :: q
!         number of particles of the current species
          integer(p_k_is):: num_par

!         rqm of the species
          real(p_k_rdoub):: rqm
!         time centered kinetic energy of the species
          real(p_k_rdoub):: dKE
!         species id
          integer(p_k_is):: sp_id

!         the physical dimensions of the simulation box
          real(p_k_rdoub), dimension(p_x_dim) :: lxmin, lxmax


          integer(p_k_is) :: n_del, n_aux, x_dim

!
!       if x_or_p eq 1 then the coord is x; if it is 2 then the coord is p
!       xp_dim is the coordinate 123 : x123 or p123

          integer(p_k_is), dimension(3) :: x_or_p, xp_dim

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before report_diag_species'
!          call file_flush(file_id_msg)


           x=>spec%x;
           p=>spec%p;
!           pE=>spec%pE;
           q=>spec%q;
!           par_xid=>spec%par_xid;
           num_par=spec%num_par;
           rqm=spec%rqm;
           dKE=spec%dKE;
           sp_id=spec%sp_id;
           lxmin=xmin(g_space);
           lxmax=xmax(g_space);



!         reports are turned off if ndump_fac is 0
          if ((this%ndump_fac_pha .ne. 0) .and. (ndump .ne. 0)) then

             !write(*,*) "In report_diag_species, ndump_fac_pha = ", this%ndump_fac_pha
             !write(*,*) "In report_diag_species, ndump = ", ndump
             !write(*,*) "lbound, ubound x", lbound(x,1), lbound(x,2), ubound(x,1), ubound(x,2) 
             !write(*,*) "lbound, ubound p", lbound(p,1), lbound(p,2), ubound(p,1), ubound(p,2) 

             n_del = ndump * this%ndump_fac_pha
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_pha
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

!               redefine n_aux to measure time in units of ndump
                n_aux = n_aux*this%ndump_fac_pha
                
!               autorange
                call acquire_range( this, num_par, p, no_co )
   
!               2D phase space diagnostics

                if ( this%if_x2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                
               ! DEBUG
               ! return
                
                if (( this%if_x3x1 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_x3x2 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif 
                if ( this%if_p2x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p1x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p2x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p3x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p2 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,lxmin,lxmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
    
!               3D phase space diagnostics
               
                if ((this%if_x3x2x1) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 1 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p1x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 1
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p2x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 2
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if ((this%if_p3x2x1) .and. (p_p_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif


!               Energy density diagnostics

                if ( this%if_x2x1_ene ) then                ! 2D energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2
                
                   call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                     this, no_co,lxmin,lxmax,n,t, &
     &                                     n_aux, xp_dim, &
     &                                     dx, coordinates )
                endif

                if (p_x_dim .eq. 3) then                    

                   if (this%if_x3x1_ene) then               ! 2D energy density (x3x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif

                   if (this%if_x3x2_ene) then               ! 2D energy density (x3x2)
                     xp_dim(1) = 2
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,lxmin,lxmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif
                   

                   if (this%if_x3x2x1_ene) then             ! 3D energy density (x3x2x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 2
                     xp_dim(3) = 3
                     
                     ! not implemented yet
                     ! call write_3D_spatial_ene
                     write(*,*) "write_3D_spatial_ene not implemented yet!"
                   endif

                endif  ! if (p_x_dim .eq. 3) 


!               Energy Binned diagnostics

                if ( this%if_x2x1_bin_ene ) then          ! 2D binned energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2

                   call write_bin_ene( x, p, q, num_par, sp_id, &
     &                                 this, no_co,lxmin,lxmax,n,t, &
     &                                 n_aux, xp_dim, &
     &                                 dx, coordinates )
                   
                endif


!                if ( this%if_nv_ene ) then
!
!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,lxmin,lxmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = p_nv_ene )
!                endif

!                if ( this%if_nv_gamma_ene ) then

!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,lxmin,lxmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = !p_nv_gamma_ene )
!                endif
                



!               gamma phase space diagnostics
               
               if (this%if_gamma) then

!                  autorange
                   call acquire_range_gamma( this, num_par, p, no_co )

                   call write_1D_gamma_phasespace( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
               endif

             endif

          endif
          if (this%if_spectrum) then
!            print*,sp_id,'report_diag_species',this%iter_spectrum,n
            if((this%iter_spectrum>0).and.(this%iter_spectrum*(n/this%iter_spectrum)==n)) then
              call write_spectrum( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
            endif
          endif
     
!         reports on particle energies are turned off if ndump_fac_ene is 0
          if ( this%ndump_fac_ene .ne. 0 )  then

             n_del = ndump * this%ndump_fac_ene
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_ene
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

                call write_particle_energies( p, q, rqm,dKE, num_par, sp_id, &
     &                                        this, no_co,n,t, dx )
             endif
          endif

!         reports on raw particle data are turned off if ndump_fac_raw is 0
          if ( this%ndump_fac_raw .ne. 0 )  then

             n_del = ndump * this%ndump_fac_raw
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_raw
             if ( ( n - n_aux*n_del ) .eq. 0 ) then
                call write_raw_new( spec, no_co, g_space, n, t, n_aux, coordinates )

             endif
          endif
  ! Tracks diagnostics
  if (spec%diag%ndump_fac_tracks > 0) then
    if ( mod( n, spec%diag%tracks%niter ) == 0 ) then
       call add_track_data( spec%diag%tracks, spec, n, t )
    endif
  endif
  
  if ( test_if_report( n, ndump, spec%diag%ndump_fac_tracks )) then
     if ( n == 0 ) then 
        ! create the diagnostics file 
        call create_file( spec%diag%tracks, spec%name, ndump, dt, &
                          x_bnd(g_space), periodic(no_co), if_move( g_space ) )
     else 
        ! write the data
        call write_tracks( spec%diag%tracks, no_co )
     endif
  endif
          

        end subroutine report_diag_species_new
!---------------------------------------------------


!---------------------------------------------------
        subroutine report_diag_species( x, p, q, num_par, rqm,dKE, sp_id, &
     &                                  this, no_co, xmin, xmax, &
     &                                  n, ndump, t, dx, coordinates )
!---------------------------------------------------
!       report on physical field - diagnostic
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

!         position and momenta of the particles
          real(p_k_rdoub), dimension(:,:), intent(inout) :: x, p
!         rqm of the individual particles
          real(p_k_rdoub), dimension(:),   intent(inout) :: q
!         number of particles of the current species
          integer(p_k_is),                 intent(in) :: num_par

!         rqm of the species
          real(p_k_rdoub)              ,   intent(in) :: rqm
!         time centered kinetic energy of the species
          real(p_k_rdoub)              ,   intent(in) :: dKE
!         species id
          integer(p_k_is),                 intent(in) :: sp_id
!         the diagnostics class of the current species
          type( t_diag_species ),          intent(inout) :: this
!         the system node configuration
          type( t_node_conf ),             intent(in) :: no_co
!         the physical dimensions of the simulation box
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
!         the iteraction number and ndump number
          integer(p_k_is),                 intent(in) :: n, ndump
!         the simulation time
          real(p_k_rdoub),                 intent(in) :: t
!         the dimensions of the simulation cell
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
!         gives information on which kind of coordinates system is used
          integer(p_k_is),                 intent(in) :: coordinates

!       local variables

          integer(p_k_is) :: n_del, n_aux, x_dim

!
!       if x_or_p eq 1 then the coord is x; if it is 2 then the coord is p
!       xp_dim is the coordinate 123 : x123 or p123

          integer(p_k_is), dimension(3) :: x_or_p, xp_dim

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before report_diag_species'
!          call file_flush(file_id_msg)

!         reports are turned off if ndump_fac is 0
          if ((this%ndump_fac_pha .ne. 0) .and. (ndump .ne. 0)) then

             !write(*,*) "In report_diag_species, ndump_fac_pha = ", this%ndump_fac_pha
             !write(*,*) "In report_diag_species, ndump = ", ndump
             !write(*,*) "lbound, ubound x", lbound(x,1), lbound(x,2), ubound(x,1), ubound(x,2) 
             !write(*,*) "lbound, ubound p", lbound(p,1), lbound(p,2), ubound(p,1), ubound(p,2) 

             n_del = ndump * this%ndump_fac_pha
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_pha
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

!               redefine n_aux to measure time in units of ndump
                n_aux = n_aux*this%ndump_fac_pha
                
!               autorange
                call acquire_range( this, num_par, p, no_co )
   
!               2D phase space diagnostics

                if ( this%if_x2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                
               ! DEBUG
               ! return
                
                if (( this%if_x3x1 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_x3x2 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif 
                if ( this%if_p2x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p1x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p2x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p3x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p2 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
    
!               3D phase space diagnostics
               
                if ((this%if_x3x2x1) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 1 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p1x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 1
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p2x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 2
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if ((this%if_p3x2x1) .and. (p_p_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif


!               Energy density diagnostics

                if ( this%if_x2x1_ene ) then                ! 2D energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2
                
                   call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                     this, no_co,xmin,xmax,n,t, &
     &                                     n_aux, xp_dim, &
     &                                     dx, coordinates )
                endif

                if (p_x_dim .eq. 3) then                    

                   if (this%if_x3x1_ene) then               ! 2D energy density (x3x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif

                   if (this%if_x3x2_ene) then               ! 2D energy density (x3x2)
                     xp_dim(1) = 2
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif
                   

                   if (this%if_x3x2x1_ene) then             ! 3D energy density (x3x2x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 2
                     xp_dim(3) = 3
                     
                     ! not implemented yet
                     ! call write_3D_spatial_ene
                     write(*,*) "write_3D_spatial_ene not implemented yet!"
                   endif

                endif  ! if (p_x_dim .eq. 3) 


!               Energy Binned diagnostics

                if ( this%if_x2x1_bin_ene ) then          ! 2D binned energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2

                   call write_bin_ene( x, p, q, num_par, sp_id, &
     &                                 this, no_co,xmin,xmax,n,t, &
     &                                 n_aux, xp_dim, &
     &                                 dx, coordinates )
                   
                endif


!                if ( this%if_nv_ene ) then
!
!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = p_nv_ene )
!                endif

!                if ( this%if_nv_gamma_ene ) then

!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = !p_nv_gamma_ene )
!                endif
                



!               gamma phase space diagnostics
               
               if (this%if_gamma) then

!                  autorange
                   call acquire_range_gamma( this, num_par, p, no_co )

                   call write_1D_gamma_phasespace( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
               endif

             endif

          endif
          if (this%if_spectrum) then
!            print*,sp_id,'report_diag_species',this%iter_spectrum,n
            if((this%iter_spectrum>0).and.(this%iter_spectrum*(n/this%iter_spectrum)==n)) then
              call write_spectrum( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
            endif
          endif
     
!         reports on particle energies are turned off if ndump_fac_ene is 0
          if ( this%ndump_fac_ene .ne. 0 )  then

             n_del = ndump * this%ndump_fac_ene
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_ene
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

                call write_particle_energies( p, q, rqm,dKE, num_par, sp_id, &
     &                                        this, no_co,n,t, dx )
             endif
          endif

!         reports on raw particle data are turned off if ndump_fac_raw is 0
          if ( this%ndump_fac_raw .ne. 0 )  then

             n_del = ndump * this%ndump_fac_raw
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_raw
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

                call  write_raw( x, p, q, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
             endif
          endif
          

        end subroutine report_diag_species
!---------------------------------------------------

! p.E
!---------------------------------------------------
        subroutine report_diag_species_pE( x, p, pE, q, num_par, rqm,dKE, sp_id, &
     &                                  this, no_co, xmin, xmax, &
     &                                  n, ndump, t, dx, coordinates )
!---------------------------------------------------
!       report on physical field - diagnostic
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

!         position and momenta of the particles
          real(p_k_rdoub), dimension(:,:), intent(inout) :: x, p
! p.E
          real(p_k_rdoub), dimension(:,:), intent(inout) :: pE
! p.E
!         rqm of the individual particles
          real(p_k_rdoub), dimension(:),   intent(inout) :: q
!         number of particles of the current species
          integer(p_k_is),                 intent(in) :: num_par

!         rqm of the species
          real(p_k_rdoub)              ,   intent(in) :: rqm
!         time centered kinetic energy of the species
          real(p_k_rdoub)              ,   intent(in) :: dKE
!         species id
          integer(p_k_is),                 intent(in) :: sp_id
!         the diagnostics class of the current species
          type( t_diag_species ),          intent(inout) :: this
!         the system node configuration
          type( t_node_conf ),             intent(in) :: no_co
!         the physical dimensions of the simulation box
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
!         the iteraction number and ndump number
          integer(p_k_is),                 intent(in) :: n, ndump
!         the simulation time
          real(p_k_rdoub),                 intent(in) :: t
!         the dimensions of the simulation cell
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
!         gives information on which kind of coordinates system is used
          integer(p_k_is),                 intent(in) :: coordinates

!       local variables

          integer(p_k_is) :: n_del, n_aux, x_dim

!
!       if x_or_p eq 1 then the coord is x; if it is 2 then the coord is p
!       xp_dim is the coordinate 123 : x123 or p123

          integer(p_k_is), dimension(3) :: x_or_p, xp_dim

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before report_diag_species'
!          call file_flush(file_id_msg)

!         reports are turned off if ndump_fac is 0
          if ((this%ndump_fac_pha .ne. 0) .and. (ndump .ne. 0)) then

             !write(*,*) "In report_diag_species, ndump_fac_pha = ", this%ndump_fac_pha
             !write(*,*) "In report_diag_species, ndump = ", ndump
             !write(*,*) "lbound, ubound x", lbound(x,1), lbound(x,2), ubound(x,1), ubound(x,2) 
             !write(*,*) "lbound, ubound p", lbound(p,1), lbound(p,2), ubound(p,1), ubound(p,2) 

             n_del = ndump * this%ndump_fac_pha
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_pha
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

!               redefine n_aux to measure time in units of ndump
                n_aux = n_aux*this%ndump_fac_pha
                
!               autorange
                call acquire_range( this, num_par, p, no_co )
   
!               2D phase space diagnostics

                if ( this%if_x2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                
               ! DEBUG
               ! return
                
                if (( this%if_x3x1 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_x3x2 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif 
                if ( this%if_p2x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p1x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p2x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p3x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p2 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
    
!               3D phase space diagnostics
               
                if ((this%if_x3x2x1) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 1 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p1x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 1
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p2x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 2
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if ((this%if_p3x2x1) .and. (p_p_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif


!               Energy density diagnostics

                if ( this%if_x2x1_ene ) then                ! 2D energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2
                
                   call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                     this, no_co,xmin,xmax,n,t, &
     &                                     n_aux, xp_dim, &
     &                                     dx, coordinates )
                endif

                if (p_x_dim .eq. 3) then                    

                   if (this%if_x3x1_ene) then               ! 2D energy density (x3x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif

                   if (this%if_x3x2_ene) then               ! 2D energy density (x3x2)
                     xp_dim(1) = 2
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif
                   

                   if (this%if_x3x2x1_ene) then             ! 3D energy density (x3x2x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 2
                     xp_dim(3) = 3
                     
                     ! not implemented yet
                     ! call write_3D_spatial_ene
                     write(*,*) "write_3D_spatial_ene not implemented yet!"
                   endif

                endif  ! if (p_x_dim .eq. 3) 


!               Energy Binned diagnostics

                if ( this%if_x2x1_bin_ene ) then          ! 2D binned energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2

                   call write_bin_ene( x, p, q, num_par, sp_id, &
     &                                 this, no_co,xmin,xmax,n,t, &
     &                                 n_aux, xp_dim, &
     &                                 dx, coordinates )
                   
                endif


!                if ( this%if_nv_ene ) then
!
!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = p_nv_ene )
!                endif

!                if ( this%if_nv_gamma_ene ) then

!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = !p_nv_gamma_ene )
!                endif
                



!               gamma phase space diagnostics
               
               if (this%if_gamma) then

!                  autorange
                   call acquire_range_gamma( this, num_par, p, no_co )

                   call write_1D_gamma_phasespace( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
               endif
            endif

          endif
          if (this%if_spectrum) then
!            print*,sp_id,'report_diag_species',this%iter_spectrum,n
            if((this%iter_spectrum>0).and.(this%iter_spectrum*(n/this%iter_spectrum)==n)) then
              call write_spectrum( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
            endif
          endif
     
!         reports on particle energies are turned off if ndump_fac_ene is 0
          if ( this%ndump_fac_ene .ne. 0 )  then

             n_del = ndump * this%ndump_fac_ene
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_ene
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

                call write_particle_energies( p, q, rqm,dKE, num_par, sp_id, &
     &                                        this, no_co,n,t, dx )
             endif
          endif

!         reports on raw particle data are turned off if ndump_fac_raw is 0
          if ( this%ndump_fac_raw .ne. 0 )  then

             n_del = ndump * this%ndump_fac_raw
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_raw
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

! p.E
                call  write_raw( x, p, pE, q, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
! p.E
             endif
          endif
          

        end subroutine report_diag_species_pE
!---------------------------------------------------

!---------------------------------------------------
        subroutine report_diag_species_pE_xid( x, p, pE, q, par_xid, num_par, &
     &                                  rqm,dKE, sp_id, &
     &                                  this, no_co, xmin, xmax, &
     &                                  n, ndump, t, dx, coordinates )
!---------------------------------------------------
!       report on physical field - diagnostic
!---------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

!         position and momenta of the particles
          real(p_k_rdoub), dimension(:,:), intent(inout) :: x, p
! p.E
          real(p_k_rdoub), dimension(:,:), intent(inout) :: pE, par_xid
! p.E
!         rqm of the individual particles
          real(p_k_rdoub), dimension(:),   intent(inout) :: q
!         number of particles of the current species
          integer(p_k_is),                 intent(in) :: num_par

!         rqm of the species
          real(p_k_rdoub)              ,   intent(in) :: rqm
!         time centered kinetic energy of the species
          real(p_k_rdoub)              ,   intent(in) :: dKE
!         species id
          integer(p_k_is),                 intent(in) :: sp_id
!         the diagnostics class of the current species
          type( t_diag_species ),          intent(inout) :: this
!         the system node configuration
          type( t_node_conf ),             intent(in) :: no_co
!         the physical dimensions of the simulation box
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
!         the iteraction number and ndump number
          integer(p_k_is),                 intent(in) :: n, ndump
!         the simulation time
          real(p_k_rdoub),                 intent(in) :: t
!         the dimensions of the simulation cell
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
!         gives information on which kind of coordinates system is used
          integer(p_k_is),                 intent(in) :: coordinates

!       local variables

          integer(p_k_is) :: n_del, n_aux, x_dim

!
!       if x_or_p eq 1 then the coord is x; if it is 2 then the coord is p
!       xp_dim is the coordinate 123 : x123 or p123

          integer(p_k_is), dimension(3) :: x_or_p, xp_dim

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before report_diag_species'
!          call file_flush(file_id_msg)

!         reports are turned off if ndump_fac is 0
          if ((this%ndump_fac_pha .ne. 0) .and. (ndump .ne. 0)) then

             !write(*,*) "In report_diag_species, ndump_fac_pha = ", this%ndump_fac_pha
             !write(*,*) "In report_diag_species, ndump = ", ndump
             !write(*,*) "lbound, ubound x", lbound(x,1), lbound(x,2), ubound(x,1), ubound(x,2) 
             !write(*,*) "lbound, ubound p", lbound(p,1), lbound(p,2), ubound(p,1), ubound(p,2) 

             n_del = ndump * this%ndump_fac_pha
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_pha
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

!               redefine n_aux to measure time in units of ndump
                n_aux = n_aux*this%ndump_fac_pha
                
!               autorange
                call acquire_range( this, num_par, p, no_co )
   
!               2D phase space diagnostics

                if ( this%if_x2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                
               ! DEBUG
               ! return
                
                if (( this%if_x3x1 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x1 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_x3x2 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 1 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p1x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif 
                if ( this%if_p2x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3x2 ) then
                   x_or_p(1)= 1 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p1x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 1
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p2x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if (( this%if_p3x3 ) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 3
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p2p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 2
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p1 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 1
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
                if ( this%if_p3p2 ) then
                   x_or_p(1)= 2 ; xp_dim(1) = 2
                   x_or_p(2)= 2 ; xp_dim(2) = 3
                   call write_phasespace( x, p, q, num_par, sp_id, &
     &                                    this, no_co,xmin,xmax,n,t, &
     &                                    n_aux,x_or_p,xp_dim, &
     &                                    dx, coordinates )
                endif
    
!               3D phase space diagnostics
               
                if ((this%if_x3x2x1) .and. (p_x_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 1 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p1x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 1
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if (this%if_p2x2x1) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 2
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif

                if ((this%if_p3x2x1) .and. (p_p_dim .eq. 3)) then
                   x_or_p(1)= 1 ; xp_dim(1) = 1
                   x_or_p(2)= 1 ; xp_dim(2) = 2    
                   x_or_p(3)= 2 ; xp_dim(3) = 3
                   call write_3D_phasespace( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux,x_or_p,xp_dim, &
     &                                       dx, coordinates )
                endif


!               Energy density diagnostics

                if ( this%if_x2x1_ene ) then                ! 2D energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2
                
                   call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                     this, no_co,xmin,xmax,n,t, &
     &                                     n_aux, xp_dim, &
     &                                     dx, coordinates )
                endif

                if (p_x_dim .eq. 3) then                    

                   if (this%if_x3x1_ene) then               ! 2D energy density (x3x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif

                   if (this%if_x3x2_ene) then               ! 2D energy density (x3x2)
                     xp_dim(1) = 2
                     xp_dim(2) = 3
                     
                     call write_spatial_ene( x, p, q, num_par, sp_id, &
     &                                       this, no_co,xmin,xmax,n,t, &
     &                                       n_aux, xp_dim, &
     &                                       dx, coordinates )
                   endif
                   

                   if (this%if_x3x2x1_ene) then             ! 3D energy density (x3x2x1)
                     xp_dim(1) = 1
                     xp_dim(2) = 2
                     xp_dim(3) = 3
                     
                     ! not implemented yet
                     ! call write_3D_spatial_ene
                     write(*,*) "write_3D_spatial_ene not implemented yet!"
                   endif

                endif  ! if (p_x_dim .eq. 3) 


!               Energy Binned diagnostics

                if ( this%if_x2x1_bin_ene ) then          ! 2D binned energy density (x2x1)
                   xp_dim(1) = 1
                   xp_dim(2) = 2

                   call write_bin_ene( x, p, q, num_par, sp_id, &
     &                                 this, no_co,xmin,xmax,n,t, &
     &                                 n_aux, xp_dim, &
     &                                 dx, coordinates )
                   
                endif


!                if ( this%if_nv_ene ) then
!
!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = p_nv_ene )
!                endif

!                if ( this%if_nv_gamma_ene ) then

!                   x_or_p(1)= 1 ; xp_dim(1) = 1
!                   x_or_p(2)= 1 ; xp_dim(2) = 2
!!                   call write_phasespace( x, p, q, num_par, sp_id, &
!!     &                                    this, no_co,xmin,xmax,n,t, &
!!     &                                    n_aux,x_or_p,xp_dim, &
!!     &                                    dx, coordinates, bin = this%ene_bin, type = !p_nv_gamma_ene )
!                endif
                



!               gamma phase space diagnostics
               
               if (this%if_gamma) then

!                  autorange
                   call acquire_range_gamma( this, num_par, p, no_co )

                   call write_1D_gamma_phasespace( p, q, num_par, sp_id, &
                                             this, no_co, n, t, n_aux, &
                                             dx, coordinates ) 
               endif

             endif

          endif
            ! very bad can cause exceptions
!          if ((this%if_spectrum).and.(this%iter_spectrum>0).and.(this%iter_spectrum*(n/this%iter_spectrum)==n)) then
!            call write_spectrum( p, q, num_par, sp_id, &
!                                        this, no_co, n, t, n_aux, &
!                                        dx, coordinates ) 
!          endif
          if (this%if_spectrum) then
            if((this%iter_spectrum>0).and.(this%iter_spectrum*(n/this%iter_spectrum)==n)) then
            call write_spectrum( p, q, num_par, sp_id, &
                                        this, no_co, n, t, n_aux, &
                                        dx, coordinates ) 
            endif
          endif

      
!         reports on particle energies are turned off if ndump_fac_ene is 0
          if ( this%ndump_fac_ene .ne. 0 )  then

             n_del = ndump * this%ndump_fac_ene
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_ene
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

                call write_particle_energies( p, q, rqm,dKE, num_par, sp_id, &
     &                                        this, no_co,n,t, dx )
             endif
          endif

!         reports on raw particle data are turned off if ndump_fac_raw is 0
          if ( this%ndump_fac_raw .ne. 0 )  then

             n_del = ndump * this%ndump_fac_raw
             n_aux = n / n_del

!            dump only when n is multiples of ndump*this%ndump_fac_raw
             if ( ( n - n_aux*n_del ) .eq. 0 ) then

! p.E
                call  write_raw( x, p, pE, q, par_xid, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
! p.E
             endif
          endif
          

        end subroutine report_diag_species_pE_xid
!---------------------------------------------------
! p.E


!---------------------------------------------------
        subroutine write_phasespace_to_hdf( x, p, q, num_par, sp_id, &
     &                                     this, no_co, xmin, xmax, &
     &                                     n, t, n_aux, &
     &                                     x_or_p, xp_dim, &
     &                                     dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

          real(p_k_rdoub), dimension(:,:), intent(in), target :: x, p
          real(p_k_rdoub), dimension(:),   intent(in) :: q
          integer(p_k_is),                 intent(in) :: num_par
          integer(p_k_is),                 intent(in) :: sp_id
          type( t_diag_species ),          intent(in) :: this
          type( t_node_conf ),             intent(in) :: no_co
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
          integer(p_k_is),                 intent(in) :: n
          real(p_k_rdoub),                 intent(in) :: t
          integer(p_k_is),                 intent(in) :: n_aux
          integer(p_k_is), dimension(:),   intent(in) :: x_or_p
          integer(p_k_is), dimension(:),   intent(in) :: xp_dim
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
          integer(p_k_is),                 intent(in) :: coordinates

!         meaning of some variable values:
!         x_or_p(i) = 1 :  use x-data for axis i of phasespace plot
!         x_or_p(i) = 2 :  use p-data for axis i of phasespace plot
!         xp_dim(i) = k :  use direction k (of x or p data) 
!                          for axis i of phasespace plot

!       local variables


!       pointer to x or p
         real(p_k_rdoub), dimension(:,:), pointer :: xp1 => NULL()
         real(p_k_rdoub), dimension(:,:), pointer :: xp2 => NULL()


!       normalised coordinates of grid point
          real(p_k_rdoub) :: nci1,nci2
!       distance to nearest grid point in each dimension (normalised)
          real(p_k_rdoub) :: xpdgp1, xpdgp2

          real(p_k_rdoub), dimension(2) :: ps_min, ps_max, ps_dx, ps_rdx
          integer(p_k_is), dimension(2) :: ps_n
          character(20)                 :: ps_name
          real(p_k_rdoub), dimension(2) :: origin

          real(p_k_rsing), dimension(:,:), allocatable :: ps_fld

          character(80)   :: path, name, name_dx, comment, unit
          character(len=80), dimension(2) :: xname, xlabel, xunits
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero
          
          integer(p_k_is) :: i, l, i1, i2, d1, d2, npt1, npt2
          real(p_k_rdoub) :: rdx1, rdx2, x1tmp, x2tmp
          real(p_k_rdoub) :: dvol, ps_norm_fac
          real(p_k_rdoub) :: dr, drh, r_cent, vol_fac_cent
          integer(p_k_is) :: r_dim

          integer :: ierr

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before write_phasespace_to_hdf'
!          call file_flush(file_id_msg)

!          write(file_id_msg, *) 'x_or_p : ', x_or_p
!          write(file_id_msg, *) 'xp_dim : ', xp_dim
!          call file_flush(file_id_msg)       

          izero =  ichar('0')

!         select the right parameters for this phasespace
          ps_name = ' '
          do i=1, 2
             if ( x_or_p(i) .eq. 1 ) then
                if (  ( this%ps_xmin(xp_dim(i)) .eq. 0.0_p_k_rdoub ) &
     &          .and. ( this%ps_xmax(xp_dim(i)) .ge. 0.0_p_k_rdoub ) ) then

!                  moving with the window of the simulation

                   ps_min(i) = xmin(xp_dim(i))
                   if (this%ps_xmax(xp_dim(i)).eq.0.0_p_k_rdoub) then
                      ps_max(i) = xmax(xp_dim(i))
                   else
                      ps_max(i) = xmin(xp_dim(i)) &
     &                          + this%ps_xmax(xp_dim(i))
                   endif
                   origin(i) = 0.0_p_k_rdoub

                else

!                  using the preset boundaries of the simulation
                   ps_min(i) = xmin(xp_dim(i))
                   ps_max(i) = xmax(xp_dim(i))
                   origin(i) = ps_min(i)

                endif

                ps_n(i) = this%ps_nx(xp_dim(i))
                ps_name = 'x'//char(izero+mod(xp_dim(i)/1,10))//trim(ps_name)
                
             else
             
                if ( this%ps_pmin(xp_dim(i)) .lt. this%ps_pmax(xp_dim(i)) ) then
                   ps_min(i) = this%ps_pmin(xp_dim(i))
                   ps_max(i) = this%ps_pmax(xp_dim(i))
                else
!                  set to some reasonable (but arbitrary ) value
!                  if not set in the input file
                   ps_min(i) = -30.0_p_k_rdoub
                   ps_max(i) =  30.0_p_k_rdoub
                endif
                ps_n(i) = this%ps_np(xp_dim(i))
                ps_name = 'p'//char(izero+mod(xp_dim(i)/1,10))//trim(ps_name)
                origin(i) = ps_min(i)
             endif
          enddo



!         allocate phasespace array
          allocate( ps_fld(ps_n(1),ps_n(2)) , stat=ierr)    
          if (ierr .ne. 0) then
            write(*,*) '[',my_aid(no_co),'] - (*error*) memory allocation failed'  
            write(*,*) '[',my_aid(no_co),'] - (*error*) in write_phasespace_to_hdf'  
            write(*,*) '[',my_aid(no_co),'] - (*error*) in file os-dspec.f90'  
            write(*,*) '[',my_aid(no_co),"] - (*error*) couldn't allocate ps_fld(",ps_n(1),',',ps_n(2),')' 
            write(*,*) '[',my_aid(no_co),'] - (*error*) bailing out!'
            
            write(file_id_msg,*) 'allocation problem in'
            write(file_id_msg,*) 'write_phasespace_to_hdf'
            write(file_id_msg,*) 'STOP'
            call file_flush(file_id_msg)
            call abort_program( -1_p_k_is )
          endif


          ps_dx  = ( ps_max - ps_min ) / ps_n
          ps_rdx = ps_n / ( ps_max - ps_min )
          ps_fld = 0.0_p_k_rsing

          rdx1  = ps_rdx(1)
          x1tmp = ps_min(1)
          d1    = xp_dim(1)
          npt1  = ps_n(1)

          rdx2  = ps_rdx(2)
          x2tmp = ps_min(2)
          d2    = xp_dim(2)
          npt2  = ps_n(2)
      
! ************** New code         
!       Selects the appropriate dimension and direction for each coordinate of
!       the phase space

		 select case (x_or_p(1))
			case (1) 
			  xname(1) = 'x'//char(izero + d1)//' axis'
			  xlabel(1) = 'x'//char(izero + d1)
			  xunits(1) = 'c / !Mw!Dp!N'
			case (2) 
			  xname(1) = 'p'//char(izero + d1)//' axis'
			  xlabel(1) = 'p'//char(izero + d1)
			  xunits(1) = 'm!De!N c'
			case default
			  write(*,*) 'variable x_or_p = ', x_or_p
			  write(*,*) 'out of range - STOP'
			  stop
		 end select 

		 ! select name and units for axis 2
	
		 select case (x_or_p(2))
			case (1) 
			  xname(2) = 'x'//char(izero + d2)//' axis'
			  xlabel(2) = 'x'//char(izero + d2)
			  xunits(2) = 'c / !Mw!Dp!N'
			case (2) 
			  xname(2) = 'p'//char(izero + d2)//' axis'
			  xlabel(2) = 'p'//char(izero + d2)
			  xunits(2) = 'm!De!N c'
			case default
			  write(*,*) 'variable x_or_p = ', x_or_p
			  write(*,*) 'out of range - STOP'
			  stop
		 end select 

		 ! select the appropriate units for this phasespace

		 unit = ' '

		 if (x_or_p(1) .ne. x_or_p(2)) then
			unit = 'e m!De!U-1!N c!U-2!N !Mw!Dp!N'
		 elseif (x_or_p(1) .eq. 1) then
			unit = 'e c!U-2!N !Mw!Dp!U2!N'
		 else
			unit = 'e m!De!U-2!N c!U-2!N'
		 endif

        ! select name and units for axis 1

          if (num_par > 0) then 
          
              ! because of a problem on the SP compiler
              ! we cannot use xp1 => x(d1,:), so we need to
              ! allocate a new array and copy the values into it
              
              !allocate(xp1(num_par), xp2(num_par))
              
              select case (x_or_p(1))
                case (1) 
                  !xp1 => x(d1,:)
                  !do l =1, num_par
                  !  xp1(l) = x(d1,l)
                  !enddo
                  xp1 => x
                case (2) 
                  !xp1 => p(d1, :)
                  !do l = 1, num_par
                  !  xp1(l) = p(d1,l)
                  !enddo
                  xp1 => p 
              end select 

              ! select name and units for axis 2
        
              select case (x_or_p(2))
                case (1) 
                  !xp2 => x(d2,:) 
                  !do l =1, num_par
                  !  xp2(l) = x(d2,l)
                  !enddo
                  xp2 => x
                case (2) 
                  !xp2 => p(d2,:) 
                  !do l =1, num_par
                  !  xp2(l) = p(d2,l)
                  !enddo
                  xp2 => p
              end select 
 
             select case (this%dep_sch)
               case (1)
				  !   Calculates phase space by linear interpolation
				  do l = 1, num_par             
					 
					 !nci1 = (xp1(l)-x1tmp)*rdx1+0.5
					 !nci2 = (xp2(l)-x2tmp)*rdx2+0.5
					 
					 nci1 = (xp1(d1,l)-x1tmp)*rdx1+0.5
					 nci2 = (xp2(d2,l)-x2tmp)*rdx2+0.5
					
					
					 i1 = nci1
					 i2 = nci2
					
					
					 xpdgp1 = nci1 - (i1)
					 xpdgp2 = nci2 - (i2)
				
					 ! i1,i2
					 if (  (i1 .ge. 1) .and. (i1 .le. npt1) .and. &
						   (i2 .ge. 1) .and. (i2 .le. npt2) ) then
								 ps_fld(i1,i2) = ps_fld(i1,i2) + &
			&                    (1.0-xpdgp1)*(1.0-xpdgp2)*q(l)
					 endif
					 
					 ! i1+1,i2
					 if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
							(i2    .ge. 1) .and.  (i2    .le. npt2) ) then
								 ps_fld((i1+1),i2) = ps_fld((i1+1),i2) + &
			&                     (xpdgp1)*(1.0-xpdgp2)*q(l)
					 endif
					 
					 ! i1,i2+1
					 if (  (    i1 .ge. 1) .and. ( i1    .le. npt1) .and. &
						   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
								 ps_fld(i1,(i2+1)) = ps_fld(i1,(i2+1)) + &
			&                     (1.0-xpdgp1)*(xpdgp2)*q(l)
					 endif
					
					 ! i1+1,i2+1
					 if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
						   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
								 ps_fld((i1+1),(i2+1)) = ps_fld((i1+1),(i2+1)) + &
			&                     (xpdgp1)*(xpdgp2)*q(l)
					endif
				  enddo      
          
               case default
                  !   Calculates phase space by nearest grid cell method
                  do l = 1, num_par
                    !i1 = (xp1(l)-x1tmp)*rdx1+1
                    !i2 = (xp2(l)-x2tmp)*rdx2+1
                    i1 = (xp1(d1,l)-x1tmp)*rdx1+1
                    i2 = (xp2(d2,l)-x2tmp)*rdx2+1
                    if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                      i2 .ge. 1 .and. i2.le.npt2 ) then
                             ps_fld(i1,i2) = ps_fld(i1,i2) + q(l)
                    endif
                  enddo
             end select

             !deallocate(xp1, xp2)
             
          endif ! num_par > 0

! ************** End New code         

!         general normalization factors
!          - these variabales are newly introduced on 05/02/99
!         needed for all simulations
          ps_norm_fac = product( ps_rdx )

!         dvol is only required to remove the normalization
!         of the charge with respect to the computational grid
          dvol = product( dx ) 

!         normalization to phase space density
          ps_fld = ps_fld * ps_norm_fac

!         collect data from all nodes on node 1
          call sum_up_array( ps_fld, no_co ) 

          me = my_aid(no_co)      

          if ( me .eq. 1 ) then

!            normalize phasespace data according to coordinate system
             select case ( coordinates )
              case ( p_cylindrical_b )
                !write(*,*) "correcting geometry, p_cylindrical_b"
                r_dim = 2
                do i=1, 2
                   if ((x_or_p(i).eq.1).and.(xp_dim(i).eq.r_dim)) then
                      dr  = ps_dx(i)
                      drh = dr / 2
                      do i2 = 1, ps_n(i)
                         r_cent = ps_min(i)+(i2-1)*dr+drh
!                         write(*,*) i2, r_cent
                         if ( abs(r_cent) .ge. drh ) then
                            vol_fac_cent = 1 / r_cent
                         else
                            vol_fac_cent = 2*dr / (abs(r_cent)+drh)**2
                            if ( r_cent .lt. 0.0_p_k_rdoub ) then
                               vol_fac_cent = - vol_fac_cent
                            endif
                         endif
!                        next code line: modification on 09/07/99 to fix a
!                        an error for phasespaces that have a radial axis
!                        -> simulations before 09/07/99 have a scaling error
                         vol_fac_cent = vol_fac_cent / ( 2 * pi )
!                         write(*,*) vol_fac_cent
                         select case (i)
                          case(1)
                            do i1 = 1, ps_n(2)
!                               write(*,*) '        ',i1,ps_fld(i2,i1)
                               ps_fld(i2,i1)=ps_fld(i2,i1)*vol_fac_cent
                            enddo
                          case(2)
                            do i1 = 1, ps_n(1)
!                               write(*,*) '        ',i1,ps_fld(i1,i2)
                               ps_fld(i1,i2)=ps_fld(i1,i2)*vol_fac_cent
                            enddo
                          case default
                            write(*,*) 'variable i = ', i
                            write(*,*) 'out of range - STOP'
                            stop
                         end select
                      enddo
                   endif
                enddo
!               factor 2 * pi is required to give 
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi 
              case ( p_cylindrical_e )
                !write(*,*) "correcting geometry, p_cylindrical_e"
                r_dim = 2
                do i=1, 2
                   if ((x_or_p(i).eq.1).and.(xp_dim(i).eq.r_dim)) then
                      dr  = ps_dx(i)
                      drh = dr / 2
                      do i2 = 1, ps_n(i)
                         r_cent = abs(ps_min(i)+(i2-1)*dr+drh)
!                         write(*,*) i2, r_cent
                         if ( r_cent .ge. drh ) then
                            vol_fac_cent = 1 / r_cent
                         else
                            vol_fac_cent = 2*dr / (r_cent+drh)**2
                         endif
!                        next code line: modification on 09/07/99 to fix a
!                        an error for phasespaces that have a radial axis
!                        -> simulations before 09/07/99 have a scaling error
                         vol_fac_cent = vol_fac_cent / ( 2 * pi )
!                         write(*,*) vol_fac_cent
                         select case (i)
                          case(1)
                            do i1 = 1, ps_n(2)
!                               write(*,*) '        ',i1,ps_fld(i2,i1)
                               ps_fld(i2,i1)=ps_fld(i2,i1)*vol_fac_cent
                            enddo
                          case(2)
                            do i1 = 1, ps_n(1)
!                               write(*,*) '        ',i1,ps_fld(i1,i2)
                               ps_fld(i1,i2)=ps_fld(i1,i2)*vol_fac_cent
                            enddo
                          case default
                            write(*,*) 'variable i = ', i
                            write(*,*) 'out of range - STOP'
                            stop
                         end select
                      enddo
                   endif
                enddo
!               factor 2 * pi is required to give
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi
              case ( p_cartesian )
                ps_fld = ps_fld * dvol
              case default
                continue
             end select           

!            prepare path and file names

             ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                // char( izero + mod( sp_id/ 1 , 10 ) ) 

!             chindx =    char( 1 + izero + mod( n_aux/1000 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/ 100 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/  10 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/   1 , 10 ) ) 

             path  = trim(path_mass) // 'PHA' // p_dir_sep &
     &            // trim(ps_name)   // p_dir_sep &
     &            // trim(ch_sp_id)  // p_dir_sep 

             comment = 'Phasespace ' // trim(ps_name)

!             name =  trim(this%file_name_pha)  &
!     &            // trim(ps_name)  // '-' &
!     &            // trim(ch_sp_id) // '-' // trim(chindx)

!             name = trim(name)//'.hdf'

            
             name = get_filename( n_aux, trim(this%file_name_pha)  &
     &                            // trim(ps_name)  // '-' // trim(ch_sp_id), &
     &                            '.hdf')           
             
             call write_hdf_sds( trim(path), trim(name), &
     &                           ps_fld, comment, 'Charge Density', unit,&
     &                           n, t, ps_n, origin, ps_dx, &
     &                           xname,xlabel,xunits)

          endif

          deallocate( ps_fld )

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now after write_phasespace_to_dx'
!          call file_flush(file_id_msg)

        end subroutine write_phasespace_to_hdf
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_spatial_ene_to_hdf( x, p, q, num_par, sp_id, &
     &                                       this, no_co, xmin, xmax, &
     &                                       n, t, n_aux, &
     &                                       x_dim, &
     &                                       dx, coordinates )
!---------------------------------------------------
!
! this will write spatially resolved average kinetic energy to file
! this implementation is a simplified version of write_phasespace
! and is not optimized since every node will allocate a global space
! array. 
!
! A more efficient implementation will only require local space and then
! merge (rather than sum) on the firt node, just like the charge deposition
!
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

          real(p_k_rdoub), dimension(:,:), intent(in) , target :: x, p
          real(p_k_rdoub), dimension(:),   intent(in) :: q
          integer(p_k_is),                 intent(in) :: num_par
          integer(p_k_is),                 intent(in) :: sp_id
          type( t_diag_species ),          intent(in) :: this
          type( t_node_conf ),             intent(in) :: no_co
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
          integer(p_k_is),                 intent(in) :: n
          real(p_k_rdoub),                 intent(in) :: t
          integer(p_k_is),                 intent(in) :: n_aux
          integer(p_k_is), dimension(:),   intent(in) :: x_dim
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
          integer(p_k_is),                 intent(in) :: coordinates

!         meaning of some variable values:
!         x_dim(i) = k :  use direction k of x data 
!                         for axis i of phasespace plot

!       local variables


!       pointer to x
          real(p_k_rdoub), dimension(:), pointer :: x1, x2

!       normalised coordinates of grid point
          real(p_k_rdoub) :: nci1,nci2
!       distance to nearest grid point in each dimension (normalised)
          real(p_k_rdoub) :: xpdgp1, xpdgp2

          real(p_k_rdoub), dimension(2) :: ps_min, ps_max, ps_dx, ps_rdx
          integer(p_k_is), dimension(2) :: ps_n
          character(20)                 :: ps_name
          real(p_k_rdoub), dimension(2) :: origin

          real(p_k_rsing), dimension(:,:), allocatable :: ps_fld, ps_den

          character(80)   :: path, name, name_dx, comment, unit
          character(len=80), dimension(2) :: xname, xlabel, xunits
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero
          
          integer(p_k_is) :: i, l, i1, i2, d1, d2, npt1, npt2
          real(p_k_rdoub) :: rdx1, rdx2, x1tmp, x2tmp
          real(p_k_rdoub) :: dvol, ps_norm_fac
          real(p_k_rdoub) :: dr, drh, r_cent, vol_fac_cent
          integer(p_k_is) :: r_dim

!       particle kinetic energy
          real(p_k_rdoub) :: kin_ene      ! abs(q) * (gamma - 1) 

          integer :: ierr

!       executable statements

          izero =  ichar('0')

!         select the right parameters for this phasespace
          ps_name = ' '
          do i=1, 2
             if (     ( this%ps_xmin(x_dim(i)) .eq. 0.0_p_k_rdoub ) &
     &          .and. ( this%ps_xmax(x_dim(i)) .ge. 0.0_p_k_rdoub ) ) then

!                  moving with the window of the simulation
                   ps_min(i) = xmin(x_dim(i))
                   if (this%ps_xmax(x_dim(i)).eq.0.0_p_k_rdoub) then
                      ps_max(i) = xmax(x_dim(i))
                   else
                      ps_max(i) = xmin(x_dim(i)) &
     &                          + this%ps_xmax(x_dim(i))
                   endif
                   origin(i) = 0.0_p_k_rdoub

             else
!                using the preset boundaries of the simulation
                 ps_min(i) = xmin(x_dim(i))
                 ps_max(i) = xmax(x_dim(i))
                 origin(i) = ps_min(i)
             endif
             ps_n(i) = this%ps_nx(x_dim(i))
             ps_name = 'x' // char(izero+mod(x_dim(i)/1,10)) &
     &                 // trim(ps_name)
          enddo

          ps_name = trim(ps_name) // "_ene"

!         allocate phasespace array
          allocate( ps_fld(ps_n(1),ps_n(2)) , stat=ierr)
     
          if (ierr .ne. 0) then
            write(file_id_msg,*) 'allocation problem in'
            write(file_id_msg,*) 'write_spatial_ene_to_hdf'
            write(file_id_msg,*) 'STOP'
            call file_flush(file_id_msg)
            call abort_program( -1_p_k_is )
          endif

          ps_dx  = ( ps_max - ps_min ) / ps_n
          ps_rdx = ps_n / ( ps_max - ps_min )
          ps_fld = 0.0_p_k_rsing

          rdx1  = ps_rdx(1)
          x1tmp = ps_min(1)
          d1    = x_dim(1)
          npt1  = ps_n(1)

          rdx2  = ps_rdx(2)
          x2tmp = ps_min(2)
          d2    = x_dim(2)
          npt2  = ps_n(2)

!       allocate density array
          allocate( ps_den(ps_n(1),ps_n(2)) , stat=ierr)
     
          if (ierr .ne. 0) then
            write(file_id_msg,*) 'allocation problem in'
            write(file_id_msg,*) 'write_spatial_ene_to_hdf'
            write(file_id_msg,*) 'STOP'
            call file_flush(file_id_msg)
            call abort_program( -1_p_k_is )
          endif
          
          ps_den = 0.0_p_k_rsing
      
        ! select name and units for axis 1

        x1 => x(d1,1:num_par)
        xname(1) = 'x'//char(izero + d1)//' axis'
        xlabel(1) = 'x'//char(izero + d1)
        xunits(1) = 'c / !Mw!Dp!N'

        ! select name and units for axis 2

        x2 => x(d2,1:num_par) 
        xname(2) = 'x'//char(izero + d2)//' axis'
        xlabel(2) = 'x'//char(izero + d2)
        xunits(2) = 'c / !Mw!Dp!N'

        ! select the appropriate units for this phasespace
    
        unit = 'm c!U2!N'
 
        
        ! deposit energy
        select case (this%dep_sch)
          case (1)
          !   Calculates phase space by linear interpolation
          do l = 1, num_par             
             
             kin_ene = abs(q(l)) * (sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1) - 1)
             
             nci1 = (x1(l)-x1tmp)*rdx1+0.5
             nci2 = (x2(l)-x2tmp)*rdx2+0.5
            
             i1 = nci1
             i2 = nci2
            
            
             xpdgp1 = nci1 - (i1)
             xpdgp2 = nci2 - (i2)
        
             ! i1,i2
             if (  (i1 .ge. 1) .and. (i1 .le. npt1) .and. &
                   (i2 .ge. 1) .and. (i2 .le. npt2) ) then
                         ps_fld(i1,i2) = ps_fld(i1,i2) + &
    &                    (1.0-xpdgp1)*(1.0-xpdgp2)*kin_ene
                         ps_den(i1,i2) = ps_den(i1,i2) + &
    &                    (1.0-xpdgp1)*(1.0-xpdgp2)*q(l)
    
             endif
             
             ! i1+1,i2
             if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
                    (i2    .ge. 1) .and.  (i2    .le. npt2) ) then
                         ps_fld((i1+1),i2) = ps_fld((i1+1),i2) + &
    &                     (xpdgp1)*(1.0-xpdgp2)*kin_ene
                         ps_den((i1+1),i2) = ps_den((i1+1),i2) + &
    &                     (xpdgp1)*(1.0-xpdgp2)*q(l)
             endif
             
             ! i1,i2+1
             if (  (    i1 .ge. 1) .and. ( i1    .le. npt1) .and. &
                   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
                         ps_fld(i1,(i2+1)) = ps_fld(i1,(i2+1)) + &
    &                     (1.0-xpdgp1)*(xpdgp2)*kin_ene
                         ps_den(i1,(i2+1)) = ps_den(i1,(i2+1)) + &
    &                     (1.0-xpdgp1)*(xpdgp2)*q(l)
             endif
            
             ! i1+1,i2+1
             if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
                   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
                         ps_fld((i1+1),(i2+1)) = ps_fld((i1+1),(i2+1)) + &
    &                     (xpdgp1)*(xpdgp2)*kin_ene
                         ps_den((i1+1),(i2+1)) = ps_den((i1+1),(i2+1)) + &
    &                     (xpdgp1)*(xpdgp2)*q(l)
            endif
          enddo      
          
          case default
             !   Calculates phase space by nearest grid cell method
            do l = 1, num_par
              i1 = (x1(l)-x1tmp)*rdx1+1
              i2 = (x2(l)-x2tmp)*rdx2+1
                if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 ) then
                     kin_ene = abs(q(l)) * (sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1) - 1)
                     ps_fld(i1,i2) = ps_fld(i1,i2) + kin_ene
                     ps_den(i1,i2) = ps_den(i1,i2) + q(l)
                endif
            enddo
        end select


!         collect data from all nodes on node 1
          call sum_up_array( ps_fld, no_co ) 
          call sum_up_array( ps_den, no_co ) 

          me = my_aid(no_co)      

          if ( me .eq. 1 ) then

!            normalize energy

             do i1 = 1, npt1
               do i2 = 1, npt2
                 
                 if (ps_den(i1,i2) .ne. 0.0) then
                   ps_fld(i1,i2) = ps_fld(i1,i2) / abs(ps_den(i1,i2))
                 endif
                 
               enddo
             enddo


!            prepare path and file names

             ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                // char( izero + mod( sp_id/ 1 , 10 ) ) 

             path  = trim(path_mass) // 'PHA' // p_dir_sep &
     &            // trim(ps_name)   // p_dir_sep &
     &            // trim(ch_sp_id)  // p_dir_sep 

             comment = 'Energy density ' // trim(ps_name)

            
             name = get_filename( n_aux, trim(this%file_name_pha)  &
     &                            // trim(ps_name)  // '-' // trim(ch_sp_id), &
     &                            '.hdf')           
             
             call write_hdf_sds( trim(path), trim(name), &
     &                           ps_fld, comment, 'Energy Density', unit,&
     &                           n, t, ps_n, origin, ps_dx, &
     &                           xname,xlabel,xunits)

          endif

          deallocate( ps_fld, ps_den )

        end subroutine write_spatial_ene_to_hdf
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_bin_ene_to_hdf( x, p, q, num_par, sp_id, &
     &                                       this, no_co, xmin, xmax, &
     &                                       n, t, n_aux, &
     &                                       x_dim, &
     &                                       dx, coordinates )
!---------------------------------------------------
!
! this will write energy binned charge density to file
! this implementation is a simplified version of write_phasespace
! and is not optimized since every node will allocate a global space
! array. 
!
! A more efficient implementation will only require local space and then
! merge (rather than sum) on the firt node, just like the charge deposition
!
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

          real(p_k_rdoub), dimension(:,:), intent(in) , target :: x, p
          real(p_k_rdoub), dimension(:),   intent(in) :: q
          integer(p_k_is),                 intent(in) :: num_par
          integer(p_k_is),                 intent(in) :: sp_id
          type( t_diag_species ),          intent(in) :: this
          type( t_node_conf ),             intent(in) :: no_co
          real(p_k_rdoub), dimension(:),   intent(in) :: xmin, xmax
          integer(p_k_is),                 intent(in) :: n
          real(p_k_rdoub),                 intent(in) :: t
          integer(p_k_is),                 intent(in) :: n_aux
          integer(p_k_is), dimension(:),   intent(in) :: x_dim
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
          integer(p_k_is),                 intent(in) :: coordinates

!         meaning of some variable values:
!         x_dim(i) = k :  use direction k of x data 
!                         for axis i of phasespace plot

!       local variables


!       pointer to x
          real(p_k_rdoub), dimension(:), pointer :: x1, x2

!       normalised coordinates of grid point
          real(p_k_rdoub) :: nci1,nci2
!       distance to nearest grid point in each dimension (normalised)
          real(p_k_rdoub) :: xpdgp1, xpdgp2

          real(p_k_rdoub), dimension(2) :: ps_min, ps_max, ps_dx, ps_rdx
          integer(p_k_is), dimension(2) :: ps_n
          character(20)                 :: ps_name
          real(p_k_rdoub), dimension(2) :: origin

          real(p_k_rsing), dimension(:,:), allocatable :: ps_fld, ps_den

          character(80)   :: path, name, name_dx, comment, unit
          character(len=80), dimension(2) :: xname, xlabel, xunits
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero
          
          integer(p_k_is) :: i, l, i1, i2, d1, d2, npt1, npt2
          real(p_k_rdoub) :: rdx1, rdx2, x1tmp, x2tmp
          real(p_k_rdoub) :: dvol, ps_norm_fac
          real(p_k_rdoub) :: dr, drh, r_cent, vol_fac_cent
          integer(p_k_is) :: r_dim

!       particle kinetic energy
          real(p_k_rdoub) :: kin_ene      ! abs(q) * (gamma - 1) 
          
          ! HDF related variables
          type(t_attr_hdf) :: sd_attr
          integer :: sdfileID
          
          ! misc
          integer :: stage, test, ierr
          integer(p_k_is) :: npartbin

!       executable statements

          if (this%n_ene_bins < 1) return
          
          
          izero =  ichar('0')

!         select the right parameters for this phasespace
          ps_name = ' '
          do i=1, 2
             if (     ( this%ps_xmin(x_dim(i)) .eq. 0.0_p_k_rdoub ) &
     &          .and. ( this%ps_xmax(x_dim(i)) .ge. 0.0_p_k_rdoub ) ) then

!                  moving with the window of the simulation
                   ps_min(i) = xmin(x_dim(i))
                   if (this%ps_xmax(x_dim(i)).eq.0.0_p_k_rdoub) then
                      ps_max(i) = xmax(x_dim(i))
                   else
                      ps_max(i) = xmin(x_dim(i)) &
     &                          + this%ps_xmax(x_dim(i))
                   endif
                   origin(i) = 0.0_p_k_rdoub

             else
!                using the preset boundaries of the simulation
                 ps_min(i) = xmin(x_dim(i))
                 ps_max(i) = xmax(x_dim(i))
                 origin(i) = ps_min(i)
             endif
             ps_n(i) = this%ps_nx(x_dim(i))
             ps_name = 'x' // char(izero+mod(x_dim(i)/1,10)) &
     &                 // trim(ps_name)
          enddo

          ps_name = trim(ps_name) // "_bin_ene"

!         allocate phasespace array
          allocate( ps_fld(ps_n(1),ps_n(2)) , stat=ierr)
     
          if (ierr .ne. 0) then
            write(file_id_msg,*) 'allocation problem in'
            write(file_id_msg,*) 'write_spatial_ene_to_hdf'
            write(file_id_msg,*) 'STOP'
            call file_flush(file_id_msg)
            call abort_program( -1_p_k_is )
          endif

          ps_dx  = ( ps_max - ps_min ) / ps_n
          ps_rdx = ps_n / ( ps_max - ps_min )

          rdx1  = ps_rdx(1)
          x1tmp = ps_min(1)
          d1    = x_dim(1)
          npt1  = ps_n(1)

          rdx2  = ps_rdx(2)
          x2tmp = ps_min(2)
          d2    = x_dim(2)
          npt2  = ps_n(2)

!         general normalization factors
          ps_norm_fac = product( ps_rdx )

!         dvol is only required to remove the normalization
!         of the charge with respect to the computational grid
          dvol = product( dx ) 


        ! select name and units for axis 1

        x1 => x(d1,1:num_par)
        xname(1) = 'x'//char(izero + d1)//' axis'
        xlabel(1) = 'x'//char(izero + d1)
        xunits(1) = 'c / !Mw!Dp!N'

        ! select name and units for axis 2

        x2 => x(d2,1:num_par) 
        xname(2) = 'x'//char(izero + d2)//' axis'
        xlabel(2) = 'x'//char(izero + d2)
        xunits(2) = 'c / !Mw!Dp!N'

        ! select the appropriate units for this phasespace
    
        unit = 'e c!U-2!N !Mw!Dp!U2!N'



        ! open file
        
        me = my_aid(no_co)      

        if (me .eq. 1) then

!            prepare path and file names

             ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                // char( izero + mod( sp_id/ 1 , 10 ) ) 

             path  = trim(path_mass) // 'PHA' // p_dir_sep &
     &            // trim(ps_name)   // p_dir_sep &
     &            // trim(ch_sp_id)  // p_dir_sep 

             comment = 'Charge Density ' // trim(ps_name)

            
             name = get_filename( n_aux, trim(this%file_name_pha)  &
     &                            // trim(ps_name)  // '-' // trim(ch_sp_id), &
     &                            '.hdf')           
             
          ! open file

          call create_hdf( trim(path), trim(name), sdfileID)

        
          ! create attributes
          
      
          call clear(sd_attr)
      
          sd_attr%rank  = 2
          sd_attr%units = unit 
          sd_attr%n     = n
          sd_attr%t     = t
          sd_attr%nx(1:2)    = ps_n(1:2)
          sd_attr%xmin(1:2)  = ps_min(1:2)
          sd_attr%dx(1:2)    = ps_dx(1:2)
      
          sd_attr%xname(1:2)  = xname(1:2)
          sd_attr%xlabel(1:2) = xlabel(1:2)
          sd_attr%xunits(1:2) = xunits(1:2)

          sd_attr%name  = 'Charge Density ' // trim(ps_name)
          sd_attr%label = 'Charge Density'
        endif

        ! Process all energy bins

        
        ! 1st stage, write particles with energies < ene_bin(1)
        ! ...
        ! nth stage, write particles with energies >= ene_bin(n-1) and < ene_bin(n)
        ! ...
        ! n_ene_bins stage, write particles with energies >= ene_bin(n_ene_bins)
        
        do stage = 1, this%n_ene_bins + 1
                
          ! set test type and label
          
          if (stage .eq. 1) then 
              test = 1
              write(sd_attr%sds_label, *) "Energy <= ", this%ene_bins(1)
              !sd_attr%sds_label = "Energy <= "//trim(sd_attr%sds_label)
          else 
            if (stage .eq. this%n_ene_bins + 1) then
              test = 3
              write(sd_attr%sds_label, *) "Energy >= ", this%ene_bins(this%n_ene_bins)
              !sd_attr%sds_label = "Energy >= "//trim(sd_attr%sds_label)

            else 
              test = 2
              write(sd_attr%sds_label, *) this%ene_bins(stage-1),&
                                      " < Energy <= ",&
                                      this%ene_bins(stage)
            endif
          endif
        
          if (me .eq. 1) write(*,*) trim(sd_attr%sds_label)
  
         
          ! initialize density
          ps_fld = 0.0_p_k_rsing

          npartbin = 0
          do l = 1, num_par             
             
             kin_ene = (sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1) - 1)
             
             select case(test)
               case(1)
                 if (kin_ene .ge. this%ene_bins(1)) cycle
               case(2)
                 if ((kin_ene .lt. this%ene_bins(stage-1)) .or. &
                     (kin_ene .ge. this%ene_bins(stage))) cycle
               case(3)
                 if (kin_ene .lt. this%ene_bins(this%n_ene_bins)) cycle
             end select
                         
             npartbin = npartbin+1
             
             nci1 = (x1(l)-x1tmp)*rdx1+0.5
             nci2 = (x2(l)-x2tmp)*rdx2+0.5
            
             i1 = nci1
             i2 = nci2
            
            
             xpdgp1 = nci1 - (i1)
             xpdgp2 = nci2 - (i2)
        
             ! i1,i2
             if (  (i1 .ge. 1) .and. (i1 .le. npt1) .and. &
                   (i2 .ge. 1) .and. (i2 .le. npt2) ) then
                         ps_fld(i1,i2) = ps_fld(i1,i2) + &
    &                    (1.0-xpdgp1)*(1.0-xpdgp2)*q(l)
    
             endif
             
             ! i1+1,i2
             if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
                    (i2    .ge. 1) .and.  (i2    .le. npt2) ) then
                         ps_fld((i1+1),i2) = ps_fld((i1+1),i2) + &
    &                     (xpdgp1)*(1.0-xpdgp2)*q(l)
             endif
             
             ! i1,i2+1
             if (  (    i1 .ge. 1) .and. ( i1    .le. npt1) .and. &
                   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
                         ps_fld(i1,(i2+1)) = ps_fld(i1,(i2+1)) + &
    &                     (1.0-xpdgp1)*(xpdgp2)*q(l)
             endif
            
             ! i1+1,i2+1
             if (  ((i1+1) .ge. 1) .and. ((i1+1) .le. npt1) .and. &
                   ((i2+1) .ge. 1) .and. ((i2+1) .le. npt2) ) then
                         ps_fld((i1+1),(i2+1)) = ps_fld((i1+1),(i2+1)) + &
    &                     (xpdgp1)*(xpdgp2)*q(l)
            endif
          enddo      

!         normalization to phase space density
          ps_fld = ps_fld * ps_norm_fac
          
         ! collect data from all nodes on node 1
          call sum_up_array( ps_fld, no_co ) 
        
         ! count total number of deposited particles
          call sum_up( no_co, npartbin )
         
         ! normalize data             
          if ( me .eq. 1 ) then

!             write(*,*) "Number of particles in bin : ", npartbin

!            normalize phasespace data according to coordinate system
             select case ( coordinates )
              case ( p_cylindrical_b )
                r_dim = 2
                do i=1, 2
                   if (x_dim(i).eq.r_dim) then
                      dr  = ps_dx(i)
                      drh = dr / 2
                      do i2 = 1, ps_n(i)
                         r_cent = ps_min(i)+(i2-1)*dr+drh
!                         write(*,*) i2, r_cent
                         if ( abs(r_cent) .ge. drh ) then
                            vol_fac_cent = 1 / r_cent
                         else
                            vol_fac_cent = 2*dr / (abs(r_cent)+drh)**2
                            if ( r_cent .lt. 0.0_p_k_rdoub ) then
                               vol_fac_cent = - vol_fac_cent
                            endif
                         endif
!                        next code line: modification on 09/07/99 to fix a
!                        an error for phasespaces that have a radial axis
!                        -> simulations before 09/07/99 have a scaling error
                         vol_fac_cent = vol_fac_cent / ( 2 * pi )
!                         write(*,*) vol_fac_cent
                         select case (i)
                          case(1)
                            do i1 = 1, ps_n(2)
!                               write(*,*) '        ',i1,ps_fld(i2,i1)
                               ps_fld(i2,i1)=ps_fld(i2,i1)*vol_fac_cent
                               
                            enddo
                          case(2)
                            do i1 = 1, ps_n(1)
!                               write(*,*) '        ',i1,ps_fld(i1,i2)
                               ps_fld(i1,i2)=ps_fld(i1,i2)*vol_fac_cent
                            enddo
                          case default
                            write(*,*) 'variable i = ', i
                            write(*,*) 'out of range - STOP'
                            stop
                         end select
                      enddo
                   endif
                enddo
!               factor 2 * pi is required to give 
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi 
                
              case ( p_cylindrical_e )
                r_dim = 2
                do i=1, 2
                   if (x_dim(i).eq.r_dim) then
                      dr  = ps_dx(i)
                      drh = dr / 2
                      do i2 = 1, ps_n(i)
                         r_cent = abs(ps_min(i)+(i2-1)*dr+drh)
!                         write(*,*) i2, r_cent
                         if ( r_cent .ge. drh ) then
                            vol_fac_cent = 1 / r_cent
                         else
                            vol_fac_cent = 2*dr / (r_cent+drh)**2
                         endif
!                        next code line: modification on 09/07/99 to fix a
!                        an error for phasespaces that have a radial axis
!                        -> simulations before 09/07/99 have a scaling error
                         vol_fac_cent = vol_fac_cent / ( 2 * pi )
!                         write(*,*) vol_fac_cent
                         select case (i)
                          case(1)
                            do i1 = 1, ps_n(2)
!                               write(*,*) '        ',i1,ps_fld(i2,i1)
                               ps_fld(i2,i1)=ps_fld(i2,i1)*vol_fac_cent
                            enddo
                          case(2)
                            do i1 = 1, ps_n(1)
!                               write(*,*) '        ',i1,ps_fld(i1,i2)
                               ps_fld(i1,i2)=ps_fld(i1,i2)*vol_fac_cent
                            enddo
                          case default
                            write(*,*) 'variable i = ', i
                            write(*,*) 'out of range - STOP'
                            stop
                         end select
                      enddo
                   endif
                enddo
!               factor 2 * pi is required to give
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi
              case ( p_cartesian )
                ps_fld = ps_fld * dvol

              case default
                continue
             end select           

             ! write dataset to file
           
             call save_sds( ps_fld, sdfileID, sd_attr )          

        
          endif ! if (me .eq. 1)
          
                    
        enddo ! stage = 1, ...

 
        !deallocate memory
        deallocate( ps_fld )


        ! close file
        if (me .eq. 1) then
          
          call close_hdf(sdfileID)
        
        endif

        end subroutine write_bin_ene_to_hdf
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_particle_energies_to_file( p, q, rqm,dKE, num_par, sp_id, &
     &                                              this, no_co, n, t , dx)
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

          real(p_k_rdoub), dimension(:,:), intent(in) , target ::  p
          real(p_k_rdoub), dimension(:),   intent(in) :: q
          real(p_k_rdoub),                 intent(in) :: rqm
          real(p_k_rdoub),                 intent(in) :: dKE
          integer(p_k_is),                 intent(in) :: num_par
          integer(p_k_is),                 intent(in) :: sp_id
          type( t_diag_species ),          intent(in) :: this
          type( t_node_conf ),             intent(in) :: no_co
          integer(p_k_is),                 intent(in) :: n
          real(p_k_rdoub),                 intent(in) :: t
!          the dimensions of the simulation cell
          real(p_k_rdoub), dimension(:),   intent(in) :: dx
          
!       local variables


          character(80)   :: path, full_name
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero, result
          real(p_k_rdoub) :: dxvol
      
          integer(p_k_is) :: i, l , total_par

          real(p_k_rdoub), dimension (1:5) :: ene_res
          
      
!       the structure of ene_res is as follows
!          ene_res(1) = Heat Flux on direction 1
!          ene_res(2) = Heat Flux on direction 2
!          ene_res(3) = Heat Flux on direction 3
!          ene_res(4) = Kinetic Energy

          real(p_k_rdoub) :: gamma_aux1, gamma_aux2, gamma
          
          integer :: ierr

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before write_particle_energies'
!          call file_flush(file_id_msg)

           izero =  ichar('0')

           ene_res = 0.0_p_k_rdoub
         
           dxvol = product( dx )

!       calculates the energy components           
           do l = 1, num_par
              gamma_aux1 = p(1,l)**2
              do i= 2, p_p_dim         
                 gamma_aux1 = gamma_aux1 + p(i,l)**2
              enddo
              gamma = sqrt(gamma_aux1 + 1.0_p_k_rdoub)
              gamma_aux1 = (gamma-1.0_p_k_rdoub)
              gamma_aux2 = gamma_aux1/gamma
!!             Kinetic Energy Flux
!              do i= 1, p_p_dim
!                ene_res(i) = ene_res(i)+q(l)*gamma_aux2*p(i,l) 
!              enddo               
!             Momentum Flux
              do i= 1, p_p_dim
                ene_res(i) = ene_res(i)+q(l)*p(i,l) 
              enddo               
!             Kinetic Energy
              ene_res(4) = ene_res(4)+q(l)*gamma_aux1 
           enddo

!        normalize kinetic energy and flux ( x mass density)

           ene_res = ene_res*rqm*dxvol
           if(dKE<0.0) then
             ene_res(5) = 0.0;
           else
             ene_res(5) = dKE*dxvol;
           endif
          
!        sums the results from all nodes
          
            call sum_up_array( ene_res, no_co ) 
            me = my_aid(no_co)      
            
            total_par = num_par
            call sum_up( no_co, total_par )

            if ( me .eq. 1 ) then
             
             
!            prepare path and file names

               ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                  // char( izero + mod( sp_id/ 1 , 10 ) ) 

               path  = trim(path_hist) 


               full_name = trim(path) // trim(this%file_name_ene)  &
     &                  // 'par' // trim(ch_sp_id) // '_ene' 
     
 
!            Open file and position at the last record
               if ( t .eq. 0.0_p_k_rdoub ) then 
                  call mkdir( path, ierr )
                  
                  open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
                  form='formatted')

!            Write Header                 
                  write(file_id_parene,'( A6, 7(1X,A14) )') &
!                       'Iter','Time    ','Total_Par ','Qx1    ','Qx2    ','Qx3    ','Kin_Energy   ','KEnergy_tc'
                       'Iter','Time    ','Total_Par ','px1    ','px2    ','px3    ','Kin_Energy   ','KEnergy_tc'
             
               else  
                  open (unit=file_id_parene, file=full_name, position = 'append', &
                  form='formatted')
               endif
! --------------------- Writes the particle energies to file
               write(file_id_parene, '( I6,1X,g14.8,1X,I14,5(1X,g14.8) )') &
                      n, t, total_par, ene_res
               close(file_id_parene)

            endif

        end subroutine write_particle_energies_to_file
!---------------------------------------------------


!---------------------------------------------------
        subroutine write_3D_phasespace_to_dx( x, p, q, num_par, sp_id, &
     &                                        this, no_co, xmin, xmax, &
     &                                        n, t, n_aux, &
     &                                        x_or_p, xp_dim, &
     &                                        dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

!       particle position and momenta x( k - direction, # part number)
          real(p_k_rdoub), dimension(:,:),intent(in), target :: x, p
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          real(p_k_rdoub), dimension(:),     intent(in) :: xmin, xmax
          integer(p_k_is),                   intent(in) :: n
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          integer(p_k_is), dimension(:),     intent(in) :: x_or_p
          integer(p_k_is), dimension(:),     intent(in) :: xp_dim
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates

!         meaning of some variable values:
!         x_or_p(i) = 1 :  use x-data for axis i of phasespace plot
!         x_or_p(i) = 2 :  use p-data for axis i of phasespace plot
!         xp_dim(i) = k :  use direction k (of x or p data) 
!                          for axis i of phasespace plot

!       local variables

!       pointer to x or p
          real(p_k_rdoub), dimension(:), pointer :: xp1, xp2, xp3 

!       normalised coordinates of grid point
          real(p_k_rdoub) :: nci1,nci2,nci3
!       distance to nearest grid point in each dimension (normalised)
          real(p_k_rdoub) :: xpdgp1, xpdgp2, xpdgp3
      

      
      
      
!       physical minimum , maximum, cells size, reciprocal of cell size
!       of phase space for each direction
          real(p_k_rdoub), dimension(3) :: ps_min, ps_max, ps_dx, ps_rdx
!       Number of cells for each direction
          integer(p_k_is), dimension(3) :: ps_n
!       Phase space name
          character(20)                 :: ps_name
!       Phase space origin
          real(p_k_rdoub), dimension(3) :: origin
!       Phase space (3D)
          real(p_k_rsing), dimension(:,:,:), allocatable :: ps_fld

!       File related variables
          character(80)   :: path, name, name_dx, comment, unit
          character(len=80), dimension(3) :: xname, xlabel, xunits
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero
!          
          integer(p_k_is) :: dep_sch
          integer(p_k_is) :: i, l, i1, i2, i3, d1, d2, d3, npt1, npt2, npt3, pax, xax
          real(p_k_rdoub) :: rdx1, rdx2, rdx3, x1tmp, x2tmp, x3tmp
          real(p_k_rdoub) :: dvol, ps_norm_fac
          real(p_k_rdoub) :: dr, drh, r_cent, vol_fac_cent
          integer(p_k_is) :: r_dim
          
          integer :: ierr

!       executable statements

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now before write_3D_phasespace_to_dx'
!          call file_flush(file_id_msg)
          
!          write(file_id_msg, *) 'x_or_p : ', x_or_p
!          write(file_id_msg, *) 'xp_dim : ', xp_dim
!          call file_flush(file_id_msg)

          izero =  ichar('0')

!         select the right parameters for this phasespace
          ps_name = ' '
          do i=1, 3
             if ( x_or_p(i) .eq. 1 ) then
                if ( (  this%ps_xmin(xp_dim(i)) &
     &          .eq. 0.0_p_k_rdoub ) &
     &          .and. ( this%ps_xmax(xp_dim(i)) &
     &          .ge. 0.0_p_k_rdoub ) ) then
!                  moving with the window of the simulation
                   ps_min(i) = xmin(xp_dim(i))
                   if (this%ps_xmax(xp_dim(i)).eq.0.0_p_k_rdoub) then
                      ps_max(i) = xmax(xp_dim(i))
                   else
                      ps_max(i) = xmin(xp_dim(i)) &
     &                          + this%ps_xmax(xp_dim(i))
                   endif
                   origin(i) = 0.0_p_k_rdoub
                else
!                  using the preset boundaries of the simulation
                   ps_min(i) = xmin(xp_dim(i))
                   ps_max(i) = xmax(xp_dim(i))
                   origin(i) = ps_min(i)
                endif
                ps_n(i) = this%ps_nx_3D(xp_dim(i))
                ps_name = 'x' // char(izero+mod(xp_dim(i)/1,10)) &
     &                 // trim(ps_name)
             else
                if ( this%ps_pmin(xp_dim(i)) &
     &          .lt. this%ps_pmax(xp_dim(i)) ) then
                   ps_min(i) = this%ps_pmin(xp_dim(i))
                   ps_max(i) = this%ps_pmax(xp_dim(i))
                else
!                  set to some reasonable (but arbitrary ) value
!                  if not set in the input file
                   ps_min(i) = -30.0_p_k_rdoub
                   ps_max(i) =  30.0_p_k_rdoub
                endif
                ps_n(i) = this%ps_np_3D(xp_dim(i))
                ps_name = 'p' // char(izero+mod(xp_dim(i)/1,10)) &
     &                 // trim(ps_name)
                origin(i) = ps_min(i)
             endif
          enddo


!         allocate phasespace array
          allocate( ps_fld(ps_n(1),ps_n(2),ps_n(3)) , stat=ierr)
     
        if (ierr .ne. 0) then
          write(file_id_msg,*) 'allocation problem in'
          write(file_id_msg,*) 'write_3D_phasespace_to_dx'
          write(file_id_msg,*) 'STOP'
          call file_flush(file_id_msg)
          call abort_program( -1_p_k_is )
        endif

          ps_dx  = ( ps_max - ps_min ) / ps_n
          ps_rdx = ps_n / ( ps_max - ps_min )
          ps_fld = 0.0_p_k_rsing

          rdx1  = ps_rdx(1)
          x1tmp = ps_min(1)
          d1    = xp_dim(1)
          npt1  = ps_n(1)

          rdx2  = ps_rdx(2)
          x2tmp = ps_min(2)
          d2    = xp_dim(2)
          npt2  = ps_n(2)

          rdx3  = ps_rdx(3)
          x3tmp = ps_min(3)
          d3    = xp_dim(3)
          npt3  = ps_n(3)
    

!       Selects the appropriate dimension and direction for each coordinate of
!       the phase space

        xname(1) = ' '
        xname(2) = ' '
        xname(3) = ' '
        xlabel(1) = ' '
        xlabel(2) = ' '
        xlabel(3) = ' '
        xunits(1) = ' '
        xunits(2) = ' '
        xunits(3) = ' '

        select case (x_or_p(1))
           case (1) 
                xp1 => x(xp_dim(1),1:num_par) 
            xname(1) = 'x'//char(izero + xp_dim(1))//' axis'
            xlabel(1) = 'x'//char(izero + xp_dim(1))
            xunits(1) = 'c / !Mw!Dp!N'
           case (2) 
                xp1 => p(xp_dim(1),1:num_par) 
                xname(1) = 'p'//char(izero + xp_dim(1))//' axis'
            xlabel(1) = 'p'//char(izero + xp_dim(1))
            xunits(1) = 'm!De!N c'
           case default
                     write(*,*) 'variable x_or_p = ', x_or_p
                     write(*,*) 'out of range - STOP'
                     stop
        end select 

        select case (x_or_p(2))
           case (1) 
                xp2 => x(xp_dim(2),1:num_par)
            xname(2) = 'x'//char(izero + xp_dim(2))//' axis'
            xlabel(2) = 'x'//char(izero + xp_dim(2))
            xunits(2) = 'c / !Mw!Dp!N'

           case (2) 
                xp2 => p(xp_dim(2),1:num_par) 
                xname(2) = 'p'//char(izero + xp_dim(2))//' axis'
            xlabel(2) = 'p'//char(izero + xp_dim(2))
            xunits(2) = 'm!De!N c'

           case default
                     write(*,*) 'variable x_or_p = ', x_or_p
                     write(*,*) 'out of range - STOP'
                     stop
        end select 

        select case (x_or_p(3))
           case (1) 
                xp3 => x(xp_dim(3),1:num_par)
            xname(3) = 'x'//char(izero + xp_dim(3))//' axis'
            xlabel(3) = 'x'//char(izero + xp_dim(3))
            xunits(3) = 'c / !Mw!Dp!N'

           case (2) 
                xp3 => p(xp_dim(3),1:num_par) 
                xname(3) = 'p'//char(izero + xp_dim(3))//' axis'
            xlabel(3) = 'p'//char(izero + xp_dim(3))
            xunits(3) = 'm!De!N c'

           case default
                     write(*,*) 'variable x_or_p = ', x_or_p
                     write(*,*) 'out of range - STOP'
                     stop
        end select 

!      Selects the appropriate units for this phasespace
    
       unit = 'e'
     
       pax = 0
       xax = 0
       
       do l=1,3
         if (x_or_p(l) .eq. 1) then 
        xax = xax + 1
     else 
         pax = pax + 1
     endif
       enddo

       if (pax .gt. 0) then
           unit = unit//' m!De!U-'//char( izero + pax )//'!N'
       endif
       
       unit= unit//' c!U-'//char( izero + pax + xax)//'!N'

       if (xax .gt. 0) then
           unit= unit//' !Mw!Dp!U-'//char(izero + xax)//'!N'
       endif
       
     
!      Calculates phasespace

    select case (this%dep_sch)
          case (1)
!            Calculates phase space by linear interpolation
         do l = 1, num_par
            nci1 = (xp1(l)-x1tmp)*rdx1+0.5
            nci2 = (xp2(l)-x2tmp)*rdx2+0.5
            nci3 = (xp3(l)-x3tmp)*rdx3+0.5
        
            i1 = nci1
            i2 = nci2
            i3 = nci3
        
            xpdgp1 = nci1 - (i1)
            xpdgp2 = nci2 - (i2)
            xpdgp3 = nci3 - (i3)    

        ! i1,i2,i3
            if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 .and. &
                  i3 .ge. 1 .and. i3.le.npt3) then
                         ps_fld(i1,i2,i3) = ps_fld(i1,i2,i3) + &
                         (1-xpdgp1)*(1-xpdgp2)*(1-xpdgp3)*q(l)
            endif
        ! i1+1,i2,i3
            if (  (i1+1) .ge. 1 .and. (i1+1).le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 .and. &
                  i3 .ge. 1 .and. i3.le.npt3) then
                         ps_fld((i1+1),i2,i3) = ps_fld((i1+1),i2,i3) + &
                         (xpdgp1)*(1-xpdgp2)*(1-xpdgp3)*q(l)
            endif
        ! i1,i2+1,i3
            if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  (i2+1) .ge. 1 .and. (i2+1).le.npt2 .and. &
                  i3 .ge. 1 .and. i3.le.npt3) then
                         ps_fld(i1,(i2+1),i3) = ps_fld(i1,(i2+1),i3) + &
                          (1-xpdgp1)*(xpdgp2)*(1-xpdgp3)*q(l)
            endif
        ! i1+1,i2+1,i3
            if (  (i1+1) .ge. 1 .and. (i1+1).le.npt1 .and. &
                  (i2+1) .ge. 1 .and. (i2+1).le.npt2 .and. &
                  i3 .ge. 1 .and. i3.le.npt3) then
                         ps_fld((i1+1),(i2+1),i3) = ps_fld((i1+1),(i2+1),i3) + &
                          (xpdgp1)*(xpdgp2)*(1-xpdgp3)*q(l)
            endif
        ! i1,i2,i3+1
            if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 .and. &
                  (i3+1) .ge. 1 .and. (i3+1).le.npt3) then
                         ps_fld(i1,i2,(i3+1)) = ps_fld(i1,i2,(i3+1)) + &
                         (1-xpdgp1)*(1-xpdgp2)*(xpdgp3)*q(l)
            endif
        ! i1+1,i2,i3+1
            if (  (i1+1) .ge. 1 .and. (i1+1).le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 .and. &
                  (i3+1) .ge. 1 .and. (i3+1).le.npt3) then
                         ps_fld((i1+1),i2,(i3+1)) = ps_fld((i1+1),i2,(i3+1)) + &
                          (xpdgp1)*(1-xpdgp2)*(xpdgp3)*q(l)
            endif
        ! i1,i2+1,i3+1
            if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  (i2+1) .ge. 1 .and. (i2+1).le.npt2 .and. &
                  (i3+1) .ge. 1 .and. (i3+1).le.npt3) then
                         ps_fld(i1,(i2+1),(i3+1)) = ps_fld(i1,(i2+1),(i3+1)) + &
                          (1-xpdgp1)*(xpdgp2)*(xpdgp3)*q(l)
            endif
        ! i1+1,i2+1,i3+1
                if (  (i1+1) .ge. 1 .and. (i1+1).le.npt1 .and. &
                      (i2+1) .ge. 1 .and. (i2+1).le.npt2 .and. &
                      (i3+1) .ge. 1 .and. (i3+1).le.npt3) then
                          ps_fld((i1+1),(i2+1),(i3+1)) = &
     &                    ps_fld((i1+1),(i2+1),(i3+1)) + &
              (xpdgp1)*(xpdgp2)*(xpdgp3)*q(l)
            endif
         enddo      
      case default
!            Calculates phase space by nearest grid cell method
         do l = 1, num_par
            i1 = (xp1(l)-x1tmp)*rdx1+1
            i2 = (xp2(l)-x2tmp)*rdx2+1
            i3 = (xp3(l)-x3tmp)*rdx3+1
                if (  i1 .ge. 1 .and. i1.le.npt1 .and. &
                  i2 .ge. 1 .and. i2.le.npt2 .and. &
              i3 .ge. 1 .and. i3.le.npt3) then
                         ps_fld(i1,i2,i3) = ps_fld(i1,i2,i3) + q(l)
                endif
         enddo
        end select
       
!         general normalization factors
!          - these variabales are newly introduced on 05/02/99
!         needed for all simulations
          ps_norm_fac = product( ps_rdx )

!         dvol is only required to remove the normalization
!         of the charge with respect to the computational grid
          dvol = product( dx )

!         normalization to phase space density
          ps_fld = ps_fld * ps_norm_fac

!         collect data from all nodes on node 1
          call sum_up_array( ps_fld, no_co ) 

          me = my_aid(no_co)      

          if ( me .eq. 1 ) then

!  ------------ Only cartesian coordinates in 3D so normalisation reduces to  
             ps_fld = ps_fld * dvol

!            prepare path and file names

             ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                // char( izero + mod( sp_id/ 1 , 10 ) ) 

!             chindx =    char( 1 + izero + mod( n_aux/1000 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/ 100 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/  10 , 10 ) ) &
!     &               //  char(     izero + mod( n_aux/   1 , 10 ) ) 

             path  = trim(path_mass) // 'PHA' // p_dir_sep &
     &            // trim(ps_name)   // p_dir_sep &
     &            // trim(ch_sp_id)  // p_dir_sep 

!             name =  trim(this%file_name_pha)  &
!     &            // trim(ps_name)  // '-' &
!     &            // trim(ch_sp_id) // '-' // trim(chindx)

         comment = 'Phasespace ' // trim(ps_name)    
!        name = trim(name)//'.hdf'
        
         name = get_filename( n_aux, trim(this%file_name_pha)  &
     &                        // trim(ps_name)  // '-' // trim(ch_sp_id),  &
     &                        '.hdf')
      
         call write_hdf_sds( trim(path), trim(name), &
     &                           ps_fld, comment, 'Charge Density', unit,&
     &                           n, t, ps_n, origin, ps_dx, &
     &                           xname,xlabel,xunits)

          endif

          deallocate( ps_fld )

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now after write_3D_phasespace_to_dx'
!         call file_flush(file_id_msg)
!          write(*, *) 'now after write_3D_phasespace_to_dx'


        end subroutine write_3D_phasespace_to_dx
!---------------------------------------------------

!---------------------------------------------------
        subroutine write_1D_gamma_phasespace_to_dx( p, q, num_par, sp_id, &
     &                                        this, no_co, n, t, n_aux, &
     &                                        dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

!       particle momenta p( k - direction, # part number)
          real(p_k_rdoub), dimension(:,:),intent(in)    :: p
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          integer(p_k_is),                   intent(in) :: n
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates


!       local variables


!       normalised coordinates of grid point
          real(p_k_rdoub) :: nci
!       distance to nearest grid point 
          real(p_k_rdoub) :: xpdgp
              
      
!       physical minimum , maximum, cells size, reciprocal of cell size
!       of phase space 
          real(p_k_rdoub) :: ps_min, ps_max, ps_rdx
          real(p_k_rdoub), dimension(1) :: ps_dx
!       Number of cells 
          integer(p_k_is), dimension(1) :: ps_n
!       Phase space name
          character(20)                 :: ps_name
!       Phase space origin
          real(p_k_rdoub),dimension(1)               :: origin
!       Phase space (1D)
          real(p_k_rsing), dimension(:), allocatable :: ps_fld

!       Total momenta
          real(p_k_rsing)               :: gamma

!       File related variables
          character(80)   :: path, name, name_dx, comment
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero
!          
          integer(p_k_is) :: dep_sch
          integer(p_k_is) :: i, l, npt
          real(p_k_rdoub) :: rdx, xtmp
          real(p_k_rdoub) :: dvol, ps_norm_fac
          real(p_k_rdoub) :: dr, drh, r_cent, vol_fac_cent
          integer(p_k_is) :: r_dim
          
          integer :: ierr

!       executable statements

!        write(file_id_msg, *) ' '
!        write(file_id_msg, *) 'now before write_1D_gamma_phasespace_to_dx'
!        call file_flush(file_id_msg)
        
        izero =  ichar('0')

!       select the right parameters for this phasespace
       
        ps_min = this%ps_gammamin
        ps_max = this%ps_gammamax

        if (this%if_ps_gamma_log) then
! SP
!         ps_min = alog10(ps_min)
!         ps_max = alog10(ps_max)
          ps_min = dlog10(ps_min)
          ps_max = dlog10(ps_max)
! SP
        endif

        ps_n(1) = this%ps_ngamma
        ps_name = 'gamma'
        origin(1) = ps_min

!       allocate phasespace array
        allocate( ps_fld(ps_n(1)) , stat=ierr)
     
        if (ierr .ne. 0) then
          write(file_id_msg,*) 'allocation problem in'
          write(file_id_msg,*) 'write_1D_gamma_phasespace_to_dx'
          write(file_id_msg,*) 'STOP'
          call file_flush(file_id_msg)
          call abort_program( -1_p_k_is )
        endif

        ps_dx(1)  = ( ps_max - ps_min ) / ps_n(1)
        ps_rdx = ps_n(1) / ( ps_max - ps_min )
        ps_fld = 0.0_p_k_rsing

        rdx  = ps_rdx
        xtmp = ps_min
        npt  = ps_n(1)

        select case (this%dep_sch)
          case (1)
!            Calculates phase space by linear interpolation
             do l = 1, num_par
               gamma = sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1)
               
               if (this%if_ps_gamma_log) gamma = alog10(gamma)
               
               nci = (gamma-xtmp)*rdx+1
               ! assigning a real to an integer effectively does
               ! rounding by truncation i.e. rounds toward zero
               ! since gamma>=0 there's no need to use the floor function
               i = nci
               xpdgp = nci - (i)
               
               ! i
               if (  i .ge. 1 .and. i.le.npt) then
                  ps_fld(i) = ps_fld(i) + (1-xpdgp)*q(l)
               endif
               
               ! i+1
               if (  (i+1) .ge. 1 .and. (i+1).le.npt) then
                  ps_fld(i+1) = ps_fld(i+1) + (xpdgp)*q(l)
               endif
             enddo
          case default
!            Calculates phase space by nearest grid cell method
             do l = 1, num_par
               gamma = sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1)
               
               if (this%if_ps_gamma_log) gamma = alog10(gamma)
               
               i = (gamma-xtmp)*rdx+1
               if (  i .ge. 1 .and. i.le.npt) then
                 ps_fld(i) = ps_fld(i) + q(l)
               endif
             enddo
        end select
        
!       general normalization factors
!        - these variabales are newly introduced on 05/02/99
!       needed for all simulations
        ps_norm_fac =  ps_rdx 

!       dvol is only required to remove the normalization
!       of the charge with respect to the computational grid
        dvol = product( dx )

!       normalization to phase space density
        ps_fld = ps_fld * ps_norm_fac

!       collect data from all nodes on node 1
        call sum_up_array( ps_fld, no_co ) 

        me = my_aid(no_co)      

        if ( me .eq. 1 ) then

!          normalize phasespace data according to coordinate system
           select case ( coordinates )
              case ( p_cylindrical_b )
!               factor 2 * pi is required to give 
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi 
              case ( p_cylindrical_e )
!               factor 2 * pi is required to give
!               the real charge of a "ring particle"
                ps_fld = ps_fld * dvol * 2 * pi
              case ( p_cartesian )
                ps_fld = ps_fld * dvol
              case default
                continue
           end select

!          prepare path and file names

           ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &                // char( izero + mod( sp_id/ 1 , 10 ) ) 

           path  = trim(path_mass) // 'PHA' // p_dir_sep &
     &            // trim(ps_name)   // p_dir_sep &
     &            // trim(ch_sp_id)  // p_dir_sep 

           comment = 'Phasespace ' // trim(ps_name)

           name = get_filename( n_aux, trim(this%file_name_pha)  &
     &                          // trim(ps_name)  // '-' // trim(ch_sp_id) , &
     &                          '.hdf' )

           call write_hdf_sds( trim(path), trim(name), &
     &                           ps_fld, comment, 'Charge Density', 'e',&
     &                           n, t, ps_n, origin, ps_dx, &
     &                           (/'gamma axis'/),(/'!Mg'/),(/''/), &
     &                           xlog = (/this%if_ps_gamma_log/))

          endif

          deallocate( ps_fld )

!          write(file_id_msg, *) ' '
!          write(file_id_msg, *) 'now after write_1D_gamma_phasespace_to_dx'
!          call file_flush(file_id_msg)

        end subroutine write_1D_gamma_phasespace_to_dx
!---------------------------------------------------

!---------------------------------------------------
        subroutine write_raw_( x, p, q, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

!       particle momenta p( k - direction, # part number)
          real(p_k_rdoub), dimension(:,:),   intent(in) :: x, p
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          real(p_k_rdoub),                   intent(in) :: rqm
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          integer(p_k_is),                   intent(in) :: n
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates


!       local variables

          integer :: sdfileID               ! HDF SD file ID number
          
!       local variables

          character(80)   :: path
          character(80)   :: full_name, file_name
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero

          real(p_k_rsing) :: gamma_limit_sqr, gamma_par_sqr
          integer :: i, l, lbuf

!        Buffer to hold particle data
          real(p_k_rsing), allocatable, dimension(:,:) :: buffer
          integer :: num_par_buffer, buffer_dim
 
!        particles written to the buffer on all nodes         
          integer, allocatable, dimension(:) :: num_par_buffer_all
 
!        local variables for HDF          
          character(20)                         :: sdname          
          integer, dimension(2)                 :: hdfstart, stride, edge, dimsize
          integer                               :: sdsID  

!        local variables for MPI
          integer :: handle
          integer, dimension(mpi_status_size):: stat
          integer :: result, count
          integer :: dest, source, tag
          integer(p_k_is) :: dim, bnd

!       Additional file data
          integer(p_k_is), dimension(p_x_dim) :: node_conf 

          integer :: ierr


          ! initialize the buffers
          buffer_dim = num_par
          num_par_buffer = 0
          allocate( buffer(p_x_dim + p_p_dim + 1,buffer_dim))

          
          ! initialize gamma limit
          
          gamma_limit_sqr = this%raw_gamma_limit**2
          
          ! loop over all particles

          l = 1
          lbuf = 0

          ! Copy all the particles to the temp buffer

          select case ( coordinates )

                  ! cylindrical coordinates not implemented yet

           case default
            
            do                                
              if (l .gt. num_par)  exit             
              gamma_par_sqr = 1. + sum( p(:,l)**2 )
              if ( gamma_par_sqr .ge. gamma_limit_sqr ) then
                 if ( random( ) .le. this%raw_fraction ) then
                       
                   lbuf = lbuf + 1
                   if (lbuf .gt. buffer_dim) then
                      write(file_id_msg,*) 'ERROR: temp buffer is full in write_raw'
                      exit
                    endif
                     
                    do i=1, p_x_dim
                       buffer(i,lbuf) = x(i,l)
                    enddo
                    do i=1, p_p_dim
                       buffer(i+p_x_dim,lbuf) = p(i,l)
                    enddo
                    buffer(1+p_x_dim+p_p_dim,lbuf) = q(l)                      
                 endif                   
              endif              
              l = l+1
            enddo
             

          end select
         
          num_par_buffer = lbuf
          dim = 0
          bnd = 0

          if (my_aid(no_co) .eq. 1) then 
            
            ! Create the file   

            ! Prepare path and file name
            izero =  ichar('0')
            ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &               // char( izero + mod( sp_id/ 1 , 10 ) ) 

          
            path  = trim(path_mass) // 'RAW' // p_dir_sep &
     &           // trim(ch_sp_id)  // p_dir_sep 
          
            
            file_name = get_filename( n_aux, ''  &
     &                          // 'RAW'  // '-' // trim(ch_sp_id) , &
     &                          '.hdf' )
     
            full_name = trim(path) // trim(file_name)
             
            ! Create directory if necessary
            call mkdir( path, ierr )
 
           
            ! Create the HDF file
            sdfileID = sfstart(full_name, DFACC_CREATE)
            ! Add simulation time to the file
            ierr = sfsnatt(sdfileID, 'TIME', DFNT_FLOAT64, 1, t) 
            ! Add simulation timestep to the file 
            ierr = sfsnatt(sdfileID, 'ITER', DFNT_INT32, 1, n)
            ! Add node number information to the file
            ierr = sfsnatt(sdfileID, 'NODE NUMBER', DFNT_INT32, 1, my_aid(no_co))
            ! Add node configuration information to the file
            node_conf = my_nx(no_co)
            dim = p_x_dim
            ierr = sfsnatt(sdfileID, 'NODE CONFIGURATION', DFNT_INT32, dim, node_conf )
            dim = 0
            ! Add gamma limit information to the file
            ierr = sfsnatt(sdfileID, 'GAMMA LIMIT', DFNT_FLOAT64, 1, this%raw_gamma_limit )
            ! Add particle fraction information to the file
            ierr = sfsnatt(sdfileID, 'PARTICLE FRACTION', DFNT_FLOAT64, 1, &
                           this%raw_fraction )          
         
            
            ! Create species dataset

            dimsize(1) = p_x_dim + p_p_dim + 1
            dimsize(2) = SD_UNLIMITED
  
 
            sdname = 'Species '// char( ichar('0') + mod( sp_id/ 10 , 10 ) ) &
     &                         // char( ichar('0') + mod( sp_id/  1 , 10 ) )

            sdsID = sfcreate(sdfileID, sdname, DFNT_FLOAT32 , 2, dimsize) 
            
            !write(*,*) "Dataset created, sdsID ", sdsID
            
            ierr = sfsblsz(sdsID, (p_x_dim + p_p_dim + 1)*1024*256)
            !write(*,*) "Block size set to ", (p_x_dim + p_p_dim + 1)*1024*256
            
            
            ! Add species ID to the dataset 
            ierr = sfsnatt(sdsID, 'SP_ID', DFNT_INT32, 1, sp_id)

            ! Add rqm to the dataset
            ierr = sfsnatt(sdsID, 'RQM', DFNT_FLOAT64, 1, rqm)

            !write(*,*) "Information added, adding data from node 0"
                     
            ! write node 0 data to file
            stride = 1
            hdfstart = 0
            edge(1) = dimsize(1)
          
            if (num_par_buffer .gt. 0) then
! DEBUG
!               write(*,*) "write_raw, Number of particles on node 0 is ", num_par_buffer
! DEBUG
               ! write data
               edge(2) = num_par_buffer
               ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
               if (ierr .ne. 0) write(*,*) "sfwdata failed for node 0 data"
               ! Update start position                   
               hdfstart(2) = hdfstart(2) + num_par_buffer
               !write(*,*) "data from node 0 added"
            endif
            
            
            
            ! clear buffer
            
            deallocate(buffer)
            
            if (no_num(no_co) .gt. 1) then            
                        
              ! get number of particles on other nodes
              allocate(num_par_buffer_all(no_num(no_co)))
              
              call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 
            
              ! loop through remaining nodes
            
              do i=2, no_num(no_co)             
               if (num_par_buffer_all(i) .gt. 0) then
! DEBUG
!                 write(*,*) "write_raw, Number of particles on node ", i-1," is ", num_par_buffer_all(i)
! DEBUG
               
                 ! prepare buffer to receive data
                 allocate( buffer(p_x_dim + p_p_dim + 1,num_par_buffer_all(i)))
              
               
                 ! get data from node i
                 tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,i), dim, bnd )
                 call send_ping( tid(no_co,i), tag )
                 
                 ! get data
                 source = tid(no_co, i)
                 tag    = tag_value( tag_code_spraw, tid(no_co, i), dim, bnd )               
                 count = (1+p_x_dim+p_p_dim) * num_par_buffer_all(i)
               
                 call MPI_RECV( buffer, count , MPI_REAL, source, &
     &                          tag, MPI_COMM_WORLD, stat, ierr )

                 ! save date from node i
               
                 edge(2) = num_par_buffer_all(i)
                 ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
                 !if (ierr .ne. 0) write(*,*) "sfwdata failed for node ", i-1, " data"
                 ! Update start position                   
                 hdfstart(2) = hdfstart(2) + num_par_buffer_all(i)
               
                 ! clear buffer
                 deallocate(buffer)
               
                 !write(*,*) "Data from node ",i-1, " added"
               endif 
               
              enddo 
               
              deallocate(num_par_buffer_all) 
            
            endif
            
            ! Close the SDS
            ierr = sfendacc(sdsID)
                           
            ! Close HDF file
            ierr = sfend(sdfileID)
         
          else          
            ! send number of particles on local node
            allocate(num_par_buffer_all(no_num(no_co)))

            call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 

            deallocate(num_par_buffer_all) 
            
            if (num_par_buffer .gt. 0) then             
              ! send particles on local node to node 1                          
              tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,my_aid(no_co)), dim, bnd )
              call recv_ping( tid(no_co, 1), tag )
           
              dest   = tid(no_co, 1)
              tag = tag_value( tag_code_spraw, tid(no_co,my_aid(no_co)), dim, bnd )
              count = (1+p_x_dim+p_p_dim) * num_par_buffer
              
              call MPI_ISEND( buffer, count, MPI_REAL, dest, tag, &
     &                        MPI_COMM_WORLD, handle, ierr )
              call MPI_WAIT( handle, stat, ierr )
            endif
            
            ! clear buffer
            deallocate(buffer)
          endif
          

          !write(file_id_msg, *) ' '
          !write(file_id_msg, *) 'now after write_raw'
          !call file_flush(file_id_msg)

        end subroutine write_raw_
!---------------------------------------------------

! p.E
!---------------------------------------------------
        subroutine write_raw_pE( x, p, pE, q, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

!       particle momenta p( k - direction, # part number)
! p.E          
          real(p_k_rdoub), dimension(:,:),   intent(in) :: x, p, pE
! p.E         
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          real(p_k_rdoub),                   intent(in) :: rqm
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          integer(p_k_is),                   intent(in) :: n
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates


!       local variables

          integer :: sdfileID               ! HDF SD file ID number
! p.E        
          integer :: p_pE_dim
! p.E
          
!       local variables

          character(80)   :: path
          character(80)   :: full_name, file_name
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero

          real(p_k_rsing) :: gamma_limit_sqr, gamma_par_sqr
          integer :: i, l, lbuf

!        Buffer to hold particle data
          real(p_k_rsing), allocatable, dimension(:,:) :: buffer
          integer :: num_par_buffer, buffer_dim
 
!        particles written to the buffer on all nodes         
          integer, allocatable, dimension(:) :: num_par_buffer_all
 
!        local variables for HDF          
          character(20)                         :: sdname          
          integer, dimension(2)                 :: hdfstart, stride, edge, dimsize
          integer                               :: sdsID  

!        local variables for MPI
          integer :: handle
          integer, dimension(mpi_status_size):: stat
          integer :: result, count
          integer :: dest, source, tag
          integer(p_k_is) :: dim, bnd

!       Additional file data
          integer(p_k_is), dimension(p_x_dim) :: node_conf 

          integer :: ierr

          ! Create species dataset
          p_pE_dim = p_p_dim
! p.E	 
          dimsize(1) = p_x_dim + p_p_dim + p_pE_dim + 1
! p.E 
          dimsize(2) = SD_UNLIMITED

          ! initialize the buffers
          buffer_dim = num_par
          num_par_buffer = 0
          allocate( buffer(dimsize(1),buffer_dim))

          
          ! initialize gamma limit
          
          gamma_limit_sqr = this%raw_gamma_limit**2
          
          ! loop over all particles

          l = 1
          lbuf = 0

          ! Copy all the particles to the temp buffer

          select case ( coordinates )

                  ! cylindrical coordinates not implemented yet

           case default
            
            do                                
              if (l .gt. num_par)  exit             
              gamma_par_sqr = 1. + sum( p(:,l)**2 )
              if ( gamma_par_sqr .ge. gamma_limit_sqr ) then
                 if ( random( ) .le. this%raw_fraction ) then
                       
                   lbuf = lbuf + 1
                   if (lbuf .gt. buffer_dim) then
                      write(file_id_msg,*) 'ERROR: temp buffer is full in write_raw'
                      exit
                    endif
                     
                    do i=1, p_x_dim
                       buffer(i,lbuf) = x(i,l)
                    enddo
                    do i=1, p_p_dim
                       buffer(i+p_x_dim,lbuf) = p(i,l)
                    enddo

! p.E                  
                    do i=1, p_pE_dim
                       buffer(i+p_x_dim+p_p_dim,lbuf) = pE(i,l)
                    enddo 
                    
                       buffer(1+p_x_dim+p_p_dim+p_pE_dim,lbuf) = q(l) 
                 endif                   
! p.E 
              endif              
              l = l+1
            enddo
             
          end select
         
          num_par_buffer = lbuf
          dim = 0
          bnd = 0

          if (my_aid(no_co) .eq. 1) then 
            
            ! Create the file   

            ! Prepare path and file name
            izero =  ichar('0')
            ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &               // char( izero + mod( sp_id/ 1 , 10 ) ) 

          
            path  = trim(path_mass) // 'RAW' // p_dir_sep &
     &           // trim(ch_sp_id)  // p_dir_sep 
          
            
            file_name = get_filename( n_aux, ''  &
     &                          // 'RAW'  // '-' // trim(ch_sp_id) , &
     &                          '.hdf' )
     
            full_name = trim(path) // trim(file_name)
             
            ! Create directory if necessary
            call mkdir( path, ierr )
 
           
            ! Create the HDF file
            sdfileID = sfstart(full_name, DFACC_CREATE)
            ! Add simulation time to the file
            ierr = sfsnatt(sdfileID, 'TIME', DFNT_FLOAT64, 1, t) 
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe:"time"',t
               write(*,*) 'abording'
               call abort_program()
            endif
            ! Add simulation timestep to the file 
            ierr = sfsnatt(sdfileID, 'ITER', DFNT_INT32, 1, n)
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe:"iter"',n
               write(*,*) 'abording'
               call abort_program()
            endif
            ! Add node number information to the file
            ierr = sfsnatt(sdfileID, 'NODE_NUMBER', DFNT_INT32, 1, my_aid(no_co))
!            write(*,*) 'No problem in write_RAW_pe:"node number"',ierr, my_aid(no_co)
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe:"node number"',my_aid(no_co)
               write(*,*) 'abording'
               call abort_program()
            endif
            ! Add node configuration information to the file
            node_conf = my_nx(no_co)
            dim = p_x_dim
            ierr = sfsnatt(sdfileID, 'NODE_CONFIGURATION', DFNT_INT32, dim, node_conf )
            dim = 0
            ! Add gamma limit information to the file
            ierr = sfsnatt(sdfileID, 'GAMMA_LIMIT', DFNT_FLOAT64, 1, this%raw_gamma_limit )
            ! Add particle fraction information to the file
            ierr = sfsnatt(sdfileID, 'PARTICLE_FRACTION', DFNT_FLOAT64, 1, &
                           this%raw_fraction )          
         
            
            ! Create species dataset
! p.E 
            p_pE_dim = p_p_dim
            dimsize(1) = p_x_dim + p_p_dim + p_pE_dim + 1
! p.E 
            dimsize(2) = SD_UNLIMITED
  
 
            sdname = 'Species'// char( ichar('0') + mod( sp_id/ 10 , 10 ) ) &
     &                         // char( ichar('0') + mod( sp_id/  1 , 10 ) )

            sdsID = sfcreate(sdfileID, sdname, DFNT_FLOAT32 , 2, dimsize) 
            
! DEBUG
!             write(*,*) "Dataset created, sdsID ", sdsID
! DEBUG
            
            ierr = sfsblsz(sdsID, (p_x_dim + p_p_dim + p_pE_dim + 1)*1024*256)
! DEBUG
!            write(*,*) "Block size set to ", (p_x_dim + p_p_dim + p_pE_dim + 1)*1024*256
! DEBUG
            
            
            ! Add species ID to the dataset 
            ierr = sfsnatt(sdsID, 'SP_ID', DFNT_INT32, 1, sp_id)

            ! Add rqm to the dataset
            ierr = sfsnatt(sdsID, 'RQM', DFNT_FLOAT64, 1, rqm)

! DEBUG
!            write(*,*) "Information added, adding data from node 0"
! DEBUG
                     
            ! write node 0 data to file
            stride = 1
            hdfstart = 0
            edge(1) = dimsize(1)
          
            if (num_par_buffer .gt. 0) then
! DEBUG
!               write(*,*) "write_raw_pE, Number of particles on node 0 is ", num_par_buffer
! DEBUG
               ! write data
               edge(2) = num_par_buffer
               ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
               if (ierr .ne. 0) write(*,*) "sfwdata failed for node 0 data"
               ! Update start position                   
               hdfstart(2) = hdfstart(2) + num_par_buffer
               !write(*,*) "data from node 0 added"
            endif
            
            
            
            ! clear buffer
            
            deallocate(buffer)
            
            if (no_num(no_co) .gt. 1) then            
                        
              ! get number of particles on other nodes
              allocate(num_par_buffer_all(no_num(no_co)))
              
              call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 
            
              ! loop through remaining nodes
            
              do i=2, no_num(no_co)             
               if (num_par_buffer_all(i) .gt. 0) then
! DEBUG
!                 write(*,*) "write_raw_pE, Number of particles on node ", i-1," is ", num_par_buffer_all(i)
! DEBUG
               
                 ! prepare buffer to receive data
! p.E             
                 allocate( buffer(p_x_dim + p_p_dim + p_pE_dim + 1,num_par_buffer_all(i)))
! p.E
              
               
                 ! get data from node i
                 tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,i), dim, bnd )
                 call send_ping( tid(no_co,i), tag )
                 
                 ! get data
                 source = tid(no_co, i)
                 tag    = tag_value( tag_code_spraw, tid(no_co, i), dim, bnd )               
! p.E
                 count = (1+p_x_dim+p_p_dim + p_pE_dim) * num_par_buffer_all(i)
! p.E         
               
                 call MPI_RECV( buffer, count , MPI_REAL, source, &
     &                          tag, MPI_COMM_WORLD, stat, ierr )

                 ! save date from node i
               
                 edge(2) = num_par_buffer_all(i)
                 ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
                 !if (ierr .ne. 0) write(*,*) "sfwdata failed for node ", i-1, " data"
                 ! Update start position                   
                 hdfstart(2) = hdfstart(2) + num_par_buffer_all(i)
               
                 ! clear buffer
                 deallocate(buffer)
               
! DEBUG
!                 write(*,*) "write_raw_pE, Data from node ",i-1, " added"
! DEBUG
               endif 
               
              enddo 
               
              deallocate(num_par_buffer_all) 
            
            endif
            
            ! Close the SDS
            ierr = sfendacc(sdsID)
                           
            ! Close HDF file
            ierr = sfend(sdfileID)
         
          else          
            ! send number of particles on local node
            allocate(num_par_buffer_all(no_num(no_co)))

            call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 

            deallocate(num_par_buffer_all) 
            
            if (num_par_buffer .gt. 0) then             
              ! send particles on local node to node 1                          
              tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,my_aid(no_co)), dim, bnd )
              call recv_ping( tid(no_co, 1), tag )
           
              dest   = tid(no_co, 1)
              tag = tag_value( tag_code_spraw, tid(no_co,my_aid(no_co)), dim, bnd )
! TZOUFRAS
              count = (1+p_x_dim+p_p_dim + p_pE_dim) * num_par_buffer
! TZOUFRAS
              
              call MPI_ISEND( buffer, count, MPI_REAL, dest, tag, &
     &                        MPI_COMM_WORLD, handle, ierr )
              call MPI_WAIT( handle, stat, ierr )
            endif
            
            ! clear buffer
            deallocate(buffer)
          endif
          

! DEBUG
!           write(*, *) ' '
!           write(*, *) 'now after write_raw'
! DEBUG
!          call file_flush(file_id_msg)

        end subroutine write_raw_pE
!---------------------------------------------------


! p.E
!---------------------------------------------------
        subroutine write_raw_pE_par_xid( x, p, pE, q, par_xid, num_par, rqm, sp_id, &
     &                        this, no_co, n, t, n_aux, &
     &                        dx, coordinates )
!---------------------------------------------------

          use m_units

          implicit none

!       dummy variables

!       particle momenta p( k - direction, # part number)
! p.E          
          real(p_k_rdoub), dimension(:,:),   intent(in) :: x, p, pE, par_xid
! p.E         
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          real(p_k_rdoub),                   intent(in) :: rqm
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          integer(p_k_is),                   intent(in) :: n
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates


!       local variables

          integer :: sdfileID               ! HDF SD file ID number
! p.E        
          integer :: p_pE_dim
! p.E
          
!       local variables

          character(80)   :: path
          character(80)   :: full_name, file_name
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero

          real(p_k_rsing) :: gamma_limit_sqr, gamma_par_sqr
          integer :: i, l, lbuf

!        Buffer to hold particle data
          real(p_k_rsing), allocatable, dimension(:,:) :: buffer
          integer :: num_par_buffer, buffer_dim
 
!        particles written to the buffer on all nodes         
          integer, allocatable, dimension(:) :: num_par_buffer_all
 
!        local variables for HDF          
          character(20)                         :: sdname          
          integer, dimension(2)                 :: hdfstart, stride, edge, dimsize
          integer                               :: sdsID  

!        local variables for MPI
          integer :: handle
          integer, dimension(mpi_status_size):: stat
          integer :: result, count
          integer :: dest, source, tag
          integer(p_k_is) :: dim, bnd

!       Additional file data
          integer(p_k_is), dimension(p_x_dim) :: node_conf 

          integer :: ierr

          ! Create species dataset
          p_pE_dim = p_p_dim
! p.E	 
          dimsize(1) = p_x_dim + p_p_dim + p_pE_dim + 1 + p_x_dim
! p.E 
          dimsize(2) = SD_UNLIMITED

          ! initialize the buffers
          buffer_dim = num_par
          num_par_buffer = 0
          allocate( buffer(dimsize(1),buffer_dim))

          
          ! initialize gamma limit
          
          gamma_limit_sqr = this%raw_gamma_limit**2
          
          ! loop over all particles

          l = 1
          lbuf = 0

          ! Copy all the particles to the temp buffer

          select case ( coordinates )

                  ! cylindrical coordinates not implemented yet

           case default
            
            do                                
              if (l .gt. num_par)  exit             
              gamma_par_sqr = 1. + sum( p(:,l)**2 )
              if ( gamma_par_sqr .ge. gamma_limit_sqr ) then
                 if ( random( ) .le. this%raw_fraction ) then
                       
                   lbuf = lbuf + 1
                   if (lbuf .gt. buffer_dim) then
                      write(file_id_msg,*) 'ERROR: temp buffer is full in write_raw'
                      exit
                    endif
                     
                    do i=1, p_x_dim
                       buffer(i,lbuf) = x(i,l)
                    enddo
                    do i=1, p_p_dim
                       buffer(i+p_x_dim,lbuf) = p(i,l)
                    enddo

! p.E                  
                    do i=1, p_pE_dim
                       buffer(i+p_x_dim+p_p_dim,lbuf) = pE(i,l)
                    enddo 
                    
                       buffer(1+p_x_dim+p_p_dim+p_pE_dim,lbuf) = q(l) 
                       
                    do i=1, p_x_dim
                       buffer(i+p_x_dim+p_p_dim+p_pE_dim+1,lbuf) = par_xid(i,l)
                    enddo 

                 endif                   
! p.E 
              endif              
              l = l+1
            enddo
             
          end select
         
          num_par_buffer = lbuf
          dim = 0
          bnd = 0

          if (my_aid(no_co) .eq. 1) then 
            
            ! Create the file   

            ! Prepare path and file name
            izero =  ichar('0')
            ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &               // char( izero + mod( sp_id/ 1 , 10 ) ) 

          
            path  = trim(path_mass) // 'RAW' // p_dir_sep &
     &           // trim(ch_sp_id)  // p_dir_sep 
          
            
            file_name = get_filename( n_aux, ''  &
     &                          // 'RAW'  // '-' // trim(ch_sp_id) , &
     &                          '.hdf' )
     
            full_name = trim(path) // trim(file_name)
             
            ! Create directory if necessary
            call mkdir( path, ierr )
 
           
            ! Create the HDF file
            sdfileID = sfstart(full_name, DFACC_CREATE)
            ! Add simulation time to the file
            ierr = sfsnatt(sdfileID, 'TIME', DFNT_FLOAT64, 1, t) 
            ! Add simulation timestep to the file 
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe_par_xid:"time"',t
               write(*,*) 'abording'
               call abort_program()
            endif
            ierr = sfsnatt(sdfileID, 'ITER', DFNT_INT32, 1, n)
            ! Add node number information to the file
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe_par_xid:"iter"',n
               write(*,*) 'abording'
               call abort_program()
            endif
            ierr = sfsnatt(sdfileID, 'NODE_NUMBER', DFNT_INT32, 1, my_aid(no_co))
!            write(*,*) 'No problem in write_RAW_pe_par_xid:"node number"',ierr, my_aid(no_co)
            if ( ierr .ne. 0 ) then
               write(*,*) 'Dumping problem in write_RAW_pe:"node number"',my_aid(no_co)
               write(*,*) 'abording'
               call abort_program()
            endif
            ! Add node configuration information to the file
            node_conf = my_nx(no_co)
            dim = p_x_dim
            ierr = sfsnatt(sdfileID, 'NODE_CONFIGURATION', DFNT_INT32, dim, node_conf )
            dim = 0
            ! Add gamma limit information to the file
            ierr = sfsnatt(sdfileID, 'GAMMA_LIMIT', DFNT_FLOAT64, 1, this%raw_gamma_limit )
            ! Add particle fraction information to the file
            ierr = sfsnatt(sdfileID, 'PARTICLE_FRACTION', DFNT_FLOAT64, 1, &
                           this%raw_fraction )          
         
            
            ! Create species dataset
! p.E 
            p_pE_dim = p_p_dim
            dimsize(1) = p_x_dim + p_p_dim + p_pE_dim + 1 + p_x_dim
! p.E 
            dimsize(2) = SD_UNLIMITED
  
 
            sdname = 'Species'// char( ichar('0') + mod( sp_id/ 10 , 10 ) ) &
     &                         // char( ichar('0') + mod( sp_id/  1 , 10 ) )

            sdsID = sfcreate(sdfileID, sdname, DFNT_FLOAT32 , 2, dimsize) 
            
! DEBUG
!             write(*,*) "Dataset created, sdsID ", sdsID
! DEBUG
            
            ierr = sfsblsz(sdsID, (p_x_dim + p_p_dim + p_pE_dim + 1 + p_x_dim)*1024*256)
! DEBUG
!            write(*,*) "Block size set to ", (p_x_dim + p_p_dim + p_pE_dim + 1 + p_x_dim)*1024*256
! DEBUG
            
            
            ! Add species ID to the dataset 
            ierr = sfsnatt(sdsID, 'SP_ID', DFNT_INT32, 1, sp_id)

            ! Add rqm to the dataset
            ierr = sfsnatt(sdsID, 'RQM', DFNT_FLOAT64, 1, rqm)

! DEBUG
!            write(*,*) "Information added, adding data from node 0"
! DEBUG
                     
            ! write node 0 data to file
            stride = 1
            hdfstart = 0
            edge(1) = dimsize(1)
          
! DEBUG
!               write(*,*) "write_raw_pE_par_xid, Number of particles on node 0 is ", num_par_buffer
! DEBUG
            if (num_par_buffer .gt. 0) then
               ! write data
               edge(2) = num_par_buffer
               ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
               if (ierr .ne. 0) write(*,*) "sfwdata failed for node 0 data"
               ! Update start position                   
               hdfstart(2) = hdfstart(2) + num_par_buffer
               !write(*,*) "data from node 0 added"
            endif
            
            
            
            ! clear buffer
            
            deallocate(buffer)
            
            if (no_num(no_co) .gt. 1) then            
                        
              ! get number of particles on other nodes
              allocate(num_par_buffer_all(no_num(no_co)))
              
              call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 
            
              ! loop through remaining nodes
            
              do i=2, no_num(no_co)             
! DEBUG
!               write(*,*) "write_raw_pE_par_xid, Number of particles on node ", i-1," is ", num_par_buffer_all(i)
! DEBUG
               if (num_par_buffer_all(i) .gt. 0) then
               
                 ! prepare buffer to receive data
! p.E             
                 allocate( buffer(p_x_dim + p_p_dim + p_pE_dim + 1 + p_x_dim,num_par_buffer_all(i)))
! p.E
              
               
                 ! get data from node i
                 tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,i), dim, bnd )
                 call send_ping( tid(no_co,i), tag )
                 
                 ! get data
                 source = tid(no_co, i)
                 tag    = tag_value( tag_code_spraw, tid(no_co, i), dim, bnd )               
! p.E
                 count = (1+p_x_dim+p_p_dim + p_pE_dim + p_x_dim) * num_par_buffer_all(i)
! p.E         
               
                 call MPI_RECV( buffer, count , MPI_REAL, source, &
     &                          tag, MPI_COMM_WORLD, stat, ierr )

                 ! save date from node i
               
                 edge(2) = num_par_buffer_all(i)
                 ierr = sfwdata(sdsID, hdfstart, stride, edge, buffer)
                 !if (ierr .ne. 0) write(*,*) "sfwdata failed for node ", i-1, " data"
                 ! Update start position                   
                 hdfstart(2) = hdfstart(2) + num_par_buffer_all(i)
               
                 ! clear buffer
                 deallocate(buffer)
               
! DEBUG
!                 write(*,*) "Data from node ",i-1, " added"
! DEBUG
               endif 
               
              enddo 
               
              deallocate(num_par_buffer_all) 
            
            endif
            
            ! Close the SDS
            ierr = sfendacc(sdsID)
                           
            ! Close HDF file
            ierr = sfend(sdfileID)
         
          else          
            ! send number of particles on local node
            allocate(num_par_buffer_all(no_num(no_co)))

            call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                         num_par_buffer_all, 1, MPI_INTEGER, &
                         tid(no_co,1), MPI_COMM_WORLD, ierr) 

            deallocate(num_par_buffer_all) 
            
            if (num_par_buffer .gt. 0) then             
              ! send particles on local node to node 1                          
              tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,my_aid(no_co)), dim, bnd )
              call recv_ping( tid(no_co, 1), tag )
           
              dest   = tid(no_co, 1)
              tag = tag_value( tag_code_spraw, tid(no_co,my_aid(no_co)), dim, bnd )
! p.E
              count = (1+p_x_dim+p_p_dim + p_pE_dim + p_x_dim) * num_par_buffer
! p.E
              
              call MPI_ISEND( buffer, count, MPI_REAL, dest, tag, &
     &                        MPI_COMM_WORLD, handle, ierr )
              call MPI_WAIT( handle, stat, ierr )
            endif
            
            ! clear buffer
            deallocate(buffer)
          endif
          

! DEBUG
!           write(*, *) ' '
!           write(*, *) 'now after write_raw'
! DEBUG
!          call file_flush(file_id_msg)

        end subroutine write_raw_pE_par_xid
!---------------------------------------------------
! p.E

!---------------------------------------------------
        subroutine acquire_range( this, num_par, p, no_co )
!---------------------------------------------------
!       acquires the ranges for autorange of p phasespace
!---------------------------------------------------

          implicit none
          
! dummy variables
           type( t_diag_species ),            intent( inout )  :: this
           integer(p_k_is),                   intent( in )     :: num_par
           real(p_k_rdoub),  dimension(:,:),  intent(in)       :: p
           type(t_node_conf),                 intent(in)       :: no_co

! local variables
!           real(p_k_rdoub),  dimension(p_p_dim)                :: p_min, p_max
           real(p_k_rdoub),  dimension(7)                      :: struc, struc_g
           integer(p_k_is)                                     :: i, l, j   !! datatype for Particle index
           logical                                             :: doit
!           real(p_k_rdoub)                                     :: mean
! SP
!          real(p_k_rdoub),  dimension(7,no_num(no_co))        :: struc_n
           real(p_k_rdoub),  dimension(:,:),allocatable        :: struc_n
!           logical,          dimension(no_num(no_co))          :: vailid_n
           real(p_k_rdoub),  dimension(3)                      :: dx
           integer                                             :: sendcount, recvcount, ierr, root 


! executable statements

         allocate(struc_n(7,no_num(no_co)))
         ! test if autorange is on
           doit = .false.
           do i=1, p_p_dim
             doit = doit .or. this%if_ps_p_auto(i)
           enddo
           if ( doit ) then 
           ! find extremas on each node 
             do i=1, p_p_dim
               if ( num_par .gt. 0 ) then
                 struc(7) = num_par
                 struc(i) = p(i,1)
                 struc(i+3) = p(i,1)
               else
                 struc(7) = -1
               endif
               do l = 1, num_par
                 if ( p(i,l) .lt. struc(i) ) then
                   struc(i) =  p(i,l)
                 elseif ( p(i,l) .gt. struc(i+3) ) then
                   struc(i+3) =  p(i,l)
                 endif
               enddo
             enddo

           ! send extremas to node 0
             sendcount = size( struc )
             recvcount = sendcount
             root = 0
             if ( no_num(no_co) .gt. 1 ) then
               call MPI_GATHER( struc, sendcount, MPI_DOUBLE_PRECISION, &
                   &            struc_n, recvcount, MPI_DOUBLE_PRECISION, &
                   &            root, MPI_COMM_WORLD, ierr) 
             else
               struc_n(:,1) = struc(:)
             endif
             
           ! on node 0:
             if ( my_aid(no_co) .eq. 1 ) then
               struc_g(7) = -1.0_p_k_rdoub
             ! find global extremas
               do j=1, no_num(no_co)
                 if ( struc_n(7,j) .gt. 0_p_k_rdoub ) then
                   if( struc_g(7) .gt. 0.0_p_k_rdoub )then
                     do i=1, 3
                       if ( struc_n(i,j) .lt. struc_g(i) ) then 
                         struc_g(i) = struc_n(i,j)
                       endif
                       if ( struc_n(i+3,j) .gt. struc_g(i+3) ) then
                         struc_g(i+3) = struc_n(i+3,j)
                       endif
                     enddo
                   else
                     struc_g(:) = struc_n(:,j)
                   endif
                 endif
               enddo
             ! compare extremas with minimal range
               do i=1, 3
                 if ( struc_g(i) .gt. this%ps_pmin_r(i) ) then
                   struc_g(i) = this%ps_pmin_r(i)
                 endif
                 if ( struc_g(i+3) .lt. this%ps_pmax_r(i) ) then
                   struc_g(i+3) = this%ps_pmax_r(i)
                 endif
                 if ( .not. this%if_ps_p_auto(i) ) then
                   struc_g(i) = this%ps_pmin(i)       !????????????
                   struc_g(i+3) = this%ps_pmax(i)
                 endif
               enddo
             ! add extra cellwith to range
               dx = ( struc_g(4:6) - struc_g(1:3) ) / ( this%ps_np -2 )
!               write(*,*) "p: ", dx(1), struc_g(1), struc_g(4)
               struc_g(4:6) = struc_g(4:6) + 2.01*dx
             endif

           ! send ranges back to all nodes
             sendcount = size( struc_g )
             if ( no_num(no_co) .gt. 1 ) then
                call MPI_BCAST( struc_g, sendcount, MPI_DOUBLE_PRECISION, root, &
                    & MPI_COMM_WORLD, ierr )
!                write(*,*) struc_g(1), struc_g(1+3)
             endif
             
           ! save range to be used by write_ps_to_dx
             do i=1, 3
!               this%ps_pmin(i) = struc_g(i) 
!               this%ps_pmax(i) = struc_g(i+3)
               this%ps_pmin(i) = struc_g(i) 
               this%ps_pmax(i) = struc_g(i+3) 
             enddo             
           endif
           deallocate(struc_n)
           end subroutine acquire_range
!---------------------------------------------------

!---------------------------------------------------
        subroutine acquire_range_gamma( this, num_par, p, no_co )
!---------------------------------------------------
!       acquires the ranges for autorange of gamma phasespace
!---------------------------------------------------

          implicit none
          
! dummy variables
           type( t_diag_species ),            intent( inout )  :: this
           integer(p_k_is),                   intent( in )     :: num_par
           real(p_k_rdoub),  dimension(:,:),  intent(in)       :: p
           type(t_node_conf),                 intent(in)       :: no_co

! local variables

           double precision                                    :: gamma
           double precision,  dimension(3)                     :: struc, struc_g
           
           ! struc is arranged as follows:
           !   struc(1) : gamma min
           !   struc(2) : gamma max
           !   struc(3) : number of particles in node (-1 if no particles in node)
           
           integer(p_k_is)                                     :: i, l, j   !! datatype for Particle index

! SP
!          double precision,  dimension(3,no_num(no_co))       :: struc_n
           double precision,  dimension(:,:),allocatable       :: struc_n
! SP
           double precision                                    :: dx
           integer                                             :: sendcount, recvcount, ierr, root 


! executable statements

         ! test if autorange is on

          allocate(struc_n(3,no_num(no_co)))
          if ( this%if_ps_gamma_auto ) then 
           
            ! write(*,*) "finding gamma auto range"
           
            ! find gamma range on each node 

            if ( num_par .gt. 0 ) then
               gamma = sqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2+1)
               struc(1) = gamma
               struc(2) = gamma
               struc(3) = num_par
            else
               struc(3) = -1
            endif
            
            do l = 2, num_par
               gamma = sqrt(p(1,l)**2+p(2,l)**2+p(3,l)**2+1)
               if ( gamma .lt. struc(1) ) then
                 struc(1) =  gamma
               elseif ( gamma .gt. struc(2) ) then
                 struc(2) =  gamma
               endif
            enddo

            ! write(*,*) my_aid(no_co), "->[",struc(1),",",struc(2),"]"

           ! send extremas to node 0
             sendcount = size( struc )
             recvcount = sendcount
             root = 0
             if ( no_num(no_co) .gt. 1 ) then
               call MPI_GATHER( struc, sendcount, MPI_DOUBLE_PRECISION, &
                   &            struc_n, recvcount, MPI_DOUBLE_PRECISION, &
                   &            root, MPI_COMM_WORLD, ierr) 
             else
               struc_n(:,1) = struc(:)
             endif
             
           ! on node 0:
             if ( my_aid(no_co) .eq. 1 ) then
               struc_g(3) = -1
               
               ! find global extremas
               
               do j=1, no_num(no_co)
                 if ( struc_n(3,j) .gt. 0_p_k_rdoub ) then
                   if( struc_g(3) .gt. 0.0_p_k_rdoub )then 
                     if ( struc_n(1,j) .lt. struc_g(1) ) then 
                        struc_g(1) = struc_n(1,j)
                     endif
                     if ( struc_n(2,j) .gt. struc_g(2) ) then
                        struc_g(2) = struc_n(2,j)
                     endif
                   else
                     struc_g(:) = struc_n(:,j)
                   endif
                 endif
               enddo
             
            ! compare extremas with minimal range
               if ( struc_g(1) .gt. this%ps_gammamin ) then
                 struc_g(1) = this%ps_gammamin
               endif
               if ( struc_g(2) .lt. this%ps_gammamax ) then
                 struc_g(2) = this%ps_gammamax
               endif

             ! add extra cellwith to range
               dx = ( struc_g(2) - struc_g(1) ) / ( this%ps_ngamma -2 )
               struc_g(2) = struc_g(2) + 2.01*dx
             endif

           ! send ranges back to all nodes
             sendcount = size( struc_g )
             if ( no_num(no_co) .gt. 1 ) then
                call MPI_BCAST( struc_g, sendcount, MPI_DOUBLE_PRECISION, root, &
                    & MPI_COMM_WORLD, ierr )
             endif
             
           ! save range to be used by write_1D_gamma_phasespace_to_dx

             this%ps_gammamin = struc_g(1) 
             this%ps_gammamax = struc_g(2) 

           endif
           
           deallocate(struc_n)
           
           end subroutine acquire_range_gamma
!---------------------------------------------------


! OS 2.0
!---------------------------------------------------
		function if_report_charge_dspec( this, n, ndump )
!---------------------------------------------------
!
!---------------------------------------------------
		   
		   implicit none
		   
		   logical :: if_report_charge_dspec
		   type( t_diag_species ), intent(in) :: this
		   integer(p_k_is),        intent(in) :: n, ndump
		   
		   if_report_charge_dspec = test_if_report( n, ndump, this%ndump_fac_charge )
		
		end function if_report_charge_dspec
!---------------------------------------------------
! OS 2.0
                   subroutine write_spectrum( p, q, num_par, sp_id, &
                                             this, no_co, nIter, t, n_aux, &
                                             dx, coordinates ) 
  implicit none

  ! passed variables
          real(p_k_rdoub), dimension(:,:),intent(in)    :: p
          real(p_k_rdoub), dimension(:),     intent(in) :: q
          integer(p_k_is),                   intent(in) :: num_par
          integer(p_k_is),                   intent(in) :: sp_id
          type( t_diag_species ),            intent(in) :: this
          type( t_node_conf ),               intent(in) :: no_co
          integer(p_k_is),                   intent(in) :: nIter
          real(p_k_rdoub),                   intent(in) :: t
          integer(p_k_is),                   intent(in) :: n_aux
          real(p_k_rdoub), dimension(:),     intent(in) :: dx
          integer(p_k_is),                   intent(in) :: coordinates

  ! local variables
          character(80)   :: psPath
          character(80)   :: full_name, file_name
          character(20)   :: chindx
          character(2)    :: ch_sp_id
          integer(p_k_is) :: me, izero


  integer,parameter::nNumBins=300;
  character*20::psIter
  integer::ios,i,ierr
  integer::mCount,indx
  logical::bpid
  integer::nSize
  integer::nNoPtsWrote
  real*8::gvsq,gami
  real::fGamma_m1_Max0,fGamma_m1_Max,fGamma_m1;
  real::finterval;
!  real(p_k_rdoub)::fSpectrum(nNumBins,11);
!  real(p_k_rdoub)::fSpectrum_out(nNumBins,11);
  real(p_k_rdoub),dimension(:,:),allocatable::fSpectrum;
  real(p_k_rdoub),dimension(:,:),allocatable::fSpectrum_out;
  integer::mycomm
  integer::nTempFile


  allocate(fSpectrum(nNumBins,11));
  allocate(fSpectrum_out(nNumBins,11));

  write(psIter,"(I7.7)") nIter

  me = my_aid(no_co)
!  print*,me,sp_id,'write_spectrum',nIter
  nTempFile=file_id_tem
  
  
            ! Prepare path and file name
            izero =  ichar('0')
            ch_sp_id =  char( izero + mod( sp_id/10 , 10 ) ) &
     &               // char( izero + mod( sp_id/ 1 , 10 ) ) 

          
!            psPath  = trim(path_hist) // 'SPECTRUM' // p_dir_sep &
!     &           // trim(ch_sp_id)  // p_dir_sep 
          
            psPath  = trim(path_hist)
            file_name = "spcs"//trim(ch_sp_id)//'Spectrum'//trim(adjustl(psIter))//'.txt'
     
            full_name = trim(psPath) // trim(file_name)
             
            ! Create directory if necessary
            call mkdir( psPath, ierr )
            


  fSpectrum=0.0;
  fGamma_m1_Max=60.;
  if(sp_id.ne.1) then
  fGamma_m1_Max0=0.;
  do mCount = 1,num_par
    gvsq=p(1,mCount)**2+p(2,mCount)**2+p(3,mCount)**2;
    gami=(1.d0/sqrt(1.d0+gvsq));                                              
    fGamma_m1=(gami+.5d0*gami*gvsq-1.d0) + .5d0*gami*gvsq;
!    fGamma=1./spcs%fGammaInv(mCount);
    if(fGamma_m1>fGamma_m1_Max0) then
      fGamma_m1_Max0=fGamma_m1;
    endif
  enddo
!now find max
  fGamma_m1_Max=fGamma_m1_Max0
  if(no_num(no_co).gt.1) then
    call MPI_ALLREDUCE(fGamma_m1_Max0,fGamma_m1_Max,1,MPI_REAL,MPI_MAX, MPI_COMM_WORLD,iErr)
    if(iErr/=MPI_SUCCESS) then
      write(file_id_msg,*)'AllReduce error finding max in writeSpectrum iErr =',iErr,',myid =',me
      return
    endif
  endif
  endif
  if(fGamma_m1_Max>0.0) then
    finterval=real(nNumBins-1)/(fGamma_m1_Max);
  else
    finterval=0.0
  endif
  if(finterval<0.0) finterval=0.0
  do mCount = 1,num_par
    gvsq=p(1,mCount)**2+p(2,mCount)**2+p(3,mCount)**2;
    gami=(1.d0/sqrt(1.d0+gvsq));                                              
    fGamma_m1=(gami+.5d0*gami*gvsq-1.d0) + .5d0*gami*gvsq;
    indx=1+int((fGamma_m1)*finterval);
!    indx=int((1./spcs%fGammaInv(mCount)-1.)*finterval);
    if(indx>nNumBins) then
      print*,me,' error in writeSpectrum (indx high) indx = ',indx,' nNumBins = ',nNumBins,'finterval',finterval,'fGamma_m1',fGamma_m1,'fGamma_m1_Max',fGamma_m1_Max,' nIter = ',nIter,' fixing'
      write(file_id_msg,*)'error in writeSpectrum (indx high) indx = ',indx,' nNumBins = ',nNumBins,' nIter = ',nIter,' fixing'
      indx = nNumBins;
    endif
    if(indx<1) then
      print*,me,' error in writeSpectrum (indx low) indx = ',indx,' nNumBins = ',nNumBins,'finterval',finterval,'fGamma_m1',fGamma_m1,'fGamma_m1_Max',fGamma_m1_Max,' nIter = ',nIter,' fixing'
      write(file_id_msg,*)'error in writeSpectrum (indx low) indx = ',indx,' nNumBins = ',nNumBins,' nIter = ',nIter,' fixing'
      indx = 1;
    endif
    fSpectrum(indx,1)=fSpectrum(indx,1)+abs(q(mCount));
    fSpectrum(indx,2)=fSpectrum(indx,2)+abs(q(mCount))*fGamma_m1;
    fSpectrum(indx,3)=fSpectrum(indx,3)+abs(q(mCount))*p(1,mCount);
    fSpectrum(indx,4)=fSpectrum(indx,4)+abs(q(mCount))*p(2,mCount);
    fSpectrum(indx,5)=fSpectrum(indx,5)+abs(q(mCount))*p(3,mCount);
    fSpectrum(indx,6)=fSpectrum(indx,6)+abs(q(mCount))*p(1,mCount)*gami;
    fSpectrum(indx,7)=fSpectrum(indx,7)+abs(q(mCount))*p(2,mCount)*gami;
    fSpectrum(indx,8)=fSpectrum(indx,8)+abs(q(mCount))*p(3,mCount)*gami;
    fSpectrum(indx,9)=fSpectrum(indx,9)+abs(q(mCount))*fGamma_m1*p(1,mCount)*gami;
    fSpectrum(indx,10)=fSpectrum(indx,10)+abs(q(mCount))*fGamma_m1*p(2,mCount)*gami;
    fSpectrum(indx,11)=fSpectrum(indx,11)+abs(q(mCount))*fGamma_m1*p(3,mCount)*gami;
  enddo
! now sum
  if(no_num(no_co).gt.1) then
    fSpectrum_out=fSpectrum
    call MPI_REDUCE(fSpectrum_out,fSpectrum,11*nNumBins,MPI_REAL,MPI_SUM,0, MPI_COMM_WORLD,iErr)
    if(iErr/=MPI_SUCCESS) then
      write(file_id_msg,*)'Reduce error summing in writeSpectrum iErr =',iErr,',myid =',me
      return
    endif
  endif
  if(me.eq.1) then
    open(nTempFile,file = full_name ,status='new',iostat=ios);
    write(nTempFile,*) ' gm1   n   nKE/Mp   np1/Mp   np2/Mp   np3/Mp  nv1   nv2   nv3   nQ1/Mp   nQ2/Mp   nQ3/Mp'
    do mCount=1,nNumBins
      write(nTempFile,*) (fGamma_m1_Max)*(real(mCount)/real(nNumBins)),fSpectrum(mCount,:)
    enddo 
    close(nTempFile);
  endif
  deallocate(fSpectrum);
  deallocate(fSpectrum_out);

end subroutine write_spectrum

!-------------------------------------------------------------------------------
! New version
!   Same memory usage, each quantity is saved in a different dataset. Before
!   getting data from other nodes we gather number of particles from all nodes
!   so we can determine exact buffer size. This means no more SD_UNLIMITED, which
!   results in significant I/O performance improvement.
!
!-------------------------------------------------------------------------------
subroutine write_raw_new( spec, no_co, g_space, n, t, n_aux, coordinates )
!-------------------------------------------------------------------------------

  use m_units
  use m_fparser

  implicit none
  
!  integer, parameter :: p_diag_prec = p_k_rsing
  integer, parameter :: p_diag_prec = p_k_rdoub
  
  
  ! dummy variables
  type( t_species ), intent(inout) :: spec
  type( t_node_conf ),               intent(in) :: no_co
  type( t_space ),      intent(in) :: g_space    ! spatial information
  integer,                   intent(in) :: n
  real(p_k_rdoub),                   intent(in) :: t
  integer,                   intent(in) :: n_aux
  integer,                   intent(in) :: coordinates


!       local variables

  integer :: sdfileID               ! HDF SD file ID number
  
!       local variables

  character(80)   :: path
  character(80)   :: full_name, file_name

  real(p_k_rdoub) :: gamma_limit_sqr, gamma_par_sqr, raw_eval
  real(p_k_rdoub), dimension(p_p_dim + p_x_dim + 2) :: raw_var
  logical :: has_raw_math_expr
  integer :: i, j, l, lbuf,iCnt

  ! Buffer to hold particle data
  real(p_diag_prec), allocatable, dimension(:,:) :: buffer
  real(p_diag_prec), allocatable, dimension(:) :: write_buffer

  integer, allocatable, dimension(:,:) :: buffer_tag
  
  integer :: num_par_buffer, buffer_dim

  real(p_k_rdoub), dimension(p_max_dim) :: box_limits
  integer, dimension(p_max_dim) :: lperiodic, lmove_c

  integer::npE,npar_xid,nbase_x,nbase_p,nbase_q,nbase_pE,nbase_par_xid,nbuffer_size
!        particles written to the buffer on all nodes         
  integer, allocatable, dimension(:) :: num_par_buffer_all

!        local variables for HDF          
  character(20) :: sdname          
  integer       :: sdsTagID, hdf_type
  integer, dimension(p_x_dim + p_p_dim + p_x_dim + p_p_dim + 1) :: sdsDataID

  integer, dimension(1) :: hdfstart, stride, edge, dimsize
  integer, dimension(2) :: hdfstartTag, strideTag, edgeTag, dimsizeTag
  
  ! local variables for MPI
  integer :: handle, mpi_type
  integer, dimension(mpi_status_size):: stat
  integer :: count
  integer :: dest, source, tag
  integer :: dim, bnd

  ! Additional file data
  integer, dimension(p_x_dim) :: node_conf 

  integer :: ierr
 
  ! set the proper types
  select case (p_diag_prec)
    case (p_k_rsing)
       mpi_type = MPI_REAL
       hdf_type = DFNT_FLOAT32
    case (p_k_rdoub)
       mpi_type = MPI_DOUBLE_PRECISION
       hdf_type = DFNT_FLOAT64
  end select
  

  ! initialize the buffers
  buffer_dim = spec%num_par
  num_par_buffer = 0
  nbuffer_size=p_x_dim + p_p_dim + 1
  nbase_x=0;
  nbase_p=nbase_x+p_x_dim;
  nbase_q=nbase_p+p_p_dim;
  nbase_pE=1+nbase_q;
  if(spec%if_pE) then
    npE=p_p_dim;
    nbuffer_size=nbuffer_size+p_p_dim
  else
    npE=0;
  endif
  nbase_par_xid=nbase_pE+npE
  if(spec%if_par_xid) then
    npar_xid=p_x_dim;
    nbuffer_size=nbuffer_size+p_x_dim
  else
    npar_xid=0;
  endif
  ! buffer must be organized like this so that each quantity is contiguos in memory
  allocate( buffer(nbuffer_size, buffer_dim), stat = ierr)
  if (ierr /= 0) then 
	write(*,*) "Unable to allocate memory for diagnostic" 
	call abort_program( p_err_alloc )
  endif

  
  if ( spec%add_tag ) then
    allocate( buffer_tag(2, buffer_dim), stat = ierr)
    if (ierr /= 0) then 
	  write(*,*) "Unable to allocate memory for diagnostic" 
	  call abort_program( p_err_alloc )
    endif
  endif
    
  ! initialize gamma limit
  
  gamma_limit_sqr = spec%diag%raw_gamma_limit**2
  
  ! loop over all particles

  l = 1
  lbuf = 0

  ! Copy all the particles to the temp buffer

  ! time for raw_math evaluation
  raw_var(p_x_dim+p_p_dim + 2) = t              

  has_raw_math_expr = (spec%diag%raw_math_expr /= '')

  select case ( coordinates )

          ! cylindrical coordinates not implemented yet

   case default
    
   do                                
      if (l > spec%num_par)  exit             
      gamma_par_sqr = 1.d0 + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2
      if ( gamma_par_sqr >= gamma_limit_sqr ) then
         
         raw_eval = 1.0_p_k_rdoub
         if (has_raw_math_expr) then
         
            ! fill evaluation variables
            raw_var(1:p_x_dim) = spec%x(1:p_x_dim,l)             ! x1-x3
            raw_var(p_x_dim+1) = spec%p(1,l)                     ! p1
            raw_var(p_x_dim+2) = spec%p(2,l)                     ! p2
            raw_var(p_x_dim+3) = spec%p(3,l)                     ! p3
            raw_var(p_x_dim+p_p_dim + 1) = sqrt(gamma_par_sqr)   ! g
            
            ! evaluate
!            raw_eval = eval( spec%diag%raw_func, real(raw_var,p_k_fparse) )
            raw_eval = eval( spec%diag%raw_func, real(raw_var,p_k_rsing) )
            
         endif
         
         if (raw_eval > 0.0_p_k_rdoub) then

!			   if ( genrand_real2() <= spec%diag%raw_fraction ) then
			   if ( random() <= spec%diag%raw_fraction ) then
					 
				  lbuf = lbuf + 1
				   
				  do i=1, p_x_dim
					 buffer(i+nbase_x, lbuf) = real(spec%x(i,l), p_diag_prec)
				  enddo
				  do i=1, p_p_dim
					 buffer(i+nbase_p, lbuf) = real(spec%p(i,l), p_diag_prec)
				  enddo
				  buffer(1+nbase_q, lbuf) = real(spec%q(l), p_diag_prec)
				  
				  if ( spec%add_tag ) then
					 buffer_tag(1,lbuf) = spec%tag(1,l)
					 buffer_tag(2,lbuf) = spec%tag(2,l)
				  endif
				  if ( spec%if_pE ) then
				    do iCnt=1,p_p_dim
				      buffer(nbase_pE+iCnt,lbuf) = spec%pE(iCnt,l)
                                    enddo
				  endif
				  if ( spec%if_par_xid ) then
				    do iCnt=1,p_x_dim
				      buffer(nbase_par_xid+iCnt,lbuf) = spec%par_xid(iCnt,l)
                                    enddo
				  endif
				  
			   endif  

	      endif
      endif              
      l = l+1
    enddo
    
  end select
 
  num_par_buffer = lbuf
  

  ! get total number of particles
  if (no_num(no_co) > 1) then            
                
      ! get number of particles on other nodes
      allocate(num_par_buffer_all(no_num(no_co)))
      
      call MPI_GATHER(num_par_buffer, 1, MPI_INTEGER, &
                 num_par_buffer_all, 1, MPI_INTEGER, &
                 tid(no_co,1), comm( no_co ), ierr) 
                 
  endif

  dim = 0
  bnd = 0

  if (my_aid(no_co) == 1) then ! root node
    
    ! Create the file   

    ! Prepare path and file name
    path  = trim(path_mass) // 'RAW' // p_dir_sep &
            // subst_spc(trim(spec%name))  // p_dir_sep 
  
    
    file_name = get_filename( n_aux, ''  &
                           // 'RAW'  // '-' // subst_spc(trim(spec%name)) , &
                           '.hdf' )

    full_name = trim(path) // trim(file_name)
     
    ! Create directory if necessary
    call mkdir( path, ierr )

   
    ! Create the HDF file
    sdfileID = sfstart(full_name, DFACC_CREATE)
    ! Add simulation time to the file
    ierr = sfsnatt(sdfileID, 'TIME', DFNT_FLOAT64, 1, t) 
    ! Add simulation timestep to the file 
    ierr = sfsnatt(sdfileID, 'ITER', DFNT_INT32, 1, n)
    
    ! Simulation box information
    box_limits(1:p_x_dim) = xmin( g_space )
    ierr = sfsnatt(sdfileID, 'XMIN', DFNT_FLOAT64, p_x_dim, box_limits) 
    box_limits(1:p_x_dim) = xmax( g_space )
    ierr = sfsnatt(sdfileID, 'XMAX', DFNT_FLOAT64, p_x_dim, box_limits) 
    
    
	! write periodic boundaries information
	do i = 1, p_x_dim
	  if (periodic( no_co, i )) then 
		 lperiodic(i) = 1
	  else
		 lperiodic(i) = 0
	  endif
	enddo
	ierr = sfsnatt(sdfileID, 'PERIODIC', DFNT_INT32, &
				   p_x_dim, lperiodic) 
	
	! write moving window information
	do i = 1, p_x_dim
	  if (if_move( g_space, i ) ) then 
		 lmove_c(i) = 1
	  else
		 lmove_c(i) = 0
	  endif
	enddo
	ierr = sfsnatt(sdfileID, 'MOVE C', DFNT_INT32, &
				   p_x_dim, lmove_c) 

    
    
    ! Add node number information to the file
    ierr = sfsnatt(sdfileID, 'NODE NUMBER', DFNT_INT32, 1, my_aid(no_co))
    ! Add node configuration information to the file
    node_conf = nx(no_co)
    dim = p_x_dim
    ierr = sfsnatt(sdfileID, 'NODE CONFIGURATION', DFNT_INT32, dim, node_conf )
    dim = 0
    ! Add gamma limit information to the file
    ierr = sfsnatt(sdfileID, 'GAMMA LIMIT', DFNT_FLOAT64, 1, spec%diag%raw_gamma_limit )
    ! Add particle fraction information to the file
    ierr = sfsnatt(sdfileID, 'PARTICLE FRACTION', DFNT_FLOAT64, 1, &
                   spec%diag%raw_fraction )          
    ! Add species ID to the file 
    ierr = sfsnatt(sdfileID, 'SP_ID', DFNT_INT32, 1, spec%sp_id)

    ! Add species Name to the file
    ierr = sfsnatt(sdfileID, 'NAME', DFNT_CHAR, len_trim(spec%name), &
                             trim(spec%name) )

    ! Add rqm to the file
    ierr = sfsnatt(sdfileID, 'RQM', DFNT_FLOAT64, 1, spec%rqm)
 
    ! get total number of particles to save
    dimsize = num_par_buffer
    do i = 2, no_num(no_co)
      dimsize = dimsize + num_par_buffer_all(i)
    enddo
    
    ! Create species datasets
    do i=1, p_x_dim
      sdname = 'x'//char(ichar('0')+i)
      sdsDataID(i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize)
!      ierr = sfsdtstr(sdsDataID(i), 'x_'//char(ichar('0')+i), 'c / !Mw!D0!N', 'F5.2', 'cartesian')
      ierr = sfsdtstr(sdsDataID(i), 'x_'//char(ichar('0')+i), 'c / \omega_0', 'F5.2', 'cartesian')
	  if (ierr /= 0) then 
		write(*,*) "Failed setting attributes" 
		call abort_program( p_err_diagfile )
	  endif
      
    enddo

    do i=1, p_p_dim
      sdname = 'p'//char(ichar('0')+i)
      sdsDataID(p_x_dim+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
      ierr = sfsdtstr(sdsDataID(p_x_dim+i), 'p_'//char(ichar('0')+i), &
                      'm_e c', 'F5.2', 'cartesian')
	  if (ierr /= 0) then 
		write(*,*) "Failed setting attributes" 
		call abort_program( p_err_diagfile )
	  endif
    enddo

    sdsDataID(p_x_dim+p_p_dim+1) = sfcreate(sdfileID, 'q', hdf_type,1, dimsize) 
    ierr = sfsdtstr( sdsDataID(p_x_dim+p_p_dim+1), 'q', 'e', 'F5.2', 'cartesian')
	if (ierr /= 0) then 
	  write(*,*) "Failed setting attributes" 
	  call abort_program( p_err_diagfile )
	endif
    
   
    ! write node 0 data to file
    stride = 1
    hdfstart = 0
    hdfstartTag = 0
  
    if (num_par_buffer > 0) then
       allocate( write_buffer(num_par_buffer) )
       ! write data
       edge = num_par_buffer
       do j = 1, p_x_dim+p_p_dim+1
		  write_buffer = buffer(j,1:num_par_buffer)
		  ierr = sfwdata(sdsDataID(j), hdfstart, stride, edge, write_buffer)
		  if (ierr /= 0) then 
			write(*,*) "Failed writing file" 
			call abort_program( p_err_diagfile )
		  endif
       enddo
       deallocate( write_buffer )
    endif
    ! write pE information if needed
    if ( spec%if_pE ) then
      do i=1, p_p_dim
        sdname = 'pE'//char(ichar('0')+i)
        sdsDataID(nbase_pE+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
        ierr = sfsdtstr(sdsDataID(nbase_pE+i), 'pE_'//char(ichar('0')+i), &
!                        'm_e c^2 !Mw!D0!N', 'F5.2', 'cartesian')
                        'm_e c^2 \omega_0', 'F5.2', 'cartesian')
	    if (ierr /= 0) then 
	      write(*,*) "Failed setting attributes" 
	      call abort_program( p_err_diagfile )
	   endif
      enddo
      ! write node 0 data to file
  
      if (num_par_buffer > 0) then
        allocate( write_buffer(num_par_buffer) )
        ! write data
        edge = num_par_buffer
        do j = 1, p_p_dim
		  write_buffer = buffer(nbase_pE+j,1:num_par_buffer)
		  ierr = sfwdata(sdsDataID(nbase_pE+j), hdfstart, stride, edge, write_buffer)
		  if (ierr /= 0) then 
			write(*,*) "Failed writing file" 
			call abort_program( p_err_diagfile )
		  endif
        enddo
        deallocate( write_buffer )
      endif
    endif
    ! write par_xid information if needed
    if ( spec%if_par_xid ) then
      do i=1, p_x_dim
        sdname = 'par_xid'//char(ichar('0')+i)
        sdsDataID(nbase_par_xid+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
        ierr = sfsdtstr(sdsDataID(nbase_par_xid+i), 'par_xid_'//char(ichar('0')+i), &
!                        'c / !Mw!D0!N', 'F5.2', 'cartesian')
                        'c / \omega_0', 'F5.2', 'cartesian')
	    if (ierr /= 0) then 
	      write(*,*) "Failed setting attributes" 
	      call abort_program( p_err_diagfile )
	   endif
      enddo
      ! write node 0 data to file
  
      if (num_par_buffer > 0) then
        allocate( write_buffer(num_par_buffer) )
        ! write data
        edge = num_par_buffer
        do j = 1, p_x_dim
		  write_buffer = buffer(nbase_par_xid+j,1:num_par_buffer)
		  ierr = sfwdata(sdsdataID(nbase_par_xid+j), hdfstart, stride, edge, write_buffer)
		  if (ierr /= 0) then 
			write(*,*) "Failed writing file" 
			call abort_program( p_err_diagfile )
		  endif
        enddo
        deallocate( write_buffer )
      endif

    endif
	   ! Update start position                   
    hdfstart = hdfstart + num_par_buffer
    
    ! clear buffer
    deallocate(buffer)

    
    ! write tag information if needed
    if ( spec%add_tag ) then
	   dimsizeTag(1) = 2
	   dimsizeTag(2) = dimsize(1)
   
	   sdname = 'tag'
	   
	   sdsTagID = sfcreate(sdfileID, sdname, DFNT_INT32 , 2, dimsizeTag) 
						   
	   ! write node 0 data to file
	   strideTag = 1
	   edgeTag(1) = dimsizeTag(1)

	   if ( num_par_buffer > 0 ) then
		  edgeTag(2) = num_par_buffer
   
		  ierr = sfwdata(sdsTagID, hdfstartTag, strideTag, edgeTag, buffer_tag)
		  if (ierr /= 0) then 
			write(*,*) "Failed writing file" 
			call abort_program( p_err_diagfile )
		  endif
	   
	   endif
	   
	   deallocate( buffer_tag )
	   hdfstartTag(2) = hdfstartTag(2) + num_par_buffer
    endif
    
    if (no_num(no_co) > 1) then            
    
      ! loop through remaining nodes
    
      do i=2, no_num(no_co)             

		 if (num_par_buffer_all(i) > 0) then
			
			! prepare buffer to receive data
			allocate( buffer(nbuffer_size, num_par_buffer_all(i)), stat = ierr)
			if (ierr /= 0) then 
			  write(*,*) "Unable to allocate memory for diagnostic" 
			  call abort_program( p_err_alloc )
			endif
		  
			! get data from node i
			tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,i), dim, bnd )
			call send_ping( no_co, tid(no_co,i), tag )
			
			! get data
			source = tid(no_co, i)
			tag    = tag_value( tag_code_spraw, tid(no_co, i), dim, bnd )               
			count = num_par_buffer_all(i) * nbuffer_size 
		  
			call MPI_RECV( buffer, count , mpi_type, source, &
							  tag, comm( no_co ), stat, ierr )
   			
			! write data
			edge = num_par_buffer_all(i)
			
			allocate( write_buffer(num_par_buffer_all(i)) )
			do j = 1, p_x_dim+p_p_dim+1
                          write_buffer = buffer(j,1:num_par_buffer_all(i))
                          ierr = sfwdata(sdsDataID(j), hdfstart, stride, edge, write_buffer)
			   if (ierr /= 0) then 
				 write(*,*) "Failed writing file" 
				 call abort_program( p_err_diagfile )
			   endif
			enddo
                        if ( spec%if_pE ) then
			  do j = 1, npE
                            write_buffer = buffer(nbase_pE+j,1:num_par_buffer_all(i))
                            ierr = sfwdata(sdsDataID(nbase_pE+j), hdfstart, stride, edge, write_buffer)
			    if (ierr /= 0) then 
			      write(*,*) "Failed writing file" 
			      call abort_program( p_err_diagfile )
			    endif
			  enddo
                        endif
                        if ( spec%if_par_xid ) then
			  do j = 1, npar_xid
                            write_buffer = buffer(nbase_par_xid+j,1:num_par_buffer_all(i))
                            ierr = sfwdata(sdsDataID(nbase_par_xid+j), hdfstart, stride, edge, write_buffer)
			    if (ierr /= 0) then 
			      write(*,*) "Failed writing file" 
			      call abort_program( p_err_diagfile )
			    endif
			  enddo
                        endif

			deallocate( write_buffer )
			
			! clear buffer
			deallocate(buffer)
			
			! Update start position        
			hdfstart = hdfstart + num_par_buffer_all(i)
			
			if ( spec%add_tag ) then
			   ! prepare buffer to receive data
			   allocate( buffer_tag(2,num_par_buffer_all(i)))
			 
			   ! get data
			   source = tid(no_co, i)
			   tag    = tag_value( tag_code_spraw, tid(no_co, i), dim, bnd )               
			   count = 2 * num_par_buffer_all(i)
			 
			   call MPI_RECV( buffer_tag, count , MPI_INTEGER, source, &
								 tag, comm( no_co ), stat, ierr )
	  
			   ! save date from node i
			   edgeTag(2) = num_par_buffer_all(i)
			   ierr = sfwdata(sdsTagID, hdfstartTag, strideTag, edgeTag, buffer_tag)
			   if (ierr /= 0) then 
				  write(*,*) "Failed writing file" 
				  call abort_program( p_err_diagfile )
			   endif
			   
			   ! clear buffer
			   deallocate(buffer_tag)
			   
			   hdfstartTag(2) = hdfstartTag(2) + num_par_buffer_all(i)
   
			endif
  
		 endif 
       
      enddo 
       
      deallocate(num_par_buffer_all) 
    
    endif
    
    ! Close the SDS
    do j = 1, nbuffer_size
      ierr = sfendacc(sdsDataID(j))
    enddo
    
    if ( spec%add_tag ) then
      ierr = sfendacc(sdsTagID)
    endif
    
    ! Close HDF file
    ierr = sfend(sdfileID)
 
  else  ! other nodes        
    
    deallocate(num_par_buffer_all) 
    
    if (num_par_buffer > 0) then             
      ! send particles on local node to node 1                          
      tag = tag_value( tag_code_spraw+tag_code_ping, tid(no_co,my_aid(no_co)), dim, bnd )
      call recv_ping( no_co, tid(no_co, 1), tag )
   
      dest   = tid(no_co, 1)
      tag = tag_value( tag_code_spraw, tid(no_co,my_aid(no_co)), dim, bnd )
      count =  num_par_buffer * nbuffer_size
      
      call MPI_ISEND( buffer, count, mpi_type, dest, tag, &
&                        comm( no_co ), handle, ierr )
      call MPI_WAIT( handle, stat, ierr )
      
      if (spec%add_tag) then
         count = 2 * num_par_buffer
         
         call MPI_ISEND( buffer_tag, count, MPI_INTEGER, dest, tag, &
   &                        comm( no_co ), handle, ierr )
         call MPI_WAIT( handle, stat, ierr )
      endif
      
    endif
    
    ! clear buffer
    deallocate(buffer)
  endif


end subroutine write_raw_new
!-------------------------------------------------------------------------------


      end module m_diag_species


