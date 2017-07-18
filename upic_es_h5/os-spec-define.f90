!-------------------------------------------------------------------------------
! Species definition module
!
! This file contains the class definition for the following classes:
!
!  t_spe_bound
!  t_piston
!  t_species
!-------------------------------------------------------------------------------

#include "os-preprocess.fpp"

module m_species_define

use m_system
use m_parameters
use m_profile
use m_fparser
use m_vdf
use m_vdf_report
use m_time_avg

! minimal block size to be used when growing particle buffers
integer, parameter :: p_spec_buf_block = 131072 ! 128k

! constant for default relmax distribution cut-off (parameter x vdist_T)
integer, parameter :: p_relmax_umax = 20

! variable needed for math func profiles
  integer, parameter :: max_math_func_len = 1024


!-------------------------------------------------------------------------------
! Particle tracks diagnostics and particle track class
!-------------------------------------------------------------------------------
integer, parameter :: p_max_tracks = 16384

type t_track
  
  integer :: npoints = 0                ! number of points in memory
  integer :: savedpoints = 0            ! number of points saved to disk
  integer, dimension(2) :: tag = 0      ! tag of particle being tracked  
  integer :: part_idx = -1              ! index of particle being tracked (-1) means
                                        ! that particle is not on local node
  
  real(p_k_rdoub), dimension(:,:), pointer :: data => null()  ! track points data
  integer, dimension(:), pointer :: n => null()               ! track points iterations

end type t_track


type t_track_set

  ! maximum number of points to hold
  integer :: maxpoints = -1
  
  ! number of iterations between track points
  integer :: niter = 1			

  ! file holding tags of particles to follow 
  character(len = p_max_filename_len) :: file_tags = ''
  
  ! filename to write tracks to
  character(len = p_max_filename_len) :: file_write = ''
  
  ! total number of tracks 
  integer :: ntracks = 0                   
  ! actual tracks
  type( t_track ), dimension(:), pointer :: tracks  => NULL()

  ! list of tags present / missing
  ! each item is the index of the track on the tracks array
  integer, dimension(:), pointer :: present, missing 

  ! list sizes (for simplicity)
  integer :: npresent, nmissing

end type t_track_set

!-------------------------------------------------------------------------------
! Phasespace diagnostics and phasespace list Class
!-------------------------------------------------------------------------------
integer, parameter :: p_max_phasespace_dims = 3
integer, parameter :: p_max_phasespace_namelen = 32
integer, parameter :: p_max_n_ene_bins = 32

! quantity to deposit
integer, parameter :: p_pstype_none   = -1     ! not defined
integer, parameter :: p_pstype_normal = 0      ! charge
integer, parameter :: p_pstype_mass   = 1      ! mass
integer, parameter :: p_pstype_ene    = 2      ! kinetic energy
integer, parameter :: p_pstype_q1     = 3      ! heat flux along x1
integer, parameter :: p_pstype_q2     = 4      ! heat flux along x2
integer, parameter :: p_pstype_q3     = 5      ! heat flux along x3
integer, parameter :: p_pstype_abs    = 6      ! abs(charge)
integer, parameter :: p_pstype_j1     = 7      ! current along x1
integer, parameter :: p_pstype_j2     = 8      ! current along x2
integer, parameter :: p_pstype_j3     = 9      ! current along x3

character(len=*), parameter :: p_psext_normal = "charge"
character(len=*), parameter :: p_psext_mass   = "m"
character(len=*), parameter :: p_psext_ene    = "ene"
character(len=*), parameter :: p_psext_q1     = "q1"
character(len=*), parameter :: p_psext_q2     = "q2"
character(len=*), parameter :: p_psext_q3     = "q3"
character(len=*), parameter :: p_psext_abs    = "|charge|"
character(len=*), parameter :: p_psext_j1     = "j1"
character(len=*), parameter :: p_psext_j2     = "j2"
character(len=*), parameter :: p_psext_j3     = "j3"

type t_phasespace
   integer           :: ndims = -1
   character(len=p_max_phasespace_namelen) :: name = "-"
   
   ! quantity to use as the axis for the phasespace
   ! 1 -> position, 2 -> momenta, 3 -> gamma  4 -> log10(gamma)
   integer, dimension( p_max_phasespace_dims ) :: x_or_p = 0
   
   ! coordinate to use as the axis for the phasespace
   ! has no meaning if x_or_p = 3, 4
   integer, dimension( p_max_phasespace_dims ) :: xp_dim = 0
   type( t_time_avg ) :: tavg
   integer :: ps_type = p_pstype_none
   
   type( t_phasespace ), pointer :: next => null()
end type t_phasespace

type t_phasespace_list
   type( t_phasespace ), pointer :: head => null()
   type( t_phasespace ), pointer :: tail => null()
end type
 
!-------------------------------------------------------------------------------
! Species Diagnostics Class
!-------------------------------------------------------------------------------
type :: t_diag_species

  ! frequency of species dianostic data dumps
  integer :: ndump_fac_pha, ndump_fac_pha_tavg, ndump_fac_ene, ndump_fac_raw, &
             ndump_fac_charge, ndump_fac_charge_lineout 

  ! number of time steps to average over for time averaged 
  ! dumps
  integer :: n_tavg

  ! physical range for phasespace data dumps
  real(p_k_rdoub), dimension(p_x_dim) :: ps_xmin, ps_xmax
  real(p_k_rdoub), dimension(p_p_dim) :: ps_pmin, ps_pmax
  real(p_k_rdoub), dimension(p_p_dim) :: ps_pmin_r, ps_pmax_r

  real(p_k_rdoub), dimension(p_p_dim) :: ps_lmin, ps_lmax
  real(p_k_rdoub), dimension(p_p_dim) :: ps_lmin_r, ps_lmax_r
  
  ! switch for autorange of p for phasespaces
  logical, dimension(p_p_dim) :: if_ps_p_auto

  ! switch for autorange of l for phasespaces
  logical, dimension(p_p_dim) :: if_ps_l_auto

  ! physical range for 1D momenta phasespace data dumps
  real(p_k_rdoub) :: ps_gammamin, ps_gammamax
  real(p_k_rdoub) :: ps_gammamin_r, ps_gammamax_r

  ! switch for autorange of ps_gamma
  logical :: if_ps_gamma_auto

  ! switch for log deposition of gamma
  logical :: if_ps_gamma_log
  

  ! resolutions for phasespace data dumps
  integer(p_k_is), dimension(p_x_dim) :: ps_nx, ps_nx_3D
  integer(p_k_is), dimension(p_p_dim) :: ps_np, ps_np_3D
  integer(p_k_is), dimension(p_p_dim) :: ps_nl, ps_nl_3D

  ! resolution for 1D momenta phasespace data dump
  integer :: ps_ngamma

  ! deposition scheme for phase space calculations
  integer :: dep_sch
  
  ! energy binned diagnostics parameters
  integer :: n_ene_bins
  real(p_k_rsing), dimension(p_max_n_ene_bins) :: ene_bins  
  
  ! parameters for raw data dump
  real(p_k_rdoub) :: raw_gamma_limit
  real(p_k_rdoub) :: raw_fraction
  ! parameters for function parser
  character(len = max_math_func_len) :: raw_math_expr = " "
  type(t_fparser) :: raw_func
  
  ! switches to turn on and off other species diagnostics 
  logical :: iforigin, iftestpar ! currently not used

  ! particle tracking data
  integer :: ndump_fac_tracks           ! frequency of track diagnostic writes
  type( t_track_set ) :: tracks

  ! new phasespace input stuff
  type( t_phasespace_list ) :: phasespace_list
  type( t_phasespace_list ) :: pha_ene_bin_list
  type( t_phasespace_list ) :: pha_cell_avg_list
  type( t_phasespace_list ) :: pha_time_avg_list
  
  ! charge lineouts
  type( t_vdf_lineout_list ) :: charge_lineouts
  
  
end type t_diag_species

!-------------------------------------------------------------------------------
! Species Boundary Condition Class
!-------------------------------------------------------------------------------
type :: t_spe_bound

  ! Boundary condition type
  integer, dimension(2,p_max_dim) :: type

  ! thermal and fluid momenta of the thermal bath BC
  real(p_k_rdoub), dimension( p_p_dim,2,p_x_dim) :: pth_bnd
  real(p_k_rdoub), dimension( p_p_dim,2,p_x_dim) :: pfl_bnd

end type t_spe_bound

!-------------------------------------------------------------------------------
! Piston Class
!-------------------------------------------------------------------------------

type :: t_piston
  ! parameters from inputdeck (processed)
  integer          :: dim        ! dimension in which piston moves
  integer          :: updown     ! edge from whitch piston starts
  real(p_k_rdoub)  :: u          ! proper velocity ( gamma v) of piston
  real(p_k_rdoub)  :: start_pos  ! position from which piston starts
  real(p_k_rdoub)  :: start_time ! time when piston launched of edge
  real(p_k_rdoub)  :: stop_time  ! time when piston disapears
  real(p_k_rdoub)  :: opacity_factor  ! uniform opacity or factor for profile
  integer          :: profile_type    ! type of profile
  type(t_fparser)  :: profile_parser  ! profile function (A, B: transverse variables

  !auxiliary variabes
  real(p_k_rdoub)  :: v          ! v
  real(p_k_rdoub)  :: squ        ! u^2 = (gamma v)^2
  real(p_k_rdoub)  :: g          ! gamma
  real(p_k_rdoub)  :: sqg        ! gamma^2
  ! real(p_k_rdoub)  :: sqv        ! v^2

  real(p_k_rdoub)  :: pos_after  ! position of piston at t
  real(p_k_rdoub)  :: pos_before ! position of piston at t-dt

  logical          :: inbox      ! true if piston is in this node
end type t_piston


!-------------------------------------------------------------------------------
! Species Class
!-------------------------------------------------------------------------------
  
! max. length of species name
integer, parameter :: p_max_spname_len = 64

type :: t_species

  ! species name
  character(len = p_max_spname_len) :: name 
  
  ! species number - internal id
  integer :: sp_id 

  ! particles per cell
  integer, dimension(p_x_dim) :: num_par_x 
  
  ! cell size
  real(p_k_rdoub), dimension(p_x_dim) :: dx

  ! particle buffer size 
  integer :: num_par_max 
  
  ! number of particles in buffer
  integer :: num_par     
  
  ! number of particles that have been created in this node
  integer :: num_created = 0
  
  !  mass to charge ratio
  real(p_k_rdoub)                     :: rqm 

  ! essential particle data position, momentum, and charge
  real(p_k_rdoub), dimension(:,:), pointer :: x => null()
  real(p_k_rdoub), dimension(:,:), pointer :: p => null()
  real(p_k_rdoub), dimension(:),   pointer :: q => null()

  ! type of current deposition / field interpolation
  integer :: interpolation 
  
  ! velocity distribution parameters
  integer :: vdist_type                           ! distribution type
  real(p_k_rdoub), dimension(p_p_dim) :: vth      ! thermo-momentum
  real(p_k_rdoub), dimension(p_p_dim) :: vfl      ! fluid-momentum
  real(p_k_rdoub)                     :: vdist_T  ! relmax temperature
  real(p_k_rdoub)                     :: umax     ! relmax max temperature
  logical :: use_class_vadd = .false.             ! force classical addition for
                                                  ! thermal + fluid velocities
  
  ! particle tags
  logical :: add_tag = .false.
  integer, dimension(:,:), pointer :: tag => null()
  
  
  !density information
  real(p_k_rdoub)   :: den_min ! approx. min. density for a particle
  type( t_profile ) :: den     ! inital density profile 

  ! Free streaming species (i.e. constant velocity). dgam can still be used
  logical :: free_stream = .false.

  ! number of acceleration steps and increase in gamma at
  ! at each step if the species is used as beam (in direction x1)
  real(p_k_rdoub) :: dgam
  integer(p_k_is) :: num_dgam ! number of acceleration time steps

  ! boundary cond. for this species
  type( t_spe_bound )    :: bnd_con
  
  ! diagnostic for this species
  type( t_diag_species ) :: diag 

  ! number of timesteps between sorting of particles
  ! n_sort = 0 turns sorting off
  integer :: n_sort 

  ! delayed push
  real(p_k_rdoub)      :: push_start_time

  ! use subcycling on this species
  logical :: subcycle

  ! local E and B fields for n_push > 1
  ! averaged for n_push times
  type(t_vdf), pointer :: E => NULL()
  type(t_vdf), pointer :: B => NULL()
            
  
  ! Numerical piston data
  integer                                 :: num_pistons
  type( t_piston ), dimension(:), pointer :: pistons  => NULL()
  
   
  ! time-centered total energy diagnostic
  logical :: if_energy = .false.
  real( p_k_rdoub ), dimension( 1 + p_p_dim ) :: energy
   
end type t_species
!-------------------------------------------------------------------------------
  


end module m_species_define