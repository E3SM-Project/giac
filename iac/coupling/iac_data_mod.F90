module iac_data_mod

  !------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains the gcam coupler variable and structure definitions
  ! I'm putting them all in one place, both input, output, and
  ! internal, to make it easier to find, but we might split them later
  !------------------------------------------------------------
  use shr_kind_mod               , only : r8 => shr_kind_r8

  implicit none
  save

  ! control structure for internal use
  type, public :: iac_ctl_type

     ! Decomposition - just in case
     integer :: begg,endg  ! Begin and end grid cells for this processor

     ! Dimensioned by grid cell - ngrid
     integer, allocatable :: gindex(:)  ! mapping grid cells to global index
     integer, allocatable :: iacmask(:) ! mask

     integer, allocatable :: ilon(:)   ! Index of lat,lon dimension,
     integer, allocatable :: jlat(:)   ! for grid cell

     ! Dimensioned (lon) or (lat) index, respectively
     ! It's easier for gcam integration if we store all our variables
     ! multidimensionally in (lon,lat,pft) order, so the only thing
     ! on the global grid is the indeces to go back and forth between
     ! (g) <-> (i,j)
     real(r8), allocatable :: lon(:)    ! longitude
     real(r8), allocatable :: lat(:)    ! latitude 

     ! Various sizes
     integer, public :: ngrid
     integer, public :: npft
     integer, public :: nlat
     integer, public :: nlon

  end type iac_ctl_type

! ! cdata datatype for GCAM and GLM internal use
  type iac_cdata_type
    character(len=640),pointer :: c(:)
    real*8 , pointer :: r(:)
    integer, pointer :: i(:)
    logical, pointer :: l(:)
  end type iac_cdata_type

  ! Making these structures, since it's easier to make sure
  ! everything is contiguous, which is important when coupling with
  ! C++ 
  type, public :: lnd2iac_type
     ! Input from lnd, dimensioned (lon,lat,pft)
     real(r8), allocatable :: npp(:,:,:)
     real(r8), allocatable :: hr(:,:,:)
     real(r8), allocatable :: pftwgt(:,:,:)

     ! These will be read once from configuration file, and are (lon,lat)
     real(r8), allocatable :: landfrac(:,:)
     real(r8), allocatable :: area(:,:)
  end type lnd2iac_type

  type, public :: iac2lnd_type
     ! Output to lnd, dimensioned (lon,lat,pft)
     real(r8), allocatable :: pct_pft(:,:,:)
 end type iac2lnd_type

  type, public :: iac2atm_type
     ! output to atm, dimensioned (lon,lat)
     real(r8), allocatable :: co2emiss(:,:)
  end type iac2atm_type

  type(lnd2iac_type) :: lnd2iac_vars
  type(iac2lnd_type) :: iac2lnd_vars
  type(iac2atm_type) :: iac2atm_vars
  type(iac_ctl_type) :: iac_ctl

  ! This may be redundant, but the downstream modules use what
  ! looks like a different EClock interface.  Rather than hunt down
  ! everything there, we'll just copy waht we need here.
  integer, pointer :: GClock(:)

  ! The gcam functions from iESM have a somewhat different "cdata"
  ! structure, so we keep cdata_z as the E3SM version and use gdata
  ! for the internal gcam version.
  type(iac_cdata_type) :: gdata

  !------------------------------------------------------
  ! Original iac_cdata parameters.  Some, most, or all might be
  ! obsolete now, and we'll figure out which later
  !-----------------------------------------------------------------

  integer, parameter, public :: iac_eclock_size  =  8
  integer, parameter, public :: iac_eclock_ymd   =  1
  integer, parameter, public :: iac_eclock_tod   =  2
  integer, parameter, public :: iac_eclock_dt    =  3
  integer, parameter, public :: iac_eclock_AclmC =  4
  integer, parameter, public :: iac_eclock_Agcam =  5
  integer, parameter, public :: iac_eclock_Aglm  =  6
  integer, parameter, public :: iac_eclock_Agcamsetden =  7

  integer, parameter, public :: iac_cdata_size            = 32
  !--- characters ---
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2base2000 = 10
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid720x360 = 11
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid288x192 = 12
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2shipbase2000 = 13
  integer, parameter, public :: iac_cdatac_gcam2emisfile_lut720x360map = 14
  integer, parameter, public :: iac_cdatac_gcam2emisfile_downscaleinfo = 15
  integer, parameter, public :: iac_cdatac_gcam2emisfile_rcp45allsteps = 16

  !--- reals ---
  real*8,  parameter, public :: iac_spval = -999.0
  integer, parameter, public :: iac_gcam_timestep =  5
  integer, parameter, public :: iac_glm_nx  = 720
  integer, parameter, public :: iac_glm_ny  = 360
  integer, parameter, public :: iac_iac_npfts  = 16

  !--- integers ---
! KVC: need to fix gcam2emiss and then can remove this
  integer, parameter, public :: iac_cdatai_logunit         =  1 
! KVC: next 3 used in glm_comp_mod but not sure if needed
  integer, parameter, public :: iac_cdatai_glm_nx          =  5 
  integer, parameter, public :: iac_cdatai_glm_ny          =  6
  integer, parameter, public :: iac_cdatai_glm_size        =  7
  integer, parameter, public :: iac_cdatai_gcam_yr1        =  8
  integer, parameter, public :: iac_cdatai_gcam_yr2        =  9

  !--- logicals ---
! KVC: below is used in restart and other logic, but isn't working right yet
  integer, parameter, public :: iac_cdatal_rest            =  1
  integer, parameter, public :: iac_cdatal_glm_present     =  4
  integer, parameter, public :: iac_cdatal_glm_prognostic  =  5
  integer, parameter, public :: iac_cdatal_initrun         =  10

  ! indices used to unpack gcamo
  integer           , public :: iac_gcamo_crop = 6
  integer           , public :: iac_gcamo_pasture = 5
  integer           , public :: iac_gcamo_woodharv = 9
  integer           , public :: iac_gcamo_forest = 4

  ! maximum glus per GCAM region
  integer           , public :: iac_max_nglu = 37
  integer           , public :: iac_first_coupled_year = 2016
  integer           , public :: iac_start_year = 2015

  ! indices used for glm input
  integer           , public :: iac_glmi_nflds = 3
  integer           , public :: iac_glmi_natveg = 1
  integer           , public :: iac_glmi_cropland = 2
  integer           , public :: iac_glmi_pasture = 3

  ! number of fields for glmo
  integer           , public :: iac_glmo_nflds = 9

end module iac_data_mod


     
     
