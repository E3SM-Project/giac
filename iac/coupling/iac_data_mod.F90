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
     integer, allocatable :: ilon(:)   ! Index of lat,lon dimension,
     integer, allocatable :: jlat(:)   ! for grid cell

     ! Dimensioned (lon) or (lat) index, respectively
     ! It's easier for gcam integration if we store all our variables
     ! multidimensionally in (lon,lat,pft) order, so the only thing
     ! on the global grid is the indeces to go back and forth between
     ! (g) <-> (i,j)
     real(r8), allocatable :: lon(:)    ! longitude
     real(r8), allocatable :: lat(:)    ! latitude 

     ! dimensioned lon,lat for use by iac
     integer, allocatable :: iacmask(:,:) ! pft land mask from gridfile
     real(r8), allocatable :: landfrac(:,:) ! landfrac from gridfile
     real(r8), allocatable :: area(:,:) ! cell area from gridfile
     real(r8), allocatable :: vegfrac(:,:) ! veg land unit frac from gridfile

     ! Various sizes
     integer, public :: ngrid
     integer, public :: npft
     integer, public :: nharvest
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
  end type lnd2iac_type

  type, public :: iac2lnd_type
     ! Output to lnd, dimensioned (lon,lat,pft)
     real(r8), allocatable :: pct_pft(:,:,:)
     real(r8), allocatable :: pct_pft_prev(:,:,:)
     ! Output to lnd, dimensioned (lon, lat, harvest frac)
     real(r8), allocatable :: harvest_frac(:,:,:)
 end type iac2lnd_type

  type, public :: iac2atm_type
     ! output to atm, dimensioned (lon,lat)
     real(r8), allocatable :: co2sfc(:,:,:)
     real(r8), allocatable :: co2airlo(:,:,:)
     real(r8), allocatable :: co2airhi(:,:,:)
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
  integer, parameter, public :: iac_iac_npfts  = 50

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
  ! Non-crop land type indices (matching luc.xml output-data order)
  integer           , public :: iac_gcamo_shrubland = 1
  integer           , public :: iac_gcamo_grassland = 2
  integer           , public :: iac_gcamo_urban = 3
  integer           , public :: iac_gcamo_forest = 4
  integer           , public :: iac_gcamo_pasture = 5
  integer           , public :: iac_gcamo_barren = 6
  integer           , public :: iac_gcamo_tundra = 7
  ! Individual GCAM crop type indices (matching luc.xml output-data order)
  integer           , public :: iac_gcamo_cornc4 = 8
  integer           , public :: iac_gcamo_fibercrop = 9
  integer           , public :: iac_gcamo_foddergrass = 10
  integer           , public :: iac_gcamo_fodderherb = 11
  integer           , public :: iac_gcamo_fodderherbc4 = 12
  integer           , public :: iac_gcamo_fruits = 13
  integer           , public :: iac_gcamo_fruitstree = 14
  integer           , public :: iac_gcamo_legumes = 15
  integer           , public :: iac_gcamo_misccrop = 16
  integer           , public :: iac_gcamo_misccropc4 = 17
  integer           , public :: iac_gcamo_misccroptree = 18
  integer           , public :: iac_gcamo_nutsseeds = 19
  integer           , public :: iac_gcamo_nutsseedstree = 20
  integer           , public :: iac_gcamo_oilcrop = 21
  integer           , public :: iac_gcamo_oilcroptree = 22
  integer           , public :: iac_gcamo_oilpalmtree = 23
  integer           , public :: iac_gcamo_otherarableland = 24
  integer           , public :: iac_gcamo_othergrain = 25
  integer           , public :: iac_gcamo_othergrainc4 = 26
  integer           , public :: iac_gcamo_rice = 27
  integer           , public :: iac_gcamo_roottuber = 28
  integer           , public :: iac_gcamo_soybean = 29
  integer           , public :: iac_gcamo_sugarcrop = 30
  integer           , public :: iac_gcamo_sugarcropc4 = 31
  integer           , public :: iac_gcamo_vegetables = 32
  integer           , public :: iac_gcamo_wheat = 33
  integer           , public :: iac_gcamo_biomassgrass = 34
  integer           , public :: iac_gcamo_biomasstree = 35
  ! Wood harvest is appended after all land types
  integer           , public :: iac_gcamo_woodharv = 36
  ! First and last crop indices for computing total crop area
  integer           , public :: iac_gcamo_crop_first = 8
  integer           , public :: iac_gcamo_crop_last = 35
  ! Number of GCAM crop types
  integer, parameter, public :: iac_num_gcam_crops = 28

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


     
     
