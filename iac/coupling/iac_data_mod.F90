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
  integer, parameter, public :: iac_cdatac_casename        =  1
  integer, parameter, public :: iac_cdatac_clm2gcam        =  2
  integer, parameter, public :: iac_cdatac_ibclmfile       =  3
  integer, parameter, public :: iac_cdatac_clmcbfndir      =  4
  integer, parameter, public :: iac_cdatac_gcam2glm_basecrop= 5
  integer, parameter, public :: iac_cdatac_gcam2glm_basepast= 6
  integer, parameter, public :: iac_cdatac_gcam2glm_baseothr= 7
  integer, parameter, public :: iac_cdatac_gcam2glm_aezmap =  8
  integer, parameter, public :: iac_cdatac_gcam2glm_basebiomass= 9
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2base2000 = 10
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid720x360 = 11
  integer, parameter, public :: iac_cdatac_gcam2emisfile_grid288x192 = 12
  integer, parameter, public :: iac_cdatac_gcam2emisfile_co2shipbase2000 = 13
  integer, parameter, public :: iac_cdatac_gcam2emisfile_lut720x360map = 14
  integer, parameter, public :: iac_cdatac_gcam2emisfile_downscaleinfo = 15
  integer, parameter, public :: iac_cdatac_gcam2emisfile_rcp45allsteps = 16


  !--- iofields ---
  integer, parameter, public :: iac_cdataio_gcami_crops     =  1
  integer, parameter, public :: iac_cdataio_gcami_flds      =  2
  integer, parameter, public :: iac_cdataio_gcamo_flds      =  3
  integer, parameter, public :: iac_cdataio_gcamo_regs      =  4
  integer, parameter, public :: iac_cdataio_glmi_flds       =  5
  integer, parameter, public :: iac_cdataio_glmo_flds       =  6

  !--- reals ---
  real*8,  parameter, public :: iac_spval = -999.0
  integer, parameter, public :: iac_gcam_nreg =  14
  integer, parameter, public :: iac_gcam_nsector =  12
  integer, parameter, public :: iac_gcamoemis_nemis=  1
  integer, parameter, public :: iac_gcam_naez =  18
  integer, parameter, public :: iac_gcam_ncrops =  27
  integer, parameter, public :: iac_gcam_timestep =  5
  integer, parameter, public :: iac_gcam_ioyears =  15
  integer, parameter, public :: iac_gcamo_ntime =  2
  integer, parameter, public :: iac_glm_nx  = 720
  integer, parameter, public :: iac_glm_ny  = 360
  integer, parameter, public :: iac_iaco_npfts =  17
  integer, parameter, public :: iac_iac_npfts  = 16
  character(len=*), parameter, public :: iac_gcam2glm_map = 'some_map_file'

  !--- integers ---
  integer, parameter, public :: iac_cdatai_logunit         =  1
  integer, parameter, public :: iac_cdatai_iac_nx          =  2
  integer, parameter, public :: iac_cdatai_iac_ny          =  3
  integer, parameter, public :: iac_cdatai_iac_size        =  4
  integer, parameter, public :: iac_cdatai_glm_nx          =  5
  integer, parameter, public :: iac_cdatai_glm_ny          =  6
  integer, parameter, public :: iac_cdatai_glm_size        =  7
  integer, parameter, public :: iac_cdatai_gcam_yr1        =  8
  integer, parameter, public :: iac_cdatai_gcam_yr2        =  9
  integer, parameter, public :: iac_cdatai_gcam_nreg       =  10
  integer, parameter, public :: iac_cdatai_gcam_naez       =  11
  integer, parameter, public :: iac_cdatai_gcam_timestep   =  12
  integer, parameter, public :: iac_cdatai_gcamo_ntime     =  13
  integer, parameter, public :: iac_cdatai_gcamo_nflds     =  14
  integer, parameter, public :: iac_cdatai_gcamo_size      =  15
  integer, parameter, public :: iac_cdatai_gcami_nflds     =  16
  integer, parameter, public :: iac_cdatai_gcami_size      =  17
  integer, parameter, public :: iac_cdatai_gcami_ncrops    =  18
  integer, parameter, public :: iac_cdatai_gcamoemis_size  =  19
  !--- logicals ---
  integer, parameter, public :: iac_cdatal_rest            =  1
  integer, parameter, public :: iac_cdatal_iac_present     =  2
  integer, parameter, public :: iac_cdatal_iac_prognostic  =  3
  integer, parameter, public :: iac_cdatal_glm_present     =  4
  integer, parameter, public :: iac_cdatal_glm_prognostic  =  5
  integer, parameter, public :: iac_cdatal_gcam_present    =  6
  integer, parameter, public :: iac_cdatal_gcam_prognostic =  7
  integer, parameter, public :: iac_cdatal_fastiac         =  8
  integer, parameter, public :: iac_cdatal_npphr           =  9
  integer, parameter, public :: iac_cdatal_initrun         =  10
  integer, parameter, public :: iac_cdatal_sneakermode     =  11
  integer, parameter, public :: iac_cdatal_nocarbonscale   =  12
  integer, parameter, public :: iac_cdatal_co2flux_coupling=  13
  integer, parameter, public :: iac_cdatal_write_rest      =  14

  integer           , public :: iac_iaci_nflds
  integer           , public :: iac_iaci_area
  integer, pointer  , public :: iac_iaci_pft(:)
  integer, pointer  , public :: iac_iaci_AGB(:)
  integer, pointer  , public :: iac_iaci_BGB(:)

  integer           , public :: iac_iaco_nflds
  integer, pointer  , public :: iac_iaco_pft(:)

  character(len=80) , pointer  , public :: iac_gcami_crop_names(:)
  character(len=80) , pointer  , public :: iac_gcami_fld_names(:)
  character(len=80) , pointer  , public :: iac_gcamo_reg_names(:)
  character(len=80) , pointer  , public :: iac_gcamo_fld_names(:)
  character(len=80) , pointer  , public :: iac_glmi_fld_names(:)
  character(len=80) , pointer  , public :: iac_glmo_fld_names(:)

  integer           , public :: iac_gcami_nflds
  integer           , public :: iac_gcami_biomass
  integer           , public :: iac_gcami_eucalyptus
  integer           , public :: iac_gcami_miscanthus
  integer           , public :: iac_gcami_willow
  integer           , public :: iac_gcami_Corn
  integer           , public :: iac_gcami_FiberCrop
  integer           , public :: iac_gcami_FodderGrass
  integer           , public :: iac_gcami_FodderHerb
  integer           , public :: iac_gcami_Forest
  integer           , public :: iac_gcami_Grassland
  integer           , public :: iac_gcami_Jatropha
  integer           , public :: iac_gcami_MiscCrop
  integer           , public :: iac_gcami_OilCrop
  integer           , public :: iac_gcami_OtherArableLand
  integer           , public :: iac_gcami_OtherGrain
  integer           , public :: iac_gcami_PalmFruit
  integer           , public :: iac_gcami_Pasture
  integer           , public :: iac_gcami_Rice
  integer           , public :: iac_gcami_RockIceDesert
  integer           , public :: iac_gcami_Root_Tuber
  integer           , public :: iac_gcami_Shrubland
  integer           , public :: iac_gcami_SugarCrop
  integer           , public :: iac_gcami_Tundra
  integer           , public :: iac_gcami_UnmanagedForest
  integer           , public :: iac_gcami_UnmanagedPasture
  integer           , public :: iac_gcami_UrbanLand
  integer           , public :: iac_gcami_Wheat
  integer,save      , public :: iac_gcami_above_ground_carbon=1
  integer,save      , public :: iac_gcami_below_ground_carbon=2

  integer, pointer  , public :: iac_gcami_cdens(:)

  integer           , public :: iac_gcamo_nflds
  integer           , public :: iac_gcamo_buildup
  integer           , public :: iac_gcamo_crop
  integer           , public :: iac_gcamo_pasture
  integer           , public :: iac_gcamo_woodharv
  integer           , public :: iac_gcamo_forest
  integer           , public :: iac_gcamo_other

  integer           , public :: iac_gcamo_reg_usa
  integer           , public :: iac_gcamo_reg_canada
  integer           , public :: iac_gcamo_reg_western_europe
  integer           , public :: iac_gcamo_reg_japan
  integer           , public :: iac_gcamo_reg_australia_nz
  integer           , public :: iac_gcamo_reg_former_soviet_union
  integer           , public :: iac_gcamo_reg_china
  integer           , public :: iac_gcamo_reg_middle_east
  integer           , public :: iac_gcamo_reg_africa
  integer           , public :: iac_gcamo_reg_latin_america
  integer           , public :: iac_gcamo_reg_southeast_asia
  integer           , public :: iac_gcamo_reg_eastern_europe
  integer           , public :: iac_gcamo_reg_korea
  integer           , public :: iac_gcamo_reg_india

  integer           , public :: iac_glmi_nflds
  integer           , public :: iac_glmi_natveg
  integer           , public :: iac_glmi_cropland
  integer           , public :: iac_glmi_pasture
  integer           , public :: iac_glmi_woodharv

  integer           , public :: iac_glmo_nflds
  integer           , public :: iac_glmo_gcrop
  integer           , public :: iac_glmo_gpast
  integer           , public :: iac_glmo_gothr
  integer           , public :: iac_glmo_gsecd
  integer           , public :: iac_glmo_gfvh1
  integer           , public :: iac_glmo_gfvh2
  integer           , public :: iac_glmo_gfsh1
  integer           , public :: iac_glmo_gfsh2
  integer           , public :: iac_glmo_gfsh3

end module iac_data_mod


     
     
