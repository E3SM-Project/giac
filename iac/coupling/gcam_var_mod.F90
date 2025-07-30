module gcam_var_mod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL
  use shr_sys_mod , only: shr_sys_abort ! cannot use endrun here due to circular dependency
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: gcam_var_set    ! Set variables
  !
  private
  save
  !
  ! !PUBLIC TYPES:
  !
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer , parameter, public ::  fname_len = SHR_KIND_CL   ! max length of file names in this module

  ! netCDF fill values as you write out
  real(r8), public, parameter :: spval    = 1.e36_r8        ! special value for real data
  integer , public, parameter :: ispval   = -9999           ! special value for int data


  !----------------------------------------------------------
  !
  ! Run control variables
  !
  ! case id
  character(len=256), public :: caseid  = ' '                            

  ! case title
  character(len=256), public :: ctitle  = ' '                            

  ! Type of run
  integer, public :: nsrest             = iundef                         

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0                          

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1                          

  ! Branch from restart files
  integer, public, parameter :: nsrBranch   = 2                          

  ! true => allow case name to remain the same for branch run
  ! by default this is not allowed
  logical, public :: brnch_retain_casename = .false.                     

  !true => no valid land points -- do NOT run
  logical, public :: noland = .false.                                    

  ! Hostname of machine running on
  character(len=256), public :: hostname = ' '                           

  ! username of user running program
  character(len=256), public :: username = ' '                           

  ! description of this source
  character(len=256), public :: source   = " " 

  ! version of program
  character(len=256), public :: version  = " "                           

  !----------------------------------------------------------
  ! Unit Numbers
  !----------------------------------------------------------
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6

  !----------------------------------------------------------
  ! To retrieve namelist
  !----------------------------------------------------------
  character(len=SHR_KIND_CL), public :: NLFilename_in ! Namelist filename
  !
  logical, private :: gcam_var_isset = .false.

  !----------------------------------------------------------
  ! instance control
  !----------------------------------------------------------

  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

  !----------------------------------------------------------
  ! Active gcam
  !----------------------------------------------------------
  logical, public :: gcam_active = .false.  ! true to turn on gcam coupling

  logical, public :: gcam_alarm =.false.    ! true in model years when gcam runs, based on the EClock iac run alarm

  ! Some sizes
  integer, public :: gcam_nlon, gcam_nlat

  ! Namelist variables for use in the iac and gcam
  ! do not initialize here anymore - defaults are set in
  !     namelist_defaults_gcam.xml
  ! Namelist variables

  ! gcam case name
  character(len=256), public ::  case_name

  ! Grid and region size parameters
  integer, public ::  num_pft             ! number of pfts in elm
  integer, public ::  num_harvest         ! number of elm harvest fields
  integer, public ::  num_lat             ! number of horizontal grid cells
  integer, public ::  num_lon             ! number of vertical grid cells
  integer, public ::  num_gcam_energy_regions   ! gcam region number
  integer, public ::  num_gcam_land_regions     ! gcam land units: regionXglu
  integer, public ::  num_iac2elm_landtypes     ! number of gcamo land fields
  integer, public ::  num_emiss_sectors         ! for emis downscaling
  integer, public ::  num_emiss_ctys            ! for emis downscaling
  integer, public ::  num_periods               ! for emis downscaling

  ! gcam config and init files
  character(len=256), public ::  gcam_config
  character(len=256), public ::  base_gcam_co2_file  ! baseline gcam out
  character(len=256), public ::  base_gcam_lu_wh_file  ! baseline gcam out
  character(len=256), public ::  base_co2_surface_file ! gridded basline
  character(len=256), public ::  base_co2_shipment_file ! gridded basline
  character(len=256), public ::  base_co2_aircraft_file ! gridded baseline
  character(len=256), public ::  base_npp_file   ! by pft, lat, lon, for land scalars
  character(len=256), public ::  base_hr_file    ! by pft, lat, lon, for land scalars
  character(len=256), public ::  base_pft_file   ! by pft, lat, lon, for land scalars
  character(len=256), public ::  gcam2elm_co2_mapping_file ! define gcam co2 out
  character(len=256), public ::  gcam2elm_luc_mapping_file ! def gcam lu out
  character(len=256), public ::  gcam2elm_woodharvest_mapping_file ! def gcam wh out
  character(len=256), public ::  gcam2elm_cdensity_mapping_file ! def gcam cdensity in

  ! Grid and mapping and initialization files
  character(len=256), public ::  gcam_gridfile ! definition of iac grid
  character(len=256), public ::  elm2gcam_mapping_file ! def gcam units to grid
  character(len=256), public :: gcam2glm_glumap
  character(len=256), public :: gcam2glm_baselu
  character(len=256), public :: gcam2glm_basebiomass

  ! Config and input files for emiss downscaling
  character(len=256), public :: country2grid_map
  character(len=256), public :: country2region_map
  character(len=256), public :: pop_iiasa_file
  character(len=256), public :: gdp_iiasa_file
  character(len=256), public :: pop_gcam_file
  character(len=256), public :: gdp_gcam_file
  character(len=256), public :: co2_gcam_file
  character(len=256), public :: surface_co2_downscaling_method
 
  ! future (>= 2015) land conversion assumptions 
  integer, public :: crop_addtreeonly = 0
  real(r8), public :: crop_setherbfracrem = 1.0
  real(r8), public :: crop_setavailtreefracrem = 0.0
  integer, public :: pasture_addtreeonly = 0
  real(r8), public :: pasture_setherbfracrem = 1.0
  real(r8), public :: pasture_setavailtreefracrem = 0.0
  
  ! Name only of the dynamic landuse timeseries file
  character(len=256), public ::  fdyndat_ehc

  ! runtime options
  logical, public :: read_scalars ! if .FALSE., scalars are calculated from npp/hr
  logical, public :: write_scalars ! scalars will be written to a file.
  ! hr, area, pft weight) are passed from e3sm.
  logical, public :: write_co2 ! gridded co2 emissions will be
  ! written to a file (in addition to passed in code).
  
  ! define coupling control variables
  ! these booleans define what is passed between gcam & e3sm.
  logical, public :: elm_ehc_agyield_scaling ! if .TRUE., changes in
  ! land productivity from elm scale ag yield in gcam.
  logical, public :: elm_ehc_carbon_scaling ! if .TRUE., changes in
  ! land productivity from elm scale carbon density in gcam.
  logical, public :: ehc_eam_co2_emissions ! if .TRUE., energy system
  ! co2 is passed from gcam to eam.

  logical, public :: gcam_spinup  ! if true do gcam spinup
  logical, public :: run_gcam     ! if true run the gcam model
    
contains

  subroutine gcam_var_set(caseid_in, ctitle_in,brnch_retain_casename_in, &
         nsrest_in, version_in, hostname_in, username_in)
    !
    ! !DESCRIPTION:
    ! Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=256), optional, intent(IN) :: caseid_in                ! case id
    character(len=256), optional, intent(IN) :: ctitle_in                ! case title
    logical,            optional, intent(IN) :: brnch_retain_casename_in ! true => allow case name to remain the 
                                                                         ! same for branch run
    integer,            optional, intent(IN) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
    character(len=256), optional, intent(IN) :: version_in               ! model version
    character(len=256), optional, intent(IN) :: hostname_in              ! hostname running on
    character(len=256), optional, intent(IN) :: username_in              ! username running job
    !-----------------------------------------------------------------------

    if ( gcam_var_isset )then
       call shr_sys_abort(' ERROR:: control variables already set, cannot call this routine')
    end if
    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(brnch_retain_casename_in) ) brnch_retain_casename = brnch_retain_casename_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in

    ! I'm no longer sure this gets set here, but for now let's turn
    ! it on if we call gcam_var_set()  
    gcam_active = .true.

    ! This is only if we pass it in to the call, which is dumb - why
    ! would we call it otherwise?  But I might have lost track of the
    ! sequence here.
    !if ( present(gcam_active     ) ) gcam_active   = gcam_active

  end subroutine gcam_var_set
end module gcam_var_mod
