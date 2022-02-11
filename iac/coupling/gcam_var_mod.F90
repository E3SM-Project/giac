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
  logical, public :: gcam_active = .true.  ! true to turn on gcam coupling

  ! Some sizes
  integer, public :: gcam_nlon, gcam_nlat

  ! Namelist variables for use in gcam
  ! Namelist variables
  character(len=256), public ::  case_name
  character(len=256), public ::  gcam_config
  character(len=256), public ::  base_co2_surface_file
  character(len=256), public ::  base_co2_aircraft_file
  character(len=256), public ::  base_npp_file
  character(len=256), public ::  base_hr_file
  character(len=256), public ::  base_pft_file
  character(len=256), public ::  gcam2elm_co2_mapping_file
  character(len=256), public ::  gcam2elm_luc_mapping_file
  character(len=256), public ::  gcam2elm_woodharvest_mapping_file
  character(len=256), public ::  elm2gcam_mapping_file
  character(len=256), public ::  gcam_gridfile

  !glm variables
  character(len=256), public :: gcam2glm_baselu
  character(len=256), public :: gcam2glm_glumap
  character(len=256), public :: gcam2glm_basebiomass

  logical, public :: read_scalars = .FALSE. ! if .FALSE., scalars are calculated from npp/hr

  logical, public :: read_elm_from_file = .TRUE. ! if .FALSE., elm data (npp,
  ! hr, area, pft weight) are passed from e3sm.

  logical, public :: write_co2 = .TRUE. ! gridded co2 emissions will be
  ! written to a file (in addition to passed in code).
  
  logical, public :: write_scalars = .TRUE. ! scalars will be written to a file.
  
  ! define coupling control variables
  ! these booleans define what is passed between gcam & e3sm.
  logical, public :: elm_iac_carbon_scaling = .TRUE.; ! if .TRUE., changes in
  ! land productivity from elm are used in gcam.
  logical, public :: iac_elm_co2_emissions = .TRUE.; ! if .TRUE., energy system
  ! co2 is passed from gcam to eam.
    
  ! define size control variables
  ! these integers define the length of the various arrays used in the coupling
  integer, public ::  num_lat = 192 ! number of horizontal grid cells
  integer, public ::  num_lon = 288 ! number of vertical grid cells
  integer, public ::  num_pft = 17  ! number of pfts in elm
  integer, public ::  num_harvest = 5 ! number of harvest fields
  integer, public ::  num_gcam_energy_regions = 32
  integer, public ::  num_gcam_land_regions = 392
  integer, public ::  num_iac2elm_landtypes = 9
  integer, public ::  num_emiss_sectors = 2
  integer, public ::  num_emiss_regions = 1

  ! define base year levels of CO2
  real*8, public :: base_co2emis_surface = 9663.0297;
  real*8, public :: base_co2emis_aircraft = 102.157;

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
