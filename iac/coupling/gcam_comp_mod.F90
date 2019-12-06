module gcam_comp_mod
  
!---------------------------------------------------------------------------
! !DESCRIPTION:
!  
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_file_mod, only : shr_file_get, shr_file_getUnit, shr_file_freeUnit
 
  use iac_data_mod         , only : cdata => gdata, EClock => GClock, &
                                    lnd2iac_type

  use iac_ctl_type         , only : nlat, nlon, npft

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: gcam_init_mod               ! gcam initialization
  public :: gcam_run_mod                ! gcam run phase
  public :: gcam_setdensity_mod         ! gcam set density phase
  public :: gcam_final_mod              ! gcam finalization/cleanup

! !PUBLIC DATA MEMBERS: None


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication


! !PRIVATE DATA MEMBERS:
! namelist values, to be used in the gcam functions
  ! probably need to be moved to the gcam_var_mod.f90 module or
  ! something, but for now i'll define them up here - they are only
  ! used within the gcam_comp_mod functions anyway.
  character(len=*) ::  case_name
  character(len=*) ::  gcam_config
  character(len=*) ::  base_co2_file
  character(len=*) ::  gcam2elm_co2_mapping_file
  character(len=*) ::  gcam2elm_luc_mapping_file
  character(len=*) ::  gcam2elm_woodharvest_mapping_file
  character(len=*) ::  elm2gcam_mapping_file

  logical :: read_scalars = .FALSE. ! if .FALSE., scalars are calculated from npp/hr

  logical :: read_elm_from_file = .TRUE. ! if .FALSE., elm data (npp,
  ! hr, area, pft weight) are passed from e3sm.

  logical :: write_co2 = .TRUE. ! gridded co2 emissions will be
  ! written to a file (in addition to passed in code).
  
  logical :: write_scalars = .TRUE. ! scalars will be written to a file.
  
  ! define coupling control variables
  ! these booleans define what is passed between gcam & e3sm.
  logical :: elm_iac_carbon_scaling = .TRUE.; ! if .TRUE., changes in
  ! land productivity from elm are used in gcam.
  logical :: iac_elm_co2_emissions = .TRUE.; ! if .TRUE., energy system
  ! co2 is passed from gcam to eam.
    
  ! define size control variables
  ! these integers define the length of the various arrays used in the coupling
  integer ::  num_lat = 180; ! number of horizontal grid cells
  integer ::  num_lon = 360; ! number of vertical grid cells
  integer ::  num_pft = 16;  ! number of pfts in elm
  integer ::  num_gcam_energy_regions = 32;
  integer ::  num_gcam_land_regions = 384;
  integer ::  num_iac2elm_landtypes = 9;

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_init_mod

! !INTERFACE:
  subroutine gcam_init_mod(gcamo, gcamoemis)

! !DESCRIPTION:
! Initialize interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(gcam_init_mod)'
    character(len=128) :: casename
    integer,save :: iulog = 6
    integer      :: i, unitn

! !REVISION HISTORY:
! Author: T Shippert - modified heavily for E3SM/GCAM coupling

!---------------------------------------------------------------------
! Namelist variables
!---------------------------------------------------------------------

    namelist /iac_inparm/ &
         case_name,gcam_config,base_co2_file, &
         gcam2elm_co2_mapping_file, gcam2elm_luc_mapping_file,&
         gcam2elm_woodharvest_mapping_file, elm2gcam_mapping_file,&
         read_scalars,read_elm_from_file, write_co2, write_scalars,&
         elm_iac_carbon_scaling, iac_elm_co2_emissions

    ! These might get renamed or removed - control values that might
    ! be duplicated elsewhere (maybe)
    namelist /iac_inparm/ &
         num_lat,num_lon,num_pft,num_gcam_energy_regions, &
         num_gcam_land_regions,num_iac2elm_landtypes
    
    ! Read in namelist
    ! Note: no concat with instance number, because we only have one
    ! instance (proc) currently
    nlfilename_iac = "gcam_in"

    inquire (file = trim(nlfilename_iac), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_iac does NOT exist:'&
            //trim(nlfilename_iac)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_rof does not exist')
    end if
  
    if (masterproc) then
       unitn = shr_file_getunit()
       write(iulog,*) 'Read in gcam_inparm namelist from: ', trim(nlfilename_iac)
       open( unitn, file=trim(nlfilename_iac), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, gcam_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on gcam_inparm read' )
          endif
       end do
       close(unitn)
       call shr_file_freeUnit(unitn)
    end if

    ! Removed a bunch of cdata() allocations, because I don't think
    ! they apply anymore - everything is handled via namelist and
    ! module variables, not configuration files.  But, if I'm wrong,
    ! see unused/giac_comp_mod.F90 and unused/cdata_frag.F90

    ! These allocations are currently wrong

    allocate(gcamo(iac_gcamo_nflds,cdata%i(iac_cdatai_gcamo_size)))
    allocate(gcamoemis(iac_gcamoemis_nemis,cdata%i(iac_cdatai_gcamoemis_size)))
    
    gcami = iac_spval
    gcamo = iac_spval
    gcamoemis = iac_spval
    
    ! create CCSM_GCAM_interface Object 
    call initCCSMInterface()
    
    ! Call initcGCAM method of CCSM/GCAM Interface 
    call initcGCAM()
    
  end subroutine gcam_init_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_run_mod

! !INTERFACE:
  subroutine gcam_run_mod(gcamo, gcamoemis)

! !DESCRIPTION:
! Run interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    real*8, pointer :: gcami(:,:)
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: iu
    integer :: i,j
    character(len=*),parameter :: subname='(gcam_run_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

  restart_now = cdata%l(iac_cdatal_rest)
  iu  = cdata%i(iac_cdatai_logunit)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iu,*) trim(subname),' date= ',ymd,tod

  !  Call runcGCAM method of E3SM Interface 
  call runcGCAM(ymd,gcamo,gcamo_emis,base_co2_file, nlon, nlat, write_co2)

  end subroutine gcam_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_setdensity_mod

! !INTERFACE:
  subroutine gcam_setdensity_mod(lnd2iac_vars)

! !DESCRIPTION:
! Setdensity interface for gcam

! !USES:
    implicit none

! !ARGUMENTS:
    type(lnd2iac_type), intent(in) :: lnd2iac_vars

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt, yyyymmdd
    integer :: iu
    integer :: i,j
    character(len=*),parameter :: subname='(gcam_setdensity_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

  restart_now = cdata%l(iac_cdatal_rest)
  iu  = cdata%i(iac_cdatai_logunit)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iu,*) trim(subname),' date= ',ymd,tod

  !  Call setdensity method of CCSM Interface 
  !call setdensitycGCAM(ymd,tod,gcami,size(gcami,dim=1),size(gcami,dim=2))
  call setdensitycGCAM(ymd, lnd2iac_vars%area, lnd2iac_vars%landfrac, &
       lnd2iac_vars%pftwgt, lnd2iac_vars%npp, lnd2iac_vars%hr, &
       nlon, nlat, npft, mapping_file, read_scalars, write_scalars)
  
  end subroutine gcam_setdensity_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_final_mod

! !INTERFACE:
  subroutine gcam_final_mod( )

! !DESCRIPTION:
! Finalize gcam model
!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:

! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!---------------------------------------------------------------------------

  !  Cleanup GCAM 
  call finalizecGCAM()

  !  Cleanup CCSM Interface Object 
  call deleteCCSMInterface()

  end subroutine gcam_final_mod

!====================================================================================

end module gcam_comp_mod
