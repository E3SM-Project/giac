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
  use shr_sys_mod , only : shr_sys_abort
  use iac_data_mod, only : iac_ctl
  use shr_kind_mod,      only: CX => SHR_KIND_CX
  use iac_spmd_mod, only : masterproc
  use gcam_var_mod

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

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_init_mod

! !INTERFACE:
  subroutine gcam_init_mod(gcamo, gcamoemis,gcamoco2sfcjan, gcamoco2sfcfeb,   &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec)

! !DESCRIPTION:
! Initialize interface for gcam

! !USES:
    !use gcam_var, only : NLFilename_in
    use iac_data_mod

    ! To null-terminate the strings we pass into C++
    use iso_c_binding

    ! For error checking
    use mct_mod

    implicit none

! !ARGUMENTS:
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)
    real*8, pointer :: gcamoco2sfcjan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcfeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcsep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcoct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcnov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcdec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlofeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlosep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlooct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlonov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlodec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhifeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhisep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhioct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhinov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhidec(:,:)  ! gcam output for eam, needs to be passed through coupler
    
! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(gcam_init_mod)'
    character(len=128) :: casename
    integer      :: i, unitn, ier, len
    logical      :: lexist
    character(len=128) :: nlfilename_in


! !REVISION HISTORY:
! Author: T Shippert - modified heavily for E3SM/GCAM coupling

!---------------------------------------------------------------------
! Namelist variables
!---------------------------------------------------------------------

    ! Here is where we put our namelist inputs into the gdata structure
    ! (remember, we call it cdata in gcam-space, but it's not the same
    ! as teh cdata_z structure used in e3sm-space).  This is how we
    ! transmit the namelist variables downstream to all the gcam and
    ! gcam/coupling functions that need them - through gdata.
    allocate(gcamo(num_iac2elm_landtypes,num_gcam_land_regions), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamo',ier)
    allocate(gcamoemis(num_emiss_sectors,num_emiss_regions), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoemis',ier)

    ! Allocate variables to store gridded CO2 emissions from GCAM
    ! KVC: Not sure if this is where this should happen in the long run or not
    allocate(gcamoco2sfcjan(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcjan',ier)
    allocate(gcamoco2sfcfeb(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcfeb',ier)
    allocate(gcamoco2sfcmar(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcmar',ier)
    allocate(gcamoco2sfcapr(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcapr',ier)
    allocate(gcamoco2sfcmay(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcmay',ier)
    allocate(gcamoco2sfcjun(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcjun',ier)
    allocate(gcamoco2sfcjul(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcjul',ier)
    allocate(gcamoco2sfcaug(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcaug',ier)
    allocate(gcamoco2sfcsep(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcsep',ier)
    allocate(gcamoco2sfcoct(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcoct',ier)
    allocate(gcamoco2sfcnov(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcnov',ier)
    allocate(gcamoco2sfcdec(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2sfcdec',ier)
    
    allocate(gcamoco2airlojan(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlojan',ier)
    allocate(gcamoco2airlofeb(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlofeb',ier)
    allocate(gcamoco2airlomar(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlomar',ier)
    allocate(gcamoco2airloapr(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airloapr',ier)
    allocate(gcamoco2airlomay(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlomay',ier)
    allocate(gcamoco2airlojun(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlojun',ier)
    allocate(gcamoco2airlojul(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlojul',ier)
    allocate(gcamoco2airloaug(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airloaug',ier)
    allocate(gcamoco2airlosep(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlosep',ier)
    allocate(gcamoco2airlooct(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlooct',ier)
    allocate(gcamoco2airlonov(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlonov',ier)
    allocate(gcamoco2airlodec(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airlodec',ier)

    allocate(gcamoco2airhijan(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhijan',ier)
    allocate(gcamoco2airhifeb(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhifeb',ier)
    allocate(gcamoco2airhimar(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhimar',ier)
    allocate(gcamoco2airhiapr(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhiapr',ier)
    allocate(gcamoco2airhimay(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhimay',ier)
    allocate(gcamoco2airhijun(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhijun',ier)
    allocate(gcamoco2airhijul(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhijul',ier)
    allocate(gcamoco2airhiaug(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhiaug',ier)
    allocate(gcamoco2airhisep(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhisep',ier)
    allocate(gcamoco2airhioct(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhioct',ier)
    allocate(gcamoco2airhinov(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhinov',ier)
    allocate(gcamoco2airhidec(num_lon,num_lat), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamoco2airhidec',ier)

    ! create CCSM_GCAM_interface Object 
    call inite3sminterface()
    
    ! Call initcGCAM method of e3sm/GCAM Interface 
    write (iulog, *) trim(case_name), &
         trim(gcam_config),&
         trim(gcam2elm_co2_mapping_file),&
         trim(gcam2elm_luc_mapping_file), & 
         trim(gcam2elm_woodharvest_mapping_file)

    ! Null terminate
    len = len_trim(case_name)
    case_name(len+1:len+1) = c_null_char
    len = len_trim(gcam_config)
    gcam_config(len+1:len+1) = c_null_char
    len = len_trim(gcam2elm_co2_mapping_file)
    gcam2elm_co2_mapping_file(len+1:len+1) = c_null_char
    len = len_trim(gcam2elm_luc_mapping_file)
    gcam2elm_luc_mapping_file(len+1:len+1) = c_null_char
    len = len_trim(gcam2elm_woodharvest_mapping_file)
    gcam2elm_woodharvest_mapping_file(len+1:len+1) = c_null_char
    len = len_trim(base_co2_surface_file)
    base_co2_surface_file(len+1:len+1) = c_null_char
    len = len_trim(base_co2_aircraft_file)
    base_co2_aircraft_file(len+1:len+1) = c_null_char

    call initcGCAM(trim(case_name), &
         trim(gcam_config),&
         trim(gcam2elm_co2_mapping_file),&
         trim(gcam2elm_luc_mapping_file),&
         trim(gcam2elm_woodharvest_mapping_file))
    
  end subroutine gcam_init_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_run_mod

! !INTERFACE:
  subroutine gcam_run_mod(gcamo, gcamoemis, gcamoco2sfcjan, gcamoco2sfcfeb,   &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec)

! !DESCRIPTION:
! Run interface for gcam

! !USES:
   use iac_data_mod, only : iac_spval, iac_cdatal_rest, iac_cdatai_logunit
   use iac_data_mod, only : iac_eclock_ymd, iac_eclock_tod, iac_eclock_dt, iac_cdatal_rest
    implicit none

! !ARGUMENTS:
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)
    real*8, pointer :: gcamoco2sfcjan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcfeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcsep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcoct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcnov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcdec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlofeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlosep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlooct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlonov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlodec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhifeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhisep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhioct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhinov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhidec(:,:)  ! gcam output for eam, needs to be passed through coupler
    
! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: i,j,w
    character(len=*),parameter :: subname='(gcam_run_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

    gcamo = iac_spval
    gcamoemis = iac_spval
    
  restart_now = cdata%l(iac_cdatal_rest)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iulog,*) trim(subname),' date= ',ymd,tod

  !  Call runcGCAM method of E3SM Interface 
  call runcGCAM(ymd,gcamo,gcamoemis)

  ! If co2 emissions need to be passed from GCAM to EAM, then call downscale CO2                                 
  if ( iac_elm_co2_emissions ) then
     ! Convert logical to int for interface with C/C++
     if ( write_co2 ) then
        w = 1
     else
        w = 0
     end if

     call downscaleemissionscgcam(gcamoemis, gcamoco2sfcjan, gcamoco2sfcfeb, &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               & 
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec, base_co2_surface_file, base_co2emis_surface,      &
          base_co2_aircraft_file, base_co2emis_aircraft, num_lon, num_lat, w, ymd)

  end if

  end subroutine gcam_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_setdensity_mod

! !INTERFACE:
  subroutine gcam_setdensity_mod()

! !DESCRIPTION:
! Setdensity interface for gcam

! !USES:
   use iac_data_mod, only : iac_spval, iac_cdatal_rest, iac_cdatai_logunit
   use iac_data_mod, only : iac_eclock_ymd, iac_eclock_tod, iac_eclock_dt, iac_cdatal_rest
   use iac_data_mod, only : iac_first_coupled_year, lnd2iac_vars
   use iso_c_binding
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt, yyyymmdd
    integer :: i,j,r,w,s,len
    character(len=*),parameter :: subname='(gcam_setdensity_mod)'


! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!-----------------------------------------------------------------------

  restart_now = cdata%l(iac_cdatal_rest)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iulog,*) trim(subname),' date= ',ymd,tod

  ! convert logical to boolean for read_scalars and write_scalars                        
  if ( read_scalars ) then
     r = 1
  else
     r = 0
  end if

  if ( write_scalars ) then
     w = 1
  else
     w = 0
  end if

  if ( elm_iac_carbon_scaling ) then
     s = 1
  else
     s = 0
  end if

  ! Null terminate                                           
  !len = len_trim(elm2gcam_mapping_file)
  !elm2gcam_mapping_file(len+1:len+1) = c_null_char
  elm2gcam_mapping_file=trim(elm2gcam_mapping_file)//c_null_char
  base_npp_file=trim(base_npp_file)//c_null_char
  base_hr_file=trim(base_hr_file)//c_null_char
  base_pft_file=trim(base_pft_file)//c_null_char

  !  Call setdensity method of GCAM-E3SM Interface 
  call setdensitycGCAM(ymd, iac_ctl%area, &
       lnd2iac_vars%pftwgt, lnd2iac_vars%npp, lnd2iac_vars%hr, &
       iac_ctl%nlon, iac_ctl%nlat, iac_ctl%npft, elm2gcam_mapping_file,&
       iac_first_coupled_year, r, w, s, base_npp_file, base_hr_file, base_pft_file) 
  
  end subroutine gcam_setdensity_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_final_mod

! !INTERFACE:
  subroutine gcam_final_mod(gcamoco2sfcjan, gcamoco2sfcfeb,   &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec )

! !DESCRIPTION:
! Finalize gcam model
!------------------------------------------------------------------------------

   implicit none
! !ARGUMENTS:
    real*8, pointer :: gcamoco2sfcjan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcfeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcmay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcjul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcsep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcoct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcnov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2sfcdec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlofeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlomay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlojul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airloaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlosep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlooct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlonov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airlodec(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijan(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhifeb(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimar(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiapr(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhimay(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijun(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhijul(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhiaug(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhisep(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhioct(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhinov(:,:)  ! gcam output for eam, needs to be passed through coupler
    real*8, pointer :: gcamoco2airhidec(:,:)  ! gcam output for eam, needs to be passed through coupler

! !REVISION HISTORY:
! Author: T Craig
! Author: JET - added interface files for cesm/gcam communication

!EOP
!---------------------------------------------------------------------------

  !  Cleanup GCAM 
  call finalizecGCAM()

  !  Cleanup CCSM Interface Object 
  call deletee3sminterface()

  !  Deallocate 
  deallocate(gcamoco2sfcjan)
  deallocate(gcamoco2sfcfeb)
  deallocate(gcamoco2sfcmar)
  deallocate(gcamoco2sfcapr)
  deallocate(gcamoco2sfcmay)
  deallocate(gcamoco2sfcjun)
  deallocate(gcamoco2sfcjul)
  deallocate(gcamoco2sfcaug)
  deallocate(gcamoco2sfcsep)
  deallocate(gcamoco2sfcoct)
  deallocate(gcamoco2sfcnov)
  deallocate(gcamoco2sfcdec)

  deallocate(gcamoco2airlojan)
  deallocate(gcamoco2airlofeb)
  deallocate(gcamoco2airlomar)
  deallocate(gcamoco2airloapr)
  deallocate(gcamoco2airlomay)
  deallocate(gcamoco2airlojun)
  deallocate(gcamoco2airlojul)
  deallocate(gcamoco2airloaug)
  deallocate(gcamoco2airlosep)
  deallocate(gcamoco2airlooct)
  deallocate(gcamoco2airlonov)
  deallocate(gcamoco2airlodec)

  deallocate(gcamoco2airhijan)
  deallocate(gcamoco2airhifeb)
  deallocate(gcamoco2airhimar)
  deallocate(gcamoco2airhiapr)
  deallocate(gcamoco2airhimay)
  deallocate(gcamoco2airhijun)
  deallocate(gcamoco2airhijul)
  deallocate(gcamoco2airhiaug)
  deallocate(gcamoco2airhisep)
  deallocate(gcamoco2airhioct)
  deallocate(gcamoco2airhinov)
  deallocate(gcamoco2airhidec)

  end subroutine gcam_final_mod

!====================================================================================

end module gcam_comp_mod
