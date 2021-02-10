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
  subroutine gcam_init_mod(gcamo, gcamoemis)

! !DESCRIPTION:
! Initialize interface for gcam

! !USES:
    !use gcam_var, only : NLFilename_in
    use iac_data_mod
    implicit none

! !ARGUMENTS:
    real*8, pointer :: gcamo(:,:)
    real*8, pointer :: gcamoemis(:,:)

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(gcam_init_mod)'
    character(len=128) :: casename
    integer      :: i, unitn, ier
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
    allocate(gcamo(iac_gcamo_nflds,gdata%i(iac_cdatai_gcamo_size)))
    allocate(gcamoemis(iac_gcamoemis_nemis,gdata&
         %i(iac_cdatai_gcamoemis_size)))

    ! create CCSM_GCAM_interface Object 
    call inite3sminterface()
    
    ! Call initcGCAM method of e3sm/GCAM Interface 
    write (iulog, *) trim(case_name), &
         trim(gcam_config),&
         trim(gcam2elm_co2_mapping_file),&
         trim(gcam2elm_luc_mapping_file), & 
         trim(gcam2elm_woodharvest_mapping_file)
         
    call initcGCAM(case_name, gcam_config, gcam2elm_co2_mapping_file,&
         gcam2elm_luc_mapping_file, gcam2elm_woodharvest_mapping_file)
    
  end subroutine gcam_init_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_run_mod

! !INTERFACE:
  subroutine gcam_run_mod(gcamo, gcamoemis)

! !DESCRIPTION:
! Run interface for gcam

! !USES:
   use iac_data_mod, only : iac_spval, iac_cdatal_rest, iac_cdatai_logunit
   use iac_data_mod, only : iac_eclock_ymd, iac_eclock_tod, iac_eclock_dt, iac_cdatal_rest
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

    gcami = iac_spval
    gcamo = iac_spval
    gcamoemis = iac_spval
    
  restart_now = cdata%l(iac_cdatal_rest)
  iu  = cdata%i(iac_cdatai_logunit)

  ymd = EClock(iac_eclock_ymd)
  tod = EClock(iac_eclock_tod)
  dt  = EClock(iac_eclock_dt)

  write(iu,*) trim(subname),' date= ',ymd,tod

  !  Call runcGCAM method of E3SM Interface 
  call runcGCAM(ymd,gcamo,gcamoemis,base_co2_file, iac_ctl%nlon, iac_ctl%nlat, write_co2)

  end subroutine gcam_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam_setdensity_mod

! !INTERFACE:
  subroutine gcam_setdensity_mod(lnd2iac_vars)

! !DESCRIPTION:
! Setdensity interface for gcam

! !USES:
   use iac_data_mod, only : iac_spval, iac_cdatal_rest, iac_cdatai_logunit
   use iac_data_mod, only : iac_eclock_ymd, iac_eclock_tod, iac_eclock_dt, iac_cdatal_rest
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
       iac_ctl%nlon, iac_ctl%nlat, iac_ctl%npft, elm2gcam_mapping_file,&
       read_scalars, write_scalars) 
  
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
  call deletee3sminterface()

  end subroutine gcam_final_mod

!====================================================================================

end module gcam_comp_mod
