
Module glm_comp_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: glm_comp_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:

  use iac_data_mod, only : cdata => gdata, EClock => GClock
  use iac_data_mod
  use gcam_var_mod

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: glm_init_mod               ! glm initialization
  public :: glm_run_mod                ! glm run phase
  public :: glm_final_mod              ! glm finalization/cleanup

! !PUBLIC DATA MEMBERS: None

  integer, parameter :: glm_data_size   = iac_glm_nx*iac_glm_ny ! should be set by glm

! !REVISION HISTORY:
! Author: T Craig


! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm_init_mod

! !INTERFACE:
  subroutine glm_init_mod(glmi, glmi_wh,glmo)

! !DESCRIPTION:
! Initialize interface for glm

! !USES:
    use iac_data_mod
    implicit none

! !ARGUMENTS:
    real*8, pointer :: glmi(:,:)
    real*8, pointer :: glmi_wh(:)
    real*8, pointer :: glmo(:,:)

! !LOCAL VARIABLES:
    integer :: numreg,numglu
    character(len=*),parameter :: subname='(glm_init_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    numreg = num_gcam_energy_regions
    numglu = num_gcam_land_regions


    cdata%l(iac_cdatal_glm_present) = .true.
    cdata%l(iac_cdatal_glm_prognostic) = .true.
    cdata%i(iac_cdatai_glm_nx) = iac_glm_nx
    cdata%i(iac_cdatai_glm_ny) = iac_glm_ny
    cdata%i(iac_cdatai_glm_size) = iac_glm_nx * iac_glm_ny

    ! Already allocated by gcam_init_mod - I think.
    allocate(glmi(iac_glmi_nflds,glm_data_size))
    allocate(glmi_wh(numglu))
    allocate(glmo(iac_glmo_nflds,glm_data_size))

    glmi = iac_spval
    glmo = iac_spval

    call initGLM()
  end subroutine glm_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm_run_mod

! !INTERFACE:
  subroutine glm_run_mod(glmi, glmi_wh, glmo)

! !DESCRIPTION:
! Run interface for glm

! !USES:
    implicit none

! !ARGUMENTS:
    real*8, pointer :: glmi(:,:)
    real*8, pointer :: glmi_wh(:)
    real*8, pointer :: glmo(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: i,j,ij,n,ni
    integer :: glmyear
    character(len=*),parameter :: subname='(glm_run_mod)'


! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    restart_now = cdata%l(iac_cdatal_rest)

    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)

    write(iulog,*) trim(subname),' date= ',ymd,tod

!-- import ---
!    do n = 1,iac_glmi_nflds
!    do ij = 1,glm_data_size
!       xxx = glmi(n,ij)
!    enddo
!    enddo

!-- advance ---
!!$    do n = 1,iac_glmo_nflds
!!$    do ij = 1,glm_data_size
!!$       ni = mod(n-1,iac_glmi_nflds) + 1
!!$       glmo(n,ij) = glmi(ni,ij) + 1000.0
!!$    enddo
!!$    enddo
    glmyear=ymd/10000
    write(6,*)'stepping glm',glmyear,ymd
!    call stepGLM(glmyear,glmi,size(glmi,dim=1),size(glmi,dim=2),glmo,size(glmo,dim=1),size(glmo,dim=2));
    call stepGLM(glmyear,            &
                 glmi,size(glmi,dim=1),size(glmi,dim=2),  &
                 glmi_wh,size(glmi_wh,dim=1),  &
                 glmo,size(glmo,dim=1),size(glmo,dim=2));
    write(6,*)'done stepping glm'

  end subroutine glm_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm_final_mod

! !INTERFACE:
  subroutine glm_final_mod( )

! !DESCRIPTION:
! Finalize glm model

!------------------------------------------------------------------------------

! !USES:
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(glm_init_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!---------------------------------------------------------------------------

     call finalizeGLM()
  end subroutine glm_final_mod

!====================================================================================

end module glm_comp_mod

