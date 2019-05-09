module gcam_cpl_indices
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !    Module containing the indices for the fields passed between gcam
  !    the driver, copied shamelessly from the clm coupler interface.
  !    This is a little overkill for just three states, but it helps
  !    to keep things organized. 
  ! !USES:
  
  use shr_sys_mod,    only : shr_sys_abort
  implicit none

  SAVE
  private                              ! By default make data private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: giac_cpl_indices_set        ! Set the coupler indices
  !
  ! !PUBLIC DATA MEMBERS:
  !

  ! iac -> drv
  ! The stuff we send back to the coupler (i.e. to lnd)
  integer, public ::index_z2x_Sz_pct_pft      ! percent pft of vegetated land unit
  integer, public ::nflds_z2x = 0

  ! drv -> iac
  ! The stuff we get from the coupler (from lnd)
  integer, public ::index_x2z_Sl_hr           ! total heterotrophic respiration
  integer, public ::index_x2z_Sl_npp          ! net primary production
  integer, public ::nflds_x2z = 0

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine giac_cpl_indices_set( )
    !
    ! !DESCRIPTION: 
    ! Set the coupler indices needed by the gcam coupler
    ! interface.
    !
    ! !USES:
    use seq_flds_mod   , only: seq_flds_x2z_fields, seq_flds_z2x_fields
    use mct_mod        , only: mct_aVect, mct_aVect_init, mct_avect_indexra
    use mct_mod        , only: mct_aVect_clean, mct_avect_nRattr
    !
    ! !ARGUMENTS:
    implicit none
    !
    ! !REVISION HISTORY:
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: z2x      ! temporary, iac to coupler
    type(mct_aVect)   :: x2z      ! temporary, coupler to iac
    integer           :: num 
    character(len= 2) :: cnum
    character(len=64) :: name
    character(len=32) :: subname = 'gcam_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2z, rList=seq_flds_x2z_fields, lsize=1)
    nflds_x2z = mct_avect_nRattr(x2z)

    call mct_aVect_init(z2x, rList=seq_flds_z2x_fields, lsize=1)
    nflds_z2x = mct_avect_nRattr(z2x)

    !-------------------------------------------------------------
    ! gcam -> drv 
    !-------------------------------------------------------------
    index_z2x_Sz_pct_pft    = mct_avect_indexra(z2x,'Sz_pct_pft')

    !-------------------------------------------------------------
    ! drv -> gcam
    !-------------------------------------------------------------
    index_x2z_Sl_hr         = mct_avect_indexra(x2z,'Sa_hr')
    index_x2z_Sl_npp        = mct_avect_indexra(x2z,'Sa_npp')

    call mct_aVect_clean(x2z)
    call mct_aVect_clean(z2x)

  end subroutine gcam_cpl_indices_set

!=======================================================================

end module gcam_cpl_indices
