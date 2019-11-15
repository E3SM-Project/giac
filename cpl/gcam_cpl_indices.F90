module gcam_cpl_indices
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !    Module containing the indices for the fields passed between gcam
  !    the driver, copied shamelessly from the clm coupler interface.
  !    Here, we need to get index arrays with one element per pft, so
  !    we need to init/allocate, too.
  ! !USES:
  
  use shr_sys_mod,    only : shr_sys_abort
  use iac_data_mod,   only : iac_ctl
  implicit none

  SAVE
  private                              ! By default make data private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: gcam_cpl_indices_init       ! Allocates the coupler index arrays
  public :: gcam_cpl_indices_set        ! Set the coupler indices
  public :: gcam_cpl_indices_finish     ! Deallocate
  ! !PUBLIC DATA MEMBERS:
  !

  ! Indeces are all arrays [pft], with one elements for each pft

  ! iac -> drv
  ! The stuff we send back to the coupler (i.e. to lnd)
  integer, public ::index_z2x_Sz_pct_pft(:)      ! percent pft of vegetated land unit
  integer, public ::index_z2x_Fazz_co2flux(:)    ! co2 from iac to atm, flux
  integer, public ::nflds_z2x = 0

  ! drv -> iac
  ! The stuff we get from the coupler (from lnd)
  integer, public ::index_x2z_Sl_hr(:)        ! total heterotrophic respiration
  integer, public ::index_x2z_Sl_npp(:)       ! net primary production
  integer, public ::index_x2z_Sl_pftwgt(:)    ! pft weights for each cell
  integer, public ::nflds_x2z = 0

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine gcam_cpl_indices_init( )
    !
    ! !DESCRIPTION:
    ! Allocate our coupler index arrays
    allocate(index_z2x_Sz_pct_pft(iac_ctl%npft))
    allocate(index_x2z_Sl_hr(iac_ctl%npft))
    allocate(index_x2z_Sl_npp(iac_ctl%npft))
    allocate(index_x2z_Sl_pftwgt(iac_ctl%npft))
  end subroutine gcam_cpl_indices_init

  !-----------------------------------------------------------------------
  subroutine gcam_cpl_indices_set( )
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
    integer           :: num, p
    character(len=4)  :: pftstr   ! Up to 1000 pfts...
    character(len=32) :: subname = 'gcam_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    ! Determine attribute vector indices

    ! create temporary attribute vectors to extract from
    call mct_aVect_init(x2z, rList=seq_flds_x2z_fields, lsize=1)
    nflds_x2z = mct_avect_nRattr(x2z)

    call mct_aVect_init(z2x, rList=seq_flds_z2x_fields, lsize=1)
    nflds_z2x = mct_avect_nRattr(z2x)

    ! Loop over pfts and get a tag to concat with
    do p=1,iac_ctl%npft

       ! We zero-offset the names, with 0 being bare ground, so tag with p-1
       write(pftstr,'(I0)') p-1
       pftstr=trim(pftstr)

       !-------------------------------------------------------------
       ! iac -> lnd
       !-------------------------------------------------------------
       index_z2x_Sz_pct_pft(p) = mct_avect_indexra(z2x,'Sz_pct_pft' // pftstr)

       !-------------------------------------------------------------
       ! lnd -> iac
       !-------------------------------------------------------------
       index_x2z_Sl_hr(p) = mct_avect_indexra(x2z,'Sl_hr_pft' // pftstr)
       index_x2z_Sl_npp(p) = mct_avect_indexra(x2z,'Sl_npp_pft' // pftstr)
       index_x2z_Sl_pftwgt(p) = mct_avect_indexra(x2z,'Sl_pftwgt_pft' // pftstr)

    end do

    ! iac -> atm
    index_z2x_Fazz_fco2_iac = mct_avect_indexra(z2x,'Fazz_fco2_iac')

    call mct_aVect_clean(x2z)
    call mct_aVect_clean(z2x)
       
  end subroutine gcam_cpl_indices_set

!=======================================================================

  subroutine gcam_cpl_indices_finish( )
    !
    ! !DESCRIPTION:
    ! Dellocate our coupler index arrays
    deallocate(index_z2x_Sz_pct_pft)
    deallocate(index_x2z_Sl_hr)
    deallocate(index_x2z_Sl_npp)
    deallocate(index_x2z_Sl_pftwgt)
  end subroutine gcam_cpl_indices_init

end module gcam_cpl_indices
