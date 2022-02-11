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
  integer, pointer, public ::index_z2x_Sz_pct_pft(:)      ! percent pft of vegetated land unit
  ! previous year percent of pft vegetated land unit (for time interpolation)
  integer, pointer, public ::index_z2x_Sz_pct_pft_prev(:)
  ! harvest fractions
  integer, pointer, public ::index_z2x_Sz_harvest_frac(:)

  integer, public ::index_z2x_Fazz_fco2_iac    ! co2 from iac to atm, flux
  integer, public ::nflds_z2x = 0

  ! drv -> iac
  ! The stuff we get from the coupler (from lnd)
  integer, pointer, public ::index_x2z_Sl_hr(:)        ! total heterotrophic respiration
  integer, pointer, public ::index_x2z_Sl_npp(:)       ! net primary production
  integer, pointer, public ::index_x2z_Sl_pftwgt(:)    ! pft weights for each cell
  integer, public ::nflds_x2z = 0

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine gcam_cpl_indices_init( )
    !
    ! !DESCRIPTION:
    ! Allocate our coupler index arrays
    !
    ! !USES:
    use mct_mod
    use shr_file_mod, only : shr_file_get, shr_file_getUnit, shr_file_freeUnit
    use iac_spmd_mod, only : masterproc
    use gcam_var_mod, only : iulog
    !
    ! !LOCAL VARIABLES:
    integer           :: ier, unitn
    character(len=32), parameter :: subname = 'gcam_cpl_indices_init'
    logical :: lexist
    character(len=32) :: nlfilename_iac
    integer :: num_pft, num_harvest

! avd - how do we get these from the namelist? they are not found until
! the comp iac_init(), which is called after this routine in iac_init_mct()

! let's try a limited namelist read-in
! this doesn't work for some reason
    namelist /gcam_inparm/ num_pft, num_harvest
    
    nlfilename_iac = "gcam_in"

    inquire (file = trim(nlfilename_iac), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_iac does NOT exist:'&
            //trim(nlfilename_iac)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_iac does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(iulog,*) 'Read in gcam_inparm num_pft, num_harvest from: ', &
          trim(nlfilename_iac)
       open( unitn, file=trim(nlfilename_iac), status='old' )
       ier = 1
       read(unitn, gcam_inparm, iostat=ier)
       if (ier < 0) then
          call shr_sys_abort( subname//' encountered end-of-file on &
             gcam_inparm read' )
       endif
       call shr_file_freeUnit( unitn )
    end if

    !write(iulog,'(A20,a,I2,a,I1)') subname , 'num_pft=' , num_pft , &
    !               'num_harvest=' , num_harvest

    iac_ctl%npft = 17
    iac_ctl%nharvest = 5

    allocate(index_z2x_Sz_pct_pft(iac_ctl%npft), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate index_z2x_Sz_pct_pft',ier)
    allocate(index_z2x_Sz_pct_pft_prev(iac_ctl%npft), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate index_z2x_Sz_pct_pft_prev',ier)
    allocate(index_z2x_Sz_harvest_frac(iac_ctl%nharvest), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate index_z2x_Sz_harvest_frac',ier)
    allocate(index_x2z_Sl_hr(iac_ctl%npft))
    if(ier/=0) call mct_die(subName,'allocate index_x2z_Sl_hr',ier)
    allocate(index_x2z_Sl_npp(iac_ctl%npft))
    if(ier/=0) call mct_die(subName,'allocate index_x2z_Sl_npp',ier)
    allocate(index_x2z_Sl_pftwgt(iac_ctl%npft))
    if(ier/=0) call mct_die(subName,'allocate index_x2z_Sl_pftwgt',ier)
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
    use gcam_var_mod, only : iulog
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

    ! KVC NOTE: iac_ctl%npft is not set at this point, setting it now
    ! Loop over pfts and get a tag to concat with
    ! avd - tried setting this in iac_init
    !iac_ctl%npft=17
    do p=1,iac_ctl%npft
       ! We zero-offset the names, with 0 being bare ground, so tag with p-1
       write(pftstr,'(I0)') p-1
       pftstr=trim(pftstr)

       !-------------------------------------------------------------
       ! iac -> lnd
       !-------------------------------------------------------------
       index_z2x_Sz_pct_pft(p) = mct_avect_indexra(z2x,trim('Sz_pct_pft' // pftstr))
       index_z2x_Sz_pct_pft_prev(p) = &
          mct_avect_indexra(z2x,trim('Sz_pct_pft_prev' // pftstr))

       if (p <= iac_ctl%nharvest) then
          index_z2x_Sz_harvest_frac(p) = &
             mct_avect_indexra(z2x,trim('Sz_harvest_frac' // pftstr))
       end if

       !-------------------------------------------------------------
       ! lnd -> iac
       !-------------------------------------------------------------
       index_x2z_Sl_hr(p) = mct_avect_indexra(x2z,trim('Sl_hr_pft' // pftstr))
       index_x2z_Sl_npp(p) = mct_avect_indexra(x2z,trim('Sl_npp_pft' // pftstr))
       index_x2z_Sl_pftwgt(p) = mct_avect_indexra(x2z,trim('Sl_pftwgt_pft' // pftstr))

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
    deallocate(index_z2x_Sz_pct_pft_prev)
    deallocate(index_z2x_Sz_harvest_frac)
    deallocate(index_x2z_Sl_hr)
    deallocate(index_x2z_Sl_npp)
    deallocate(index_x2z_Sl_pftwgt)
  end subroutine gcam_cpl_indices_finish

end module gcam_cpl_indices
