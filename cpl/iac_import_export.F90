module iac_import_export
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Module to deal with importing and exporting coupled vars with
  !  internal iac variables.  
  ! !uses:
  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type
  use lnd2iacType  , only: lnd2iac_type
  use iac2lndType  , only: iac2lnd_type
  use iac2atmType  , only: iac2atm_type
  use GridcellType , only: grc_pp          ! for access to gridcell topology
  use gcam_cpl_indices
  use mct_mod
  !
  implicit none
  public :: iac_import    ! import lnd vars from coupler
  public :: iac_export    ! export iac vars to lnd and atm via coupler

contains
  subroutine lnd_import(bounds, x2z, lnd2iac_vars)
    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to iac
    !
    ! !USES:
    use gcam_ctl
    use controlMod       , only: NLFilename
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use domainMod        , only: ldomain
    use shr_kind_mod     , only: r8 => shr_kind_r8, CL => shr_kind_CL
    use fileutils        , only: getavu, relavu
    use spmdmod          , only: masterproc, mpicom, MPI_REAL8
    use netcdf
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2z(:,:) ! driver import state to iacmodel
    type(iac2lnd_type) , intent(inout) :: iac2lnd_vars    ! gcam internal input data type
    
  end subroutine lnd_import
   !===============================================================================

  subroutine iac_export( bounds, iac2lnd_vars, iac2atm_vars, z2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the clm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use iac_varctl         
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(iac2lnd_type), intent(inout) :: iac2lnd_vars ! gcam to land output
    type(iac2atm_type), intent(inout) :: iac2atm_vars ! gcam to atm output
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid

  end subroutine iac_export
end module iac_import_export
