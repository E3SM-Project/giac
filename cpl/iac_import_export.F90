module iac_import_export
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Module to deal with importing and exporting coupled vars with
  !  internal iac variables.  
  ! !uses:
  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
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
  subroutine iac_import(x2z, lnd2iac_vars)
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
    real(r8)           , intent(in)    :: x2z(:,:) ! driver import state to iacmodel
    type(lnd2iac_type) , intent(inout) :: lnd2iac_vars    ! gcam
    !
    ! LOCAL VARIABLES
    integer :: n,n1
    integer :: begr, endr
    character(len=32), parameter :: sub = 'iac_import'

    ! Gcam expects things in npp_m[lon][lat][pft] format, so we need
    ! to extract from the flattened column representation.
    
    ! Once again, domain decomp isn't meaningful with one proc, but keep it for
    ! consistency with everybody else
    begg = iac_ctl%begg
    endg = iac_ctl%endg

    ! This won't work, if we are of different ranks.  We could  use
    ! reshape, but that gets complicated if we do have some domain
    ! decomposition.  What I don't currently know is how to go from
    ! grid cell index to lat/lon in a domain decomposition.  So this
    ! means we want probably want to keep our vars in a grid
    ! decomposition and extract lats and lons later on.
 
    ! Also, *this* won't work either, because we have pfts in our
    ! dimensioning.  I need to review how attribute vectors work
    ! again.  Flattened (g,pft) or (lon,lat,pft)
    lnd2iac_vars%npp(begg:endg,:) = x2z(index_x2z_Sl_npp,begg:endg)
    lnd2iac_vars%hr(begg:endg,:) = x2z(index_x2z_Sl_hr,:)

  end subroutine lnd_import
   !===============================================================================

  subroutine iac_export(iac2lnd_vars, iac2atm_vars, z2x)

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
    type(iac2lnd_type), intent(inout) :: iac2lnd_vars ! gcam to land output
    type(iac2atm_type), intent(inout) :: iac2atm_vars ! gcam to atm output
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid

  end subroutine iac_export
end module iac_import_export
