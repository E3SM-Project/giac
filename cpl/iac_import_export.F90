module iac_import_export
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Module to deal with importing and exporting coupled vars with
  !  internal iac variables.  
  ! !uses:
  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use iac_data_mod
  use gcam_cpl_indices
  use mct_mod
  use iac_data_mod,   only : iac_ctl
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
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use shr_kind_mod     , only: r8 => shr_kind_r8, CL => shr_kind_CL
    use netcdf
    !
    ! !ARGUMENTS:
    real(r8), intent(in)    :: x2z(:,:) ! driver import state to iacmodel
    type(lnd2iac_type) , intent(inout) :: lnd2iac_vars    ! gcam
    !
    ! LOCAL VARIABLES
    integer :: n,n1,p,g,i,j
    integer :: begg, endg
    character(len=32), parameter :: sub = '(iac_import)'

    ! Gcam expects things in npp_m[lon][lat][pft] format, so we need
    ! to extract from the flattened column representation.
    
    ! Once again, domain decomp isn't meaningful with one proc, but keep it for
    ! consistency with everybody else
    begg = iac_ctl%begg
    endg = iac_ctl%endg

    ! So, the idea here is to take the 17 (npft) fields for each of
    ! hr, npp, and pftwgtg and map them into a single 1D array of
    ! nlat*nlon*npft for each of those.  So we first loop over npft,
    ! then extract the global-indexed data, then pack them into the
    ! lnd2iac_vars(:,:) structure.

    ! pftwgt is in actual fraction of grid cell and is properly
    !    addressed in the gcam setdensity function

    do p=1,iac_ctl%npft

       ! This for logging purposes...
       !write(pftstr,'(I0)') p
       !pftstr=trim(pftstr)

       ! Now loop over our global index...
       ! g=1,ngrid - one domain.
       do g=iac_ctl%begg,iac_ctl%endg
          i=iac_ctl%ilon(g)
          j=iac_ctl%jlat(g)  ! extract lat, lon indeces from g

          ! i (lon) varies fastest, p slowest:
          ! n=i+nlon*(j-1)+nlat*nlon*(p-1)
          lnd2iac_vars%hr(i,j,p) = x2z(index_x2z_Sl_hr(p),g)
          lnd2iac_vars%npp(i,j,p) = x2z(index_x2z_Sl_npp(p),g)
          lnd2iac_vars%pftwgt(i,j,p) = x2z(index_x2z_Sl_pftwgt(p),g)
       end do ! global index g
    end do ! pft index p

  end subroutine iac_import
   !===============================================================================

  subroutine iac_export(iac2lnd_vars, iac2atm_vars, z2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the clm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use gcam_var_mod       , only : iulog
    !
    ! !ARGUMENTS:
    implicit none
    type(iac2lnd_type), intent(inout) :: iac2lnd_vars ! gcam to land output
    type(iac2atm_type), intent(inout) :: iac2atm_vars ! gcam to atm output
    real(r8),    intent(out)   :: z2x(:,:)! land to coupler export state on land grid

    ! LOCAL VARIABLES
    integer :: n,n1,g,i,j,p,m
    integer :: begg, endg
    character(len=32), parameter :: sub = '(iac_export)'

    ! Once again, domain decomp isn't meaningful with one proc, but keep it for
    ! consistency with everybody else
    begg = iac_ctl%begg
    endg = iac_ctl%endg

! avd
write(iulog,*) trim(sub),' export pft and harvest data'

    ! g=1,ngrid - one domain.
    do g=iac_ctl%begg,iac_ctl%endg
       i=iac_ctl%ilon(g)
       j=iac_ctl%jlat(g)  ! extract lat, lon indeces from g

       ! Co2flux to atm
       ! Convention has fluxes negative from lnd to atm, so we
       ! assume for iac to atm as well
       !z2x(index_z2x_Fazz_fco2_iac,g) = -iac2atm_vars%co2emiss(i,j)
       ! Monthly sfc, low alt air, and high alt air values
       do m=1,12
          z2x(index_z2x_Fazz_co2sfc_iac(m),g) = -iac2atm_vars%co2sfc(m,i,j)
          z2x(index_z2x_Fazz_co2airlo_iac(m),g) = -iac2atm_vars%co2airlo(m,i,j)
          z2x(index_z2x_Fazz_co2airhi_iac(m),g) = -iac2atm_vars%co2airhi(m,i,j)
       end do

       ! these iac values are percent/fraction of veg land unit, but for proper
       !   coupling they need to be converted to fraction of actual grid cell

       ! Now the 17 iac->lnd pfts
       ! need the new pfts and the previous pfts
       do p=1,iac_ctl%npft
          z2x(index_z2x_Sz_pct_pft(p),g) = iac2lnd_vars%pct_pft(i,j,p) / 100.0_R8 * &
             iac_ctl%iacmask(i,j) * iac_ctl%vegfrac(i,j) * iac_ctl%landfrac(i,j)
          z2x(index_z2x_Sz_pct_pft_prev(p),g) = &
             iac2lnd_vars%pct_pft_prev(i,j,p) / 100.0_R8 * iac_ctl%iacmask(i,j) * &
             iac_ctl%vegfrac(i,j) * iac_ctl%landfrac(i,j)

       end do ! pft index p

       ! Now the 5 harvest fields
       do p=1,iac_ctl%nharvest
          z2x(index_z2x_Sz_harvest_frac(p),g) = &
             iac2lnd_vars%harvest_frac(i,j,p) * iac_ctl%iacmask(i,j) * &
             iac_ctl%vegfrac(i,j) * iac_ctl%landfrac(i,j) 

       end do ! harvest index p

    end do ! global index g

  end subroutine iac_export
end module iac_import_export
