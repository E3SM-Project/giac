module iac_data_mod

  !------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module contains the gcam coupler variable and structure definitions
  ! I'm putting them all in one place, both input, output, and
  ! internal, to make it easier to find, but we might split them later
  !------------------------------------------------------------
  use shr_kind_mod               , only : r8 => shr_kind_r8

  implicit none
  save

  ! control structure for internal use
  type, public :: iac_ctl_type

     ! Decomposition - just in case
     integer :: begg,endg  ! Begin and end grid cells for this processor
     integer, allocatable :: gindex(:)  ! mapping grid cells to global index
     integer, allocatable :: long(:)    ! longitude of grid cell
     integer, allocatable :: latg(:)    ! latitude of grid cell
     integer, allocatable :: area(:)    ! area of grid cell (?)
  end type iac_ctl_type

  ! These might need to be pointers rather than allocatable, with => null().
  type, public :: lnd2iac_type
     ! The npp and hr are the inputs from clm
     real(r8), allocatable :: npp(:,:)
     real(r8), allocatable :: hr(:,:)
  end type lnd2iac_type

  type, public :: iac2lnd_type
     real(r8), allocatable :: foo(:,:)
  end type iac2lnd_type

  type, public :: iac2atm_type
     real(r8), allocatable :: co2emiss(:,:)
  end type iac2atm_type

  public

  type(lnd2iac_type) :: lnd2iac_vars
  type(iac2lnd_type) :: iac2lnd_vars
  type(iac2atm_type) :: iac2atm_vars
  type(iac_ctl_type) :: iac_ctl

end module iac_data_mod


     
     
