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

     ! Dimensioned by grid cell - ngrid
     integer, allocatable :: gindex(:)  ! mapping grid cells to global index
     integer, allocatable :: iacmask(:) ! mask

     integer, allocatable :: ilon(:)   ! Index of lat,lon dimension,
     integer, allocatable :: jlat(:)   ! for grid cell

     ! Dimensioned (lon) or (lat) index, respectively
     ! It's easier for gcam integration if we store all our variables
     ! multidimensionally in (lon,lat,pft) order, so the only thing
     ! on the global grid is the indeces to go back and forth between
     ! (g) <-> (i,j)
     real(r8), allocatable :: lon(:)    ! longitude
     real(r8), allocatable :: lat(:)    ! latitude 

     ! Various sizes
     integer, public :: ngrid
     integer, public :: npft
     integer, public :: nlat
     integer, public :: nlon

  end type iac_ctl_type

  ! Making these structures, since it's easier to make sure
  ! everything is contiguous, which is important when coupling with
  ! C++ 
  type, public :: lnd2iac_type
     ! Input from lnd, dimensioned (lon,lat,pft)
     real(r8), allocatable :: npp(:,:,:)
     real(r8), allocatable :: hr(:,:,:)
     real(r8), allocatable :: pftwgt(:,:,:)

     ! These will be read once from configuration file, and are (lon,lat)
     real(r8), allocatable :: landfrac(:,:)
     real(r8), allocatable :: area(:,:)
  end type lnd2iac_type

  type, public :: iac2lnd_type
     ! Output to lnd, dimensioned (lon,lat,pft)
     real(r8), allocatable :: pct_pft(:,:,:)
  end type iac2lnd_type

  type, public :: iac2atm_type
     ! output to atm, dimensioned (lon,lat)
     real(r8), allocatable :: co2emiss(:,:)
  end type iac2atm_type

  type(lnd2iac_type) :: lnd2iac_vars
  type(iac2lnd_type) :: iac2lnd_vars
  type(iac2atm_type) :: iac2atm_vars
  type(iac_ctl_type) :: iac_ctl

end module iac_data_mod


     
     
