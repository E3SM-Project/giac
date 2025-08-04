!#define DEBUG
Module gcam2glm_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: gcam2glm_mod
!
!  Interface of the human component in E3SM
!
! !DESCRIPTION:
!
! !USES:

  use iac_data_mod, only : cdata => gdata, EClock => GClock
  use iac_data_mod
  use shr_file_mod, only: shr_file_getunit, shr_file_freeunit
  use gcam_var_mod
  use shr_cal_mod
  use shr_sys_mod
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4
  use netcdf
  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: gcam2glm_init_mod               ! clm initialization
  public :: gcam2glm_run_mod                ! clm run phase
  public :: gcam2glm_final_mod              ! clm finalization/cleanup
  public :: handle_err
  public :: fround
  private :: D_mrgrnk

! !PUBLIC DATA MEMBERS: 
    integer, parameter :: kdp = selected_real_kind(15)
    real(r8), dimension(:, :), allocatable :: hydeGCROP2015,&
         hydeGPAST2015,   &
         hydeGWH2015,     &
         hydeGOTHR2015,   &
         hydeGSECD2015,   &
         cellarea,        &
         cellarea_forest, &
         cellarea_nonforest,&
         glm_crop_ann,    &
         glm_past_ann,    &
         glm_othr_ann,    &
         fnfforest,       &
         fnfnonforest,    &
         pot_veg,         &
         pot_veg_rev,         &
         crop_area,       &
         avail_land0,     &
         avail_landA,     &
         gcam_past,       &
         gcam_forest_area,&
         gcam_crop,       &
         gcam_farea

    real(r8), dimension(:), allocatable :: unmet_neg_past,unmet_neg_crop,unmet_farea, &
         cumsum_sorted_farea

    real(r8), dimension(:, :), allocatable,save :: gcam_wh, &
         pctland_in2015,glu_weights_rev,sortsitesup,sortsitesdn

    integer, dimension(:,:), allocatable      ::  rglus
    real(r8), dimension(:, :, :), allocatable,save :: glm_crop, &
         glm_past

    real(r8), dimension(:,:,:), allocatable, save :: glu_weights

    integer, dimension(:), allocatable :: rgmin, rgmax

    integer,save :: year1,year2

    real(r8), dimension(:), allocatable :: datearr, gcam_wh_ann, glm_wh_ann
    real(r8), PUBLIC, dimension(:), allocatable :: lon,lat
    real(r8), dimension(:), pointer :: array1d
    real(r8)  :: forested_past, forested_past_percent
    real(r8)  :: nonforested_past, nonforested_past_percent
    real(r8)  :: forested_crop, forested_crop_percent
    real(r8)  :: nonforested_crop, nonforested_crop_percent

    real(r4)  :: miss_val = 1.0e36

    integer :: n,np1,nflds,gcamsize,nglu,nregions,io,r,g1
    integer :: lonx,latx,max_nglu
    real(r8) :: x,y,weight
    character(256) :: region_name, glu_name ! dummy vars to read csv
    integer, dimension(3) :: start3,count3
    integer, dimension(2) :: start2,count2
    integer :: ncid,tmp(1), &
         lonDimID, latDimId, timeDimId, &
         numTimes,    &
         status,GCROPVarId,timevarid,varid
    integer, PUBLIC :: numLats, numLons
    integer, dimension(nf90_max_var_dims) :: dimIDs
    character(len=*),parameter :: gcam2glm_restfile = 'gcam2glm_restart.'
    character(256) :: filename

! !REVISION HISTORY:
! Author: T Craig
! Author: A Di Vittorio: major revision to correct downscaling errors
! Note that there is still a problem with cells containing multiple land units
!    tracking only the fraction of the whole cell, with just the fixed land unit
!    weights, and not tracking the changing land unit weights of
!    fores/nonforest, causes errors
! Note that these erros are small (so far < 0.5%) for global crop/pasture area
!   and also for the area changes, that we are considering this acceptable
!   this is because the fix requires serious structural changes

! !PRIVATE DATA MEMBERS:
    real(r8), allocatable :: gcamo_base(:,:)

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam2glm_init_mod

! !INTERFACE:
  subroutine gcam2glm_init_mod()

! !DESCRIPTION:
! Initialize interface for glm

! !USES:
    use iac_data_mod
    use mct_mod
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    logical :: lexist, restart_run
    integer :: iun,tmpyears(2),ier,t,yy,g
    integer :: ymd, year, mon, day
    real(r8) :: v
    character(len=*),parameter :: subname='(gcam2glm_init_mod)'

    character(len=512) :: dum

! !REVISION HISTORY:
! Author: T Craig
! Author: JET          ! rewrite of matlat preprocessing script new_grids_matchforest4.m

!EOP
!-----------------------------------------------------------------------

#ifdef DEBUG
     write(iulog,*) subname,' starting subroutine '
     write(iulog,*) subname,'gcam_var_mod nsrest is', nsrest
#endif
    if(nsrest == nsrContinue .or. nsrest == nsrBranch) then
       restart_run = .true.
    else
       restart_run = .false.
    end if

    nregions=num_gcam_energy_regions
    nglu=num_gcam_land_regions

    gcamsize=nglu

! get current ehc clock time and extract year
! this is a year ahead of the e3sm model year
!    except during the ehc run timestep
ymd = EClock(iac_EClock_ymd)
call shr_cal_date2ymd(ymd,year,mon,day)

! initialize two level time indexes 

    n=1
    np1=2

! Variable to keep track of when to calculate the next GCAM time step.
! Initialize this with the year 2015 data - from the single glm initial file

    status= nf90_open(gcam2glm_baselu,nf90_nowrite,ncid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inq_varid(ncid, "gcrop", GCROPVarId)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_variable(ncid, GCROPVarId, dimids = dimIDs)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(1), len = numLons)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(2), len = numLats)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inquire_dimension(ncid, dimIDs(3), len = numTimes)
    if(status /= nf90_NoErr) call handle_err(status)
    allocate(lon(numLons), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate lon',ier)
    allocate(lat(numLats))
    if(ier/=0) call mct_die(subName,'allocate lat',ier)
    
    status = nf90_inq_varid(ncid, "lon", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,lon)
    if(status /= nf90_NoErr) call handle_err(status)
    
    status = nf90_inq_varid(ncid, "lat", varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,lat)
    if(status /= nf90_NoErr) call handle_err(status)
    
    allocate(hydeGCROP2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate hydeGCROP2015',ier)
    allocate(hydeGPAST2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate hydeGPAST2015',ier)
    allocate(hydeGOTHR2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate hydeGOTHR2015',ier)
    allocate(hydeGSECD2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate hydeGSECD2015',ier)
    allocate(hydeGWH2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate hydeGWH2015',ier)
    allocate(cellarea(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cellarea',ier)
    allocate(cellarea_forest(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cellarea_forest',ier)
    allocate(cellarea_nonforest(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cellarea_nonforest',ier)
    allocate(glm_crop_ann(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_crop_ann',ier)
    allocate(glm_past_ann(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_past_ann',ier)
    allocate(glm_othr_ann(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_othr_ann',ier)
    allocate(cumsum_sorted_farea(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cumsum_sorted_farea',ier)
    allocate(glm_wh_ann(nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_wh_ann',ier)
    allocate(fnfforest(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate fnfforest',ier)
    allocate(fnfnonforest(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate fnfnonforest',ier)
    allocate(pot_veg(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate pot_veg',ier)
    allocate(pot_veg_rev(numLats, numLons), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate pot_veg_rev',ier)
    allocate(crop_area(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate crop_area',ier)
    allocate(pctland_in2015(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate pctland_in2015',ier)
    allocate(sortsitesup(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortsitesup',ier)
    allocate(sortsitesdn(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortsitesdn',ier)
    allocate(datearr(numTimes), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate datearr',ier)
    allocate(glm_crop(numLons, numLats, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_crop',ier)
    allocate(glm_past(numLons, numLats, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glm_past',ier)
    allocate(gcam_crop(gcamsize, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcam_crop',ier)
    allocate(gcam_wh(gcamsize, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcam_wh',ier)
    allocate(gcam_past(gcamsize, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcam_past',ier)
    allocate(gcam_forest_area(gcamsize, 2), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcam_forest_area',ier)
    allocate(unmet_neg_past(gcamsize), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate unmet_neg_past',ier)
    allocate(unmet_neg_crop(gcamsize), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate unmet_neg_crop',ier)
    allocate(unmet_farea(gcamsize), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate unmet_farea',ier)
    allocate(avail_land0(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate avail_land0',ier)
    allocate(avail_landA(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate avail_landA',ier)

    allocate(glu_weights(gcamsize, numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glu_weights',ier)
    allocate(glu_weights_rev(numLats, numLons), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glu_weights_rev',ier)
    allocate(rglus(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rglus',ier)

    allocate(rgmin(nregions), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rgmin',ier)
    allocate(rgmax(nregions), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate rgmax',ier)

    allocate(gcamo_base(num_iac2elm_landtypes,nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gcamo_base',ier)

    glm_crop=iac_spval
    glm_past=iac_spval
    avail_land0=iac_spval
    avail_landA=iac_spval
    hydeGCROP2015=iac_spval
    hydeGPAST2015=iac_spval
    hydeGOTHR2015=iac_spval
    hydeGSECD2015=iac_spval
    hydeGWH2015=iac_spval
    cellarea=iac_spval
    cellarea_forest=iac_spval
    cellarea_nonforest=iac_spval
    pctland_in2015=iac_spval
    datearr=iac_spval
    gcamo_base=0.


    status = nf90_inq_varid(ncid,'time',timeVarId)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,timeVarId,datearr)
    if(status /= nf90_NoErr) call handle_err(status)
    start3(1)=1
    count3(1)=numLons
    start3(2)=1
    count3(2)=numLats
    ! there should be only one 2015 record in this file
    tmp(1)=0
    do t=1,size(datearr)
       if(datearr(t) == 2015) then
          tmp(1) = t
          exit
       endif
    enddo
    if(tmp(1)==0) then
       write(iulog,*) subname,' ERROR: start year 2015 not found in file ', gcam2glm_baselu
    endif

    start3(3)=tmp(1)
    count3(3)=1
    start2=1
    count2(1)=numLons
    count2(2)=numLats
    ! read in luh2 data - now from the single glm init file

    status = nf90_inq_varid(ncid,'gcrop',varid)
    status = nf90_get_var(ncid,varid,hydeGCROP2015,start3,count3)
    status = nf90_inq_varid(ncid,'gpast',varid)
    status = nf90_get_var(ncid,varid,hydeGPAST2015,start3,count3)
    status = nf90_inq_varid(ncid,'gothr',varid)
    status = nf90_get_var(ncid,varid,hydeGOTHR2015,start3,count3)
    status = nf90_inq_varid(ncid,'gsecd',varid)
    status = nf90_get_var(ncid,varid,hydeGSECD2015,start3,count3)
    status = nf90_inq_varid(ncid,'cell_area',varid)
    status = nf90_get_var(ncid,varid,cellarea,start3,count3)
    status = nf90_close(ncid)
    cellarea=cellarea/1.e6

    status = nf90_open(gcam2glm_basebiomass,nf90_nowrite,ncid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inq_varid(ncid,'biomass',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_get_var(ncid,varid,pot_veg,start2,count2)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Read in glumap csv file
    ! the glu_weights are the proportions of the cell covered by each land unit 
    ! the weights sum to 1 within each cell
    ! there were calculated based on the proportions of land covered by each land unit

    glu_weights=0.
    rgmin = 1000
    rgmax = 0

    open(5,file=gcam2glm_glumap)
    ! Read header
    read(5,*) dum
    do
       read(5,*,iostat=io) r,g,x,y,region_name,glu_name,weight
       if (io < 0) then
          exit
       endif

       ! Need to find the lon and lat index for these x,y (which are lon,lat values at center of cell)
       ! 70s way, cuz findloc and minloc and etc are messed up
       do lonx=1,numLons
          if (x .eq. lon(lonx)) then 
             exit
          endif
       end do
       
       do latx=1,numLats
          if (y .eq. lat(latx)) then 
             exit
          endif
       end do
     
       ! Basic checking
       if (glu_weights(g,lonx,latx) .ne. 0.) then
          write(iulog,*) subname,' ERROR: multiple weights for ',g,lonx,latx
          call shr_sys_abort(subname//' ERROR: weights file')
       endif
       glu_weights(g,lonx,latx) = weight
       where ( glu_weights(g,:,:) < 0. )
          glu_weights(g,:,:) = 0.
       end where

       ! Find region of glus for each r
       if (g .lt. rgmin(r)) rgmin(r) = g
       if (g .gt. rgmax(r)) rgmax(r) = g

    end do

    close(5)
    
    cellarea_nonforest(:,:)=0.
    cellarea_forest(:,:)=0.
    pot_veg=pot_veg*0.75
    pctland_in2015(:,:) = 0.
    pctland_in2015=hydeGCROP2015 + hydeGPAST2015 + hydeGOTHR2015 + hydeGSECD2015


    if (.not. restart_run) then

        ! if gcam_spinup == .true. this file does not exist yet
        if (.not. gcam_spinup) then
           ! Read in gcam base file 
           open(5,file=base_gcam_lu_wh_file)
           ! Read header
           read(5,*) dum
           do
             read(5,*,iostat=io) g,t,yy,v
             if (io < 0) then
               exit
             endif
             gcamo_base(t,g) = v
          end do
          close(5)
       end if

       glm_crop(:,:,n)=hydeGCROP2015
       glm_past(:,:,n)=hydeGPAST2015
       glm_crop(:,:,np1)=hydeGCROP2015
       glm_past(:,:,np1)=hydeGPAST2015

       ! Set initial gcam land
       ! current baseline has been converted from billion m^3 of biomass
       !   to the needed MgC here
       !      conv fact is 0.250 tonnes C per m^3 (MgC per m^3); same as in gcam
       ! but do it only if not gcam spinup
       if (.not. gcam_spinup) then
          gcam_crop(:,n) = gcamo_base(iac_gcamo_crop,:)
          gcam_past(:,n) = gcamo_base(iac_gcamo_pasture,:)
          gcam_wh(:,n) = gcamo_base(iac_gcamo_woodharv,:)
          gcam_forest_area(:,n) = gcamo_base(iac_gcamo_forest,:)
       end if

       ! Set years 
       cdata%i(iac_cdatai_gcam_yr1)=iac_start_year
       cdata%i(iac_cdatai_gcam_yr2)=iac_start_year+iac_gcam_timestep
    else
       ! read restart and set crop and past
     
       ! the correct restart file is associated with the current ehc year
       !    which matches e3sm model year only during the ehc run timestep
       !    otherwise it is a year ahead because of the ehc clock advance
       ! this is because the ehc clock is advanced at the end of its run call 
       ! this ensures the correct data are loaded for the next ehc run call
       ! note that init is called no matter the start/restart timestep

       write(filename,'(a,i4.4,a,i2.2,a)') trim(gcam2glm_restfile)//'r.',year,'.nc'

       inquire(file=trim(filename),exist=lexist)
       if (lexist) then

           write(iulog,*) subname,' read restart file ',trim(filename)          
          
          status= nf90_open(filename,nf90_nowrite,ncid)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_inq_varid(ncid,'gcam_crop',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,gcam_crop)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_inq_varid(ncid,'gcam_past',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,gcam_past)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_inq_varid(ncid,'gcam_wh',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,gcam_wh)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_inq_varid(ncid,'gcam_forest_area',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,gcam_forest_area)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_inq_varid(ncid,'gcam_years',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,tmpyears)
          if(status /= nf90_NoErr) call handle_err(status)
	  cdata%i(iac_cdatai_gcam_yr1)=tmpyears(1)
	  cdata%i(iac_cdatai_gcam_yr2)=tmpyears(2)
          
          status = nf90_inq_varid(ncid,'glm_crop',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,glm_crop)
          if(status /= nf90_NoErr) call handle_err(status)
          
          status = nf90_inq_varid(ncid,'glm_past',varid)
          if(status /= nf90_NoErr) call handle_err(status)
          status = nf90_get_var(ncid,varid,glm_past)
          if(status /= nf90_NoErr) call handle_err(status)

          status = nf90_close(ncid)
       else
          write(iulog,*) subname,' ERROR: restart file NOT found ',trim(filename)
          call shr_sys_abort(subname//' ERROR: missing file')
       end if ! restart file exist
    end if ! restart run
  end subroutine gcam2glm_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: gcam2glm_run_mod

! !INTERFACE:
  subroutine gcam2glm_run_mod(gcamo, glmi, glmi_wh)

! !DESCRIPTION:
! Run interface for glm

! !USES:
    use iac_data_mod
    use mct_mod
    implicit none

! !ARGUMENTS:
    real(r8), pointer :: gcamo(:,:)
    real(r8), pointer :: glmi(:,:)
    real(r8), pointer :: glmi_wh(:)

! !PARAMETERS:

    character(len=*),parameter :: subname='(gcam2glm_run_mod) '

    real(r8), parameter :: crop_forest_abandon_percent = 0.9
    real(r8), parameter :: past_forest_abandon_percent = 0.9

! !LOCAL VARIABLES:

    character(len=512) :: dum
    character*4 :: yearc
    character(256) :: filename
    integer :: i,j,ij,r,i1,j1,aez,h,z
    integer :: row,g,t,y,yy
    integer :: iun,iyr,ier
    integer :: ymd, tod, dt,naez,nreg,ii,year,mon,day
    logical :: restart_run
    real(r8)  :: crop_d,past_d,crop_neg,crop_pos,past_neg,past_pos,farea_d,v
    real(r8)  :: gcam_crop_tmp(2,18,14),gcam_past_tmp(2,18,14),gcam_forest_area_tmp(2,18,14)
    real(r8)  :: fact1,fact2,eclockyr,fact1yrp1, fact2yrp1,delyr,eclockyrp1
    real(r8)  :: tmp0
    integer,save :: ncid,varid,dimid,dimid3(3), dimid_gcam(2)
    integer :: totrglus, nrglu
    integer, allocatable  ::indxdn(:),indxup(:),sortlatsup(:),sortlonsup(:),sortlatsdn(:),sortlonsdn(:),indxa(:),indxadn(:),v1u(:),v2u(:),v1d(:),v2d(:)
    real(r8), allocatable   :: tmparr(:)
    integer :: sortind(1),sortind2(1),max_aez_ind(1),regional_unmet_reassign
    real(r8)  :: crop_pos_f,crop_pos_nf,pot_forest_to_crop,pot_forest_from_crop,reassign_crop, &
               crop_neg_f,crop_neg_nf,crop_nfarea,crop_farea, &
               past_pos_f,past_pos_nf,pot_forest_to_past,pot_forest_from_past,reassign_past, &
               past_neg_f,past_neg_nf,past_nfarea,past_farea,GLM_nfarea,GLM_farea,regional_farea_needed, &
               crop_after_decrease,past_after_decrease,crop_before_decrease,past_before_decrease, &
               total_ag_decrease,crop_decrease_ratio,past_decrease_ratio,ag_area_avail,crop_area_avail, &
               final_area_needed,f_diff,reassign_ag_at_max_aez_ind,sumavail_land0,sumavail_landA

   real(r8) :: pot_nonforest_from_crop, pot_nonforest_from_past, ag_farea_new, &
               pot_nonforest_to_crop, pot_nonforest_to_past, GLM_nfarea_new

   real(r8), allocatable    ::  avail_farea(:),avail_nfarea(:),avail_ag_farea(:),reassign_ag(:), &
                              unmet_aez_farea(:),cumsum_sorted_reassign_ag(:),unmet_regional_farea(:)
   integer ntimes,nntimes,zz

   integer :: nmode, ierr
   character(len=128) :: hfile

   real(r8) :: curr_crop, curr_past, gcam_glb_forest_change, gcam_glb_past_change
   real(r8) :: curr_crop_g, curr_past_g, temp_crop, temp_past
   real(r8) :: start_crop, start_crop_g, start_past, start_past_g

   integer,allocatable :: glatinds(:), gloninds(:)
   integer :: count_gcells
   logical :: found
   real(r8), allocatable :: temp_past_grid(:,:)
   real(r8), parameter :: diff_tol = 1.e-8

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------
    regional_unmet_reassign=1
    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)
    eclockyr=ymd/10000
    ! this is .true. if gcam is running this model year
    write(iulog,*) trim(subname),'model year is ', eclockyr, 'gcam alarm is ', gcam_alarm
#ifdef DEBUG
    write(iulog,*) trim(subname),' date= ',ymd,tod
    write(iulog,*) trim(subname), 'gcam_var_mod nsrest is', nsrest
#endif

    if(nsrest == nsrContinue .or. nsrest == nsrBranch) then
       restart_run  = .true.
    else
       restart_run = .false.
    end if
    year1=cdata%i(iac_cdatai_gcam_yr1)
    year2=cdata%i(iac_cdatai_gcam_yr2)
    
    nglu=num_gcam_land_regions
    nreg=num_gcam_energy_regions

    ! Maxiumum number of glus in a region
    max_nglu=iac_max_nglu

    allocate(indxup(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate indxup',ier)
    allocate(indxdn(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate indxdn',ier)
    allocate(sortlatsup(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortlatsup',ier)
    allocate(sortlatsdn(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortlatsdn',ier)
    allocate(sortlonsup(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortlonsup',ier)
    allocate(sortlonsdn(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate sortlonsdn',ier)
    allocate(tmparr(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate tmparr',ier)
    allocate(indxa(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate indxa',ier)
    allocate(indxadn(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate indxadn',ier)
    allocate(avail_farea(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate avail_farea',ier)
    allocate(avail_nfarea(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate avail_nfarea',ier)
    allocate(avail_ag_farea(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate avail_ag_farea',ier)
    allocate(reassign_ag(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate reassign_ag',ier)
    allocate(unmet_aez_farea(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate unmet_aez_farea',ier)
    allocate(cumsum_sorted_reassign_ag(max_nglu), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate cumsum_sorted_reassign_ag',ier)
    allocate(unmet_regional_farea(nreg), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate unmet_regional_farea',ier)

    allocate(glatinds(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate glatinds',ier)
    allocate(gloninds(numLons*numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate gloninds',ier)

    allocate(temp_past_grid(numLons, numLats), stat=ier)
    if(ier/=0) call mct_die(subName,'allocate temp_past_grid',ier)

    avail_farea=0.
    avail_nfarea=0.
    avail_ag_farea=0.
    reassign_ag=0.
    unmet_aez_farea=0.
    cumsum_sorted_reassign_ag=0.
    unmet_regional_farea=0.
    unmet_farea=0.

    ! this should only do the processing when gcam runs
    ! this is because nothing changes until gcam runs
    ! the annual interpolation is after this conditional if
    if (gcam_alarm) then

    ! the glm np1 values are initially set to the n state, as needed
    ! this should always be the case in a gcam model run year
    ! for intermediate years (gcam_alarm==.false.) the n and np1 values do not change,
    ! and if there is a restart the appropriate values are read from the restart file

#ifdef DEBUG
    write(iulog,*) subname, 'these n and np1 should initially be the same'
    write(iulog,*) 'sum start n glm_crop='
    write(iulog,fmt="(1ES25.15)") sum(glm_crop(:,:,n))
    write(iulog,*) 'sum start n glm_past='
    write(iulog,fmt="(1ES25.15)") sum(glm_past(:,:,n))
    write(iulog,*) 'sum start np1 glm_crop='
    write(iulog,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(iulog,*) 'sum start np1 glm_past='
    write(iulog,fmt="(1ES25.15)") sum(glm_past(:,:,np1))
#endif     
    ! Unpack gcamo field
    ! the previous field (n) is initialized by file or by previous year, so read gcamo into next (np1) field
   
    ! if gcam_spinup == .true. the base data need to be read now to initialize
    !    but only if it is the first model year 
    !    and only if this is not a restart run, in which the data are filled
    !    from restart file on initialization
    if (gcam_spinup .and. (.not. restart_run) .and. eclockyr == iac_start_year) then
       ! Read in gcam base file 
       open(5,file=base_gcam_lu_wh_file)
       ! Read header
       read(5,*) dum
       do
          read(5,*,iostat=io) g,t,yy,v
          if (io < 0) then
             exit
          endif
          gcamo_base(t,g) = v
       end do
       close(5)

       gcam_crop(:,n) = gcamo_base(iac_gcamo_crop,:)
       gcam_past(:,n) = gcamo_base(iac_gcamo_pasture,:)
       gcam_wh(:,n) = gcamo_base(iac_gcamo_woodharv,:)
       gcam_forest_area(:,n) = gcamo_base(iac_gcamo_forest,:)
    end if
 
    ! avd - write the harvest data to a diag file
    !write(hfile,'(a)') 'gcam2glm_harvest.nc'
    !ierr = nf90_create(trim(hfile),nf90_clobber,ncid)
    !ierr = nf90_def_dim(ncid,'num_gcam_units',gcamsize,dimid)
    !ierr = nf90_def_var(ncid,'gcam_harvest',NF90_DOUBLE,dimid,varid)
    !ierr = nf90_enddef(ncid)
    !ierr = nf90_put_var(ncid,varid,gcamo(iac_gcamo_woodharv,:))
    !if(ierr /= nf90_NoErr) call handle_err(ierr)
    !ierr = nf90_close(ncid)

    ! avd - write harvest data to the log
    !do g=1,gcamsize
    !  write(iulog,*) subname,' gunit ', g, ' gcamo wood harvest =', &
    !                 gcamo(iac_gcamo_woodharv,g)
    !end do

    ! note that the GCAM coupling now converts wh to MgC

    gcam_crop(:,np1) = gcamo(iac_gcamo_crop,:)
    gcam_past(:,np1) = gcamo(iac_gcamo_pasture,:)
    gcam_wh(:,np1) = gcamo(iac_gcamo_woodharv,:)
    gcam_forest_area(:,np1) = gcamo(iac_gcamo_forest,:)

    write(iulog,*) subname, 'initial global change in GCAM crop area km^2 = ', &
       real(nint(sum(1000*gcam_crop(:,np1))*100))/100. - real(nint(sum(1000*gcam_crop(:,n))*100))/100.
    write(iulog,*) subname, 'initial global change in GCAM pasture area km^2 = ', &
       sum(real(nint(1000*gcam_past(:,np1)*100))/100.) - sum(real(nint(1000*gcam_past(:,n)*100))/100.)
    write(iulog,*) subname, 'initial global change in GCAM forest area km^2 = ', &
       sum(real(nint(1000*gcam_forest_area(:,np1)*100))/100.) - sum(real(nint(1000*gcam_forest_area(:,n)*100))/100.)

    write(iulog,*) subname, 'initial n global GCAM crop area km^2 = ', &
       real(nint(sum(1000*gcam_crop(:,n))*100))/100.
    write(iulog,*) subname, 'initial np1 global GCAM crop area km^2 = ', &
       real(nint(sum(1000*gcam_crop(:,np1))*100))/100.

    write(iulog,*) subname, 'initial n global GCAM pasture area km^2 = ', &
       sum(real(nint(1000*gcam_past(:,n)*100))/100.)
    write(iulog,*) subname, 'initial np1 global GCAM pasture area km^2 = ', &
       sum(real(nint(1000*gcam_past(:,np1)*100))/100.)

    write(iulog,*) subname, 'initial n global GCAM forest area km^2 = ', &
       sum(real(nint(1000*gcam_forest_area(:,n)*100))/100.)
    write(iulog,*) subname, 'initial np1 global GCAM forest area km^2 = ', &
       sum(real(nint(1000*gcam_forest_area(:,np1)*100))/100.)


!avd - write some data to the log
!do g=1,gcamsize
!   write(iulog,*) subname,' gunit ', g, ' gcam_wh np1 =', &
!                  gcam_wh(g,np1)
!   write(iulog,*) subname,' gunit ', g, ' gcam_wh n =', &
!                  gcam_wh(g,n)
!end do

! avd - test this by setting no change
!    gcam_crop(:,np1) = gcamo_base(iac_gcamo_crop,:)
!    gcam_past(:,np1) = gcamo_base(iac_gcamo_pasture,:)
!    gcam_wh(:,np1) = gcamo_base(iac_gcamo_woodharv,:)
!    gcam_forest_area(:,np1) = gcamo_base(iac_gcamo_forest,:)

#ifdef DEBUG
    write(iulog,*) subname,' year1 year2 ',year1,year2
#endif
    call shr_sys_flush(iulog)

    gcam_glb_forest_change = 0.
    gcam_glb_past_change = 0.

    do g=1,nglu

#ifdef DEBUG       
       write(iulog,*) subname, '***** land unit index for transitions = ', g

       write(iulog,*) subname, 'initial n g GCAM crop area km^2 = ', &
          real(nint(1000*gcam_crop(g,n)*100))/100.
       write(iulog,*) subname, 'initial np1 g GCAM crop area km^2 = ', &
          real(nint(1000*gcam_crop(g,np1)*100))/100.

       write(iulog,*) subname, 'initial n g GCAM pasture area km^2 = ', &
          real(nint(1000*gcam_past(g,n)*100))/100.
       write(iulog,*) subname, 'initial np1 g GCAM pasture area km^2 = ', &
          real(nint(1000*gcam_past(g,np1)*100))/100.

       write(iulog,*) subname, 'initial n g GCAM forest area km^2 = ', &
          real(nint(1000*gcam_forest_area(g,n)*100))/100.
       write(iulog,*) subname, 'initial np1 g GCAM forest area km^2 = ', &
          real(nint(1000*gcam_forest_area(g,np1)*100))/100.
#endif

       ! convert to km^2 and Round GCAM results to 2 decimal places before using as LUC
       crop_d = real(nint(1000*gcam_crop(g,np1)*100))/100.-real(nint(1000*gcam_crop(g,n)*100))/100.
       past_d = real(nint(1000*gcam_past(g,np1)*100))/100.-real(nint(1000*gcam_past(g,n)*100))/100.
       farea_d = real(nint(1000*gcam_forest_area(g,np1)*100))/100.-real(nint(1000*gcam_forest_area(g,n)*100))/100.

       gcam_glb_forest_change = gcam_glb_forest_change + farea_d
       gcam_glb_past_change = gcam_glb_past_change + past_d

       ! convert area changes into positive and negative changes

       if (crop_d <= 0.) then
          crop_neg = -1.*crop_d
          crop_pos = 0.
       else
          crop_neg = 0.
          crop_pos = crop_d
       end if
       if (past_d <= 0.) then
          past_neg = -1.*past_d
          past_pos = 0.
       else
          past_neg = 0.
          past_pos = past_d
       end if

       ! the glu_weights are the proportions of each land unit in each cell (they sum to 1 in a cell)
       ! these were calculated based on their respective proportions of land area in each cell
       ! using these and the land area fraction we can determine the actual
       ! areas of different components in a cell and apply the appropriate changes

       ! TRS - we are now looping by glu.  So instead of finding this
       ! specific combo of aez and region, each glu implies a region.

       ! Find lon/lat where we have any weights at all
       rglus=0
       where ( glu_weights(g,:,:) .gt. 0. )
          rglus=1
       end where
       totrglus=sum(rglus)      
       !write(iulog,*), subname, 'totrglus = ', totrglus 
       glatinds = 0
       gloninds = 0
       count_gcells = 0
       do j=1,numLats
          do i=1,numLons
             if (rglus(i,j) == 1) then
                count_gcells = count_gcells+1
                glatinds(count_gcells) = j
                gloninds(count_gcells) = i
             endif
          enddo
       enddo
       if (count_gcells /= totrglus) then
          write(iulog,*) subname, 'Error: count_gcells ', count_gcells, ' != totrglus ', totrglus
       endif

       ! note that a cell is potentially forest or non-forest; it is never mixed
       ! fnf values are the glu weights for this land unit, masked by forest or non-forest
       ! cellarea_forest and cellarea_nonforest are the fractions of cell area
       ! in this land unit, masked by forest or non-forest
       cellarea_forest(:,:)=0.
       cellarea_nonforest(:,:)=0.
       fnfnonforest=0.
       fnfforest=0.
       where ( pot_veg > 1.)
          cellarea_forest(:,:)=cellarea(:,:)
          fnfforest(:,:)=1.
       elsewhere
          fnfnonforest(:,:)=1.
          cellarea_nonforest(:,:)=cellarea(:,:)
       end where
       fnfforest=fnfforest*glu_weights(g,:,:)
       fnfnonforest=fnfnonforest*glu_weights(g,:,:)
       cellarea_nonforest=cellarea_nonforest*fnfnonforest
       cellarea_forest=cellarea_forest*fnfforest

       ! these always need to start with the current np1 state
       ! otherwise cells with multiple land units will not be updated correctly
       ! np1 values initially equal n
       crop_farea =sum(glm_crop(:,:,np1)*cellarea_forest)
       crop_nfarea =sum(glm_crop(:,:,np1)*cellarea_nonforest)
       past_farea =sum(glm_past(:,:,np1)*cellarea_forest)
       past_nfarea =sum(glm_past(:,:,np1)*cellarea_nonforest)
       GLM_nfarea = sum((pctland_in2015 - glm_crop(:,:,np1) - glm_past(:,:,np1))*cellarea_nonforest)
       GLM_farea = sum((pctland_in2015 - glm_crop(:,:,np1) - glm_past(:,:,np1))*cellarea_forest)
       if (allocated(v1u)) deallocate(v1u)
       if (allocated(v2u)) deallocate(v2u)
       if (allocated(v1d)) deallocate(v1d)
       if (allocated(v2d)) deallocate(v2d)

       allocate(v1u(totrglus), stat=ier)
       if(ier/=0) call mct_die(subName,'allocate v1u',ier)
       allocate(v2u(totrglus), stat=ier)
       if(ier/=0) call mct_die(subName,'allocate v2u',ier)
       allocate(v1d(totrglus), stat=ier)
       if(ier/=0) call mct_die(subName,'allocate v1d',ier)
       allocate(v2d(totrglus), stat=ier)
       if(ier/=0) call mct_die(subName,'allocate v2d',ier)
       indxup(:)=(/(i,i=1,numlons*numlats)/)
       indxdn=indxup

       ! force the opposite sign to ensure the non-glu cells are not sorted into
       ! the glu space
       pot_veg_rev=transpose(pot_veg)
       glu_weights_rev=transpose(glu_weights(g,:,:))

       where (transpose(rglus) == 0)
          glu_weights_rev = -1.
          pot_veg_rev = 1.
       end where

       call D_mrgrnk(pot_veg_rev*glu_weights_rev,indxup,numlons*numlats)
       call D_mrgrnk(pot_veg_rev*glu_weights_rev*-1.,indxdn,numlons*numlats)
       !	  
       !  The sortxxxup and sortxxxdn arrays are only good 1:totrglus these arrays are also based on arrays lat first (360,720)
       !  sorting to match original matlab scripts.  Thats why I transposed the arrays above.  Also didn't use qsort because it
       !  isn't a stable sort.
       !  These calcs below reflect that the sort was done with transposed matrices
       sortlonsdn=(indxdn-1)/numlats+1
       sortlatsdn=mod(indxdn-1,numlats)+1
       sortlonsup=(indxup-1)/numlats+1
       sortlatsup=mod(indxup-1,numlats)+1
       ! these are correct now because non-glu cells are beyond the totrglus range
       v1u=sortlonsup(numlons*numlats-totrglus+1:numlons*numlats)
       v2u=sortlatsup(numlons*numlats-totrglus+1:numlons*numlats)
       v1d=(sortlonsdn(:totrglus))
       v2d=(sortlatsdn(:totrglus))

#ifdef DEBUG
       ! check that sorted cells are all in the land unit
       z = 0
       do ij=1,totrglus
          ! start with lon
          found = .false.
          do i1=1,count_gcells
            if (v1d(ij) == gloninds(i1)) then
               ! check for matching lat
               if(v2d(ij) == glatinds(i1)) then
                  z=z+1
                  found = .true.
               endif
            endif
          enddo
          if ( .not. found ) then
             write(iulog,*) subname, 'Error: Did not find sorted cell ', ij, ' lonind=', v1d(ij), ' latind=', v2d(ij)
          endif
       enddo
       if (z /= totrglus) then
          write(iulog,*) subname, 'Error: found ', z, ' sorted cells out of ', totrglus
       endif
#endif
       ! limit abandonment to existing glm crop and pasture
       if ((past_neg >= (past_farea + past_nfarea)).and.(past_neg>0.)) then
          unmet_neg_past(g) = past_neg - (past_farea + past_nfarea)
          past_neg = past_farea + past_nfarea
#ifdef DEBUG
             write(iulog,*) subname, 'not enough pasture: reducing pasture loss  by km^2', unmet_neg_past(g)
#endif
       end if
       if ((crop_neg >= (crop_farea + crop_nfarea)).and.(crop_neg>0.)) then
          unmet_neg_crop(g) = crop_neg - (crop_farea + crop_nfarea)
          crop_neg = crop_farea + crop_nfarea
#ifdef DEBUG
             write(iulog,*) subname, 'not enough crop: reducing crop loss by km^2', unmet_neg_crop(g)
#endif
       end if

       ! compute the potential glm forest that could be created from cropland
       ! abandonment, assuming crop2forest is abandoned first
       ! also calculate the potential glm nonforest that could be created
       ! if there is more crop abandonemnt than there is crop2forest
       if (crop_neg <= crop_farea) then
          pot_forest_from_crop = crop_neg
          pot_nonforest_from_crop = 0.0
       else
          pot_forest_from_crop = crop_farea
          pot_nonforest_from_crop = crop_neg - crop_farea
       end if

       ! compute the potential glm forest that could be created from pasture
       ! abandonment, assuming pasture2forest is abandoned first
       ! also calculate the potential glm nonforest that could be created
       ! if there is more pasture abandonemnt than there is pasture2forest
       if (past_neg <= past_farea) then
          pot_forest_from_past = past_neg
          pot_nonforest_from_past = 0.0
       else
          pot_forest_from_past = past_farea
          pot_nonforest_from_past = past_neg - past_farea
       end if

       ! compute the potential glm forest and nonforest that could be lost from cropland and
       ! pasture expansion
       ! include abandoned pasture and cropland, and limit to total available area
       ! the assumption here is that ag expansion occurs on nonforest land first
       ! note that pot_### values are already limited to glm limits

       ! cropland
       if ((crop_pos>0.).and.( crop_pos > (GLM_nfarea+pot_nonforest_from_past) )) then
          pot_forest_to_crop = crop_pos - (GLM_nfarea+pot_nonforest_from_past)
          pot_nonforest_to_crop = GLM_nfarea+pot_nonforest_from_past
          if (pot_forest_to_crop >= GLM_farea) then
#ifdef DEBUG
             write(iulog,*) subname, 'not enough land area: reducing crop expansion by km^2', pot_forest_to_crop - GLM_farea
#endif
             pot_forest_to_crop = GLM_farea
          end if
       else
          pot_forest_to_crop = 0.
          pot_nonforest_to_crop = crop_pos
       end if

       ! pasture
       if ((past_pos>0.).and.( past_pos > (GLM_nfarea - pot_nonforest_to_crop + pot_nonforest_from_crop) )) then
          pot_forest_to_past = past_pos - (GLM_nfarea - pot_nonforest_to_crop + pot_nonforest_from_crop)
          pot_nonforest_to_past = GLM_nfarea - pot_nonforest_to_crop + &
                                  pot_nonforest_from_crop
          if (pot_forest_to_past >= GLM_farea - pot_forest_to_crop) then
#ifdef DEBUG
             write(iulog,*) subname, 'not enough land area: reducing pasture expansion by km^2', pot_forest_to_past - (GLM_farea - pot_forest_to_crop)
#endif
             pot_forest_to_past = GLM_farea - pot_forest_to_crop
          end if
       else
          pot_forest_to_past = 0.
          pot_nonforest_to_past = past_pos
       end if

#ifdef DEBUG
          write(iulog,*) subname, 'all units are km^2'
          write(iulog,*) subname, 'current global crop = ' ,sum(glm_crop(:,:,np1)*cellarea),' pasture = ', sum(glm_past(:,:,np1)*cellarea)
          write(iulog,*) subname, 'current g crop = ' , sum(glm_crop(:,:,np1)*cellarea*glu_weights(g,:,:)), &
             ' pasture = ', sum(glm_past(:,:,np1)*cellarea*glu_weights(g,:,:))
          write(iulog,*) subname, 'crop_farea=', crop_farea, ' crop_nfarea=', &
                              crop_nfarea
          write(iulog,*) subname, 'past_farea=', past_farea, ' past_nfarea=', &
                              past_nfarea
          write(iulog,*) subname, 'GLM_farea=', GLM_farea, ' GLM_nfarea=', &
                              GLM_nfarea
          write(iulog,*) subname, 'crop_pos=', crop_pos, ' crop_neg=', crop_neg
          write(iulog,*) subname, 'past_pos=', past_pos, ' past_neg=', past_neg
          write(iulog,*) subname, 'potential transition areas'
          write(iulog,*) subname, 'pot_forest_to_crop=', pot_forest_to_crop
          write(iulog,*) subname, 'pot_nonforest_to_crop=', pot_nonforest_to_crop
          write(iulog,*) subname, 'pot_forest_from_crop=', pot_forest_from_crop
          write(iulog,*) subname, 'pot_nonforest_from_crop=', pot_nonforest_from_crop
          write(iulog,*) subname, 'pot_forest_to_past=', pot_forest_to_past
          write(iulog,*) subname, 'pot_nonforest_to_past=', pot_nonforest_to_past
          write(iulog,*) subname, 'pot_forest_from_past=', pot_forest_from_past
          write(iulog,*) subname, 'pot_nonforest_from_past=', &
                              pot_nonforest_from_past
#endif

       ! the problem here is that GCAM forest area changes are not necessarily
       !    tied to agricultural expansion/contraction
       ! however, in order to translate GCAM forest expansion we need to try to
       !    match positive farea_d
       ! GCAM forest contraction is not attempted to be matched here

       ! all cases start with the general transitions
       ! calculate generic crop and pasture transitions from
       ! forest/nonforest
       ! don't worry about matching forest loss for now
       ! note that this is same behavior as for meeting positive farea_d
       crop_pos_f = pot_forest_to_crop
       crop_pos_nf = pot_nonforest_to_crop
       crop_neg_f = pot_forest_from_crop
       ! no !note that crop_neg_nf is a fraction
       !crop_neg_nf = pot_nonforest_from_crop/(crop_nfarea+1e-12)
       crop_neg_nf = pot_nonforest_from_crop
       past_pos_f = pot_forest_to_past
       past_pos_nf = pot_nonforest_to_past
       past_neg_f = pot_forest_from_past
       ! note that past_neg_nf is a fraction
       !past_neg_nf = pot_nonforest_from_past/(past_nfarea+1e-12)
       past_neg_nf = pot_nonforest_from_past

       ! only if the GCAM forest expansion is not met is there an adjustment
       if (farea_d > 0. .and. (pot_forest_from_crop + pot_forest_from_past - &
           pot_forest_to_crop - pot_forest_to_past) < farea_d) then

          f_diff = farea_d - (pot_forest_from_crop + pot_forest_from_past -&
                              pot_forest_to_crop - pot_forest_to_past)

          ! create net transitions from ag forest to ag nonforest
          !    make sure all cases are covered for estimating the
          !    crop/pasture breakdown
          ! ensure that there is enough area available
          !    if not, store unmet area for creating transitions later to
          !    meet regional area
          ag_farea_new = crop_farea + past_farea -  &
             pot_forest_from_crop - pot_forest_from_past + &
             pot_forest_to_crop + pot_forest_to_past

          ! skip if no ag forest area to convert
          if (ag_farea_new > 0. ) then
             reassign_crop = f_diff * (crop_farea - pot_forest_from_crop + &
                pot_forest_to_crop) / ag_farea_new
             GLM_nfarea_new = GLM_nfarea - pot_nonforest_to_crop - &
                pot_nonforest_to_past + pot_nonforest_from_crop + &
                pot_nonforest_from_past

#ifdef DEBUG
             write(iulog,*) subname, 'matching farea_d'
             write(iulog,*) subname, 'farea_d=', farea_d, ' f_diff=', f_diff
             write(iulog,*) subname, 'reassign_crop1=', reassign_crop, &
                              ' GLM_nfarea_new1=', GLM_nfarea_new
#endif

             ! do crop transitions first
             if ((crop_farea - pot_forest_from_crop + pot_forest_to_crop) < reassign_crop .or. GLM_nfarea_new < reassign_crop) then
                reassign_crop = min((crop_farea - pot_forest_from_crop + pot_forest_to_crop), GLM_nfarea_new)
             end if
             ! adjust nf crop trans
             if(crop_pos_nf > 0.) then
                crop_pos_nf = crop_pos_nf + reassign_crop
             else if (crop_neg_nf > 0. .and. crop_neg_nf >= reassign_crop) then
                crop_neg_nf = crop_neg_nf - reassign_crop
             else if (crop_neg_nf > 0. .and. crop_neg_nf < reassign_crop) then
                crop_neg_nf = 0.
                crop_pos_nf = reassign_crop - crop_neg_nf
             else
                ! both crop_neg_nf and crop_pos_nf are zero; they can never both be positive
                crop_pos_nf = crop_pos_nf + reassign_crop
             endif
             ! adjust f crop trans
             if(crop_neg_f > 0.) then
                crop_neg_f = crop_neg_f + reassign_crop
             else if (crop_pos_f > 0. .and. crop_pos_f >= reassign_crop) then
                crop_pos_f = crop_pos_f - reassign_crop
             else if (crop_pos_f > 0. .and. crop_pos_f < reassign_crop) then
                crop_pos_f = 0.
                crop_neg_f = reassign_crop - crop_pos_f
             else
                ! both crop_neg_f and crop_pos_f are zero; they can never both be positive
                crop_neg_f = crop_neg_f + reassign_crop
             endif
             ! determine the pasture reassignment
             GLM_nfarea_new = GLM_nfarea_new - reassign_crop
             reassign_past = f_diff - reassign_crop

#ifdef DEBUG
             write(iulog,*) subname, 'reassign_crop2=', reassign_crop
             write(iulog,*) subname, 'reassign_past1=', reassign_past, &
                              ' GLM_nfarea_new2=', GLM_nfarea_new
#endif

             ! do pasture transitions second
             if ((past_farea - pot_forest_from_past + pot_forest_to_past) < reassign_past .or. &
                 GLM_nfarea_new < reassign_past) then
                reassign_past = min((past_farea - pot_forest_from_past + pot_forest_to_past), GLM_nfarea_new)
             end if
             ! adjust nf pasture trans
             if(past_pos_nf > 0.) then
                past_pos_nf = past_pos_nf + reassign_past
             else if (past_neg_nf > 0. .and. past_neg_nf >= reassign_past) then
                past_neg_nf = past_neg_nf - reassign_past
             else if (past_neg_nf > 0. .and. past_neg_nf < reassign_past) then
                past_neg_nf = 0.
                past_pos_nf = reassign_past - past_neg_nf
             else
                ! both past_neg_nf and past_pos_nf are zero; they can never both be positive
                past_pos_nf = past_pos_nf + reassign_past
             endif
             ! adjust f pasture trans
             if(past_neg_f > 0.) then
                past_neg_f = past_neg_f + reassign_past
             else if (past_pos_f > 0. .and. past_pos_f >= reassign_past) then
                past_pos_f = past_pos_f - reassign_past
             else if (past_pos_f > 0. .and. past_pos_f < reassign_past) then
                past_pos_f = 0.
                past_neg_f = reassign_past - past_pos_f
             else
                ! both past_neg_f and past_pos_f are zero; they can never both be positive
                past_neg_f = past_neg_f + reassign_past
             endif

          else
             ! no ag area to reassign
             reassign_crop = 0.
             reassign_past = 0.
          end if
 
          ! now set unmet_farea
          unmet_farea(g) = f_diff - (reassign_crop + reassign_past)

#ifdef DEBUG
          write(iulog,*) subname, 'reassign_past2=', reassign_past, &
                              ' unmet_farea(',g,')=', unmet_farea(g)
#endif

       end if

#ifdef DEBUG
          write(iulog,*) subname, 'crop_pos_f=', crop_pos_f, ' crop_pos_nf=',crop_pos_nf
          write(iulog,*) subname, 'crop_neg_f=', crop_neg_f, ' crop_neg_nf =',crop_neg_nf
          write(iulog,*) subname, 'past_pos_f=', past_pos_f, ' past_pos_nf=',past_pos_nf
          write(iulog,*) subname, 'past_neg_f=', past_neg_f, ' past_neg_nf =',past_neg_nf
#endif

       ! apply the cropland and pature changes (both forested and non-forested)
       ! to individual gridcells in each region/AEZ

       ! the approach is to expand onto nonforest first,
       ! and to abandon on forest land first

       ! np1 values initially = n values, and so always use np1 for calcs
       ! otherwise cells with multiple land units will not be updated correctly

       !!!!!!!!!!!!!!!
       ! the problem is that the glu_weights for forest, nonforest, crop, and
       ! pasture change with each transition, while the their land weights
       ! remain the same
       ! so to do this properly the type-specific weights need to be updated and
       ! tracked, and also stored and restored via restart
       ! for now, just use rglus mask to determine whether the appropriate cells
       ! have changed properly

       start_crop = sum(glm_crop(:,:,np1)*cellarea)
       start_crop_g = sum(glm_crop(:,:,np1)*cellarea*rglus)
       start_past = sum(glm_past(:,:,np1)*cellarea)
       start_past_g = sum(glm_past(:,:,np1)*cellarea*rglus)

       ! crop_neg_nf and crop_neg_f !!!!!!!!!!!!!!!!!!!!!!!
       curr_crop = sum(glm_crop(:,:,np1)*cellarea)
       curr_crop_g = sum(glm_crop(:,:,np1)*cellarea*rglus)
#ifdef DEBUG
       write(iulog,*) subname, 'before crop loss, global crop ', curr_crop
       write(iulog,*) subname, 'before crop loss, g crop ', curr_crop_g
       !avail_landA = 0.
       !where ( glu_weights(g,:,:) > 0. )
       !   avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       !end where
       !sumavail_landA=sum(avail_landA)
       !write(iulog,*) subname, 'sumavail_landA  = ', sumavail_landA
       write(iulog,*) subname, '*crop to nonforest ', sum((glm_crop(:,:,np1)*cellarea_nonforest)*crop_neg_nf/(crop_nfarea+1e-12))
#endif
       if (crop_nfarea > 0.) then
#ifdef DEBUG
          if (abs(sum((glm_crop(:,:,np1)*cellarea_nonforest)*crop_neg_nf/(crop_nfarea)) - crop_neg_nf) > diff_tol) then
             write(iulog,*) subname, 'Warning: crop to nonforest != crop_neg_nf ', crop_neg_nf 
          end if
#endif
          where (glu_weights(g,:,:) > 0. .and. glm_crop(:,:,np1) > 0.) 
             glm_crop(:,:,np1) = glm_crop(:,:,np1)   -  ((glm_crop(:,:,np1)*cellarea_nonforest)*crop_neg_nf/crop_nfarea)/cellarea
          endwhere
       endif

#ifdef DEBUG
       !avail_landA = 0.
       !where ( glu_weights(g,:,:) > 0. )
       !   avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       !end where
       !sumavail_landA=sum(avail_landA)
       !write(iulog,*) subname, 'sumavail_landA  = ', sumavail_landA
       write(iulog,*) subname, 'global crop after crop2nf ', sum(glm_crop(:,:,np1)*cellarea)
       write(iulog,*) subname, 'g crop after crop2nf ',sum(glm_crop(:,:,np1)*cellarea*rglus)
#endif
       if (crop_neg_f>0.) then
          temp_crop = 0. 
          cumsum_sorted_farea=0.
          call cumsum(glm_crop(:,:,np1)*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)

          ! for some reason minloc returns 1 now if the array is all zeros or the mask is all false
          ! so don't use minloc; it isn't faster than a loop
          sortind(1) = 0
          do i1=1,totrglus
             if(cumsum_sorted_farea(i1)>crop_neg_f) then
                sortind(1) = i1
                exit
             endif
          enddo
#ifdef DEBUG
          write(iulog,*) subname, 'initial sortind(1) = ', sortind(1)
          write(iulog,*) subname, 'total cumsum crop_fa ', cumsum_sorted_farea(totrglus)
#endif
          if (sortind(1)==0) then
             if ( abs(crop_neg_f-cumsum_sorted_farea(totrglus)) .lt. 1e-10) then 
                crop_neg_f=cumsum_sorted_farea(totrglus)
#ifdef DEBUG
                write(iulog,*) subname, 'setting crop_neg_f to total cumsum ', cumsum_sorted_farea(totrglus)
#endif
                sortind(1) = 0
                do i1=1,totrglus
                   if(cumsum_sorted_farea(i1) >= crop_neg_f) then
                      sortind(1) = i1
                      exit
                   endif
                enddo
                if(sortind(1)==0) then
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum pasture farea ', cumsum_sorted_farea(totrglus), ' to meet past_neg_f ', past_neg_f
#endif
                endif
             else
                sortind(1) = totrglus
#ifdef DEBUG
                write(iulog,*) subname, 'Warning: not enough cum crop farea ', cumsum_sorted_farea(totrglus), ' to meet crop_neg_f ', crop_neg_f
#endif
             end if
          end if
#ifdef DEBUG
          write(iulog,*) subname, 'final sortind(1) = ', sortind(1)
#endif
          if (sortind(1)>1) then
             z=0
             sortsitesdn=0
             do ij=1,sortind(1)-1
                found = .false.
                ! get lon first
                do i1=1,count_gcells
                   if (v1d(ij) == gloninds(i1)) then
                      ! check for matching lat
                      if(v2d(ij) == glatinds(i1)) then
                         z=z+1
                         found = .true.
                      endif
                   endif
                enddo
#ifdef DEBUG
                if ( .not. found ) then
                   write(iulog,*) subname, 'Error crop loss: Did not find sorted cell ', ij, ' lonind=', v1d(ij), ' latind=', v2d(ij)
                endif
#endif
                sortsitesdn(v1d(ij),v2d(ij))=1
             end do
#ifdef DEBUG
             if (z /= sortind(1)-1) then
                write(iulog,*) subname, 'Error: found ', z, ' sorted cells out of ', sortind(1)-1
             endif
#endif
             temp_crop = sum(glm_crop(:,:,np1) * fnfforest * cellarea * sortsitesdn)
             where(sortsitesdn>0)
                glm_crop(:,:,np1)= glm_crop(:,:,np1) - glm_crop(:,:,np1) * fnfforest
             end where

             final_area_needed =  crop_neg_f - cumsum_sorted_farea(sortind(1)-1)
          else
             final_area_needed =  crop_neg_f
          end if
#ifdef DEBUG
          write(iulog,*) subname, 'temp_crop ', temp_crop
          write(iulog,*) subname, 'final_area_needed ', final_area_needed
          write(iulog,*) subname, 'last cell crop_fa ', &
             glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))*cellarea(v1d(sortind(1)),v2d(sortind(1)))
          write(iulog,*) subname, '*crop to forest ', temp_crop + &
          min(final_area_needed,&
                   glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))* cellarea(v1d(sortind(1)),v2d(sortind(1))))
          if( abs(temp_crop + &
              min(final_area_needed,&
                   glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))* cellarea(v1d(sortind(1)),v2d(sortind(1)))) - &
              crop_neg_f) > diff_tol) then
              write(iulog,*) subname, 'Error: crop to forest != crop_neg_f ', crop_neg_f
          end if
#endif
          glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) = &
               glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) - &
               min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
               glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
       end if

#ifdef DEBUG
       write(iulog,*) subname, 'after crop loss, global crop ', sum(glm_crop(:,:,np1)*cellarea)
       write(iulog,*) subname, 'after crop loss, g crop ', sum(glm_crop(:,:,np1)*cellarea*rglus)
       write(iulog,*) subname, '*crop change (global) ', sum(glm_crop(:,:,np1)*cellarea) - curr_crop
       write(iulog,*) subname, '*crop change (g) ', sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g
       if ( abs( (sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g) -  -(crop_neg_nf + crop_neg_f)) >  diff_tol ) then
          write(iulog,*) subname, 'Error: crop change (g) != -(crop_neg_f+crop_neg_nf) ', -(crop_neg_f+crop_neg_nf)
       endif
       write(iulog,*) subname, 'before pasture loss,  global pasture ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'before pasture loss, g pasture ', sum(glm_past(:,:,np1)*cellarea*rglus)
#endif


       ! past_neg_f and past_neg_nf !!!!!!!!!!!!!!!!!!!!!!!!!!!
       curr_past = sum(glm_past(:,:,np1)*cellarea)
       curr_past_g = sum(glm_past(:,:,np1)*cellarea*rglus)
#ifdef DEBUG
!       avail_landA = 0.
!       where ( glu_weights(g,:,:) > 0. )
!          avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
!       end where
!       sumavail_landA=sum(avail_landA)
       !write(iulog,*) subname, 'sumavail_landA  = ', sumavail_landA
       write(iulog,*) subname, '*pasture to nonforest ', sum((glm_past(:,:,np1)*cellarea_nonforest)*past_neg_nf/(past_nfarea+1e-12))
 !      temp_past_grid = 0.
 !      where (glu_weights(g,:,:) > 0. .and. glm_past(:,:,np1) > 0.)
 !         temp_past_grid = (glm_past(:,:,np1)*cellarea_nonforest)*past_neg_nf/(past_nfarea+1e-12)
 !      endwhere 
 !      write(iulog,*) subname, '*pasture to nonforest filtered ',sum(temp_past_grid)
#endif
       if (past_nfarea > 0.) then
#ifdef DEBUG
          if ( abs( sum((glm_past(:,:,np1)*cellarea_nonforest)*past_neg_nf/(past_nfarea)) - past_neg_nf) > diff_tol ) then
             write(iulog,*) subname, 'Warning: pasture to nonforest != past_neg_nf ', past_neg_nf
          endif
#endif
          where(glu_weights(g,:,:) > 0. .and. glm_past(:,:,np1) > 0.) 
             glm_past(:,:,np1) = glm_past(:,:,np1) - ((glm_past(:,:,np1)*cellarea_nonforest)*past_neg_nf/past_nfarea)/cellarea
          end where

       endif
#ifdef DEBUG
!       avail_landA = 0.
!       where ( glu_weights(g,:,:) > 0. )
!          avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
!       end where
!       sumavail_landA=sum(avail_landA)
       !write(iulog,*) subname, 'sumavail_landA  = ', sumavail_landA
       write(iulog,*) subname, 'global pasture after past2nf ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'g pasture after past2nf  ',sum(glm_past(:,:,np1)*cellarea*rglus)
!       avail_landA = 0.
!       where ( glu_weights(g,:,:) > 0. )
!          avail_landA=(pctland_in2015*cellarea_nonforest - glm_crop(:,:,np1)*cellarea_nonforest-glm_past(:,:,np1)*cellarea_nonforest)
!       end where
!       sumavail_landA=sum(avail_landA)
       !write(iulog,*) subname, 'alternative sumavail_landA  = ', sumavail_landA
#endif

       if (past_neg_f>0) then
          temp_past = 0.
          cumsum_sorted_farea=0.
          call cumsum(glm_past(:,:,np1)*cellarea_forest(:,:),v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)
          sortind(1) = 0
          do i1=1,totrglus
             if(cumsum_sorted_farea(i1) > past_neg_f) then
                sortind(1) = i1
                exit
             endif
          enddo
#ifdef DEBUG
          write(iulog,*) subname, 'initial sortind(1) = ', sortind(1)
          write(iulog,*) subname, 'total cumsum past_fa ', cumsum_sorted_farea(totrglus)
#endif
          if (sortind(1)==0) then
             if ( abs(past_neg_f-cumsum_sorted_farea(totrglus)) .lt. 1e-10) then 
                past_neg_f=cumsum_sorted_farea(totrglus)
#ifdef DEBUG
                write(iulog,*) subname, 'setting past_neg_f to total cumsum ', cumsum_sorted_farea(totrglus)
#endif
                sortind(1) = 0
                do i1=1,totrglus
                   if(cumsum_sorted_farea(i1) >= past_neg_f) then
                      sortind(1) = i1
                      exit
                   endif
                enddo
                if(sortind(1)==0) then
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum pasture farea ', cumsum_sorted_farea(totrglus), ' to meet past_neg_f ', past_neg_f
#endif
                endif
             else
                sortind(1) = totrglus
#ifdef DEBUG
                write(iulog,*) subname, 'Warning: not enough cum pasture farea ', cumsum_sorted_farea(totrglus), ' to meet past_neg_f ', past_neg_f
#endif
             end if
          end if
#ifdef DEBUG
          write(iulog,*) subname, 'final sortind(1) = ', sortind(1)
#endif
          if (sortind(1)>1) then
             sortsitesdn=0
             z=0
             do i=1,sortind(1)-1
                z=z+1
                sortsitesdn(v1d(i),v2d(i))=1
             end do
             temp_past = sum(glm_past(:,:,np1) * fnfforest * cellarea * sortsitesdn)
#ifdef DEBUG
             write(iulog,*) subname, 'z = ', z, 'sorting(1) = ', sortind(1), ' totrglus = ', totrglus
             write(iulog,*) subname, 'temp_past = ', temp_past
#endif
             where(sortsitesdn>0)
                glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,np1) * fnfforest
             end where
             final_area_needed =  past_neg_f - cumsum_sorted_farea(sortind(1)-1)
          else
             final_area_needed =  past_neg_f
          end if
#ifdef DEBUG
          write(iulog,*) subname, 'temp_past ', temp_past
          write(iulog,*) subname, 'final_area_needed ', final_area_needed
          write(iulog,*) subname, 'last cell past_fa ', &
             glm_past(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))*cellarea(v1d(sortind(1)),v2d(sortind(1)))
          write(iulog,*) subname, '*pasture to forest ', temp_past + &
          min(final_area_needed,&
                   glm_past(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))*cellarea(v1d(sortind(1)),v2d(sortind(1))))

          if (  abs( temp_past + &
                min(final_area_needed,&
                   glm_past(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1)))*cellarea(v1d(sortind(1)),v2d(sortind(1)))) - &
                past_neg_f) > diff_tol ) then
             write(iulog,*) subname, 'Error: pasture to forest != past_neg_f ', past_neg_f
          endif

#endif
          glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) = &
               glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) - &
               min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
               glm_past(v1d(sortind(1)),v2d(sortind(1)),np1)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
       end if

#ifdef DEBUG
          write(iulog,*) subname, 'after pasture loss,  global pasture ', sum(glm_past(:,:,np1)*cellarea)
          write(iulog,*) subname, 'after pasture loss,  g pasture ',sum(glm_past(:,:,np1)*cellarea*rglus)
          write(iulog,*) subname, '*pasture change (global)', sum(glm_past(:,:,np1)*cellarea) - curr_past
          write(iulog,*) subname, '*pasture change (g) ', sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g
          if ( abs( (sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g) - -(past_neg_nf + past_neg_f) ) > diff_tol ) then
             write(iulog,*) subname, 'Error: pasture change (g) != -(past_neg_nf + past_neg_f) ', -(past_neg_nf + past_neg_f)
          endif
          write(iulog,*) subname, 'before crop nf gain,  global crop ', sum(glm_crop(:,:,np1)*cellarea)
          write(iulog,*) subname, 'before crop nf gain,  g crop ', sum(glm_crop(:,:,np1)*cellarea*rglus)
#endif


       ! crop_pos_nf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       curr_crop = sum(glm_crop(:,:,np1)*cellarea)
       curr_crop_g = sum(glm_crop(:,:,np1)*cellarea*rglus)
       avail_land0 = 0.
       where ( glu_weights(g,:,:)>0. .and. glm_crop(:,:,np1)>0. )
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       endwhere
       sumavail_land0=sum(avail_land0)

       if ( sumavail_land0 >= crop_pos_nf .or. abs(crop_pos_nf)<=1e-6) then
          if (abs(crop_pos_nf)<=1e-6) then
             crop_pos_nf=0.
          end if

         if (sumavail_land0 > 0.) then

#ifdef DEBUG
            write(iulog,*) subname, 'cropland increase on non-forested land - land available'
            !write(iulog,*) subname, 'sumavail_land0  = ', sumavail_land0
            write(iulog,*) subname, '*nonforest to crop 0 ', sum((avail_land0/(sumavail_land0)*crop_pos_nf))
            if ( abs( sum((avail_land0/(sumavail_land0)*crop_pos_nf)) - crop_pos_nf) > diff_tol ) then
               write(iulog,*) subname, 'Warning: nonforest to crop 0 != crop_post_nf ', crop_pos_nf
            endif
#endif
            where(glu_weights(g,:,:) > 0. .and. glm_crop(:,:,np1)>0. ) 
               glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_land0/(sumavail_land0)*crop_pos_nf)/cellarea
            end where

          endif

#ifdef DEBUG
          !avail_land0 = 0.
          !where ( glu_weights(g,:,:) > 0. .and. glm_crop(:,:,np1)>0.)
          !   avail_land0=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          !end where
          !sumavail_land0=sum(avail_land0)
          !write(iulog,*) subname, 'sumavail_land0  = ', sumavail_land0
#endif

       else
          avail_landA = 0.
          where ( glu_weights(g,:,:) > 0. )
             avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          end where
          sumavail_landA=sum(avail_landA)
          if (sumavail_landA >= crop_pos_nf .and. (sumavail_landA-sumavail_land0) > 0.) then
#ifdef DEBUG
             write(iulog,*) subname, 'cropland increase on non-forest - land available'
             !write(iulog,*) subname, 'sumavail_landA  = ', sumavail_landA
             write(iulog,*) subname, '*nonforest to crop A ', sum(avail_land0) + &
                sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0)*(crop_pos_nf-sumavail_land0)))
             if ( abs( sum(avail_land0) + &
                sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0)*(crop_pos_nf-sumavail_land0))) - crop_pos_nf) > diff_tol ) then
                write(iulog,*) subname, 'Warning: nonforest to crop A != crop_pos_nf ', crop_pos_nf
             endif
#endif 
             where(glu_weights(g,:,:) > 0.) 
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0/cellarea
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                     ((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0)*(crop_pos_nf-sumavail_land0))/cellarea
             end where
          else

#ifdef DEBUG
             write(iulog,*) subname, 'Warning: crop increase on non-forest - land not available'
             write(iulog,*) subname, 'reducing crop increase to sumavail_landA'
             write(iulog,*) subname, 'crop_pos_nf=', crop_pos_nf, &
                                     ' sumavail_landA=', sumavail_landA, &
                                     ' sumavail_land0=', sumavail_land0
             write(iulog,*) subname, '*nonforest to crop ', sum(avail_land0) + &
                sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0+1e-12)*(crop_pos_nf-sumavail_land0)))
# endif

             crop_pos_nf = sumavail_landA
             where(glu_weights(g,:,:) > 0.)
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                   avail_land0/cellarea
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                     ((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0+1e-12)*(crop_pos_nf-sumavail_land0))/cellarea
             end where

          end if
       end if
#ifdef DEBUG
       write(iulog,*) subname, 'after crop nf gain, global crop ', sum(glm_crop(:,:,np1)*cellarea)
       write(iulog,*) subname, 'after crop nf gain, g crop ', sum(glm_crop(:,:,np1)*cellarea*rglus)
       write(iulog,*) subname, '*crop change (global) ', sum(glm_crop(:,:,np1)*cellarea) - curr_crop
       write(iulog,*) subname, '*crop change (g) ', sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g
       if ( abs(sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g - crop_pos_nf) > diff_tol  ) then
          write(iulog,*) subname, 'Warning: crop change (g) != crop_pos_nf ', crop_pos_nf
       endif
       write(iulog,*) subname, 'before pasture nf gain, global pasture ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'before pasture nf gain, g pasture ', sum(glm_past(:,:,np1)*cellarea*rglus)
#endif


       ! past_pos_nf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       curr_past = sum(glm_past(:,:,np1)*cellarea)
       curr_past_g = sum(glm_past(:,:,np1)*cellarea*rglus)
       avail_land0 = 0.
       where ( glu_weights(g,:,:) > 0. .and.glm_past(:,:,np1)>0.)
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       end where
       sumavail_land0=sum(avail_land0)
       if ( sumavail_land0 >= past_pos_nf .or. abs(past_pos_nf)<=1e-6) then
          if (abs(past_pos_nf)<=1e-6) then
             past_pos_nf=0.
          end if

          if (sumavail_land0 > 0.) then
#ifdef DEBUG
             write(iulog,*) subname, 'pasture increase on non-forested land - land available'
             write(iulog,*) subname, '*nonforest to pasture 0 ', sum((avail_land0/(sumavail_land0)*past_pos_nf))
             if ( abs( sum((avail_land0/(sumavail_land0)*past_pos_nf)) - past_pos_nf) > diff_tol  ) then
                write(iulog,*) subname, 'Warning: nonforest to pasture 0 != past_pos_nf ', past_pos_nf
             endif
#endif
             where(glu_weights(g,:,:) > 0. .and. glm_past(:,:,np1)>0.) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0/(sumavail_land0)*past_pos_nf)/cellarea
             end where
          endif
       else
          avail_landA = 0.
          where ( glu_weights(g,:,:) > 0. )
             avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          end where
          sumavail_landA=sum(avail_landA)
          if ( sumavail_landA >= past_pos_nf .and. sumavail_landA-sumavail_land0 > 0. ) then
#ifdef DEBUG
             write(iulog,*) subname, 'pasture increase on non-forested land - land available'
             write(iulog,*) subname, '*nonforest to pasture A ', sum(avail_land0) + &
                sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0)*(past_pos_nf-sumavail_land0)))
             if ( abs( sum(avail_land0) + &
                  sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0)*(past_pos_nf-sumavail_land0))) - past_pos_nf) > diff_tol  ) then
                write(iulog,*) subname, 'Warning: nonforest to pasture A != past_pos_nf ', past_pos_nf
             end if
#endif
             where(glu_weights(g,:,:) > 0.) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0/cellarea
                glm_past(:,:,np1) = glm_past(:,:,np1) + &
                     ((avail_landA-avail_land0)/(sumavail_landA -sumavail_land0) * &
                     (past_pos_nf-sumavail_land0))/cellarea
             end where
          else
#ifdef DEBUG
             write(iulog,*) subname, 'Warning: pasture increase on non-forest - land NOT available, reducing pasture increase to accomodate'
             write(iulog,*) subname, 'past_pos_nf=', past_pos_nf, &
                                     'sumavail_landA=', sumavail_landA, &
                                     'sumavail_land0=', sumavail_land0
             write(iulog,*) subname, '*nonforest to pasture ', sum(avail_land0) + &
                sum(((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0+1e-12)*(past_pos_nf-sumavail_land0)))
#endif
             past_pos_nf = sumavail_landA
             where(glu_weights(g,:,:) > 0.) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0/cellarea)
                glm_past(:,:,np1) = glm_past(:,:,np1) + ((avail_landA - avail_land0) / &
                     (sumavail_landA - sumavail_land0 + 1e-12)*(past_pos_nf-sumavail_land0))/cellarea
             end where
          end if
       end if

#ifdef DEBUG
       write(iulog,*) subname, 'after pasture nf gain, global pasture ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'after pasture nf gain, g pasture ', sum(glm_past(:,:,np1)*cellarea*rglus)
       write(iulog,*) subname, '*pasture change (global)', sum(glm_past(:,:,np1)*cellarea) - curr_past
       write(iulog,*) subname, '*pasture change (g) ', sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g
       if ( abs(sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g - past_pos_nf) > diff_tol ) then
          write(iulog,*) subname, 'Warning: pasture change (g) != past_pos_nf ', past_pos_nf
       endif
       write(iulog,*) subname, 'before crop f gain, global crop ', sum(glm_crop(:,:,np1)*cellarea)
       write(iulog,*) subname, 'before crop f gain, g crop ', sum(glm_crop(:,:,np1)*cellarea*rglus)
#endif


       ! crop_pos_f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       curr_crop = sum(glm_crop(:,:,np1)*cellarea)
       curr_crop_g = sum(glm_crop(:,:,np1)*cellarea*rglus)
       if (abs(crop_pos_f)>1e-6) then
          if (sumavail_land0>=crop_pos_f) then
#ifdef DEBUG
             write(iulog,*) subname, 'cropland increase on forested land - land available'
#endif
             cumsum_sorted_farea=0.
             call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
             sortind(1) = 0
             do i1=1,totrglus
                if(cumsum_sorted_farea(i1) >= crop_pos_f) then
                   sortind(1) = i1
                   exit
                endif
             enddo
             if (sortind(1)==0) then
                ! may need another check here for precision tolerance?
                sortind(1) = totrglus
                write(iulog,*) subname, 'Warning: not enough cum avail forest 0 ',cumsum_sorted_farea(totrglus) , ' to meet crop_pos_f ', crop_pos_f
             end if
             if (sortind(1)>1) then
                sortsitesup=0
                do i=1,sortind(1)-1
                   sortsitesup(v1u(i),v2u(i))=1
                end do
                where(sortsitesup>0 .and. glm_crop(:,:,np1)>0.)
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0 / cellarea
                end where
                final_area_needed =  crop_pos_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  crop_pos_f
             end if
             glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                  glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                  min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                  avail_land0(v1u(sortind(1)),v2u(sortind(1))) / &
                  cellarea(v1u(sortind(1)),v2u(sortind(1))))
          else
             avail_landA = 0.
             where ( glu_weights(g,:,:) > 0 )
                avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
             end where
             sumavail_landA=sum(avail_landA)
             if (sumavail_landA >=crop_pos_f) then
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0 .and. glm_crop(:,:,np1)>0.)
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0/cellarea
                end where
#ifdef DEBUG
                write(iulog,*) subname, 'cropland increase on forested land - land available'
#endif                   
                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                sortind(1) = 0
                do i1=1,totrglus
                   if( cumsum_sorted_farea(i1) >= (crop_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)) ) then
                      sortind(1) = i1
                      exit
                   endif
                enddo
                if (sortind(1)==0) then
                   ! may need another precision check here?
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum avail forest A ',cumsum_sorted_farea(totrglus) ,&
                      ' to meet remaining crop_pos_f ', (crop_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0))
#endif
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_landA - avail_land0) / cellarea
                   end where
                   final_area_needed =  crop_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  crop_pos_f - sumavail_land0
                end if
                glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) +&
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                     (avail_landA(v1u(sortind(1)),v2u(sortind(1)))) / &
                     cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
#ifdef DEBUG
                write(iulog,*) subname, 'Warning: crop increase on forest - land not available'
                write(iulog,*) subname, 'reducing crop increase to sumavail_landA'
                write(iulog,*) subname, 'crop_pos_f=', crop_pos_f, &
                           'sumavail_landA=', sumavail_landA, &
                           'sumavail_land0=', sumavail_land0
#endif
                crop_pos_f = sumavail_landA
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,&
                   cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0 .and. glm_crop(:,:,np1)>0.)
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) +&
                      avail_land0/cellarea
                end where

                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,&
                   cumsum_sorted_farea(:totrglus),totrglus)
                sortind(1) = 0
                do i1=1,totrglus
                   if( cumsum_sorted_farea(i1) >= (crop_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)) ) then
                      sortind(1) = i1
                      exit
                   endif
                enddo

                if (sortind(1)==0) then
                   ! may need another precision check here?
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum avail forest A ',cumsum_sorted_farea(totrglus) ,& 
                      ' to meet remaining crop_pos_f ', (crop_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0))
#endif
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                         (avail_landA - avail_land0) / cellarea
                   end where
                   final_area_needed =  crop_pos_f - &
                      sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  crop_pos_f - sumavail_land0
                end if
                glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                   glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) +&
                   min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))),&
                   (avail_landA(v1u(sortind(1)),v2u(sortind(1)))) / &
                   cellarea(v1u(sortind(1)),v2u(sortind(1))))   

             end if
          end if

       end if
#ifdef DEBUG
       write(iulog,*) subname, 'after crop f gain, global crop ', sum(glm_crop(:,:,np1)*cellarea)
       write(iulog,*) subname, 'after crop f gain, g crop ', sum(glm_crop(:,:,np1)*cellarea*rglus)
       write(iulog,*) subname, '*crop change (global) ', sum(glm_crop(:,:,np1)*cellarea) - curr_crop
       write(iulog,*) subname, '*crop change (g) ', sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g
       if ( abs(sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g - crop_pos_f) > diff_tol ) then
          write(iulog,*) subname, 'Warning: crop change (g) != crop_pos_f ', crop_pos_f
       endif
       write(iulog,*) subname, 'before pasture f gain, global pasture ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'before pasture f gain, g pasture ', sum(glm_past(:,:,np1)*cellarea*rglus)
#endif


       ! past_pos_f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       curr_past = sum(glm_past(:,:,np1)*cellarea)
       curr_past_g = sum(glm_past(:,:,np1)*cellarea*rglus)
       avail_land0 = 0.
       where ( glu_weights(g,:,:) > 0 .and.glm_past(:,:,np1)>0.)
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
       end where
       sumavail_land0=sum(avail_land0)
       if (abs(past_pos_f)>1e-6) then 
          if (sumavail_land0>=past_pos_f) then
#ifdef DEBUG
             write(iulog,*) subname, 'pasture increase on forest - land available'
#endif                
             cumsum_sorted_farea=0.
             call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
             sortind(1) = 0
             do i1=1,totrglus
                if (cumsum_sorted_farea(i1) >= past_pos_f) then
                   sortind(1) = i1
                   exit
                endif
             enddo

             if (sortind(1)==0) then
                ! may need another precision check here?
                sortind(1) = totrglus
#ifdef DEBUG
                write(iulog,*) subname, 'Warning: not enough cum avail forest 0 ',cumsum_sorted_farea(totrglus) , ' to meet past_pos_f ', past_pos_f
#endif
             end if
             if (sortind(1)>1) then 
                sortsitesup=0
                do i=1,sortind(1)-1
                   sortsitesup(v1u(i),v2u(i))=1
                end do
                where(sortsitesup>0 .and.glm_past(:,:,np1)>0.)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 / cellarea
                end where
                final_area_needed =  past_pos_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  past_pos_f
             end if
             glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                  glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                  min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))) , &
                  avail_land0(v1u(sortind(1)),v2u(sortind(1))) / &
                  cellarea(v1u(sortind(1)),v2u(sortind(1))))
          else
             avail_landA = 0.
             where ( glu_weights(g,:,:) > 0 )
                avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
             end where
             sumavail_landA=sum(avail_landA)
             if (sumavail_landA >= past_pos_f) then 
#ifdef DEBUG
                write(iulog,*) subname, 'pasture increase on forest - land available'
#endif
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0 .and.glm_past(:,:,np1)>0.)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 / cellarea
                end where

                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                sortind(1) = 0
                do i1=1,totrglus
                   if ( cumsum_sorted_farea(i1) >= (past_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)) ) then
                      sortind(1) = i1
                      exit
                   endif
                enddo

                if (sortind(1)==0) then
                   ! may need another precision check here?
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum avail forest A ',cumsum_sorted_farea(totrglus) ,&
                      ' to meet remaining past_pos_f ', (past_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0))
#endif
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1)+(avail_landA-avail_land0)/cellarea
                   end where
                   final_area_needed =  past_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  past_pos_f - sumavail_land0
                end if
                glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                     (avail_landA(v1u(sortind(1)),v2u(sortind(1))) - &
                     avail_landA(v1u(sortind(1)),v2u(sortind(1)))) / &
                     cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
#ifdef DEBUG
                write(iulog,*) subname, 'Warning: pasture increase on forest - land not available'
                write(iulog,*) subname, 'reducing pasture increase to sumavail_landA'
                write(iulog,*) subname, 'past_pos_f=', past_pos_f, &
                        'sumavail_landA=', sumavail_landA, &
                        'sumavail_land0=', sumavail_land0
#endif
                past_pos_f = sumavail_landA
                cumsum_sorted_farea=0.

                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0 .and.glm_past(:,:,np1)>0.)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 / cellarea
                end where

                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                sortind(1) = 0
                do i1=1,totrglus
                   if ( cumsum_sorted_farea(i1) >= (past_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)) ) then
                      sortind(1) = i1
                      exit
                   endif
                enddo

                if (sortind(1)==0) then
                   ! may need another precision check here?
                   sortind(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum avail forest A ',cumsum_sorted_farea(totrglus) ,&
                      ' to meet remaining past_pos_f ', (past_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0))
#endif
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1)+(avail_landA-avail_land0)/cellarea
                   end where
                   final_area_needed =  past_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  past_pos_f - sumavail_land0
                end if
                glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                   glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                   min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))),&
                   (avail_landA(v1u(sortind(1)),v2u(sortind(1))) - &
                   avail_landA(v1u(sortind(1)),v2u(sortind(1)))) / &
                   cellarea(v1u(sortind(1)),v2u(sortind(1))))

             end if
          end if

       end if

#ifdef DEBUG
       write(iulog,*) subname, 'units = km^2'
       write(iulog,*) subname, 'after pasture f gain, global pasture ', sum(glm_past(:,:,np1)*cellarea)
       write(iulog,*) subname, 'after pasture f gain, g pasture ', sum(glm_past(:,:,np1)*cellarea*rglus)
       write(iulog,*) subname, '*pasture change (global)', sum(glm_past(:,:,np1)*cellarea) - curr_past
       write(iulog,*) subname, '*pasture change (g) ', sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g
       if ( abs(sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g - past_pos_f) > diff_tol ) then
          write(iulog,*) subname, 'Error: pasture change (g) != past_pos_f ', past_pos_f
       endif
       write(iulog,*) subname, '***** at end of initial land unit g loop, g cell values: ', g
       write(iulog,*) subname, 'crop before = ', start_crop_g 
       write(iulog,*) subname, 'pasture before = ', start_past_g 
       write(iulog,*) subname, 'prescribed crop change = ', crop_d 
       write(iulog,*) subname, 'prescribed pasture change = ', past_d
       write(iulog,*) subname, 'calced crop change = ', sum(glm_crop(:,:,np1)*cellarea*rglus) - start_crop_g
       write(iulog,*) subname, 'calced pasture change = ', sum(glm_past(:,:,np1)*cellarea*rglus) - start_past_g 
       if ( abs(sum(glm_crop(:,:,np1)*cellarea*rglus) - start_crop_g - crop_d) > diff_tol ) then
          write(iulog,*) subname, 'Warning: g cell prescribed crop change != prescribed crop change '
       endif
       if ( abs(sum(glm_past(:,:,np1)*cellarea*rglus) - start_past_g - past_d) > diff_tol ) then
          write(iulog,*) subname, 'Warning: g cell prescribed pasture change != prescribed pasture change '
       endif

#endif

    end do ! end glu loop


#ifdef DEBUG
    cellarea_forest(:,:)=0.
    where ( pot_veg > 1.)
       cellarea_forest(:,:)=cellarea(:,:)
    end where

    write(iulog,*) subname, '***** after initial land unit loop, model year = ', eclockyr
    write(iulog,*) subname, 'global GCAM forest area km^2 = ',sum(1000*gcam_forest_area(:,n)),'model year'
    write(iulog,*) subname, 'global g2g forest area km^2 = ',sum((pctland_in2015-glm_crop(:,:,n)-glm_past(:,:,n))*cellarea_forest),'model year'
    write(iulog,*) subname, 'global change in GCAM forest area km^2 = ', &
       sum(real(nint(1000*gcam_forest_area(:,np1)*100))/100.) - sum(real(nint(1000*gcam_forest_area(:,n)*100))/100.)
    write(iulog,*) subname, 'global change in g2g forest area km^2 = ', sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest) - &
       sum((pctland_in2015-glm_crop(:,:,n)-glm_past(:,:,n))*cellarea_forest)
    write(iulog,*) subname, 'global g2g crop on forest area km^2 =',sum(glm_crop(:,:,n)*cellarea_forest),'model year'
    write(iulog,*) subname, 'global g2g crop n area km^2 = ', sum(glm_crop(:,:,n)*cellarea)
    write(iulog,*) subname, 'global g2g crop np1 area km^2 = ', sum(glm_crop(:,:,np1)*cellarea)
    write(iulog,*) subname, 'global g2g crop change km^2 = ', sum(glm_crop(:,:,np1)*cellarea)-sum(glm_crop(:,:,n)*cellarea) 
    write(iulog,*) subname, 'global GCAM crop change km^2 = ', &
       real(nint(sum(1000*gcam_crop(:,np1))*100))/100. - real(nint(sum(1000*gcam_crop(:,n))*100))/100.
    write(iulog,*) subname, 'global g2g pasture n area km^2 = ',sum(glm_past(:,:,n)*cellarea)
    write(iulog,*) subname, 'global g2g pasture np1 area km^2 = ',sum(glm_past(:,:,np1)*cellarea)
    write(iulog,*) subname, 'global g2g pasture change km^2 = ',sum(glm_past(:,:,np1)*cellarea)-sum(glm_past(:,:,n)*cellarea)
    write(iulog,*) subname, 'global GCAM past change km^2 = ', &
       sum(real(nint(1000*gcam_past(:,np1)*100))/100.) - sum(real(nint(1000*gcam_past(:,n)*100))/100.)
#endif


    if (regional_unmet_reassign==1) then
       ! Here we are looping over regions, using the rgmin(r) and
       ! rgmax(r) to find the g's inside this region
       do r = 1,nreg

          regional_farea_needed = sum(unmet_farea(rgmin(r):rgmax(r)))
          
          ! number of glus in this region
          ! nrglu = rgmax(r)-rgmin(r)-1
          nrglu = rgmax(r)-rgmin(r)+1

          ! Just in case
          avail_farea = 0.
          avail_nfarea = 0.
          avail_ag_farea = 0.
          reassign_ag = 0.
          unmet_aez_farea = 0.

          do g = rgmin(r),rgmax(r)
             g1 = g-rgmin(r)+1

             cellarea_forest(:,:)=0.
             cellarea_nonforest(:,:)=0.
             fnfnonforest=0.
             fnfforest=0.
             where ( pot_veg > 1.)
                cellarea_forest(:,:)=cellarea(:,:)
                fnfforest(:,:)=1.
             elsewhere
                fnfnonforest(:,:)=1.
                cellarea_nonforest(:,:)=cellarea(:,:)
             end where
             cellarea_nonforest=cellarea_nonforest*glu_weights(g,:,:)
             cellarea_forest=cellarea_forest*glu_weights(g,:,:)
             fnfforest=fnfforest*glu_weights(g,:,:)
             fnfnonforest=fnfnonforest*glu_weights(g,:,:)

             ! above these are zero-init before each loop, becasue we
             ! now have variable nrglu per region.
             avail_farea(g1) = sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest)
             avail_nfarea(g1) = sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest)
             avail_ag_farea(g1) = sum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest)
             reassign_ag(g1) = min(avail_ag_farea(g1), avail_nfarea(g1), regional_farea_needed)
             unmet_aez_farea(g1) = regional_farea_needed - reassign_ag(g1)
          end do
          
          indxa=(/(i,i=1,nrglu)/)
          indxadn=(/(i,i=1,nrglu)/)

          call D_mrgrnk(reassign_ag,indxa,nrglu)
          call D_mrgrnk(reassign_ag*-1.,indxadn,nrglu)

          cumsum_sorted_reassign_ag(1)=reassign_ag(indxadn(1))
          do i=2,nrglu
             cumsum_sorted_reassign_ag(i)=cumsum_sorted_reassign_ag(i-1)+reassign_ag(indxadn(i))
          end do
          sortind(1) = 0
          do i1=1,nrglu
             if(cumsum_sorted_reassign_ag(i1) >= regional_farea_needed) then
                sortind(1) = i1
                exit
             endif
          enddo

          if (sortind(1).eq.0) then
             sortind(1) = nrglu
             reassign_ag_at_max_aez_ind = reassign_ag(indxadn(sortind(1)))
             unmet_regional_farea(r) = regional_farea_needed - sum(reassign_ag(indxadn))
          elseif (sortind(1)>1) then 
             reassign_ag_at_max_aez_ind = regional_farea_needed - cumsum_sorted_reassign_ag(sortind(1)-1)
             unmet_regional_farea(r) = 0
          else
             reassign_ag_at_max_aez_ind = regional_farea_needed
             unmet_regional_farea(r) = 0
          end if
          
          ! TRS - Here's my take on what's happening here.  We loop
          ! over the sortind(1) most important glus in this region,
          ! and then reassign forest in them.
#ifdef DEBUG
          write(iulog,*) subname, '***** land unit r ', r, 'model year ', eclockyr
          write(iulog,*) subname, 'regional_farea_needed = ', regional_farea_needed
          write(iulog,*) subname, 'unmet_regional_area = ', unmet_regional_farea(r)
          write(iulog,*) subname, 'sortind(1) = ', sortind(1)
#endif
          do zz=1,sortind(1)
              z=indxadn(zz)
              ! Remember, z = 1,nrglu(r).  So convert to a global g.
              g=z+rgmin(r)-1

#ifdef DEBUG
              write(iulog,*) subname, '** land unit in forest reassign loop r = ', r, 'g = ', g
              write(iulog,*) subname, 'global crop before reassign = ', sum(glm_crop(:,:,np1)*cellarea)
              write(iulog,*) subname, 'global pasture before reassign = ', sum(glm_past(:,:,np1)*cellarea)
              write(iulog,*) subname, 'g cell crop before reassign = ', sum(glm_crop(:,:,np1)*cellarea*rglus)
              write(iulog,*) subname, 'g cell pasture before reassign = ', sum(glm_past(:,:,np1)*cellarea*rglus)
#endif

              if (reassign_ag(z)>0) then


#ifdef DEBUG
              write(iulog,*) subname, 'reassign_ag = ', reassign_ag(z)
#endif
                crop_before_decrease = sum(glm_crop(:,:,np1)*cellarea)
                past_before_decrease = sum(glm_past(:,:,np1)*cellarea)
                rglus=0
                where ( glu_weights(g,:,:) .gt. 0)
                   rglus=1
                end where
                totrglus=sum(rglus)

                curr_crop = sum(glm_crop(:,:,np1)*cellarea)
                curr_crop_g = sum(glm_crop(:,:,np1)*cellarea*rglus)
                curr_past = sum(glm_past(:,:,np1)*cellarea)
                curr_past_g = sum(glm_past(:,:,np1)*cellarea*rglus)
               
                glatinds = 0
                gloninds = 0
                count_gcells = 0
                do j=1,numLats
                   do i=1,numLons
                      if (rglus(i,j) == 1) then
                         count_gcells = count_gcells+1
                         glatinds(count_gcells) = j
                         gloninds(count_gcells) = i
                      endif
                   enddo
                enddo
                if (count_gcells /= totrglus) then
                   write(iulog,*) subname, 'Error: count_gcells ', count_gcells, ' != totrglus ', totrglus
                endif
 
                cellarea_forest(:,:)=0.
                cellarea_nonforest(:,:)=0.
                fnfnonforest=0.
                fnfforest=0.
                where ( pot_veg > 1.)
                   cellarea_forest(:,:)=cellarea(:,:)
                   fnfforest(:,:)=1.
                elsewhere
                   fnfnonforest(:,:)=1.
                   cellarea_nonforest(:,:)=cellarea(:,:)
                end where
                cellarea_nonforest=cellarea_nonforest*glu_weights(g,:,:)
                cellarea_forest=cellarea_forest*glu_weights(g,:,:)
                fnfforest=fnfforest*glu_weights(g,:,:)
                fnfnonforest=fnfnonforest*glu_weights(g,:,:)
                indxup(:)=(/(i,i=1,numlons*numlats)/)
                indxdn=indxup
                pot_veg_rev=transpose(pot_veg)
                glu_weights_rev=transpose(glu_weights(g,:,:))

                ! forcing the opposite sign for non-glu cells
                ! to make sure they are outside the totrglus range
                where (transpose(rglus) == 0)
                   glu_weights_rev = -1.
                   pot_veg_rev = 1.
                end where

                call D_mrgrnk(pot_veg_rev*glu_weights_rev, indxup,numlons*numlats)
                call D_mrgrnk(pot_veg_rev*glu_weights_rev*-1., indxdn,numlons*numlats)
                     
                !jt  The sortxxxup and sortxxxdn arrays are only good 1:totrglus these arrays are also based on row major
                !jt  sorting to match original matlab scripts.
                ! these calcs reflect the transposed sort 
                sortlonsdn=(indxdn-1)/numlats+1
                sortlatsdn=mod(indxdn-1,numlats)+1
                sortlonsup=(indxup-1)/numlats+1
                sortlatsup=mod(indxup-1,numlats)+1
                
                if (allocated(v1u)) deallocate(v1u)
                if (allocated(v2u)) deallocate(v2u)
                if (allocated(v1d)) deallocate(v1d)
                if (allocated(v2d)) deallocate(v2d)

                allocate(v1u(totrglus), stat=ier)
                if(ier/=0) call mct_die(subName,'allocate v1u',ier)
                allocate(v2u(totrglus), stat=ier)
                if(ier/=0) call mct_die(subName,'allocate v2u',ier)
                allocate(v1d(totrglus), stat=ier)
                if(ier/=0) call mct_die(subName,'allocate v1d',ier)
                allocate(v2d(totrglus), stat=ier)
                if(ier/=0) call mct_die(subName,'allocate v2d',ier)
                ! these are set correctly because the cells are sorted correctly
                v1u=sortlonsup(numlons*numlats-totrglus+1:numlons*numlats)
                v2u=sortlatsup(numlons*numlats-totrglus+1:numlons*numlats)
                v1d=(sortlonsdn(:totrglus))
                v2d=(sortlatsdn(:totrglus))
               

#ifdef DEBUG
       ! check that sorted cells are all in the land unit
       j1 = 0
       do ij=1,totrglus
          ! start with lon
          found = .false.
          do i1=1,count_gcells
            if (v1d(ij) == gloninds(i1)) then
               ! check for matching lat
               if(v2d(ij) == glatinds(i1)) then
                  j1=j1+1
                  found = .true.
               endif
            endif
          enddo
          if ( .not. found ) then
             write(iulog,*) subname, 'Error: Did not find sorted cell ', ij, ' lonind=', v1d(ij), ' latind=', v2d(ij)
          endif
       enddo
       if (j1 /= totrglus) then
          write(iulog,*) subname, 'Error: found ', j1, ' sorted cells out of ', totrglus
       endif
#endif
 
                cumsum_sorted_farea=0.
                call cumsum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)
                sortind2(1) = 0
                do i1=1,totrglus
                   if(cumsum_sorted_farea(i1) >= reassign_ag(z)) then
                      sortind2(1) = i1
                      exit
                   endif
                enddo

#ifdef DEBUG
                write(iulog,*) subname, 'sortind2(1) = ' ,sortind2(1), ' totrglus(totcells) = ', totrglus
#endif

                if (sortind2(1)==0) then
                   ! may need another precision check here?
                   sortind2(1) = totrglus
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: not enough cum ag forest', cumsum_sorted_farea(totrglus) ,' to meet reassign_ag ', reassign_ag(z)
#endif
                end if
                if (sortind2(1)>1) then
                   sortsitesdn=0
                   do i=1,sortind2(1)-1
                      sortsitesdn(v1d(i),v2d(i))=1
                   end do
                   where(sortsitesdn>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) - glm_crop(:,:,np1) * fnfforest
                   end where

                   where(sortsitesdn>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,np1) * fnfforest
                   end where

                   final_area_needed =  reassign_ag(z) - cumsum_sorted_farea(sortind2(1)-1)
                else
                   final_area_needed =  reassign_ag(z)
                end if

                ag_area_avail = (glm_crop(v1d(sortind2(1)),v2d(sortind2(1)),np1) + glm_past(v1d(sortind2(1)),v2d(sortind2(1)),np1))*&
                   cellarea_forest(v1d(sortind2(1)),v2d(sortind2(1)))
                crop_area_avail = glm_crop(v1d(sortind2(1)),v2d(sortind2(1)),np1)*cellarea_forest(v1d(sortind2(1)),v2d(sortind2(1)))

#ifdef DEBUG
                write(iulog,*) subname, 'final_area_needed before crop decrease = ' , final_area_needed
                write(iulog,*) subname, 'current cell crop forest area = ', crop_area_avail
#endif

                glm_crop(v1d(sortind2(1)),v2d(sortind2(1)),np1) = &
                     glm_crop(v1d(sortind2(1)),v2d(sortind2(1)),np1) - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed / &
                     cellarea(v1d(sortind2(1)),v2d(sortind2(1))), &
                     glm_crop(v1d(sortind2(1)),v2d(sortind2(1)),np1)*fnfforest(v1d(sortind2(1)),v2d(sortind2(1))))

                final_area_needed = final_area_needed - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed,crop_area_avail)

                temp_past = glm_past(v1d(sortind2(1)),v2d(sortind2(1)),np1)*cellarea_forest(v1d(sortind2(1)),v2d(sortind2(1)))

#ifdef DEBUG
                write(iulog,*) subname, 'final_area_needed before pasture decrease =' , final_area_needed
                write(iulog,*) subname, 'current cell pasture forest area = ', temp_past
#endif

                glm_past(v1d(sortind2(1)),v2d(sortind2(1)),np1) = &
                     glm_past(v1d(sortind2(1)),v2d(sortind2(1)),np1) - &
                     min(final_area_needed / cellarea(v1d(sortind2(1)),v2d(sortind2(1))), &
                     glm_past(v1d(sortind2(1)),v2d(sortind2(1)),np1)*fnfforest(v1d(sortind2(1)),v2d(sortind2(1))))

                ! this should be zero after this calc
                final_area_needed = final_area_needed - &
                     min(final_area_needed,temp_past)
                if (final_area_needed < 0) then
                   final_area_needed = 0.
                endif
                
                crop_after_decrease = sum(glm_crop(:,:,np1)*cellarea)
                past_after_decrease = sum(glm_past(:,:,np1)*cellarea)
                total_ag_decrease = crop_before_decrease+past_before_decrease - crop_after_decrease - past_after_decrease
                crop_decrease_ratio = (crop_before_decrease - crop_after_decrease)/(total_ag_decrease+1e-12)
                past_decrease_ratio = 1-crop_decrease_ratio

#ifdef DEBUG
                write(iulog,*) subname, '*reassign to meet forest demand (km^2), based on global values'               
                write(iulog,*) subname, 'crop decrease = ', (crop_before_decrease - crop_after_decrease)
                write(iulog,*) subname, 'past decrease = ', (past_before_decrease - past_after_decrease)
                write(iulog,*) subname, 'final_area_needed (should be zero) = ', final_area_needed
                if (abs(final_area_needed) > diff_tol) then
                   write(iulog,*) subname, 'Warning: final_area_needed not zero!'
                endif

                write(iulog,*) subname, 'global crop after decrease = ', sum(glm_crop(:,:,np1)*cellarea)
                write(iulog,*) subname, 'global pasture after decrease = ', sum(glm_past(:,:,np1)*cellarea)
                write(iulog,*) subname, 'g cell crop after decrease = ', sum(glm_crop(:,:,np1)*cellarea*rglus)
                write(iulog,*) subname, 'g cell pasture after decrease = ', sum(glm_past(:,:,np1)*cellarea*rglus)
#endif
                reassign_ag(z) = reassign_ag(z) - final_area_needed

                avail_land0 = 0.
                
                where ( glu_weights(g,:,:) > 0 .and.(glm_crop(:,:,np1)>0..or.glm_past(:,:,np1)>0.))
                   avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                end where
                sumavail_land0=sum(avail_land0)
                if (sumavail_land0 >=reassign_ag(z)) then 

                   where(glu_weights(g,:,:) > 0 .and. (glm_crop(:,:,np1)>0. .or. glm_past(:,:,np1)>0.))
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z)/ cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z) /cellarea
                   end where

#ifdef DEBUG
                   write(iulog,*) subname, 'crop and pasture increase on nonforest - land available'
                   write(iulog,*) subname, 'actual crop increase = ', sum(crop_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z))
                   write(iulog,*) subname, 'actual pasture increase = ', sum(past_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z))
#endif

                else
                   avail_landA = 0.
                   where ( glu_weights(g,:,:) > 0 )
                      avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                   end where
                   sumavail_landA=sum(avail_landA)
                   if (reassign_ag(z) > sumavail_landA) then 
                      reassign_ag(z) = sum(avail_landA)
#ifdef DEBUG
                      write(iulog,*) subname, 'Warning: unavailable nonforest'
                      write(iulog,*) subname, '--when trying to meet reg forest'
                      write(iulog,*) subname, 'reducing conversion to sumavail_landA'
                      write(iulog,*) subname, 'reassign_ag(',z,')=',reassign_ag(z)&
                                 ,'sumavail_landA=', sumavail_landA, &
                                 'sumavail_land0=', sumavail_land0
#endif
                   end if
#ifdef DEBUG
                   write(iulog,*) subname, 'Warning: cropland and pasture increase on non-forest - land not available'
#endif
                   where(glu_weights(g,:,:) > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*(avail_land0)/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*(avail_land0)/cellarea
                   end where
                      
                   if (sumavail_landA-sumavail_land0 > 0.) then
                      where(glu_weights(g,:,:) > 0)
                         glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*((avail_landA-avail_land0) / &
                              (sumavail_landA-sumavail_land0)*(reassign_ag(z)-sumavail_land0))/cellarea
                         glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*((avail_landA-avail_land0)/ &
                              (sumavail_landA-sumavail_land0)*(reassign_ag(z)-sumavail_land0))/cellarea
                   end where
                   endif
                end if

#ifdef DEBUG
                write(iulog,*) subname, 'global km^2: crop before decrease ', crop_before_decrease, ' crop after increase ', sum(glm_crop(:,:,np1)*cellarea)
                write(iulog,*) subname, 'global km^2: pasture before decrease ', past_before_decrease, ' pasture after increase ', sum(glm_past(:,:,np1)*cellarea)
#endif

             end if

#ifdef DEBUG
             write(iulog,*) subname, 'global crop after reassign km^2 = ', sum(glm_crop(:,:,np1)*cellarea)
             write(iulog,*) subname, 'global pasture after reassign km^2 = ', sum(glm_past(:,:,np1)*cellarea)             
             write(iulog,*) subname, 'g cell crop after reassign km^2 = ', sum(glm_crop(:,:,np1)*cellarea*rglus)
             write(iulog,*) subname, 'g cell pasture after reassign km^2 = ', sum(glm_past(:,:,np1)*cellarea*rglus)

             write(iulog,*) subname, '*reassign crop change (global) ',sum(glm_crop(:,:,np1)*cellarea) - curr_crop
             write(iulog,*) subname, '*reassign crop change (g) ', sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g
             write(iulog,*) subname, '*reassign pasture change (global) ',sum(glm_past(:,:,np1)*cellarea) - curr_past
             write(iulog,*) subname, '*reassign pasture change (g) ', sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g

             if (abs(sum(glm_crop(:,:,np1)*cellarea*rglus) - curr_crop_g) > diff_tol) then
                write(iulog,*) subname, 'Warning: crop change (g) not zero!'
             endif
             if (abs(sum(glm_past(:,:,np1)*cellarea*rglus) - curr_past_g) > diff_tol) then
                write(iulog,*) subname, 'Warning: pasture change (g) not zero!'
             endif
#endif
          end do ! end z loop
       end do !  end r loop
    end if   ! ! if regional_unmet_reassign


#ifdef DEBUG
    write(iulog,*) subname, '***** after reassign loop'
    write(iulog,*) subname, 'sum final gcrop/gpast n'
    write(iulog,fmt="(1ES25.15)") sum(glm_crop(:,:,n))
    write(iulog,fmt="(1ES25.15)") sum(glm_past(:,:,n))
    write(iulog,*) 'sum final gcrop/gpast np1'
    write(iulog,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(iulog,fmt="(1ES25.15)") sum(glm_past(:,:,np1))

    cellarea_forest(:,:)=0.
    where ( pot_veg > 1.)
    cellarea_forest(:,:)=cellarea(:,:)
    end where
    
    write(iulog,*) subname, '*** these are the final values'
    write(iulog,*) 'final np1 g2g global forest area km^2', sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest)
    write(iulog,*) 'final change in g2g global forest area km^2', sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest) - &
       sum((pctland_in2015-glm_crop(:,:,n)-glm_past(:,:,n))*cellarea_forest)
    write(iulog,*) subname, 'global change in GCAM forest area km^2 = ', &
       sum(real(nint(1000*gcam_forest_area(:,np1)*100))/100.) -sum(real(nint(1000*gcam_forest_area(:,n)*100))/100.) 
    write(iulog,*) 'final np1 g2g global crop on forest area km^2', sum(glm_crop(:,:,np1)*cellarea_forest)
    write(iulog,*) subname, 'global g2g crop n area km^2 = ', sum(glm_crop(:,:,n)*cellarea)
    write(iulog,*) subname, 'global g2g crop np1 area km^2 = ', sum(glm_crop(:,:,np1)*cellarea)
    write(iulog,*) subname, 'global g2g crop change km^2 = ', sum(glm_crop(:,:,np1)*cellarea)-sum(glm_crop(:,:,n)*cellarea)
    write(iulog,*) subname, 'global GCAM crop change km^2 = ', &
       real(nint(1000*sum(gcam_crop(:,np1))*100))/100. - real(nint(1000*sum(gcam_crop(:,n))*100))/100.
    write(iulog,*) subname, 'global unmet forest increase km^2 = ', sum(unmet_regional_farea)

    write(iulog,*) subname, 'global g2g pasture n area km^2 = ',sum(glm_past(:,:,n)*cellarea)
    write(iulog,*) subname, 'global g2g pasture np1 area km^2 = ',sum(glm_past(:,:,np1)*cellarea)
    write(iulog,*) subname, 'global g2g pasture change km^2 = ',sum(glm_past(:,:,np1)*cellarea)-sum(glm_past(:,:,n)*cellarea)
    write(iulog,*) subname, 'global GCAM pasture change km^2 = ', &
       sum(real(nint(1000*gcam_past(:,np1)*100))/100.) - sum(real(nint(1000*gcam_past(:,n)*100))/100.)
    write(iulog,*) subname, 'global unmet forest increase km^2 = ', sum(unmet_regional_farea(:))
#endif 

 end if ! end if gcam runs this model year

 if (allocated(v1u)) deallocate(v1u)
 if (allocated(v2u)) deallocate(v2u)
 if (allocated(v1d)) deallocate(v1d)
 if (allocated(v2d)) deallocate(v2d)

!
! Interpolate data to years needed by GLM
!
! glm calculates year+1 using 
! gcam constructed states at year+1
! and harvest transitions at year
!
! note that eclock is the model year
! if the years are correct the factors should be correct
! if GCAM is at annual time step, these should be 0 and 1

    eclockyr=ymd/10000
    eclockyrp1=eclockyr+1
    delyr= year2-year1
    fact1yrp1=(year2-eclockyrp1)/delyr
    fact2yrp1=(eclockyrp1-year1)/delyr
    fact1=(year2-eclockyr)/delyr
    fact2=(eclockyr-year1)/delyr

#ifdef DEBUG
     write(iulog,*) subname, 'crop/pasture interpolation factors fact1yrp1,fact2yrp1,year1,year2=',fact1yrp1,fact2yrp1,year1,year2
     write(iulog,*) subname, 'wood harvest interpolation factors fact1,fact2,year1,year2=',fact1,fact2,year1,year2
#endif

! use eclock year year for interpolating crop past and othr
    glm_crop_ann(:,:)=glm_crop(:,:,n)*fact1yrp1+glm_crop(:,:,np1)*fact2yrp1
    glm_past_ann(:,:)=glm_past(:,:,n)*fact1yrp1+glm_past(:,:,np1)*fact2yrp1
    glm_othr_ann(:,:)=pctland_in2015-glm_past_ann-glm_crop_ann

! use previous year for fractions fact1 and fact2 for woodharvest
    glm_wh_ann(:)=gcam_wh(:,n)*fact1+gcam_wh(:,np1)*fact2

!    do j=1,360
!       do i=1,720
!          glm_crop_ann(i,j)=fround(glm_crop_ann(i,j),6)
!          glm_past_ann(i,j)=fround(glm_past_ann(i,j),6)
!          glm_othr_ann(i,j)=fround(glm_othr_ann(i,j),6)
!       end do
!    end do


!  Conversion factor of 1.3 to account for a slash fraction 
!  of 30% for non-harvested carbon lost as a result of wood harvest

    glm_wh_ann(:)=glm_wh_ann(:)*1.3

    write(yearc,fmt="(I4)") ymd/10000
    open (unit=55,file='gcrop'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_crop_ann(:,j)
    end do
    close(55)

    open (unit=55,file='gpast'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_past_ann(:,j)
    end do
    close(55)

    open (unit=55,file='gothr'//yearc//'.txt',action="write",status="unknown")
    do j=1,360
       write(55,fmt="(720F9.6)") glm_othr_ann(:,j)
    end do
    close(55)

    ij=0
    do j=1,numLats
       do i=1,numLons
          ij=ij+1
          glmi(iac_glmi_cropland,ij)=glm_crop_ann(i,j)
          glmi(iac_glmi_pasture,ij)=glm_past_ann(i,j)
          glmi(iac_glmi_natveg,ij)=glm_othr_ann(i,j)
       end do
    end do
    glmi_wh(:)=glm_wh_ann(:)

!
! If we are at the boundary year of the gcam data. Need to move np1 info
! into timelevel n to calculate the next set of years.
! also advance the gcam years
! this should always be true if gcam is at annual time step
    if (eclockyr==(year2-1)) then 
       glm_crop(:,:,n)=glm_crop(:,:,np1)
       glm_past(:,:,n)=glm_past(:,:,np1)
       gcam_crop(:,n)=gcam_crop(:,np1)
       gcam_past(:,n)=gcam_past(:,np1)
       gcam_wh(:,n)=gcam_wh(:,np1)
       gcam_forest_area(:,n)=gcam_forest_area(:,np1)
       cdata%i(iac_cdatai_gcam_yr1) = year2
       cdata%i(iac_cdatai_gcam_yr2) = cdata%i(iac_cdatai_gcam_yr1) + &
          iac_gcam_timestep
    end if

    ! write a restart file each year,

    call shr_cal_date2ymd(ymd,year,mon,day)
    write(filename,'(a,i4.4,a,i2.2,a)') trim(gcam2glm_restfile)//'r.',year+1,'.nc'

    write(iulog,*) subname,' write_restart file     ',trim(filename)

    status= nf90_create(filename,nf90_clobber,ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_dim(ncid,'lon',numlons,dimid3(1))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'lat',numlats,dimid3(2))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'time',2,dimid3(3))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'lu_gcam',num_gcam_land_regions,dimid_gcam(1))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_inq_dimid(ncid,'time',dimid_gcam(2))
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(1)
    status = nf90_def_var(ncid,'lon',NF90_DOUBLE,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","degrees_east")
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(2)
    status = nf90_def_var(ncid,'lat',NF90_DOUBLE,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","degrees_north")
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'gcam_crop',NF90_DOUBLE,dimid_gcam,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","thous_sqkm")
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'gcam_past',NF90_DOUBLE,dimid_gcam,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","thous_sqkm")
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'gcam_wh',NF90_DOUBLE,dimid_gcam,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","MgC")
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'gcam_forest_area',NF90_DOUBLE,dimid_gcam,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"units","thous_sqkm")
    if(status /= nf90_NoErr) call handle_err(status)

    dimid = dimid3(3)
    status = nf90_def_var(ncid,'gcam_years',NF90_INT,dimid,varid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'glm_crop',NF90_DOUBLE,dimid3,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"fraction","frac")
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"missing_value",miss_val)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_var(ncid,'glm_past',NF90_DOUBLE,dimid3,varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"fraction","frac")
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_att(ncid,varid,"missing_value",miss_val)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_enddef(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'lon',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,lon)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'lat',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,lat)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_crop',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,gcam_crop)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_past',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,gcam_past)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_wh',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,gcam_wh)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_forest_area',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,gcam_forest_area)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'gcam_years',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,(/cdata%i(iac_cdatai_gcam_yr1),cdata%i(iac_cdatai_gcam_yr2)/))
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'glm_crop',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,glm_crop)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_inq_varid(ncid,'glm_past',varid)
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_put_var(ncid,varid,glm_past)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    deallocate(indxup)
    deallocate(indxdn)
    deallocate(sortlatsup)
    deallocate(sortlatsdn)
    deallocate(sortlonsup)
    deallocate(sortlonsdn)
    deallocate(tmparr)
    deallocate(indxa)
    deallocate(indxadn)
    deallocate(avail_farea)
    deallocate(avail_nfarea)
    deallocate(avail_ag_farea)
    deallocate(reassign_ag)
    deallocate(unmet_aez_farea)
    deallocate(cumsum_sorted_reassign_ag)
    deallocate(unmet_regional_farea)

    deallocate(glatinds)
    deallocate(gloninds)
    deallocate(temp_past_grid)

  end subroutine gcam2glm_run_mod
  
  
!---------------------------------------------------------------------------
!BOP
  
! !IROUTINE: gcam2glm_final_mod

! !INTERFACE:
  subroutine gcam2glm_final_mod( )

! !DESCRIPTION:
! Finalize glm model
! !USES:
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(gcam2glm_final_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP

!---------------------------------------------------------------------------

    deallocate(hydeGCROP2015)
    deallocate(hydeGPAST2015)
    deallocate(hydeGOTHR2015)
    deallocate(hydeGSECD2015)
    deallocate(hydeGWH2015)
    deallocate(cellarea)
    deallocate(cellarea_forest)
    deallocate(cellarea_nonforest)
    deallocate(glm_crop_ann)
    deallocate(glm_past_ann)
    deallocate(glm_othr_ann)
    deallocate(cumsum_sorted_farea)
    deallocate(glm_wh_ann)
    deallocate(fnfforest)
    deallocate(fnfnonforest)
    deallocate(pot_veg)
    deallocate(pot_veg_rev)
    deallocate(crop_area)
    deallocate(pctland_in2015)
    deallocate(sortsitesup)
    deallocate(sortsitesdn)
    deallocate(rglus)
    deallocate(datearr)
    deallocate(glm_crop)
    deallocate(glm_past)
    deallocate(gcam_crop)
    deallocate(gcam_wh)
    deallocate(gcam_past)
    deallocate(gcam_forest_area)
    deallocate(unmet_neg_past)
    deallocate(unmet_neg_crop)
    deallocate(unmet_farea)
    deallocate(avail_land0)
    deallocate(avail_landA)
    deallocate(lon)
    deallocate(lat)

    deallocate(glu_weights)
    deallocate(glu_weights_rev)

    deallocate(rgmin)
    deallocate(rgmax)

    deallocate(gcamo_base)

  end subroutine gcam2glm_final_mod
!====================================================================================
real(r8) function fround(n, d)
real(r8) n
integer d
  fround= floor(n * 10.**d + .5) / 10.**d
end function fround
!====================================================================================
 subroutine handle_err (status)
    implicit none
    integer, intent (in) :: status
    print *, nf90_strerror(status)
    stop
 end subroutine handle_err
 subroutine cumsum (arrayin,sortcol,sortrow,cumsumvecout,n)
   implicit none
   real(r8), intent (in) :: arrayin(:,:)
   integer, intent (in) :: sortcol(:),sortrow(:)
   integer, intent (in) :: n
   real(r8), intent (out):: cumsumvecout(:)
   integer             ::j

   cumsumvecout(1)=arrayin(sortcol(1),sortrow(1))
   do j=2,n
      cumsumvecout(j)=cumsumvecout(j-1)+arrayin(sortcol(j),sortrow(j))
   end do
 end subroutine cumsum
!====================================================================================
Subroutine D_mrgrnk (XDONT, IRNGT,N)
!                                                                                                                                                      
!   MRGRNK = Merge-sort ranking of an array                                                                                                            
!   For performance reasons, the first 2 passes are taken                                                                                              
!   out of the standard loop, and use dedicated coding.                                                                                                
!                                                                                                                                                      
!   The routine is part of ORDERPACK 2.0 -- Unconditional, Unique, and Partial Ranking,                                                                
!   Sorting, and Permutation Downloadable Fortran 90 source code                                                                                       
!   Author: Michel Olagnon                                                                                                                             
!   http://www.fortran-2000.com/rank/                                                                                                                  
!   Users can freely download ORDERPACK 2.0 from this site.                                                                                            
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (N), Intent (In) :: XDONT
      Integer, Dimension (N), Intent (Out) :: IRNGT
      Integer, Intent (In) :: N
! __________________________________________________________
      Real (kind=kdp) :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine D_mrgrnk


!====================================================================================
end module gcam2glm_mod

