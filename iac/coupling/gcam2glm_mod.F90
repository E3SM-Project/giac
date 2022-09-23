#define DEBUG
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
		 hydeGSECD2015,	  &
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

    integer :: n,np1,nflds,gcamsize,nglu,nregions,io,r,g,g1
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
    character(len=*),parameter :: gcam2glm_rpointer = 'rpointer.gcam2glm'
    character(256) :: filename

! !REVISION HISTORY:
! Author: T Craig


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
    logical :: restart_run,lexist
    logical :: initial_run
    integer :: iun,tmpyears(2),ier,t,yy
    real(r8) :: v
    character(len=*),parameter :: subname='(gcam2glm_init_mod)'

    !character(len=512) :: gcam2glm_baselu
    !character(len=512) :: gcam2glm_glumap
    !character(len=512) :: gcam2glm_basebiomass

    character(len=512) :: dum

! !REVISION HISTORY:
! Author: T Craig
! Author: JET          ! rewrite of matlat preprocessing script new_grids_matchforest4.m

!EOP
!-----------------------------------------------------------------------

#ifdef DEBUG
     write(iulog,*) subname,' starting subroutine '
#endif
    restart_run  = cdata%l(iac_cdatal_rest)

    nregions=num_gcam_energy_regions
    nglu=num_gcam_land_regions

    gcamsize=nglu
    initial_run = cdata%l(iac_cdatal_initrun)

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
    if(ier/=0) call mct_die(subName,'allocate gcam_forest_ara',ier)
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

!avd
write(iulog,*) subname,'glu_weights dims: ',gcamsize,numLons,numLats

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
    tmp=MAXLOC(datearr,mask=datearr.EQ.2015)
    start3(3)=tmp(1)
    count3(3)=1
    start2=1
    count2(1)=numLons
    count2(2)=numLats
! read in hyde data - now from the single glm init file

    status = nf90_inq_varid(ncid,'gcrop',varid)
    status = nf90_get_var(ncid,varid,hydeGCROP2015,start3,count3)
    status = nf90_inq_varid(ncid,'gpast',varid)
    status = nf90_get_var(ncid,varid,hydeGPAST2015,start3,count3)
! need to add together the primary (gothr) and secondary (gsecd) from the initial file to get total other
    status = nf90_inq_varid(ncid,'gothr',varid)
    status = nf90_get_var(ncid,varid,hydeGOTHR2015,start3,count3)
    status = nf90_inq_varid(ncid,'gsecd',varid)
    status = nf90_get_var(ncid,varid,hydeGSECD2015,start3,count3)
    status = nf90_inq_varid(ncid,'cell_area',varid)
    status = nf90_get_var(ncid,varid,cellarea,start3,count3)
    status = nf90_close(ncid)
    hydeGOTHR2015 = hydeGOTHR2015 + hydeGSECD2015
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
    glu_weights=0.
    rgmin = 1000
    rgmax = 0

!avd
write(iulog,*) subname,'nglu=',nglu

    open(5,file=gcam2glm_glumap)
    ! Read header
    read(5,*) dum
    do
       read(5,*,iostat=io) r,g,x,y,region_name,glu_name,weight
       if (io < 0) then
          exit
       endif

       ! Need to find the lon and lat index for these x,y (which are lon,lat values at center of cell)
       ! Probably should do some idiot checking to see if 
       !lonx=findloc(lon,x)
       !laty=findloc(lat,y)

       ! If findloc isn't available (it's f2008 standard), do it the
       ! 70s way
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

       ! Find region of glus for each r
       if (g .lt. rgmin(r)) rgmin(r) = g
       if (g .gt. rgmax(r)) rgmax(r) = g

    end do

    close(5)
    
    cellarea_nonforest(:,:)=0.
    cellarea_forest(:,:)=0.
    pot_veg=pot_veg*0.75
    pctland_in2015(:,:) = 0.
    pctland_in2015=hydeGCROP2015+hydeGPAST2015+hydeGOTHR2015

    ! Read in gcam base file 
    open(5,file='./gcamo_base.csv')
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


! KVC: Temp -- set this to false. 
restart_run = .false.
    if (.not.restart_run) then
       glm_crop(:,:,1)=hydeGCROP2015;
       glm_past(:,:,1)=hydeGPAST2015;

       ! Set initial gcam land
       ! current baseline is still in billion m^3 of biomass
       !   need MgC here; conv fact is 0.288 tonnes C per m^3 (MgC per m^3)
       gcam_crop(:,1) = gcamo_base(iac_gcamo_crop,:)
       gcam_past(:,1) = gcamo_base(iac_gcamo_pasture,:)
       gcam_wh(:,1) = gcamo_base(iac_gcamo_woodharv,:) * 288000000.
       gcam_forest_area(:,1) = gcamo_base(iac_gcamo_forest,:)

       ! Set years 
       cdata%i(iac_cdatai_gcam_yr1)=iac_start_year
       cdata%i(iac_cdatai_gcam_yr2)=iac_start_year+iac_gcam_timestep
    else
       ! read restart and set crop and past
      
!!!! avd - need to add the gcam crop, past and forest area and wh to the restart file!
!!!! and use them here!
 
       inquire(file=trim(gcam2glm_rpointer),exist=lexist)
       if (lexist) then
          
#ifdef DEBUG
           write(iulog,*) subname,' read_restart rpointer ',trim(gcam2glm_rpointer)
#endif
          
          iun = shr_file_getunit()
          open(iun,file=trim(gcam2glm_rpointer),form='formatted')
          read(iun,'(a)') filename
          close(iun)
          call shr_file_freeunit(iun)
          
#ifdef DEBUG
           write(iulog,*) subname,' read_restart file ',trim(filename)
#endif
          
          inquire(file=trim(filename),exist=lexist)
          if (.not.lexist) then
             write(iulog,*) subname,' ERROR: missing file ',trim(filename)
             call shr_sys_abort(subname//' ERROR: missing file')
          endif
          
          status= nf90_open(filename,nf90_nowrite,ncid)
          if(status /= nf90_NoErr) call handle_err(status)

	  ! Need to set cdata year1 and year 2 for restart

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
       else
          write(iulog,*) subname,' read_restart rpointer NOT found ',trim(gcam2glm_rpointer)
          call shr_sys_abort(subname//' ERROR: missing file')
       end if ! rpointer exist
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

    character(len=*),parameter :: subname='(gcam2glm_run_mod)'

    real(r8), parameter :: crop_forest_abandon_percent = 0.9
    real(r8), parameter :: past_forest_abandon_percent = 0.9

! !LOCAL VARIABLES:

    character(len=512) :: dum
    character*4 :: yearc
    character(256) :: filename
    integer :: i,j,ij,r,i1,j1,aez,ind,h,z
    integer :: row,g,t,y
    integer :: iun,iyr,ier
    integer :: ymd, tod, dt,naez,nreg,ii,year,mon,day
    logical :: restart_now,gcam_alarm
    real(r8)  :: crop_d,past_d,crop_neg,crop_pos,past_neg,past_pos,farea_d,v
    real(r8)  :: gcam_crop_tmp(2,18,14),gcam_past_tmp(2,18,14),gcam_forest_area_tmp(2,18,14)
    real(r8)  :: fact1,fact2,eclockyr,fact1yrp1, fact2yrp1,delyr,eclockyrp1
    real(r8)  :: tmp0
    integer,save :: ncid,varid,dimid,dimid3(3)
    integer :: totrglus
    integer, allocatable  ::indxdn(:),indxup(:),sortlatsup(:),sortlonsup(:),sortlatsdn(:),sortlonsdn(:),indxa(:),indxadn(:),v1u(:),v2u(:),v1d(:),v2d(:)
    real(r8), allocatable   :: tmparr(:)
    integer :: sortind(1),max_aez_ind(1),regional_unmet_reassign
    real(r8)  :: crop_pos_f,crop_pos_nf,pot_forest_to_crop,pot_forest_from_crop,reassign_crop, &
               crop_neg_f,crop_neg_nf,crop_nfarea,crop_farea, &
               past_pos_f,past_pos_nf,pot_forest_to_past,pot_forest_from_past,reassign_past, &
               past_neg_f,past_neg_nf,past_nfarea,past_farea,GLM_nfarea,GLM_farea,regional_farea_needed, &
               crop_after_decrease,past_after_decrease,crop_before_decrease,past_before_decrease, &
               total_ag_decrease,crop_decrease_ratio,past_decrease_ratio,ag_area_avail,crop_area_avail, &
               final_area_needed,f_diff,reassign_ag_at_max_aez_ind,sumavail_land0,sumavail_landA

   real(r8), allocatable    ::  avail_farea(:),avail_nfarea(:),avail_ag_farea(:),reassign_ag(:), &
                              unmet_aez_farea(:),cumsum_sorted_reassign_ag(:),unmet_regional_farea(:)
   integer ntimes,nntimes,zz

! avd
integer :: nmode, ierr
character(len=128) :: hfile

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------
    regional_unmet_reassign=1
    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)
    gcam_alarm=(EClock(iac_EClock_Agcam)==1)
#ifdef DEBUG
    write(iulog,*) trim(subname),' date= ',ymd,tod
#endif

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

    avail_farea=0.
    avail_nfarea=0.
    avail_ag_farea=0.
    reassign_ag=0.
    unmet_aez_farea=0.
    cumsum_sorted_reassign_ag=0.
    unmet_regional_farea=0.
    unmet_farea=0.

    ! KVC: This does not appear to be set. Manually setting for now
    gcam_alarm = .true.
    if (gcam_alarm) then

    ! TRS refactor
    ! This stuff attempts to set to default anywhere not covered by an
    ! aez.  I'm not sure if that's still a concern here, but for now
    ! we'll have to implement it as a loop, because we don't have a
    ! one-to-one mapping of lat,lon to regions anymore
    !where (  aez_regions  >= 1 .and.   aez_regions <= 14)
    !where ( glu_val >= 1 .and. glu_val <= nglu )
    !   glm_crop(:,:,np1)=0
    !   glm_past(:,:,np1)=0
    !elsewhere
    !   glm_crop(:,:,np1)=hydeGCROP2015;
    !   glm_past(:,:,np1)=hydeGPAST2015;
    !end where

    glm_crop(:,:,np1)=hydeGCROP2015;
    glm_past(:,:,np1)=hydeGPAST2015;

    do g=1,nglu
       where (glm_crop(:,:,np1) == 0 .or. glu_weights(g,:,:) > 0)
          glm_crop(:,:,np1)=0
          glm_past(:,:,np1)=0
       end where
    enddo

#ifdef notdef
    ! Ugly triple loop
    do g=1,nglu
       do lonx=1,numLons
          do latx,1,numLats
             if (glu_weights(g,lonx,latx) > 0) then 
                glm_crop(lonx,latx,np1)=0
                glm_past(lonx,latx,np1)=0
             endif
          enddo
       enddo
    enddo
#endif 

#ifdef DEBUG
    write(6,*) 'sum 0 glm_crop='
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(6,*) 'sum 0 glm_past='
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,np1))
#endif     
    ! Unpack gcamo field
    ! the previous field (n) is initialized by file or by previous year, so read into next (np1) field
    
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
    ind=1

    !do r = 1,nreg
    !   do aez = 1,naez
    do g=1,nglu
       ! Round GCAM results to two decimal places before using as LUC
       crop_d = (real(nint(gcam_crop(g,np1)*100)/100)-real(nint(gcam_crop(g,n)*100)/100))*1000
       past_d = (real(nint(gcam_past(g,np1)*100)/100)-real(nint(gcam_past(g,n)*100)/100))*1000
       farea_d = (real(nint(gcam_forest_area(g,np1)*100)/100)-real(nint(gcam_forest_area(g,n)*100)/100))*1000

       ! convert area changes into positive and negative changes

       if (crop_d <= 0) then
          crop_neg = -1.*crop_d
          crop_pos = 0.
       else
          crop_neg = 0.
          crop_pos = crop_d
       end if
       if (past_d <= 0) then
          past_neg = -1.*past_d
          past_pos = 0.
       else
          past_neg = 0.
          past_pos = past_d
       end if

       ! compute the potentially forested and potentially non-forested
       ! area currently occupied by cropland or pasture  in each
       ! region/AEZ. Also compute the *actual* forested and non-forested area
       ! currently available in each region/AEZ        

       ! TRS - we are now looping by glu.  So instead of finding this
       ! specific combo of aez and region, each glu implies a region.

       ! Find lon/lat where we have any weights at all
       rglus=0
       where ( glu_weights(g,:,:) .gt. 0 )
          rglus=1
       end where

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

       crop_farea =sum(glm_crop(:,:,n)*cellarea_forest)
       crop_nfarea =sum(glm_crop(:,:,n)*cellarea_nonforest)
       past_farea =sum(glm_past(:,:,n)*cellarea_forest)
       past_nfarea =sum(glm_past(:,:,n)*cellarea_nonforest)
       GLM_nfarea = sum((pctland_in2015 - glm_crop(:,:,n) - glm_past(:,:,n))*cellarea_nonforest)
       GLM_farea = sum((pctland_in2015 - glm_crop(:,:,n) - glm_past(:,:,n))*cellarea_forest)
       totrglus=sum(rglus)
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
       pot_veg_rev=transpose(pot_veg)
       glu_weights_rev=transpose(glu_weights(g,:,:))
       call D_mrgrnk(pot_veg_rev*glu_weights_rev,indxup,numlons*numlats)
       call D_mrgrnk(pot_veg_rev*glu_weights_rev*-1.,indxdn,numlons*numlats)
       !	  
       !  The sortxxxup and sortxxxdn arrays are only good 1:totrglus these arrays are also based on arrays lat first (360,720)
       !  sorting to match original matlab scripts.  Thats why I transposed the arrays above.  Also didn't use qsort because it
       !  isn't a stable sort.
       !
       sortlonsdn=(indxdn-1)/numlats+1
       sortlatsdn=mod(indxdn-1,numlats)+1
       sortlonsup=(indxup-1)/numlats+1
       sortlatsup=mod(indxup-1,numlats)+1
       v1u=sortlonsup(numlons*numlats-totrglus+1:numlons*numlats)
       v2u=sortlatsup(numlons*numlats-totrglus+1:numlons*numlats)
       v1d=(sortlonsdn(:totrglus))
       v2d=(sortlatsdn(:totrglus))

       if ((past_neg >= (past_farea + past_nfarea)).and.(past_neg>0)) then
          unmet_neg_past(ind) = past_neg - (past_farea + past_nfarea)
          past_neg = past_farea + past_nfarea
       end if
       if ((crop_neg >= (crop_farea + crop_nfarea)).and.(crop_neg>0)) then
          unmet_neg_crop(ind) = crop_neg - (crop_farea + crop_nfarea)
          crop_neg = crop_farea + crop_nfarea
       end if

       ! compute the potential forest that could be created from cropland
       ! and pasture abandonment
       if (crop_neg <= crop_farea) then
          pot_forest_from_crop = crop_neg
       else
          pot_forest_from_crop = crop_farea
       end if

       if (past_neg <= past_farea) then
          pot_forest_from_past = past_neg
       else
          pot_forest_from_past = past_farea
       end if

       ! compute the potential forest that could be lost from cropland and
       ! pasture expansion
       if ((crop_pos>0).and.(crop_pos > GLM_nfarea)) then
          pot_forest_to_crop = crop_pos - GLM_nfarea
          ! QUESTION: what about using abandoned pasture if needed?
          ! (MODIFY IN FUTURE?)
       else
          pot_forest_to_crop = 0
       end if

       if ((past_pos>0).and.(past_pos > (GLM_nfarea - (crop_pos - pot_forest_to_crop)))) then
          pot_forest_to_past = past_pos - (GLM_nfarea - (crop_pos - pot_forest_to_crop))
          ! QUESTION: what about using abandoned cropland if needed?
          ! (MODIFY IN FUTURE)
       else
          pot_forest_to_past = 0
       end if

       ! if GCAM forest area not expanding, compute the crop and pasture
       ! changes that need to come from forest and non-forest
       if (farea_d <= 0) then
          ! farea not expanding
          ! MODIFY THIS IN FUTURE? SO THAT WE DON'T "OVERSHOOT" FAREA_D
          crop_pos_f = pot_forest_to_crop
          crop_pos_nf = crop_pos - pot_forest_to_crop
          crop_neg_f = pot_forest_from_crop
          crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
          past_pos_f = pot_forest_to_past
          past_pos_nf = past_pos - pot_forest_to_past
          past_neg_f = pot_forest_from_past
          past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)

       elseif ((pot_forest_from_crop + pot_forest_from_past - pot_forest_to_crop - pot_forest_to_past) >= farea_d) then
          ! compute the areas of crop and past needed to move off forest

          ! enough forested area available just by abandoning
          ! cropland and pasture on pot. forest
          ! compute forest percentages etc

          ! MODIFY THIS IN FUTURE? ... currently doing max forest abaondon ... might need less
          ! WHAT ABOUT CROPLAND ABANDONMENT THAT ALLOWS PASTURE
          ! EXPANSION?
          crop_pos_f = pot_forest_to_crop
          crop_pos_nf = crop_pos - pot_forest_to_crop
          crop_neg_f = pot_forest_from_crop
          crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
          past_pos_f = pot_forest_to_past
          past_pos_nf = past_pos - pot_forest_to_past
          past_neg_f = pot_forest_from_past
          past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)

       else
          ! in addition to abandoning crop and past on pot. forest,
          ! also reassign some crop and past to non-forest
          ! compute forest percentages and also reassignment
          ! percentages
          f_diff = farea_d - (pot_forest_from_crop + pot_forest_from_past - pot_forest_to_crop - pot_forest_to_past)
          if (((GLM_nfarea - (crop_pos - pot_forest_to_crop) -(past_pos - pot_forest_to_past))>=f_diff) .and. ((crop_farea+past_farea - pot_forest_from_crop - pot_forest_from_past)>=f_diff)) then
             reassign_crop = (crop_farea - pot_forest_from_crop)/(crop_farea + past_farea - pot_forest_from_crop - pot_forest_from_past)*f_diff
             reassign_past = (past_farea - pot_forest_from_past)/(crop_farea + past_farea - pot_forest_from_crop - pot_forest_from_past)*f_diff
             ! what if availability is mostly in crop or past?
          else
             reassign_crop = min(crop_farea - pot_forest_from_crop, GLM_nfarea - (crop_pos - pot_forest_to_crop) - (past_pos - pot_forest_to_past))
             reassign_past = min(past_farea - pot_forest_from_past, GLM_nfarea - (past_pos - pot_forest_to_past) - (crop_pos - pot_forest_to_crop) - reassign_crop)
             unmet_farea(ind) = f_diff - (reassign_crop + reassign_past)
          end if

          crop_pos_f = pot_forest_to_crop
          crop_pos_nf = crop_pos - pot_forest_to_crop + reassign_crop
          crop_neg_f = (pot_forest_from_crop + reassign_crop)
          crop_neg_nf = (crop_neg - pot_forest_from_crop)/(crop_nfarea+1e-12)
          past_pos_f = pot_forest_to_past
          past_pos_nf = past_pos - pot_forest_to_past + reassign_past
          past_neg_f = (pot_forest_from_past + reassign_past)
          past_neg_nf = (past_neg - pot_forest_from_past)/(past_nfarea+1e-12)

       end if

       ! apply the cropland and pature changes (both forested and non-forested)
       ! to individual gridcells in each region/AEZ

       ! crop_neg_nf and crop_neg_f
       where (glu_weights(g,:,:) > 0) 
          glm_crop(:,:,np1) = glm_crop(:,:,n)   -  glm_crop(:,:,n)*fnfnonforest*crop_neg_nf
       endwhere
       !jt             [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
       if (crop_neg_f>0) then 
          cumsum_sorted_farea=0.
          call cumsum(glm_crop(:,:,n)*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)
          sortind=MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>crop_neg_f)
          if (sortind(1)==0) then
             if ( abs(crop_neg_f-cumsum_sorted_farea(totrglus)) .lt. 1e-10) then 
                crop_neg_f=cumsum_sorted_farea(totrglus)
                sortind=MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>=crop_neg_f)
             else
                sortind(1) = totrglus
             end if
          end if
          if (sortind(1)>1) then
             sortsitesdn=0
             do i=1,sortind(1)-1
                sortsitesdn(v1d(i),v2d(i))=1
             end do
             where(sortsitesdn>0)
                glm_crop(:,:,np1)= glm_crop(:,:,np1) - glm_crop(:,:,n) * fnfforest
             end where
             final_area_needed =  crop_neg_f - cumsum_sorted_farea(sortind(1)-1)
          else
             final_area_needed =  crop_neg_f
          end if
          glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) = &
               glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) - &
               min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
               glm_crop(v1d(sortind(1)),v2d(sortind(1)),n)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
       end if

       ! past_neg_f and past_neg_nf

       where(glu_weights(g,:,:) > 0) 
          glm_past(:,:,np1) = glm_past(:,:,n) - glm_past(:,:,n)*fnfnonforest*past_neg_nf
       end where
       if (past_neg_f>0) then
          !jt             [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
          cumsum_sorted_farea=0.
          call cumsum(glm_past(:,:,n)*cellarea_forest(:,:),v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)
          sortind=MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>past_neg_f)
          if (sortind(1)==0) then
             if ( abs(past_neg_f-cumsum_sorted_farea(totrglus)) .lt. 1e-10) then 
                past_neg_f=cumsum_sorted_farea(totrglus)
                sortind=MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>=past_neg_f)
             else
                sortind(1) = totrglus
             end if
          end if
          if (sortind(1)>1) then
             sortsitesdn=0
             do i=1,sortind(1)-1
                sortsitesdn(v1d(i),v2d(i))=1
             end do
             where(sortsitesdn>0)
                glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,n) * fnfforest(:,:)
             end where
             final_area_needed =  past_neg_f - cumsum_sorted_farea(sortind(1)-1)
          else
             final_area_needed =  past_neg_f
          end if
          glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) = &
               glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) - &
               min(final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))), &
               glm_past(v1d(sortind(1)),v2d(sortind(1)),n)*fnfforest(v1d(sortind(1)),v2d(sortind(1))))
       end if

       ! crop_pos_nf
       avail_land0 = 0
       where ( glu_weights(g,:,:)>0. .and. glm_crop(:,:,np1)>0. )
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       endwhere
       sumavail_land0=sum(avail_land0)

       if ( sumavail_land0 >= crop_pos_nf .or. abs(crop_pos_nf)<=1e-6) then
          if (abs(crop_pos_nf)<=1e-6) then
             crop_pos_nf=0.
          end if
#ifdef DEBUG
          write(6,*)'cropland increase on non-forested land - land available'
#endif
          where(glu_weights(g,:,:) > 0) 
             glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_land0/(sumavail_land0+1e-12)*crop_pos_nf*fnfnonforest)/cellarea
          end where
       else
          avail_landA = 0.
          where ( glu_weights(g,:,:) > 0 )
             avail_landA=(pctland_in2015 - glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          end where
          sumavail_landA=sum(avail_landA)
          if (crop_pos_nf - sumavail_landA < 1e-6) then
#ifdef DEBUG
             write(6,*)'cropland increase - land available'
#endif 
             where(glu_weights(g,:,:) > 0) 
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0*fnfnonforest/cellarea
                glm_crop(:,:,np1) = glm_crop(:,:,np1) + &
                     ((avail_landA-avail_land0)/(sumavail_landA-sumavail_land0+1e-12)*(crop_pos_nf-sumavail_land0)*fnfnonforest)/cellarea
             end where
          else
             write(6,*)'crop increase on non-forest - land not available'
             call abort
          end if
       end if

       ! past_pos_nf
       avail_land0 = 0
       where ( glu_weights(g,:,:) > 0 .and.glm_past(:,:,np1)>0.)
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
       end where
       sumavail_land0=sum(avail_land0)
       if ( sumavail_land0 >= past_pos_nf .or. abs(past_pos_nf)<=1e-6) then
          if (abs(past_pos_nf)<=1e-6) then
             past_pos_nf=0
          end if
#ifdef DEBUG
          write(6,*)'pasture increase on non-forested land - land available'
#endif
          where(glu_weights(g,:,:) > 0) 
             glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0/(sumavail_land0+1e-12)*past_pos_nf*fnfnonforest)/cellarea
          end where
       else
          avail_landA = 0
          where ( glu_weights(g,:,:) > 0 )
             avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
          end where
          sumavail_landA=sum(avail_landA)
          if ( sumavail_landA >= past_pos_nf) then
#ifdef DEBUG
             write(6,*)'pasture increase on non-forested land - land available'
#endif
             where(glu_weights(g,:,:) > 0) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0*fnfnonforest/cellarea
                glm_past(:,:,np1) = glm_past(:,:,np1) + &
                     ((avail_landA-avail_land0)/(sumavail_landA -sumavail_land0+1e-12) * &
                     (past_pos_nf-sumavail_land0) * fnfnonforest)/cellarea
             end where
          else
#ifdef DEBUG
             write(6,*)'pasture increase on non-forest - land NOT available, reducing pasture increase to accomodate'
#endif
             past_pos_nf = sumavail_landA
             where(glu_weights(g,:,:) > 0) 
                glm_past(:,:,np1) = glm_past(:,:,np1) + (avail_land0*fnfnonforest/cellarea)
                glm_past(:,:,np1) = glm_past(:,:,np1) + ((avail_landA - avail_land0) / &
                     (sumavail_landA - sumavail_land0 + 1e-12)*(past_pos_nf-sumavail_land0)*fnfnonforest)/cellarea
             end where
          end if
       end if

       ! crop_pos_f
       avail_land0 = 0
       where ( glu_weights(g,:,:) > 0 .and.glm_crop(:,:,np1)>0.)
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
       end where
       sumavail_land0=sum(avail_land0)
       if (abs(crop_pos_f)>1e-6) then
          if (sumavail_land0>=crop_pos_f) then
#ifdef DEBUG
             write(6,*)'cropland increase on forested land - land available'
#endif
             !jt              [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
             cumsum_sorted_farea=0.
             call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
             sortind=MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>crop_pos_f)
             if (sortind(1)==0) then
                sortind(1) = totrglus
             end if
             if (sortind(1)>1) then
                sortsitesup=0
                do i=1,sortind(1)-1
                   sortsitesup(v1u(i),v2u(i))=1
                end do
                where(sortsitesup>0)
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0 * fnfforest / cellarea
                end where
                final_area_needed =  crop_pos_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  crop_pos_f
             end if
             glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                  glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                  min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                  avail_land0(v1u(sortind(1)),v2u(sortind(1))) * &
                  fnfforest(v1u(sortind(1)),v2u(sortind(1)))/cellarea(v1u(sortind(1)),v2u(sortind(1))))
          else
             avail_landA = 0.
             where ( glu_weights(g,:,:) > 0 )
                avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
             end where
             sumavail_landA=sum(avail_landA)
             if (sumavail_landA >=crop_pos_f) then
                !jt  [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0)
                   glm_crop(:,:,np1) = glm_crop(:,:,np1) + avail_land0*fnfforest/cellarea
                end where
#ifdef DEBUG
                write(6,*)'cropland increase on forested land - land available'
#endif                   
                !jt  [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                sortind = MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>(crop_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)))
                if (sortind(1)==0) then
                   sortind(1) = totrglus
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + (avail_landA - avail_land0) * fnfforest / cellarea
                   end where
                   final_area_needed =  crop_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  crop_pos_f - sumavail_land0
                end if
                glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_crop(v1u(sortind(1)),v2u(sortind(1)),np1) +&
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                     (avail_landA(v1u(sortind(1)),v2u(sortind(1)))) * &
                     fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                     cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
                write(6,*)'crop increase on forest - land not available'
                call abort
             end if
          end if

       end if

       ! past_pos_f
       avail_land0 = 0.
       where ( glu_weights(g,:,:) > 0 .and.glm_past(:,:,np1)>0.)
          avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
       end where
       sumavail_land0=sum(avail_land0)
       if (abs(past_pos_f)>1e-6) then 
          if (sumavail_land0>=past_pos_f) then
#ifdef DEBUG
             write(6,*)'pasture increase on forest - land available'
#endif                
             !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
             cumsum_sorted_farea=0.
             call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
             sortind = MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>past_pos_f)
             if (sortind(1)==0) then
                sortind(1) = totrglus
             end if
             if (sortind(1)>1) then 
                sortsitesup=0
                do i=1,sortind(1)-1
                   sortsitesup(v1u(i),v2u(i))=1
                end do
                where(sortsitesup>0)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 * fnfforest / cellarea
                end where
                final_area_needed =  past_pos_f - cumsum_sorted_farea(sortind(1)-1)
             else
                final_area_needed =  past_pos_f
             end if
             glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                  glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                  min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))) , &
                  avail_land0(v1u(sortind(1)),v2u(sortind(1))) * &
                  fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                  cellarea(v1u(sortind(1)),v2u(sortind(1))))
          else
             avail_landA = 0.
             where ( glu_weights(g,:,:) > 0 )
                avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest
             end where
             sumavail_landA=sum(avail_landA)
             if (sumavail_landA >= past_pos_f) then 
#ifdef DEBUG
                write(6,*)'pasture increase on forest - land available'
#endif
                !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')
                cumsum_sorted_farea=0.
                call cumsum(avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                where(glu_weights(g,:,:) > 0)
                   glm_past(:,:,np1) = glm_past(:,:,np1) + avail_land0 * fnfforest / cellarea
                end where

                !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'ascend')

                call cumsum(avail_landA(:,:)-avail_land0(:,:),v1u,v2u,cumsum_sorted_farea(:totrglus),totrglus)
                sortind = MINLOC(cumsum_sorted_farea(:totrglus),mask=cumsum_sorted_farea(:totrglus)>(past_pos_f-sum(avail_land0(:,:),mask=glu_weights(g,:,:)>0)))
                if (sortind(1)==0) then
                   sortind(1) = totrglus
                end if
                if (sortind(1)>1) then
                   sortsitesup=0
                   do i=1,sortind(1)-1
                      sortsitesup(v1u(i),v2u(i))=1
                   end do
                   where(sortsitesup>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1)+(avail_landA-avail_land0)*fnfforest/cellarea
                   end where
                   final_area_needed =  past_pos_f - sumavail_land0-cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  past_pos_f - sumavail_land0
                end if
                glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) = &
                     glm_past(v1u(sortind(1)),v2u(sortind(1)),np1) + &
                     min(final_area_needed/cellarea(v1u(sortind(1)),v2u(sortind(1))), &
                     (avail_landA(v1u(sortind(1)),v2u(sortind(1))) - &
                     avail_landA(v1u(sortind(1)),v2u(sortind(1)))) * &
                     fnfforest(v1u(sortind(1)),v2u(sortind(1))) / &
                     cellarea(v1u(sortind(1)),v2u(sortind(1))))
             else
                write(6,*)'pasture increase on forest - land not available'
                call abort
             end if
          end if

       end if
    end do ! end glu loop

!       end do ! end aez loop
!    end do ! end nreg loop

    if (regional_unmet_reassign==1) then
       ! Here we are looping over regions, using the rgmin(r) and
       ! rgmax(r) to find the g's inside this region
       do r = 1,nreg

          !regional_farea_needed = sum(unmet_farea(((r-1)*18+1):r*18))
          regional_farea_needed = sum(unmet_farea(rgmin(r):rgmax(r)))
          
          ! number of glus in this region
          ! nglu = rgmax(r)-rgmin(r)-1
          ! TRS - zounds! I'm pretty sure -1 is wrong, just from arithmetic
          nglu = rgmax(r)-rgmin(r)+1

          ! Just in case
          avail_farea = 0.
          avail_nfarea = 0.
          avail_ag_farea = 0.
          reassign_ag = 0.
          unmet_aez_farea = 0.

          do g = rgmin(r),rgmax(r)
             ! g1 =  1:nglu
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

             ! May need to zero-init all these each loop, becasue we
             ! now have variable nglu per region.
             avail_farea(g1) = sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_forest)
             avail_nfarea(g1) = sum((pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest)
             avail_ag_farea(g1) = sum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest)
             reassign_ag(g1) = min(avail_ag_farea(g1), avail_nfarea(g1), regional_farea_needed)
             unmet_aez_farea(g1) = regional_farea_needed - reassign_ag(g1)
          end do
          
          !jt [sorted_reassign_ag,sort_aez] = sort(reassign_ag,'descend')
          ! TRS - not sure this will work with variable nglu
          indxa=(/(i,i=1,nglu)/)
          indxadn=(/(i,i=1,nglu)/)

          call D_mrgrnk(reassign_ag,indxa,nglu)
          call D_mrgrnk(reassign_ag*-1.,indxadn,nglu)

          cumsum_sorted_reassign_ag(1)=reassign_ag(indxadn(1))
          do i=2,nglu
             cumsum_sorted_reassign_ag(i)=cumsum_sorted_reassign_ag(i-1)+reassign_ag(indxadn(i))
          end do
          sortind = MINLOC(cumsum_sorted_reassign_ag(:),mask=cumsum_sorted_reassign_ag(:) >= regional_farea_needed)
          if (sortind(1).eq.0) then
             sortind(1) = nglu
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
          ! and then reassign forest in them.  I'm not sure what is
          ! going to happen when we only have one glu, or one
          ! important glu, but whatever, let's go with it.
          do zz=1,sortind(1)
              z=indxadn(zz)
              ! Remember, z = 1,nglu(r).  So convert to a global g.
              g=z+rgmin(r)-1

              if (reassign_ag(z)>0) then
                 crop_before_decrease = sum(glm_crop(:,:,np1)*cellarea)
                 past_before_decrease = sum(glm_past(:,:,np1)*cellarea)
                rglus=0
                where ( glu_weights(g,:,:) .gt. 0)
                   rglus=1
                end where
                totrglus=sum(rglus)
                
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
                call D_mrgrnk(pot_veg_rev*glu_weights_rev, indxup,numlons*numlats)
                call D_mrgrnk(pot_veg_rev*glu_weights_rev*-1., indxdn,numlons*numlats)
                
                !jt  The sortxxxup and sortxxxdn arrays are only good 1:totrglus these arrays are also based on row major
                !jt  sorting to match original matlab scripts.  We can get rid of this after validation.
                
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
                v1u=sortlonsup(numlons*numlats-totrglus+1:numlons*numlats)
                v2u=sortlatsup(numlons*numlats-totrglus+1:numlons*numlats)
                v1d=(sortlonsdn(:totrglus))
                v2d=(sortlatsdn(:totrglus))
                
                !jt [sorted_pot_veg,sort_ind] = sort(pot_veg(rAEZ_sites),'descend')
                cumsum_sorted_farea=0.
                call cumsum((glm_crop(:,:,np1)+glm_past(:,:,np1))*cellarea_forest,v1d,v2d,cumsum_sorted_farea(:totrglus),totrglus)
                sortind = MINLOC(cumsum_sorted_farea(:totrglus),cumsum_sorted_farea(:totrglus) > reassign_ag(z))
                if (sortind(1)==0) then
                      sortind(1) = totrglus
                end if
                if (sortind(1)>1) then
                   sortsitesdn=0
                   do i=1,sortind(1)-1
                      sortsitesdn(v1d(i),v2d(i))=1
                   end do
                   where(sortsitesdn>0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) - glm_crop(:,:,np1) * fnfforest
                   end where

                   where(sortsitesdn>0)
                      glm_past(:,:,np1) = glm_past(:,:,np1) - glm_past(:,:,np1) * fnfforest
                   end where

                   final_area_needed =  reassign_ag(z) - cumsum_sorted_farea(sortind(1)-1)
                else
                   final_area_needed =  reassign_ag(z)
                end if
                ag_area_avail = (glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) + glm_past(v1d(sortind(1)),v2d(sortind(1)),np1))*cellarea(v1d(sortind(1)),v2d(sortind(1)))
                crop_area_avail = glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1)*cellarea(v1d(sortind(1)),v2d(sortind(1)))
                glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                     glm_crop(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1))), &
                     crop_area_avail/cellarea(v1d(sortind(1)),v2d(sortind(1))))

                glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) = &
                     glm_past(v1d(sortind(1)),v2d(sortind(1)),np1) - &
                     (final_area_needed/cellarea(v1d(sortind(1)),v2d(sortind(1))) - &
                     min(crop_area_avail/(ag_area_avail+1e-12)*final_area_needed / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1))),crop_area_avail / &
                     cellarea(v1d(sortind(1)),v2d(sortind(1)))))
                
                crop_after_decrease = sum(glm_crop(:,:,np1)*cellarea)
                past_after_decrease = sum(glm_past(:,:,np1)*cellarea)
                total_ag_decrease = crop_before_decrease+past_before_decrease - crop_after_decrease - past_after_decrease
                crop_decrease_ratio = (crop_before_decrease - crop_after_decrease)/(total_ag_decrease+1e-12)
                past_decrease_ratio = 1-crop_decrease_ratio
                
                avail_land0 = 0.
                
                where ( glu_weights(g,:,:) > 0 .and.(glm_crop(:,:,np1)>0..or.glm_past(:,:,np1)>0.))
                   avail_land0=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                end where
                sumavail_land0=sum(avail_land0)
                if (sumavail_land0 >=reassign_ag(z)) then 
#ifdef DEBUG
                   write(6,*)'crop and pasture increase on nonforest - land available'
#endif
                   where(glu_weights(g,:,:) > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z)*fnfnonforest/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*avail_land0 / &
                           (sumavail_land0+1e-12) * reassign_ag(z)*fnfnonforest/cellarea
                   end where
                else
                   avail_landA = 0.
                   where ( glu_weights(g,:,:) > 0 )
                      avail_landA=(pctland_in2015-glm_crop(:,:,np1)-glm_past(:,:,np1))*cellarea_nonforest
                   end where
                   sumavail_landA=sum(avail_landA)
                   if (reassign_ag(z) > sumavail_landA) then 
                      reassign_ag(z) = sum(avail_landA)
                      call abort
                   end if
#ifdef DEBUG
                   write(6,*)'cropland and pasture increase on non-forest - land available'
#endif
                   where(glu_weights(g,:,:) > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*(avail_land0*fnfnonforest)/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*(avail_land0*fnfnonforest)/cellarea
                   end where
                      
                   where(glu_weights(g,:,:) > 0)
                      glm_crop(:,:,np1) = glm_crop(:,:,np1) + crop_decrease_ratio*((avail_landA-avail_land0) / &
                           (sumavail_landA-sumavail_land0+1e-12)*(reassign_ag(z)-sumavail_land0)*fnfnonforest)/cellarea
                      glm_past(:,:,np1) = glm_past(:,:,np1) + past_decrease_ratio*((avail_landA-avail_land0)/ &
                           (sumavail_landA-sumavail_land0+1e-12)*(reassign_ag(z)-sumavail_land0)*fnfnonforest)/cellarea
                   end where
                end if
             end if
          end do ! end z loop
       end do !  end r loop
    end if   ! ! if regional_unmet_reassign
 end if

 if (allocated(v1u)) deallocate(v1u)
 if (allocated(v2u)) deallocate(v2u)
 if (allocated(v1d)) deallocate(v1d)
 if (allocated(v2d)) deallocate(v2d)

#ifdef DEBUG
    write(6,*) 'sum final gcrop/gpast n'
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,n))
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,n))
    write(6,*) 'sum final gcrop/gpast np1'
    write(6,fmt="(1ES25.15)") sum(glm_crop(:,:,np1))
    write(6,fmt="(1ES25.15)") sum(glm_past(:,:,np1))
#endif 
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
     write(iulog,*)'crop interpolation factors fact1,fact2,year1,year2=',fact1,fact2,year1,year2
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
! this should always be true if gcam is at annual time step
    if (eclockyr==(year2-1)) then 
       glm_crop(:,:,n)=glm_crop(:,:,np1)
       glm_past(:,:,n)=glm_past(:,:,np1)
       gcam_crop(:,n)=gcam_crop(:,np1)
       gcam_past(:,n)=gcam_past(:,np1)
       gcam_wh(:,n)=gcam_wh(:,np1)
       gcam_forest_area(:,n)=gcam_forest_area(:,np1)
    end if

    ! lets write a restart every time this routine is called

    call shr_cal_date2ymd(ymd,year,mon,day)
    write(filename,'(a,i4.4,a,i2.2,a)') trim(gcam2glm_restfile)//'r.',year,'.nc'

    iun = shr_file_getunit()
    open(iun,file=trim(gcam2glm_rpointer),form='formatted')
    write(iun,'(a)') trim(filename)
    close(iun)
    call shr_file_freeunit(iun)

    write(iulog,*) subname,' write_restart rpointer ',trim(gcam2glm_rpointer)
    write(iulog,*) subname,' write_restart file     ',trim(filename)

    status= nf90_create(filename,nf90_clobber,ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    status = nf90_def_dim(ncid,'lon',numlons,dimid3(1))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'lat',numlats,dimid3(2))
    if(status /= nf90_NoErr) call handle_err(status)
    status = nf90_def_dim(ncid,'time',2,dimid3(3))
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

