#define DEBUG
Module glm2iac_mod
  
!---------------------------------------------------------------------------
!BOP
!
! !MODULE: glm2iac_mod
!
!  Interface of the integrated assessment component in CCSM
!
! !DESCRIPTION:
!
! !USES:

  use iac_data_mod, only : cdata => gdata, EClock => GClock, iac2lnd_vars, &
iac_ctl
  use gcam_var_mod
  use shr_cal_mod
  use netcdf
  use gcam2glm_mod, only : lon, lat, numLats, handle_err
  use mksurfdat, only : fdyndat
  use abortutils, only : endrun
  use shr_log_mod, only : errMsg => shr_log_errMsg

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: glm2iac_init_mod               ! clm initialization
  public :: glm2iac_run_mod                ! clm run phase
  public :: glm2iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS:

  character(len=64), allocatable :: harvest_names(:)

! !REVISION HISTORY:
! Author: T Craig

  real*8, pointer, save :: plodata(:,:)

! !PRIVATE DATA MEMBERS:

!EOP
!===============================================================
contains
!===============================================================

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_init_mod

! !INTERFACE:
  subroutine glm2iac_init_mod( glmo )

! !DESCRIPTION:
! Initialize interface for glm

! !USES:
    use iac_data_mod
    implicit none

! !ARGUMENTS:
    real*8, pointer :: glmo(:,:)

! !LOCAL VARIABLES:

    integer :: nflds, nsize
    character(len=*),parameter :: subname='(glm2iac_init_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    nflds = iac_iac_npfts + 7
    nsize = cdata%i(iac_cdatai_glm_size)
    ! npfts + extra pft + vh1,vh2,sh1,sh2,sh3,grazing
    allocate(plodata(nflds,nsize))
    plodata = 0.0

    allocate(harvest_names(iac_ctl%nharvest))
    harvest_names = ['HARVEST_VH1', 'HARVEST_VH2', 'HARVEST_SH1', &
                     'HARVEST_SH2', 'HARVEST_SH3']

#ifdef DEBUG
    write(iulog,*) trim(subname),' allocate plodata ',nflds,nsize
#endif
  end subroutine glm2iac_init_mod

!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_run_mod

! !INTERFACE:
  subroutine glm2iac_run_mod( glmo )

! !DESCRIPTION:
! Run interface for glm

! !USES:
    use iac_data_mod
    use mksurfdat, only: mksurfdat_run
    use mct_mod, only : mct_die
    implicit none

! !ARGUMENTS:
    real*8, pointer :: glmo(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt, indnew, indprev
    integer :: i,j,n,ij,ierr,nmode
    integer :: dimid(3),varid,ncid, lat_id, lon_id, rawdims(2), rid
    real*8, pointer :: array3d(:,:,:)
    character(len=128) :: fname,casename,hfile
    integer :: myear, mon, day, e3smyear
    character(len=*),parameter :: subname='(glm2iac_run_mod)'
    ! these are extracted from the dynamic land file
    integer :: nlon, nlat, ntime, npft
    ! these are for reading the dynamic land file
    integer, dimension(3) :: start3, count3
    integer, dimension(4) :: start4, count4
    integer, allocatable :: lsf_years(:)

    real*8, allocatable :: lat_rev(:)

! !REVISION HISTORY:
! Author: T Craig

!EOP
!-----------------------------------------------------------------------

    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)

#ifdef DEBUG
    do j = 1,iac_glmo_nflds
       write(iulog,*) trim(subname),' glmo minmax ',j,minval(glmo(j,:)),maxval(glmo(j,:))
    enddo
#endif

    allocate(lat_rev(numLats))
    lat_rev = lat(numLats:1:-1)

    call shr_cal_date2ymd(ymd,e3smyear,mon,day)
    
    ! the iac component runs ahead of E3SM by one year since it is defining                          
    ! boundary conditions for the other E3SM components. Increment the clock by a year               
    myear=e3smyear+1
    write(iulog,*) 'e3smyear is ',e3smyear,' iac year is ',myear

    casename = trim(case_name)
    write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a)') 'iac.hglmo.',myear,'-',mon,'-',day,'.nc'

#ifdef DEBUG
    write(iulog,*) trim(subname),' writing history file ',trim(hfile)
#endif

    nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
    ierr = nf90_create(trim(hfile),nmode,ncid)
    ierr = nf90_def_dim(ncid,'glmo_nx' ,iac_glm_nx,dimid(1))
    ierr = nf90_def_dim(ncid,'glmo_ny' ,iac_glm_ny,dimid(2))
    ierr = nf90_def_dim(ncid,'glmo_nf' ,size(glmo,dim=1),dimid(3))
    ierr = nf90_def_var(ncid,'glmodata',NF90_DOUBLE,dimid,varid)

    ierr = nf90_def_var(ncid,'lat',NF90_DOUBLE,dimid(2),lat_id)
    ierr = nf90_def_var(ncid,'lon',NF90_DOUBLE,dimid(1),lon_id)
    rawdims(1)=size(glmo,dim=1)
    rawdims(2) = size(glmo,dim=2)
    ierr = nf90_def_var(ncid,'raw_glmo',NF90_DOUBLE,rawdims,rid)

    ierr = nf90_enddef(ncid)
    allocate(array3d(iac_glm_nx,iac_glm_ny,size(glmo,dim=1)))

    ! flip the latitude into this diagnostic array and also in glmo
    ! GLM's origin is +90,-180, while the LUT origin is -90, -180

    do n = 1,size(glmo,dim=1)
    ij = 0
    do j = iac_glm_ny,1,-1
    do i = 1,iac_glm_nx
       ij = ij + 1
       array3d(i,j,n) = glmo(n,ij)
    enddo
    enddo
    enddo
    ierr = nf90_put_var(ncid,varid,array3d)

    ierr = nf90_put_var(ncid,lat_id,lat_rev)
    ierr = nf90_put_var(ncid,lon_id,lon) 

    ! avd - this isn't writing for some reason
    ierr = nf90_put_var(ncid,rid,glmo)

    ! now refill glmo
    do n = 1,size(glmo,dim=1)
       ij = 0
       do j = 1,iac_glm_ny
          do i = 1,iac_glm_nx
             ij = ij + 1
             glmo(n,ij) = array3d(i,j,n)
          enddo
       enddo
    enddo
    
    deallocate(array3d)
    ierr = nf90_close(ncid)

#ifdef DEBUG
    write(iulog,*) trim(subname),' date= ',ymd,tod
    write(iulog,*) trim(subname),' myear = ',myear
#endif

    write(iulog,*) trim(subname),' running LUT  '

    ! now using the double precision code
    call updateannuallanduse(glmo,plodata,myear)

#ifdef DEBUG
    do j = 1,size(plodata,dim=1)
       write(iulog,*) trim(subname),' plodata minmax ',j,minval(plodata(j,:)),maxval(plodata(j,:))
    enddo
#endif

    casename = trim(case_name)
    write(hfile,'(a,i4.4,a,i2.2,a,i2.2,a)') 'iac.hplo.',myear,'-',mon,'-',day,'.nc'
#ifdef DEBUG
    write(iulog,*) trim(subname),' writing history file ',trim(hfile)
#endif
    nmode = ior(NF90_CLOBBER,NF90_64BIT_OFFSET)
    ierr = nf90_create(trim(hfile),nmode,ncid)
    ierr = nf90_def_dim(ncid,'plo_nx' ,iac_glm_nx,dimid(1))
    ierr = nf90_def_dim(ncid,'plo_ny' ,iac_glm_ny,dimid(2))
    ierr = nf90_def_dim(ncid,'plo_nf' ,size(plodata,dim=1),dimid(3))
    ierr = nf90_def_var(ncid,'plodata',NF90_DOUBLE,dimid,varid)

    ierr = nf90_def_var(ncid,'lat',NF90_DOUBLE,dimid(2),lat_id)
    ierr = nf90_def_var(ncid,'lon',NF90_DOUBLE,dimid(1),lon_id)

    ierr = nf90_enddef(ncid)
    allocate(array3d(iac_glm_nx,iac_glm_ny,size(plodata,dim=1)))
    do n = 1,size(plodata,dim=1)
       ij = 0
       do j = 1,iac_glm_ny
          do i = 1,iac_glm_nx
             ij = ij + 1
             array3d(i,j,n) = plodata(n,ij)
          enddo
       enddo
    enddo
    ierr = nf90_put_var(ncid,varid,array3d)
    deallocate(array3d)

    ierr = nf90_put_var(ncid,lat_id,lat_rev)
    ierr = nf90_put_var(ncid,lon_id,lon)

    ierr = nf90_close(ncid)

    call mksurfdat_run(myear,plodata)

    deallocate(lat_rev)

    ! read in the dynamic land surface file fdyndat
    ! this was just updated by the mksurfdat_run call
    ! and fill the iac2lnd_vars
    ! myear is already advanced to model year +1
    ! myear in the file contains myear start pfts and myear-1 harvest

    ierr = nf90_open(fdyndat,nf90_nowrite,ncid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! first get the dimensions
    ierr= nf90_inq_dimid(ncid, "lsmlon", dimid(1))
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid, dimid(1), len=nlon)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid, "lsmlat", dimid(1))
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid, dimid(1), len=nlat)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid, "time", dimid(1))
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid, dimid(1), len=ntime)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid, "natpft", dimid(1))
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid, dimid(1), len=npft)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
   
    ! get the years in the dynamic land surface file
    ierr= nf90_inq_varid(ncid, "YEAR", varid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    allocate(lsf_years(ntime), stat=ierr)
    if(ierr/=0) call endrun(msg='ERROR reading dyn land file years'// &
                             errMsg(__FILE__, __LINE__))
    ierr= nf90_get_var(ncid, varid, lsf_years)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! get the indices for the new year and the previous year
    ! these should be ntime and ntime-1, respectively
    ! but find them in case a year is duplicated
    ! find the last matching element

    indnew = findloc(lsf_years, myear, dim=1, back=.true.)
    indprev = findloc(lsf_years, myear - 1, dim=1, back=.true.)

    ! pct pft
    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) =  indnew
    count4(1) = nlon
    count4(2) = nlat
    count4(3) = npft
    count4(4) = 1
    ierr= nf90_inq_varid(ncid, "PCT_NAT_PFT", varid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_get_var(ncid, varid, iac2lnd_vars%pct_pft, start=start4, count=count4)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! prev pct pft
    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = indprev
    count4(1) = nlon
    count4(2) = nlat
    count4(3) = npft
    count4(4) = 1
    ierr= nf90_inq_varid(ncid, "PCT_NAT_PFT", varid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_get_var(ncid, varid, iac2lnd_vars%pct_pft_prev, start=start4, count=count4)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! harvest frac
    start3(1) = 1
    start3(2) = 1
    start3(3) = indnew
    count3(1) = nlon
    count3(2) = nlat
    count3(3) = 1
    do n = 1,iac_ctl%nharvest
       ierr= nf90_inq_varid(ncid, harvest_names(n), varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, iac2lnd_vars%harvest_frac(:,:,n), start=start3, count=count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
    enddo  

    ierr= nf90_close(ncid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

  end subroutine glm2iac_run_mod


!---------------------------------------------------------------------------
!BOP

! !IROUTINE: glm2iac_final_mod

! !INTERFACE:
  subroutine glm2iac_final_mod( )

! !DESCRIPTION:
! Finalize glm model
! !USES:
    implicit none

! !ARGUMENTS:

! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='(glm2iac_final_mod)'

! !REVISION HISTORY:
! Author: T Craig

!EOP

!---------------------------------------------------------------------------

!    iu  = cdata%i(iac_cdatai_logunit)
!    write(iu,*) trim(subname)
    deallocate(plodata)

  end subroutine glm2iac_final_mod

!====================================================================================

end module glm2iac_mod

