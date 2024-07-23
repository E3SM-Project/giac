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
    nsize = iac_glm_nx*iac_glm_ny 
    ! npfts + extra pft + vh1,vh2,sh1,sh2,sh3,grazing
    allocate(plodata(nflds,iac_glm_nx*iac_glm_ny))
    plodata(:,:) = 0.0

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
    integer :: ymd, tod, dt, indnew, indprev
    integer :: i,j,n,ij,ierr,nmode
    integer :: dimid(3),varid,ncid, lat_id, lon_id, rawdims(2), rid
    real*8, pointer :: array3d(:,:,:)
    integer, dimension(166) :: array1d
    character(len=128) :: fname,casename,hfile
    integer :: myear, mon, day, e3smyear
    character(len=*),parameter :: subname='(glm2iac_run_mod)'
    integer, dimension(3) :: start3, count3

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

    ! need this for both cases to reduce looping
    allocate(array3d(iac_glm_nx,iac_glm_ny,size(glmo,dim=1)))

    ! if running gcam (default for the future), flip the glmo data 
    !   otherwise read in the historical data
    ! no need for glm history file when reading historical data
    if ( run_gcam ) then

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

       ierr = nf90_def_dim(ncid,'glmo_nf_raw', size(glmo,dim=1), rawdims(1))
       ierr = nf90_def_dim(ncid,'glmo_nxy', size(glmo,dim=2), rawdims(2))
       ierr = nf90_def_var(ncid,'raw_glmo',NF90_DOUBLE,rawdims,rid)
       ierr = nf90_enddef(ncid)

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

       ierr = nf90_put_var(ncid,rid,glmo)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

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
    
       ierr = nf90_close(ncid)
    else
       ! read in the historical data
       ! the data order is the same as expected by the LUT below
       ! hardcode these file names for now

       ! land type data is first
       ierr = nf90_open('LUH2_HIST_LUH1f_c07182019.nc',&
          nf90_nowrite,ncid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ! first find the correct time index
       ! these data can be used for harvest below
       !    because the harvest time array matches but
       !    is one element shorter at the end
       ! the desired year for land cover is myear
       ! the desired year for harvest is e3smyear, which is myear-1
       ierr= nf90_inq_varid(ncid, "TIME", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array1d(:))
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       start3(1) = 1
       start3(2) = 1
       start3(3) = findloc(array1d, myear, dim=1)
       count3(1) = iac_glm_nx
       count3(2) = iac_glm_ny
       count3(3) = 1      
 
       ! need to read into a multi dim array then transfer
       
       ! order is crop, pasture, primary, secondary
       ierr= nf90_inq_varid(ncid, "GCROP", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,1), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GPAST", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,2), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GOTHR", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,3), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GSECD", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,4), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_close(ncid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ! wood harvest data is second
       ! order is primary forest, primary non-forest,
       !   secondary mature forest, secondary young forest, secondary non-forest
       ierr = nf90_open('LUH2_HIST_LUH1f_c07182019_harvest.nc',&
          nf90_nowrite,ncid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ! need the previous year data for harvest
       start3(3) = findloc(array1d, e3smyear, dim=1)

       ierr= nf90_inq_varid(ncid, "GFVH1", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,5), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GFVH2", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,6), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GFSH1", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,7), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GFSH2", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,8), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_inq_varid(ncid, "GFSH3", varid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)
       ierr= nf90_get_var(ncid, varid, array3d(:,:,9), start3, count3)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ierr= nf90_close(ncid)
       if(ierr /= nf90_NoErr) call handle_err(ierr)

       ! now loop over the cells to transfer the data because glmo is a pointer
       do n = 1,size(glmo,dim=1)
          ij = 0
          do j = 1,iac_glm_ny
             do i = 1,iac_glm_nx
                ij = ij + 1
                glmo(n,ij) = array3d(i,j,n)
             enddo
          enddo
       enddo

    end if

    deallocate(array3d)

#ifdef DEBUG
    do j = 1,iac_glmo_nflds
       write(iulog,*) trim(subname),' glmo minmax after flip ',&
                       j,minval(glmo(j,:)),maxval(glmo(j,:))
    enddo
#endif

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
    deallocate(lat_rev)

    ! move the stored pft percents to previous before filling with
    !    newly calculated values for myear
    ! recall that myear is model year + 1
    !    and that pfts are for start of myear and harvest is for model year
    ! mksurfdat fills iac2lnd_vars%pct_pft and iac2lnd_vars%harvest_frac

    do n=1,iac_ctl%npft
       do j=1,iac_ctl%nlat
          do i=1,iac_ctl%nlon
             iac2lnd_vars%pct_pft_prev(i,j,n) = iac2lnd_vars%pct_pft(i,j,n)
          enddo
       enddo
    enddo

    call mksurfdat_run(myear,plodata)


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

