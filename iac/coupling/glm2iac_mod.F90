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

  use iac_data_mod, only : cdata => gdata, EClock => GClock
  use gcam_var_mod
  use shr_cal_mod
  use netcdf
  use gcam2glm_mod, only : lon, lat, numLats

  implicit none
  SAVE
  private                              ! By default make data private

! !PUBLIC MEMBER FUNCTIONS:

  public :: glm2iac_init_mod               ! clm initialization
  public :: glm2iac_run_mod                ! clm run phase
  public :: glm2iac_final_mod              ! clm finalization/cleanup

! !PUBLIC DATA MEMBERS: None


! !REVISION HISTORY:
! Author: T Craig

  real*8, pointer, save :: plodata(:,:)
  real, pointer, save :: plodataf(:,:)
  real, pointer, save :: glmof(:,:)

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
    allocate(plodataf(nflds,nsize))
    allocate(glmof(size(glmo,dim=1),size(glmo,dim=2)))
    plodata = 0.0

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
    implicit none

! !ARGUMENTS:
    real*8, pointer :: glmo(:,:)

! !LOCAL VARIABLES:
    logical :: restart_now
    integer :: ymd, tod, dt
    integer :: i,j,n,ij,ierr,nmode
    integer :: dimid(3),varid,ncid, lat_id, lon_id, rawdims(2), rid
    real*8, pointer :: array3d(:,:,:)
    character(len=128) :: fname,casename,hfile
    integer :: myear, mon, day, e3smyear
    character(len=*),parameter :: subname='(glm2iac_run_mod)'

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

!this isn't writing for some reason
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

! now using the doule precision code

!    plodataf=plodata
!    glmof=glmo
!    call updateannuallanduse(glmof,plodataf,myear)
    call updateannuallanduse(glmo,plodata,myear)
!    plodata=plodataf
!    glmo=glmof

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
!    call mksurfdata(plodata)

deallocate(lat_rev)

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
    deallocate(plodataf)
    deallocate(glmof)

  end subroutine glm2iac_final_mod

!====================================================================================

end module glm2iac_mod

