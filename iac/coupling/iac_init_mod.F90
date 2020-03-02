module iac_init_mod
  !------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module to contain initialization routine(s?) for gcam/iac
  !------------------------------------------------------------
  use shr_sys_mod , only : shr_sys_abort
  use shr_file_mod, only : shr_file_get, shr_file_getUnit, shr_file_freeUnit
  use iac_spmd_mod, only : masterproc
  !use iac_io_mod
  use gcam_var_mod, only : iulog
  use iac_data_mod
  use iac_io_mod, only : ncd_pio_openfile, ncd_pio_closefile, ncd_io
  use iac_io_mod, only : file_desc_t
  use iac_io_mod

  implicit none
  save

contains
  subroutine iac_init()
    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Read in namelist, define grid, allocate variables
    !---------------------------------------------------------------------------
    use shr_kind_mod, only: CL => SHR_KIND_CL, CXX => SHR_KIND_CXX

    character(len=CXX) :: gcam_gridfile
    character(len=CL) :: nlfilename_iac
    integer :: iac_npft, unitn, ier, dimid, i, j, n
    type(file_desc_t) :: ncid
    logical :: lexist, found
    real(r8), pointer :: tempr(:,:)
    integer, pointer :: itempr(:,:)
    character(len=32) :: subname = 'iac_init'

    ! GCAM namelist
    namelist /gcam_inparm/ gcam_gridfile, iac_npft

    nlfilename_iac = "gcam_in"

    inquire (file = trim(nlfilename_iac), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_iac does NOT exist:'&
            //trim(nlfilename_iac)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_iac does not exist')
    end if

    if (masterproc) then
       !unitn = sr_file_getunit()
       write(iulog,*) 'Read in gcam_inparm namelist from: ', trim(nlfilename_iac)
       open( unitn, file=trim(nlfilename_iac), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, gcam_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on gcam_inparm read' )
          endif
       end do
       call shr_file_freeUnit( unitn )
    end if

    ! Don't need to mpi_bcast the namelist, because we don't have other procs
    ! to communitate with.  But if we did, we'd mpi_bcast them around
    ! mpicom_iac here.  

    if (masterproc) then
       write(iulog,*) 'define GCAM run:'
       ! write out namelist stuff here

       !if (nsrest == nsrStartup .and. finidat_rtm /= ' ') then
       !   write(iulog,*) '   MOSART initial data   = ',trim(finidat_rtm)
       !end if
    endif

    ! Initialize some iac_data stuff
    allocate(GClock(iac_eclock_size))
    allocate(gdata%c(iac_cdata_size))
    allocate(gdata%r(iac_cdata_size))
    allocate(gdata%i(iac_cdata_size))
    allocate(gdata%l(iac_cdata_size))

    ! Read grid file
    call ncd_pio_init()

    inquire(file=gcam_gridfile,exist=lexist)
    if (.NOT. lexist) then 
       call shr_sys_abort ('IAC INIT: FAILED to find '//trim(gcam_gridfile))
    endif

    ! Read in dimension lengths, so we know how to allocate stuff
    call ncd_pio_openfile (ncid, trim(gcam_gridfile), 0)
    call ncd_inqdid(ncid,'lon',dimid)
    call ncd_inqdlen(ncid,dimid,iac_ctl%nlon)
    call ncd_inqdid(ncid,'lat',dimid)
    call ncd_inqdlen(ncid,dimid,iac_ctl%nlat)

    iac_ctl%ngrid=iac_ctl%nlon*iac_ctl%nlat

    ! No decomposition
    iac_ctl%begg = 1
    iac_ctl%endg = iac_ctl%ngrid

    ! From the namelist/config 
    iac_ctl%npft = iac_npft

    ! Start allocating vars and control vectors
    allocate(iac_ctl%ilon(iac_ctl%nlon))
    allocate(iac_ctl%jlat(iac_ctl%nlat))

    allocate(iac_ctl%gindex(iac_ctl%ngrid))
    allocate(iac_ctl%iacmask(iac_ctl%ngrid))

    ! Dimension order: (lon,lat,pft)
    allocate(lnd2iac_vars%npp(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(lnd2iac_vars%hr(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(lnd2iac_vars%pftwgt(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))

    allocate(lnd2iac_vars%landfrac(iac_ctl%nlon,iac_ctl%nlat))
    allocate(lnd2iac_vars%area(iac_ctl%nlon,iac_ctl%nlat))

    allocate(iac2lnd_vars%pct_pft(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))

    allocate(iac2atm_vars%co2emiss(iac_ctl%nlon,iac_ctl%nlat))

    ! Assign our global index, ilon, and jlat, to go from (g)
    ! dimension to (lon,lat) dimensions.
    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          n = (j-1)*iac_ctl%nlat + i
          iac_ctl%gindex(n)=n  ! Dumb, but allows decomp someday, maybe
          iac_ctl%ilon(n)=i
          iac_ctl%jlat(n)=j
       end do
    end do

    ! Now we start reading our grid variables.  I'm assuming all the
    ! lat,lon,area, and landfrac variables are in this netcdf
    ! configuration file - extracted (or used straight up?) from a
    ! clm2.h1 file?

    ! Apparently this kind of netcdf interface is all the rage
    ! I'm going to copy existing interfaces (rtm, specifically) in
    ! these calls, because there's probably better error checking
    ! than if I rewrote them for the simple tasks of doing our reads. 
    ! But it *is* annoying not to have any simple examples of netcdf
    ! interaction 

    allocate(tempr(iac_ctl%nlon,iac_ctl%nlat))  
    allocate(itempr(iac_ctl%nlon,iac_ctl%nlat))  

    ! lon
    call ncd_io(ncid=ncid, varname='lon', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read GCAM longitudes')
    if (masterproc) write(iulog,*) 'Read lon ',minval(tempr),maxval(tempr)
    do i=1,iac_ctl%nlon
       iac_ctl%lon(i) = tempr(i,1)
    enddo

    ! lat
    call ncd_io(ncid=ncid, varname='lat', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read GCAM latitudes')
    if (masterproc) write(iulog,*) 'Read lat ',minval(tempr),maxval(tempr)
    do j=1,iac_ctl%nlat
       iac_ctl%lat(j) = tempr(1,j)
    enddo

    ! area
    call ncd_io(ncid=ncid, varname='area', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read GCAM area')
    if (masterproc) write(iulog,*) 'Read area ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          ! area needed by lnd2iac_var, input to gcam
          lnd2iac_vars%area(i,j) = tempr(i,j)
       enddo
    enddo

    ! landfrac
    call ncd_io(ncid=ncid, varname='landfrac', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read GCAM landfrac')
    if (masterproc) write(iulog,*) 'Read landfrac ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          ! area needed by lnd2iac_var, input to gcam
          lnd2iac_vars%landfrac(i,j) = tempr(i,j)
       enddo
    enddo

    ! landmask - same as iacmask?  Needed for gsl map
    call ncd_io(ncid=ncid, varname='landmask', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read GCAM landmask')
    if (masterproc) write(iulog,*) 'Read landmask ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          ! For global seg mapping
          n = (j-1)*iac_ctl%nlon + i
          iac_ctl%iacmask(n) = itempr(i,j)
       enddo
    enddo

    deallocate(tempr)
    deallocate(itempr)             

    call ncd_pio_closefile(ncid)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize Restart
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize history handler
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine iac_init
end module iac_init_mod
