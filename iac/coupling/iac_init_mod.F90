module iac_init_mod
  !------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module to contain initialization routine(s?) for gcam/iac
  !------------------------------------------------------------
  use shr_sys_mod , only : shr_sys_abort
  use shr_file_mod, only : shr_file_get, shr_file_getUnit, shr_file_freeUnit
  use iac_spmd_mod, only : masterproc
  !use iac_io_mod
  !use gcam_var_mod, only : iulog
  use gcam_var_mod
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

    character(len=CL) :: nlfilename_iac
    integer :: unitn, ier, dimid, i, j, n
    type(file_desc_t) :: ncid
    logical :: lexist, found
    real(r8), pointer :: tempr(:,:)
    integer, pointer :: itempr(:,:)
    character(len=32) :: subname = 'iac_init'

    ! GCAM namelist
    namelist /gcam_inparm/ &
         gcam_gridfile, gcam_config,&
         case_name, gcam2elm_co2_mapping_file, gcam2elm_luc_mapping_file,&
         gcam2elm_woodharvest_mapping_file, elm2gcam_mapping_file,&
         base_co2_surface_file, base_co2_aircraft_file,&
         base_npp_file, base_hr_file, base_pft_file, &
         read_scalars,read_elm_from_file, write_co2, write_scalars, &
         elm_iac_carbon_scaling, iac_elm_co2_emissions, num_lat,&
         num_lon, num_pft, num_harvest, num_gcam_energy_regions,&
         num_gcam_land_regions, num_iac2elm_landtypes,&
         num_emiss_sectors, num_emiss_regions,&
         gcam2glm_baselu, gcam2glm_basebiomass, gcam2glm_glumap,&
         base_co2emis_surface, base_co2emis_aircraft

    nlfilename_iac = "gcam_in"

    inquire (file = trim(nlfilename_iac), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_iac does NOT exist:'&
            //trim(nlfilename_iac)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_iac does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(iulog,*) 'Read in gcam_inparm namelist from: ', trim(nlfilename_iac)
       open( unitn, file=trim(nlfilename_iac), status='old' )
       ier = 1
       !do while ( ier /= 0 )
          read(unitn, gcam_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on gcam_inparm read' )
          endif
       !end do
       call shr_file_freeUnit( unitn )
    end if

    ! Just make sure that worked

    ! Don't need to mpi_bcast the namelist, because we don't have other procs
    ! to communitate with.  But if we did, we'd mpi_bcast them around
    ! mpicom_iac here.  

    if (masterproc) then
       write(iulog,*) 'define GCAM run:'
       ! write out namelist stuff here
       write(iulog, '(A,A)') "base_co2_surface_file = ", trim(base_co2_surface_file )
       write(iulog, '(A,A)') "base_co2_aircraft_file = ", trim(base_co2_aircraft_file )
       write(iulog, '(A,A)') "base_hr_file = ", trim(base_hr_file)
       write(iulog, '(A,A)') "base_npp_file = ", trim(base_npp_file )
       write(iulog, '(A,A)') "base_pft_file = ", trim(base_pft_file )
       write(iulog, '(A,A)') "case_name = ", trim(case_name)
       write(iulog, '(A,A)') "elm2gcam_mapping_file = ", trim(elm2gcam_mapping_file)
       write(iulog, '(A,L)') "elm_iac_carbon_scaling = ", elm_iac_carbon_scaling
       write(iulog, '(A,A)') "gcam2elm_co2_mapping_file = ", trim(gcam2elm_co2_mapping_file )
       write(iulog, '(A,A)') "gcam2elm_luc_mapping_file = ", trim(gcam2elm_luc_mapping_file)
       write(iulog, '(A,A)') "gcam2elm_woodharvest_mapping_file = ", trim(gcam2elm_woodharvest_mapping_file)
       write(iulog, '(A,A)') "gcam_config = ", trim(gcam_config)
       write(iulog, '(A,A)') "gcam_gridfile = ", trim(gcam_gridfile)
       write(iulog, '(A,L)') "iac_elm_co2_emissions = ", iac_elm_co2_emissions
       write(iulog, '(A,I)') "num_emiss_regions = ",num_emiss_regions
       write(iulog, '(A,I)') "num_emiss_sectors = ",num_emiss_sectors
       write(iulog, '(A,I)') "num_gcam_energy_regions = ",num_gcam_energy_regions
       write(iulog, '(A,I)') "num_gcam_land_regions = ",num_gcam_land_regions
       write(iulog, '(A,I)') "num_iac2elm_landtypes = ",num_iac2elm_landtypes
       write(iulog, '(A,I)') "num_lat = ",num_lat
       write(iulog, '(A,I)') "num_lon = ",num_lon
       write(iulog, '(A,I)') "num_pft = ",num_pft
       write(iulog, '(A,L)') "read_elm_from_file = ",read_elm_from_file
       write(iulog, '(A,L)') "read_scalars = ",read_scalars
       write(iulog, '(A,L)') "write_co2 = ",write_co2
       write(iulog, '(A,L)') "write_scalars = ",write_scalars
       write(iulog, '(A,F)') "base_co2emis_surface = ",base_co2emis_surface
       write(iulog, '(A,F)') "base_co2emis_aircraft = ",base_co2emis_aircraft

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
    call ncd_inqdid(ncid,'lsmlon',dimid)
    call ncd_inqdlen(ncid,dimid,iac_ctl%nlon)
    call ncd_inqdid(ncid,'lsmlat',dimid)
    call ncd_inqdlen(ncid,dimid,iac_ctl%nlat)
    iac_ctl%ngrid=iac_ctl%nlon*iac_ctl%nlat

    ! No decomposition
    iac_ctl%begg = 1
    iac_ctl%endg = iac_ctl%ngrid

    ! From the namelist/config 
    iac_ctl%npft = num_pft
    iac_ctl%nharvest = num_harvest

    ! Start allocating vars and control vectors

    ! These are the 1D coordinate vectors in each direction
    allocate(iac_ctl%lon(iac_ctl%nlon))
    allocate(iac_ctl%lat(iac_ctl%nlat))

    ! These are the grid-based coordinate values
    allocate(iac_ctl%ilon(iac_ctl%ngrid))
    allocate(iac_ctl%jlat(iac_ctl%ngrid))

    allocate(iac_ctl%gindex(iac_ctl%ngrid))

    ! dimensions: lon,lat
    allocate(iac_ctl%iacmask(iac_ctl%nlon,iac_ctl%nlat))
    allocate(iac_ctl%landfrac(iac_ctl%nlon,iac_ctl%nlat))
    allocate(iac_ctl%area(iac_ctl%nlon,iac_ctl%nlat))
    allocate(iac_ctl%vegfrac(iac_ctl%nlon,iac_ctl%nlat))

    ! Dimension order: (lon,lat,pft)
    allocate(lnd2iac_vars%npp(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(lnd2iac_vars%hr(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(lnd2iac_vars%pftwgt(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))

    allocate(iac2lnd_vars%pct_pft(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(iac2lnd_vars%pct_pft_prev(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%npft))
    allocate(iac2lnd_vars%harvest_frac(iac_ctl%nlon,iac_ctl%nlat,iac_ctl%nharvest))

    allocate(iac2atm_vars%co2emiss(iac_ctl%nlon,iac_ctl%nlat))

    ! Assign our global index, ilon, and jlat, to go from (g)
    ! dimension to (lon,lat) dimensions.
    do j=1,iac_ctl%nlat
       do i=1,iac_ctl%nlon
          n = (j-1)*iac_ctl%nlon + i
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

    ! Netcdf variables stored as (lat,lon) in netcdf (row-major), which
    ! is (lon,lat) in fortran (column-major)

    ! lon
    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read IAC longitudes')
    if (masterproc) write(iulog,*) 'Read LONGXY ',minval(tempr),maxval(tempr)
    do i=1,iac_ctl%nlon
       ! lon is first index in column major indexing
       iac_ctl%lon(i) = tempr(i,1)
    enddo

    ! lat
    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read IAC latitudes')
    if (masterproc) write(iulog,*) 'Read LATIXY ',minval(tempr),maxval(tempr)
    do j=1,iac_ctl%nlat
       ! lat is second index in column major indexing
       iac_ctl%lat(j) = tempr(1,j)
    enddo

    ! grid cell area
    call ncd_io(ncid=ncid, varname='AREA', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read IAC area')
    if (masterproc) write(iulog,*) 'Read AREA ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          iac_ctl%area(i,j) = tempr(i,j)
       enddo
    enddo

    ! landfrac
    call ncd_io(ncid=ncid, varname='LANDFRAC_PFT', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read IAC landfrac')
    if (masterproc) write(iulog,*) 'Read LANDFRAC_PFT ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          iac_ctl%landfrac(i,j) = tempr(i,j)
       enddo
    enddo

    ! veg land unit frac 
    call ncd_io(ncid=ncid, varname='PCT_NATVEG', flag='read', data=tempr, &
                readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)// &
          ' ERROR: read IAC vegfrac')
    if (masterproc) write(iulog,*) &
          'Read PCT_NATVEG ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          iac_ctl%vegfrac(i,j) = tempr(i,j) / 100.0_R8
       enddo
    enddo

    ! pft land mask
    call ncd_io(ncid=ncid, varname='PFTDATA_MASK', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read IAC landmask')
    if (masterproc) write(iulog,*) 'Read PFTDATA_MASK ',minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          iac_ctl%iacmask(i,j) = itempr(i,j)
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
