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
  use seq_timemgr_mod
  use netcdf
  use gcam2glm_mod, only : handle_err
  use ESMF

  implicit none
  save


contains
  subroutine iac_init(EClock)
    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Read in namelist, define grid, allocate variables
    !---------------------------------------------------------------------------

    ! !Uses:
    use shr_kind_mod, only: CL => SHR_KIND_CL, CXX => SHR_KIND_CXX

    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock     ! Input synchronization clock

    ! !LOCAL VARIABLES:
    character(len=CL) :: nlfilename_iac
    integer :: unitn, ier, dimid, i, j, n
    type(file_desc_t) :: ncid
    logical :: lexist, found
    real(r8), pointer :: tempr(:,:)
    integer, pointer :: itempr(:,:)
    character(len=32) :: subname = 'iac_init'
    ! these are for reading the dynamic land file
    integer :: curr_yr, ncid_int                               ! current model year
    integer :: nlon, nlat, ntime, npft, indprev, varid, ierr
    integer, dimension(4) :: start4, count4
    integer, allocatable :: lsf_years(:)

    ! GCAM namelist
    namelist /gcam_inparm/ &
         case_name, &
         num_pft, num_harvest, num_lat, num_lon, &
         num_gcam_energy_regions, num_gcam_land_regions, &
         num_iac2elm_landtypes, num_emiss_sectors, &
         num_emiss_ctys, num_periods, &
         gcam_config, base_gcam_co2_file, base_gcam_lu_wh_file, &
         base_co2_surface_file, base_co2_shipment_file, base_co2_aircraft_file, &
         base_npp_file, base_hr_file, base_pft_file, &
         gcam2elm_co2_mapping_file, gcam2elm_luc_mapping_file, &
         gcam2elm_woodharvest_mapping_file, gcam2elm_cdensity_mapping_file, &
         gcam_gridfile, elm2gcam_mapping_file, &
         gcam2glm_glumap, gcam2glm_baselu, gcam2glm_basebiomass, &
         country2grid_map, country2region_map, pop_iiasa_file, gdp_iiasa_file, &
         pop_gcam_file, gdp_gcam_file, co2_gcam_file, &
         surface_co2_downscaling_method, &
         crop_addtreeonly, crop_setherbfracrem, crop_setavailtreefracrem, &
         pasture_addtreeonly, pasture_setherbfracrem, pasture_setavailtreefracrem, &         
         fdyndat_ehc, &
         read_scalars, write_scalars, write_co2, &
         elm_ehc_agyield_scaling, elm_ehc_carbon_scaling, ehc_eam_co2_emissions, &
         gcam_spinup, run_gcam 

    nlfilename_iac = "gcam_in"

    inquire (file = trim(nlfilename_iac), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_iac does NOT exist:'&
            //trim(nlfilename_iac)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_iac does not exist')
    end if

    if (masterproc) then
       unitn = shr_file_getunit()
       write(iulog,*) '('//trim(subname)//')', 'Read in gcam_inparm namelist from: ', trim(nlfilename_iac)
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

    ! get the current model year for readin in the start-of-year pft and harvest data
    call seq_timemgr_EClockGetData(EClock, curr_yr=curr_yr)
    write(iulog,*) "(", subname, ") Current model year is ", curr_yr

    ! Don't need to mpi_bcast the namelist, because we don't have other procs
    ! to communitate with.  But if we did, we'd mpi_bcast them around
    ! mpicom_iac here.  

    ! write the namelist values to the log to make sure they are read in
    if (masterproc) then
       write(iulog,*) '('//trim(subname)//')', 'Namelist for IAC/GCAM run:'
       write(iulog,*) 'GCAM case name:', trim(case_name)
       write(iulog,*) 'grid and region size parameters:'
       write(iulog, '(A,I20)') "num_pft = ",num_pft
       write(iulog, '(A,I20)') "num_harvest = ",num_harvest
       write(iulog, '(A,I20)') "num_lat = ",num_lat
       write(iulog, '(A,I20)') "num_lon = ",num_lon
       write(iulog, '(A,I20)') "num_gcam_energy_regions = ", &
                             num_gcam_energy_regions
       write(iulog, '(A,I20)') "num_gcam_land_regions = ",num_gcam_land_regions
       write(iulog, '(A,I20)') "num_iac2elm_landtypes = ",num_iac2elm_landtypes
       write(iulog, '(A,I20)') "num_emiss_sectors = ",num_emiss_sectors
       write(iulog, '(A,I20)') "num_emiss_ctys = ", num_emiss_ctys
       write(iulog, '(A,I20)') "num_periods = ", num_periods

       write(iulog,*) 'GCAM config and init files:'
       write(iulog, '(A,A)') "gcam_config = ", trim(gcam_config)
       write(iulog, '(A,A)') "base_gcam_co2_file = ", trim(base_gcam_co2_file)
       write(iulog, '(A,A)') "base_gcam_lu_wh_file = ", trim(base_gcam_lu_wh_file)
       write(iulog, '(A,A)') "base_co2_surface_file = ", &
                             trim(base_co2_surface_file )
       write(iulog, '(A,A)') "base_co2_shipment_file = ", trim(base_co2_shipment_file )
       write(iulog, '(A,A)') "base_co2_aircraft_file = ", trim(base_co2_aircraft_file )
       write(iulog, '(A,A)') "base_npp_file = ", trim(base_npp_file )
       write(iulog, '(A,A)') "base_hr_file = ", trim(base_hr_file)
       write(iulog, '(A,A)') "base_pft_file = ", trim(base_pft_file )
       write(iulog, '(A,A)') "gcam2elm_co2_mapping_file = ", trim(gcam2elm_co2_mapping_file )
       write(iulog, '(A,A)') "gcam2elm_luc_mapping_file = ", trim(gcam2elm_luc_mapping_file)
       write(iulog, '(A,A)') "gcam2elm_woodharvest_mapping_file = ", trim(gcam2elm_woodharvest_mapping_file)
       write(iulog, '(A,A)') "gcam2elm_cdensity_mapping_file = ", trim(gcam2elm_cdensity_mapping_file)

       write(iulog,*) 'grid mapping and initialization files:'
       write(iulog, '(A,A)') "gcam_gridfile = ", trim(gcam_gridfile)
       write(iulog, '(A,A)') "elm2gcam_mapping_file = ", &
                             trim(elm2gcam_mapping_file)
       write(iulog, '(A,A)') "gcam2glm_glumap = ", trim(gcam2glm_glumap)
       write(iulog, '(A,A)') "gcam2glm_baselu = ", trim(gcam2glm_baselu)
       write(iulog, '(A,A)') "gcam2glm_basebiomass = ", trim(gcam2glm_basebiomass)

       write(iulog,*) 'config and input files for emiss downscaling:'
       write(iulog, '(A,A)') "country2grid_map = ", trim(country2grid_map)
       write(iulog, '(A,A)') "country2region_map = ", trim(country2region_map)
       write(iulog, '(A,A)') "pop_iiasa_file = ", trim(pop_iiasa_file)
       write(iulog, '(A,A)') "gdp_iiasa_file = ", trim(gdp_iiasa_file)
       write(iulog, '(A,A)') "pop_gcam_file = ", trim(pop_gcam_file)
       write(iulog, '(A,A)') "gdp_gcam_file = ", trim(gdp_gcam_file) 
       write(iulog, '(A,A)') "co2_gcam_file = ", trim(co2_gcam_file)
       write(iulog, '(A,A)') "surface_co2_downscaling_method = ", trim(surface_co2_downscaling_method)

       write(iulog,*) 'future land conversion assumptions:'
       write(iulog, '(A,I)') "crop_addtreeonly = ",crop_addtreeonly
       write(iulog, '(A,F)') "crop_setherbfracrem = ",crop_setherbfracrem
       write(iulog, '(A,F)') "crop_setavailtreefracrem = ",crop_setavailtreefracrem
       write(iulog, '(A,I)') "pasture_addtreeonly = ",pasture_addtreeonly
       write(iulog, '(A,F)') "pasture_setherbfracrem = ",pasture_setherbfracrem
       write(iulog, '(A,F)') "pasture_setavailtreefracrem = ",pasture_setavailtreefracrem

       write(iulog,*) 'name of dynamic landuse timeseries file:'
       write(iulog, '(A,A)') "fdyndat_ehc = ", trim(fdyndat_ehc)

       write(iulog,*) 'rumtime options:'
       write(iulog, '(A,L)') "read_scalars = ",read_scalars
       write(iulog, '(A,L)') "write_scalars = ",write_scalars
       write(iulog, '(A,L)') "write_co2 = ",write_co2
       write(iulog, '(A,L)') "elm_ehc_agyield_scaling = ", elm_ehc_agyield_scaling
       write(iulog, '(A,L)') "elm_ehc_carbon_scaling = ", elm_ehc_carbon_scaling
       write(iulog, '(A,L)') "ehc_eam_co2_emissions = ", ehc_eam_co2_emissions
       write(iulog, '(A,L)') "gcam_spinup = ",gcam_spinup
       write(iulog, '(A,L)') "run_gcam = ",run_gcam

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
       call shr_sys_abort ('IAC_INIT: FAILED to find '//trim(gcam_gridfile))
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

    ! Dimension: month,lon,lat
    allocate(iac2atm_vars%co2sfc(12,iac_ctl%nlon,iac_ctl%nlat))
    allocate(iac2atm_vars%co2airlo(12,iac_ctl%nlon,iac_ctl%nlat))
    allocate(iac2atm_vars%co2airhi(12,iac_ctl%nlon,iac_ctl%nlat))

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
    ! configuration file - this is a landuse.timeseries file

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
    call ncd_io(ncid=ncid, varname='PFTDATA_MASK', flag='read', data=itempr,&
                readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//&
                             ' ERROR: read IAC landmask')
    if (masterproc) write(iulog,*) 'Read PFTDATA_MASK',&
                                minval(tempr),maxval(tempr)

    do i=1,iac_ctl%nlon
       do j=1,iac_ctl%nlat
          iac_ctl%iacmask(i,j) = itempr(i,j)
       enddo
    enddo
 
    deallocate(tempr)
    deallocate(itempr)             

    call ncd_pio_closefile(ncid)

    ! read in the dynamic land surface file fdyndat_ehc
    !   and fill the iac2lnd_vars%pct_pft with the start-of-current-year data
    ! this will be copied into iac2lnd_vars%pct_pft_prev as needed
    !    so this will work as-is for restarts also
    ! the iac2lnd_vars%pct_pft and iac2lnd_vars%harvest_frac are filled as needed (after the copy)
    ! this is so that these data no longer have to be read from the land surface
    !    file at runtime; they are now set in memory by mksurfdata,
    !    (and the current moved to previous before the new current are set)

    ierr = nf90_open(fdyndat_ehc,nf90_nowrite,ncid_int)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! first get the dimensions
    ierr= nf90_inq_dimid(ncid_int, "lsmlon", dimid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid_int, dimid, len=nlon)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid_int, "lsmlat", dimid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid_int, dimid, len=nlat)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid_int, "time", dimid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid_int, dimid, len=ntime)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inq_dimid(ncid_int, "natpft", dimid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_inquire_dimension(ncid_int, dimid, len=npft)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
   
    ! get the years in the dynamic land surface file
    ierr= nf90_inq_varid(ncid_int, "YEAR", varid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    allocate(lsf_years(ntime), stat=ierr)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_get_var(ncid_int, varid, lsf_years)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ! find the last matching year
    indprev = findloc(lsf_years, curr_yr, dim=1, back=.true.)

    ! prev pct pft, but put it in pct_pft for now
    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = indprev
    count4(1) = nlon
    count4(2) = nlat
    count4(3) = npft
    count4(4) = 1
    ierr= nf90_inq_varid(ncid_int, "PCT_NAT_PFT", varid)
    if(ierr /= nf90_NoErr) call handle_err(ierr)
    ierr= nf90_get_var(ncid_int, varid, iac2lnd_vars%pct_pft, start=start4, count=count4)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    ierr= nf90_close(ncid_int)
    if(ierr /= nf90_NoErr) call handle_err(ierr)

    deallocate(lsf_years)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize Restart
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ToDo: adivi: so far I think this is the same for a restart

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize history handler
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (masterproc) write(iulog,*) '('//trim(subname)//') finished'

  end subroutine iac_init
end module iac_init_mod
