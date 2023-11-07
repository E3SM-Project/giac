module iac_comp_mct
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Interface of the active gcam model to the MCT drive
  !
  ! !uses:
  use seq_flds_mod
  use shr_kind_mod     , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use shr_taskmap_mod  , only : shr_taskmap_write
  use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
  use seq_comm_mct     , only : info_taskmap_comp
  use seq_timemgr_mod  
  use seq_infodata_mod 
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name

  use gcam_cpl_indices   , only : gcam_cpl_indices_set
  use iac_data_mod        ! including gdata, GClock
  use iac_init_mod
  use iac_import_export
  use gcam_comp_mod
  use glm_comp_mod
  use gcam2glm_mod
  use glm2iac_mod
  use gcam_var_mod        , only : run_gcam, gcam_nlon, gcam_nlat,  iulog, &
                                nsrStartup, nsrContinue, nsrBranch, & 
                                inst_index, inst_suffix, inst_name, &
                                gcam_active, gcam_var_set, num_lon, num_lat, &
                                gcam_alarm
  use iac_spmd_mod
  !use iac_fields_mod
  use ESMF
  use mct_mod

! Stub in other moduals for running iac
! use IacMod
! use IacVar

  implicit none
  SAVE
  private                              ! By default make data private

! Internal variable structures - probably move to a mod or something later
! We probably need to allocate these during init, if I can figure out
! the sizes 

! I'm naming these with a _data suffix to differentiate between these
  ! and the lnd2iac_vars, which are AVects.

  real*8, pointer :: gcam2glm_data(:,:)  ! gcam output for glm, set by gcam, passed to gcam2glm
  real*8, pointer :: gcam_emis_data(:,:) ! co2 flux output for atm
  real*8, pointer :: glmi_data(:,:)      ! glm input, converted from gcam2glm_data
  real*8, pointer :: glm2lnd_data(:,:)   ! This is used for glm output, but it is called glmo in subroutines that use it
  real*8, pointer :: glm_wh_data(:)
  real*8, pointer :: gcamoco2sfcjan(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcfeb(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcmar(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcapr(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcmay(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcjun(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcjul(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcaug(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcsep(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcoct(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcnov(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2sfcdec(:,:)  ! gcam output for eam, needs to be passed through coupler                                        
  real*8, pointer :: gcamoco2airlojan(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlofeb(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlomar(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airloapr(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlomay(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlojun(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlojul(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airloaug(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlosep(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlooct(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlonov(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airlodec(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhijan(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhifeb(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhimar(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhiapr(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhimay(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhijun(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhijul(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhiaug(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhisep(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhioct(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhinov(:,:)  ! gcam output for eam, needs to be passed through coupler                                      
  real*8, pointer :: gcamoco2airhidec(:,:)  ! gcam output for eam, needs to be passed through coupler


  !
  ! PUBLIC MEMBER FUNCTIONS:

  public :: iac_init_mct               ! iac initialization
  public :: iac_run_mct                ! iac run phase
  public :: iac_final_mct              ! iac finalization/cleanup
  !
  ! PRIVATE MEMBER FUNCTIONS:
  
  ! PRIVATE MEMBER FUNCTIONS:

  private :: iac_SetgsMap_mct         ! Set the iac model MCT GS map
  private :: iac_domain_mct           ! Set the iac model domain information
  private :: gcam2iac_copy            ! Copy z->a gcam output to coupler variables


! REVISION HISTORY:
! Author: Tim Shippert
!===============================================================
contains
!===============================================================

!========================================================================

  subroutine iac_init_mct( EClock, cdata_z, x2z_z, z2x_z, NLFilename)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initialize iac model and hook in arrays from lnd module
    !
    ! !ARGUMENTS:
    type(ESMF_Clock),           intent(inout) :: EClock     ! Input synchronization clock
    type(seq_cdata),            intent(inout) :: cdata_z    ! Input iac driver data
    type(mct_aVect) ,           intent(inout) :: x2z_z      ! Iac import state
    type(mct_aVect),            intent(inout) :: z2x_z      ! Iac export state
    character(len=*), optional, intent(in)    :: NLFilename ! Namelist filename to read
    !
    ! !LOCAL VARIABLES:
    ! 
    integer :: IACID	                             ! iac identifyer
    integer :: mpicom_iac                            ! mpi communicator
    type(mct_gsMap),         pointer :: gsMap_iac    ! iac model MCT GS map
    type(mct_gGrid),         pointer :: dom_z        ! iac model domain
    type(seq_infodata_type), pointer :: infodata     ! CESM driver level info data
    integer :: lsize                                 ! size of attribute vector
    integer :: g,i,j,n                               ! indices
    logical :: exists                                ! true if file exists
    integer :: nsrest                                ! restart type
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level

    ! Possibly not useful, as we expect gcam to run with just one proc
    integer :: begg, endg                            ! Region indeces

    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=32), parameter :: sub = 'iac_init_mct'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    !---------------------------------------------------------------------------
   
    ! Pull out our cdata and communication stuff, from wherever it is set
    ! these days
    call seq_cdata_setptrs(cdata_z, ID=IACID, mpicom=mpicom_iac, &
         gsMap=gsMap_iac, dom=dom_z, infodata=infodata)

    ! Initialize gcam MPI communicator 
    call spmd_init(mpicom_iac, IACID)

    if (masterproc) then
       write(iulog,*) "masterproc before gcam_cpl_indices_init is ", masterproc
    end if

    ! Determine attriute vector indices
    call gcam_cpl_indices_init()
    call gcam_cpl_indices_set()

    ! I see this kind of thing everywhere, so why not
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','iac_init_mct:start::',lbnum)
    endif
#endif                      

    ! Initialize io logs
    inst_name   = seq_comm_name(IACID)
    inst_index  = seq_comm_inst(IACID)
    inst_suffix = seq_comm_suffix(IACID)

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       ! I'm copying from other comps, but we may want to change teh
       ! namelist file here
       inquire(file='iac_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('iac_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "IAC initialization"
    else
       iulog = shrlogunit
    end if
    
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Generic timer and infodata stuff
    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )
    call seq_infodata_GetData(infodata, case_name=caseid,                  &
                              case_desc=ctitle, start_type=starttype,      &
                              brnch_retain_casename=brnch_retain_casename, &
                              model_version=version,                       &
                              hostname=hostname, username=username)

    ! Startup type - we'll be generic here too for now
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( sub//' ERROR: unknown starttype' )
    end if

    ! Modeled on RTM, we'll modify later - passing stuff back to gcam
    ! from our namelists and configs
    call gcam_var_set(caseid_in=caseid, ctitle_in=ctitle,             &
         brnch_retain_casename_in=brnch_retain_casename, &
         nsrest_in=nsrest, version_in=version,           &
         hostname_in=hostname, username_in=username)

    ! Do whatever init gcam needs
    ! I.e. read namelist and grid, etc.
    call iac_init(EClock)

    ! Now that we have the grid, we can allocate some of our working
    ! arrays. 
    call gcam_init_mod(gcam2glm_data,gcam_emis_data,                          &
          gcamoco2sfcjan, gcamoco2sfcfeb,                                     &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec)
    call gcam2glm_init_mod()
    call glm_init_mod(glmi_data, glm_wh_data, glm2lnd_data)
    call glm2iac_init_mod(glm2lnd_data)

    ! Here it gets tricky, as I copy from rof and other components.
    ! I'm going to go ahead and use the whole begg,endg domain
    ! indeces for now - I'm not sure if there is any meaningful
    ! overhead if we do this with only one proc, but if we ever go
    ! multiproc it will be better to do this now, and it lets me cut
    ! and paste for now.
    if ( .not. gcam_active) then
       call seq_infodata_PutData( infodata, iac_present   =.false.)
       call seq_infodata_PutData( infodata, iac_prognostic=.false.)
    else
       ! Initialize memory for input state
       begg = iac_ctl%begg
       endg = iac_ctl%endg
       !allocate (totrunin(begg:endg,nt_gcam))
       
       ! Initialize iac gsMap 
       call iac_SetgsMap_mct( mpicom_iac, IACID, gsMap_iac)

       ! Initialize iac domain - I feel like this should be mostly
       ! null, but I think the global seg mapping stuff needs it
       lsize = mct_gsMap_lsize(gsMap_iac, mpicom_iac)
       call iac_domain_mct( lsize, gsMap_iac, dom_z )
       
       ! Initialize input attribute vectors
       call mct_aVect_init(x2z_z, rList=seq_flds_x2z_fields, lsize=lsize)
       call mct_aVect_zero(x2z_z)
       
       ! Initialize output attribute vectors
       call mct_aVect_init(z2x_z, rList=seq_flds_z2x_fields, lsize=lsize)
       call mct_aVect_zero(z2x_z) 

       ! Create mct iac export state - try to review what this is.
       call iac_export(iac2lnd_vars, iac2atm_vars, z2x_z%rattr)
    end if

    ! Fill in infodata - of course, have to review all this
    ! this happens before lnd init, so the elm flag is set later 
    call seq_infodata_PutData( infodata, iac_present=gcam_active, &
         iac_prognostic=gcam_active, &
         iac_nx = gcam_nlon, iac_ny = gcam_nlat)

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','iac_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif


  end subroutine iac_init_mct

!---------------------------------------------------------------------------

  subroutine iac_run_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-------------------------------------------------------
    ! DESCRIPTION:
    ! Run IAC model

    ! !USES:

    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for iac model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from iac model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from iac model

    ! LOCAL VARIABLES:
    integer      :: ymd_sync             ! Sync date (YYYYMMDD)
    integer      :: yr_sync              ! Sync current year
    integer      :: mon_sync             ! Sync current month
    integer      :: day_sync             ! Sync current day
    integer      :: tod_sync             ! Sync current time of day (sec)
    integer      :: ymd                  ! gcam current date (YYYYMMDD)
    integer      :: yr                   ! gcam current year
    integer      :: mon                  ! gcam current month
    integer      :: day                  ! gcam current day
    integer      :: tod                  ! gcam current time of day (sec)
    integer      :: dtime                ! time step increment (sec)
    integer,save :: nstep                ! time step index
!    logical      :: nlend_sync           ! Flag signaling last time-step
!    logical      :: nlend                ! .true. ==> last time-step
    logical      :: dosend               ! true => send data back to driver
    real(r8)     :: nextsw_cday          ! calday from clock of next radiation computation
    real(r8)     :: caldayp1             ! gcam calday plus dtime offset
    integer      :: shrlogunit,shrloglev ! old values for share log unit and log level
    integer      :: lbnum                ! input to memory diagnostic
    integer      :: g,i,lsz              ! counters
    real(r8)     :: calday               ! calendar day for nstep
    real(r8)     :: recip                ! reciprical
    logical,save :: first_call = .true.  ! first call work
    integer      :: rc                   ! return value from Clock functions.
    type(seq_infodata_type),pointer :: infodata             ! CESM information from the driver
    type(mct_gGrid),        pointer :: dom_z                ! iac model domain data
    !type(bounds_type)               :: bounds               ! bounds
    character(len=32)               :: rdate                ! date char string for restart file names
    character(len=32), parameter    :: sub = "(iac_run_mct)"


    ! I feel this is probably important
    if (.not.gcam_active) return

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    if (first_call) then
       nstep = 1
    else
       nstep = nstep+1
    endif

    ! Get the model time - I'm not certain what should be sync and
    ! what shouldn't, but for now we'll just follow the pattern
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod,  &
         curr_yr=yr, curr_mon=mon, curr_day=day,&
         dtime=dtime)

    ! We are going to do our time sync check up top, because iac runs
    ! first at the top of the year.

    ! Check that internal clock is in sync with master clock
    ! Again, may need time manager module with these kind of functions
    ! call get_curr_date( yr, mon, day, tod )
    !ymd = yr*10000 + mon*100 + day
    !tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*) trim(sub), ' gcam ymd=',ymd     ,'  gcam tod= ',tod
       write(iulog,*)trim(sub), 'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: GCAM clock is not in sync with Master Sync clock" )
    end if

    write(iulog,*) trim(sub), "Nstep: ", nstep
    write(iulog,*) trim(sub), "EClock: ", yr,mon,day,tod, ymd, dtime

    ! Assign the date to GClock
    GClock(iac_eclock_ymd) = ymd
    GClock(iac_eclock_tod) = tod
    GClock(iac_eclock_dt) = dtime
 
    ! store whether GCAM runs this model year or not
    ! use the hardcoded gcam time step to determine this
    gcam_alarm = .false.
    if (modulo(yr, iac_gcam_timestep)==0) then
       gcam_alarm = .true.
    endif 

    ! Import from land model
    !call t_startf('iac_import')
    call iac_import(x2z_z%rattr, lnd2iac_vars)
    !call iac_import(x2z_z, lnd2iac_vars)
    !call t_stopf('iac_import')
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,tod
!    nlend = seq_timemgr_StopAlarmIsOn( EClock )

    ! Run gcam and glm, unless it is the historical period
    ! TODO: make the iac res match the land res
    ! For now, running the historical period will get the appropriate
    ! land surface files, until we make the iac res match the land
    ! Without gcam/glm, glm2iac_run_mod reads in the reformatted luh2 data

    ! note that gcam handles its own restarts
    ! note that the iac does not run in a timestep that handles model restart
    !    writing
    ! so gcam2glm and glm will always write restart files each year,
    ! then the land model code (elm_driver) will write the two appropriate
    !    restart pointer files when it is time to write the land restart file

    if ( run_gcam ) then  

       ! This calcs and pushes the carbon density scalars to gcam
       !   the push happens only if elm_iac_carbon_scaling is true
       ! this is called regardless in order to calculate and write scalars
       !    and also to interpolate the co2 outputs to annual values and
       !    to downscale them

       ! Run GCAM, for this timestep.  
       ! Inputs taken care of in gcam_comp_mod, so these
       ! are outputs, for input to glm and for export to atm (emis)
       call gcam_run_mod(gcam2glm_data, gcam_emis_data,                          &
             gcamoco2sfcjan, gcamoco2sfcfeb,                                     &
             gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
             gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
             gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
             gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
             gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
             gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
             gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
             gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
             gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
             gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
             gcamoco2airhidec)

       ! Set up to run glm
       call gcam2glm_run_mod(gcam2glm_data, glmi_data, glm_wh_data)

       ! Run glm, which produces the coupled vars to lnd (z->l)
       call glm_run_mod(glmi_data, glm_wh_data, glm2lnd_data)

       ! Copy gcamo atm variables to iac2atm_vars
       call gcam2iac_copy()
    
    end if

    ! Run glm2iac, which runs mksurfdat and produces land input files for ELM
    call glm2iac_run_mod(glm2lnd_data)

    ! Now export the land and atmosphere data
    call iac_export(iac2lnd_vars, iac2atm_vars, z2x_z%rattr)

    ! Advance our timestep in the gcam clock to next year
    call ESMF_ClockAdvance( EClock, rc=rc )
    if (rc .ne. ESMF_SUCCESS) then
       call shr_sys_abort( sub//":: GCAM clock did not advance correctly" )
    endif

    ! You know what - let's just check that
    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod,  &
         curr_yr=yr, curr_mon=mon, curr_day=day,&
         dtime=dtime)

    write(iulog,*) trim(sub), "EClock post-iac run: ", yr,mon,day,tod, ymd, dtime

    first_call = .false.

    ! Reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine iac_run_mct

!===============================================================================

  subroutine iac_final_mct( EClock, cdata_z, x2z_z, z2x_z)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Finalize iac surface model
    !
    ! ARGUMENTS:
    implicit none
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(seq_cdata)  , intent(inout) :: cdata_z   ! Input driver data for iac model
    type(mct_aVect)  , intent(inout) :: x2z_z     ! Import state from iac model
    type(mct_aVect)  , intent(inout) :: z2x_z     ! Export state from iac model
    !-----------------------------------------------------

    ! KVC: do we need to call more finalize subroutines?
    call gcam_final_mod(gcamoco2sfcjan, gcamoco2sfcfeb,   &
          gcamoco2sfcmar, gcamoco2sfcapr, gcamoco2sfcmay, gcamoco2sfcjun,     &
          gcamoco2sfcjul, gcamoco2sfcaug, gcamoco2sfcsep, gcamoco2sfcoct,     &
          gcamoco2sfcnov, gcamoco2sfcdec, gcamoco2airlojan, gcamoco2airlofeb, &
          gcamoco2airlomar, gcamoco2airloapr, gcamoco2airlomay,               &
          gcamoco2airlojun, gcamoco2airlojul, gcamoco2airloaug,               &
          gcamoco2airlosep, gcamoco2airlooct, gcamoco2airlonov,               &
          gcamoco2airlodec, gcamoco2airhijan, gcamoco2airhifeb,               &
          gcamoco2airhimar, gcamoco2airhiapr, gcamoco2airhimay,               &
          gcamoco2airhijun, gcamoco2airhijul, gcamoco2airhiaug,               &
          gcamoco2airhisep, gcamoco2airhioct, gcamoco2airhinov,               &
          gcamoco2airhidec)
    call gcam2glm_final_mod()
  end subroutine iac_final_mct

!================================================================================
! Internal functions to deal with gsmap and domain

  subroutine iac_SetgsMap_mct( mpicom_z, IACID, gsMap_iac)

    !-----------------------------------------------------
    ! DESCRIPTION:
    ! Set the MCT GS map for the iac model
    !
    ! ARGUMENTS:
    implicit none
    integer        , intent(in)    :: mpicom_z      ! MPI communicator for iac model
    integer        , intent(in)    :: IACID         ! iac model identifier
    type(mct_gsMap), intent(inout) :: gsMap_iac     ! MCT gsmap for iac -> land data
    !
    ! LOCAL VARIABLES
    integer,allocatable :: gindex(:)         ! indexing for iac grid cells
    integer :: n, ni                         ! indices
    integer :: lsize,gsize                   ! size of iac data and number of grid cells
    integer :: begg, endg                   ! beg, end iac indices
    integer :: ier                           ! error code
    character(len=32), parameter :: sub = '(iac_SetgsMap_mct)'
    !-----------------------------------------------------

    begg  = iac_ctl%begg
    endg  = iac_ctl%endg

    gcam_nlon = iac_ctl%nlon
    gcam_nlat = iac_ctl%nlat

    lsize = iac_ctl%endg - iac_ctl%begg + 1
    gsize = gcam_nlon*gcam_nlat

    ! Check 
    ni = 0
    do n = begg,endg
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR iac count',n,ni,lsize
          call shr_sys_abort( sub//' ERROR: iac > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR iac total count',ni,lsize
       call shr_sys_abort( sub//' ERROR: iac not equal to expected' )
    endif

    ! Determine gsmap_iac
    allocate(gindex(lsize),stat=ier)
    if(ier/=0) call mct_die(sub,'allocate gindex',ier)
    ni = 0
    do n = begg,endg
       ni = ni + 1
       gindex(ni) = iac_ctl%gindex(n)
    end do
    call mct_gsMap_init( gsMap_iac, gindex, mpicom_z, IACID, lsize, gsize )
    deallocate(gindex)

  end subroutine iac_SetgsMap_mct

!===============================================================================

  subroutine iac_domain_mct( lsize, gsMap_z, dom_z )

    !-----------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Send the iac model domain information to the coupler
    ! 
    ! Note: GCAM is currently only single proc, so domain information may not be
    ! necessary or useful.  But I want the infrastructure in case
    ! that changes, and so that the calling structure is the same for
    ! this componant is the same as for the others.  Also, its
    ! completely possible that I don't really understand what
    ! "domain" means in the context of E3SM.
    !
    ! !ARGUMENTS:
    implicit none
    integer        , intent(in)    :: lsize       ! Size of iac domain information
    type(mct_gsMap), intent(inout) :: gsMap_z     ! Output MCT GS map for iac model
    type(mct_ggrid), intent(out)   :: dom_z       ! Domain information from the iac model
    !
    ! LOCAL VARIABLES
    integer :: n, ni, i, j, ier   ! index
    integer , pointer :: idata(:) ! temporary
    real(r8), pointer :: data(:)  ! temporary
    real(r8) :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)
    character(len=32), parameter :: sub = '(iac_domain_mct)'
    !-----------------------------------------------------

    ! the domain info is for iac-atm only
    ! lat/lon in degrees,  area in radians^2, mask is 1 (atm everywhere) 
    ! the iac frac need to be 1 everywhere, as the data are for grid cell 
    ! TRS - this is something I need to review
    call mct_gGrid_init( GGrid=dom_z, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Allocate memory
    allocate(data(lsize), stat=ier)
    if(ier/=0) call mct_die(sub,'allocate data',ier)

    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    call mct_gsMap_orderedPoints(gsMap_z, iam, idata)
    call mct_gGrid_importIAttr(dom_z,'GlobGridNum',idata,lsize)

    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    data(:) = -9999.0_r8 
    call mct_gGrid_importRAttr(dom_z,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_z,"aream",data,lsize) 
    data(:) = 0.0_r8
    call mct_gGrid_importRAttr(dom_z,"mask" ,data,lsize)

    ! Determine bounds numbering consistency
    ni = 0
    do n = iac_ctl%begg,iac_ctl%endg
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR iac count',n,ni,lsize
          call shr_sys_abort( sub//' ERROR: iac > expected' )
       end if
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR iac total count',ni,lsize
       call shr_sys_abort( sub//' ERROR: iac not equal to expected' )
    endif

    ! lon 
    ni = 0
    do n = iac_ctl%begg,iac_ctl%endg
       ni = ni + 1
       i = iac_ctl%ilon(n)
       data(ni) = iac_ctl%lon(i)
    end do
    call mct_gGrid_importRattr(dom_z,"lon",data,lsize) 

    ! lat
    ni = 0
    do n = iac_ctl%begg,iac_ctl%endg
       ni = ni + 1
       j = iac_ctl%jlat(n)
       data(ni) = iac_ctl%lat(j)
    end do
    call mct_gGrid_importRattr(dom_z,"lat",data,lsize) 

    ! cell Area - pulled from iac_ctl array in km^2
    ni = 0
    do n = iac_ctl%begg,iac_ctl%endg
       ni = ni + 1
       i = iac_ctl%ilon(n)
       j = iac_ctl%jlat(n)
       ! scale area to radians
       data(ni) = iac_ctl%area(i,j)/(re*re)
    end do
    call mct_gGrid_importRattr(dom_z,"area",data,lsize) 

    ! frac - this is one because atm covers whole grid cell
    data(:) = 1.0_r8
    call mct_gGrid_importRattr(dom_z,"frac",data,lsize) 

    ! iacmask = this is one because global coverage
    data(:) = 1.0_r8
    call mct_gGrid_importRAttr(dom_z,"mask" ,data,lsize)

    deallocate(data)
    deallocate(idata)

  end subroutine iac_domain_mct

!====================================================================================

  subroutine gcam2iac_copy()
    !-----------------------------------------------------
    !
    ! !DESCRIPTION:
    ! Copy the gcamoco2* variables into iac2atm_vars, which is what the
    ! coupler code needs.

    ! !ARGUMENTS:
    implicit none
    ! LOCAL VARIABLES
    !-----------------------------------------------------

    ! This is tedious, but seems to be the only way to convert these
    ! text based variables to index based, which is what we'll need
    ! inside the export function.  The alternative is to simply declare
    ! all 36 of these expclitly as iac2atm_var elements, but then we
    ! couldn't loop in iac export and clm import, so we'd be paying the
    ! tedium forward.  So we do it just once, in this subroutine.

    ! Surface co2 flux
    iac2atm_vars%co2sfc(1,:,:) = gcamoco2sfcjan(:,:)
    iac2atm_vars%co2sfc(2,:,:) = gcamoco2sfcfeb(:,:)
    iac2atm_vars%co2sfc(3,:,:) = gcamoco2sfcmar(:,:)
    iac2atm_vars%co2sfc(4,:,:) = gcamoco2sfcapr(:,:)
    iac2atm_vars%co2sfc(5,:,:) = gcamoco2sfcmay(:,:)
    iac2atm_vars%co2sfc(6,:,:) = gcamoco2sfcjun(:,:)
    iac2atm_vars%co2sfc(7,:,:) = gcamoco2sfcjul(:,:)
    iac2atm_vars%co2sfc(8,:,:) = gcamoco2sfcaug(:,:)
    iac2atm_vars%co2sfc(9,:,:) = gcamoco2sfcsep(:,:)
    iac2atm_vars%co2sfc(10,:,:) = gcamoco2sfcoct(:,:)
    iac2atm_vars%co2sfc(11,:,:) = gcamoco2sfcnov(:,:)
    iac2atm_vars%co2sfc(12,:,:) = gcamoco2sfcdec(:,:)

    ! Low alt air co2
    iac2atm_vars%co2airlo(1,:,:) = gcamoco2airlojan(:,:)
    iac2atm_vars%co2airlo(2,:,:) = gcamoco2airlofeb(:,:)
    iac2atm_vars%co2airlo(3,:,:) = gcamoco2airlomar(:,:)
    iac2atm_vars%co2airlo(4,:,:) = gcamoco2airloapr(:,:)
    iac2atm_vars%co2airlo(5,:,:) = gcamoco2airlomay(:,:)
    iac2atm_vars%co2airlo(6,:,:) = gcamoco2airlojun(:,:)
    iac2atm_vars%co2airlo(7,:,:) = gcamoco2airlojul(:,:)
    iac2atm_vars%co2airlo(8,:,:) = gcamoco2airloaug(:,:)
    iac2atm_vars%co2airlo(9,:,:) = gcamoco2airlosep(:,:)
    iac2atm_vars%co2airlo(10,:,:) = gcamoco2airlooct(:,:)
    iac2atm_vars%co2airlo(11,:,:) = gcamoco2airlonov(:,:)
    iac2atm_vars%co2airlo(12,:,:) = gcamoco2airlodec(:,:)

    ! High alt air co2
    iac2atm_vars%co2airhi(1,:,:) = gcamoco2airhijan(:,:)
    iac2atm_vars%co2airhi(2,:,:) = gcamoco2airhifeb(:,:)
    iac2atm_vars%co2airhi(3,:,:) = gcamoco2airhimar(:,:)
    iac2atm_vars%co2airhi(4,:,:) = gcamoco2airhiapr(:,:)
    iac2atm_vars%co2airhi(5,:,:) = gcamoco2airhimay(:,:)
    iac2atm_vars%co2airhi(6,:,:) = gcamoco2airhijun(:,:)
    iac2atm_vars%co2airhi(7,:,:) = gcamoco2airhijul(:,:)
    iac2atm_vars%co2airhi(8,:,:) = gcamoco2airhiaug(:,:)
    iac2atm_vars%co2airhi(9,:,:) = gcamoco2airhisep(:,:)
    iac2atm_vars%co2airhi(10,:,:) = gcamoco2airhioct(:,:)
    iac2atm_vars%co2airhi(11,:,:) = gcamoco2airhinov(:,:)
    iac2atm_vars%co2airhi(12,:,:) = gcamoco2airhidec(:,:)
end subroutine 

!===============================================================================

!==============================================================
! Local functions
!==============================================================

end module iac_comp_mct
