module patmo_volc
  use patmo_commons, only: cellsNumber, chemSpeciesNumber, reactionsNumber, maxNameLength
  implicit none
  private

  type volcano_event
     real*8::startTime
     real*8::duration
     real*8::center
     real*8::sigma
     real*8::so2ColumnFlux
     real*8::ashTau550
     real*8::ashLifetime
     real*8::ashHorizontalLifetime
     real*8::ashSettling
     real*8::ashWavelengthExp
     real*8::ashRiseStart
     real*8::ashRiseTime
     real*8::ashVerticalDiffusion
     real*8::ashParticleRadius
     real*8::ashParticleDensity
     real*8::ashAirViscosity
     real*8::ashGravity
     real*8::ashCunninghamFactor
     real*8::sulfateTau550
     real*8::sulfateLifetime
     real*8::sulfateHorizontalLifetime
     real*8::sulfateSettling
     real*8::sulfateWavelengthExp
     real*8::sulfateRiseStart
     real*8::sulfateRiseTime
     real*8::sulfateVerticalDiffusion
     real*8::sulfateParticleRadius
     real*8::sulfateParticleDensity
     real*8::sulfateAirViscosity
     real*8::sulfateGravity
     real*8::sulfateCunninghamFactor
     real*8::sulfateStartDelay
     real*8::sulfateFormationTime
  end type volcano_event

  logical::volcanoEnabled = .false.
  logical::volcanoClockRunning = .true.
  logical::ashStateReady = .false.
  logical::warnedMissingSO2 = .false.
  logical::hasSO2Source = .false.
  logical::hasAshOpacity = .false.
  logical::hasSulfateOpacity = .false.
  integer::volcanoEventsNumber = 0
  integer::so2IndexCache = -2
  integer,parameter::volcComponentAsh = 1
  integer,parameter::volcComponentSulfate = 2
  integer,parameter::volcComponentTotal = 3
  real*8::volcanoElapsedTime = 0d0
  real*8::ashStateTime = 0d0
  real*8::earlyOutputStep = 0d0, earlyOutputUntil = 0d0
  type(volcano_event),allocatable::events(:)
  real*8,allocatable::ashLayerTau550(:,:)
  real*8,allocatable::sulfateLayerTau550(:,:)
  real*8,allocatable::sulfatePrecursor(:,:)

  public::patmo_volc_loadEvents
  public::patmo_volc_reset
  public::patmo_volc_setEnabled
  public::patmo_volc_setClockRunning
  public::patmo_volc_setTime
  public::patmo_volc_getTime
  public::patmo_volc_advanceTime
  public::patmo_volc_holdUntilEquilibrium
  public::patmo_volc_startAfterEquilibrium
  public::patmo_volc_injectSO2Column
  public::patmo_volc_addSources
  public::patmo_volc_applyAshOpacity
  public::patmo_volc_dumpOpticalDepth
  public::patmo_volc_dumpAshProfile
  public::patmo_volc_dumpState
  public::patmo_volc_limitStep
  public::patmo_volc_setEarlyOutput
  public::patmo_volc_isActive

  !======================================================================
  ! Run settings and equilibrium monitor (volcano_run.in)
  !======================================================================
  real*8 :: spinup_max_years=40d0, volcano_duration_day=730d0
  real*8 :: equilibrium_relative_limit=1d-4, equilibrium_mixing_floor=1d-15
  integer :: stable_days_required=30
  namelist /volcano_run/ spinup_max_years, volcano_duration_day, &
       equilibrium_relative_limit, equilibrium_mixing_floor, stable_days_required
  real*8::background_previous(cellsNumber,chemSpeciesNumber)
  integer::stable_days=0

  !======================================================================
  ! Formal diagnostic state (volcano_history.in)
  !======================================================================
  ! Optional build/volcano_history.in overrides these defaults using a namelist:
  ! &history_output
  ! species_names='SO2','SO3','H2SO4', species_every_s=600, species_end_s=259200,
  ! reaction_ids=33,34,56, reaction_every_s=600, reaction_end_s=259200,
  ! solar_every_s=21600, solar_end_s=259200,
  ! history_prefix='volcano_history'
  ! /
  ! An interval <= 0 disables that output. End times are seconds AFTER the
  ! first event, independent of the full simulation duration. Final endpoints
  ! are included even when they are not multiples of the interval.
  character(len=maxNameLength) :: species_names(chemSpeciesNumber)=''
  integer :: reaction_ids(reactionsNumber)=0
  real*8 :: species_every_s=3600d0, species_end_s=259200d0
  real*8 :: reaction_every_s=3600d0, reaction_end_s=259200d0
  real*8 :: solar_every_s=21600d0, solar_end_s=259200d0
  character(len=128) :: history_prefix='volcano_history'
  namelist /history_output/ species_names,reaction_ids,species_every_s,species_end_s, &
       reaction_every_s,reaction_end_s,solar_every_s,solar_end_s,history_prefix
  real*8 :: eruption_start, every(3), finish(3), next_sample(3)
  integer :: units(3), ns, nr, species_indices(chemSpeciesNumber)
  integer(kind=8) :: sample_index(3)=0
  logical :: enabled(3), opened=.false.

  !======================================================================
  ! Standalone pre-run settings (volcano_prerun.in)
  !======================================================================

  type,public :: volcano_prerun_settings
     character(len=256)::eventFile="volcano_events.dat",profileFile="profile.dat"
     character(len=256)::photoMetricFile="xsecs/photoMetric.dat"
     character(len=256)::outputFile="volcano_optical_depth.dat"
     character(len=256)::ashProfileFile="volcano_ash_profile.dat"
     character(len=16)::outputTimeUnit="hour"
     real*8::outputTimeStep=3600d0,outputTimeUnitSeconds=3600d0
     real*8::wavelengthStepNm=5d0,tauFloor=1d-6,endAfterStart=-1d0
     real*8::earlyStep=300d0,earlyUntil=3600d0
  end type volcano_prerun_settings

  public::volc_transport,volc_form_sulfate
  public::patmo_volc_readRunSettings,patmo_volc_getRunLimits
  public::patmo_volc_beginEquilibriumCheck,patmo_volc_checkEquilibrium
  public::patmo_volc_requireEquilibrium
  public::patmo_volc_historyConfigure,patmo_volc_historyBegin,patmo_volc_historyLimit
  public::patmo_volc_historySample,patmo_volc_historyClose
  public::patmo_volc_configurePreRun,patmo_volc_runPreRun

contains

  !======================================================================
  ! Run configuration and equilibrium gate
  !======================================================================
  subroutine patmo_volc_readRunSettings()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    logical::exists
    integer::u,ios
    character(len=512)::message
    inquire(file='volcano_run.in',exist=exists)
    if(exists) then
       open(newunit=u,file='volcano_run.in',status='old',action='read')
       read(u,nml=volcano_run,iostat=ios,iomsg=message)
       close(u)
       if(ios/=0) then
          print *,trim(message)
          error stop 'Invalid volcano_run.in namelist'
       end if
    end if
    if(.not.all(ieee_is_finite([spinup_max_years,volcano_duration_day, &
         equilibrium_relative_limit,equilibrium_mixing_floor]))) error stop 'Nonfinite run setting'
    if(min(spinup_max_years,volcano_duration_day,equilibrium_relative_limit, &
         equilibrium_mixing_floor)<=0d0.or.stable_days_required<1) error stop 'Invalid run setting'
    if(spinup_max_years*365d0<stable_days_required) error stop 'Spin-up shorter than stability window'
    print '(A,F10.3)', 'Maximum background spin-up [year]: ',spinup_max_years
    print '(A,F10.3)', 'Volcanic simulation duration [day]: ',volcano_duration_day
  end subroutine patmo_volc_readRunSettings

  subroutine patmo_volc_getRunLimits(spinupEnd,volcanoEnd)
    use patmo_constants, only: secondsPerDay
    real*8,intent(out)::spinupEnd,volcanoEnd
    spinupEnd=secondsPerDay*365d0*spinup_max_years
    volcanoEnd=secondsPerDay*volcano_duration_day
  end subroutine patmo_volc_getRunLimits

  subroutine patmo_volc_beginEquilibriumCheck()
    use patmo_parameters, only: nall
    background_previous(:,:)=nall(:,1:chemSpeciesNumber)
    stable_days=0
  end subroutine patmo_volc_beginEquilibriumCheck

  function patmo_volc_checkEquilibrium(change) result(equilibrated)
    use patmo_parameters, only: nall
    real*8,intent(out)::change
    real*8::air_scale(cellsNumber,chemSpeciesNumber)
    logical::equilibrated
    ! Preserve the daily relative-change metric, including the trace floor.
    air_scale=spread(max(0.5d0*sum(nall(:,1:chemSpeciesNumber),2),1d-99), &
         2,chemSpeciesNumber)*equilibrium_mixing_floor
    change=maxval(abs(nall(:,1:chemSpeciesNumber)-background_previous) &
         /max(abs(nall(:,1:chemSpeciesNumber)),abs(background_previous),air_scale))
    background_previous(:,:)=nall(:,1:chemSpeciesNumber)
    stable_days=merge(stable_days+1,0,change<equilibrium_relative_limit)
    equilibrated=stable_days>=stable_days_required
  end function patmo_volc_checkEquilibrium

  subroutine patmo_volc_requireEquilibrium(equilibrated)
    logical,intent(in)::equilibrated
    if(equilibrated) return
    print *,"ERROR: background atmosphere did not reach steady state."
    print *,"       Volcano forcing was not started."
    print *,"       Inspect unconverged species and extend spinup_tend before proceeding."
    error stop 1
  end subroutine patmo_volc_requireEquilibrium

  !======================================================================
  ! Volcano clock and event state
  !======================================================================

  function patmo_volc_isActive() result(active)
    logical::active
    active=volcanoEnabled.and.volcanoClockRunning
  end function patmo_volc_isActive

  subroutine patmo_volc_setEarlyOutput(stepSeconds,untilSeconds)
    real*8,intent(in)::stepSeconds,untilSeconds
    earlyOutputStep=max(stepSeconds,0d0)
    earlyOutputUntil=max(untilSeconds,0d0)
  end subroutine patmo_volc_setEarlyOutput

  function patmo_volc_nextOutput(t,firstTime,dt) result(nextTime)
    real*8,intent(in)::t,firstTime,dt
    real*8::nextTime
    nextTime=t+dt
    if(earlyOutputStep>0d0.and.t<firstTime+earlyOutputUntil) &
         nextTime=min(nextTime,t+earlyOutputStep,firstTime+earlyOutputUntil)
  end function patmo_volc_nextOutput

  !***************
  subroutine patmo_volc_reset()
    implicit none

    if(allocated(events)) deallocate(events)
    if(allocated(ashLayerTau550)) deallocate(ashLayerTau550)
    if(allocated(sulfateLayerTau550)) deallocate(sulfateLayerTau550)
    if(allocated(sulfatePrecursor)) deallocate(sulfatePrecursor)
    volcanoEnabled = .false.
    volcanoClockRunning = .true.
    ashStateReady = .false.
    warnedMissingSO2 = .false.
    hasSO2Source = .false.
    hasAshOpacity = .false.
    hasSulfateOpacity = .false.
    volcanoEventsNumber = 0
    so2IndexCache = -2
    volcanoElapsedTime = 0d0
    ashStateTime = 0d0

  end subroutine patmo_volc_reset

  !***************
  subroutine patmo_volc_setEnabled(enabled)
    implicit none
    logical,intent(in)::enabled

    volcanoEnabled = enabled

  end subroutine patmo_volc_setEnabled

  !***************
  subroutine patmo_volc_setClockRunning(running)
    implicit none
    logical,intent(in)::running

    volcanoClockRunning = running

  end subroutine patmo_volc_setClockRunning

  !***************
  subroutine patmo_volc_setTime(timeSeconds)
    implicit none
    real*8,intent(in)::timeSeconds
    real*8::newTime

    newTime = max(timeSeconds,0d0)
    if(ashStateReady.and.newTime<ashStateTime) call patmo_volc_resetAshState()
    volcanoElapsedTime = newTime

  end subroutine patmo_volc_setTime

  !***************
  function patmo_volc_getTime()
    implicit none
    real*8::patmo_volc_getTime

    patmo_volc_getTime = volcanoElapsedTime

  end function patmo_volc_getTime

  !***************
  subroutine patmo_volc_advanceTime(dt)
    implicit none
    real*8,intent(in)::dt
    real*8::newTime

    if(.not.volcanoClockRunning) return
    newTime = max(volcanoElapsedTime + dt,0d0)
    if(newTime<volcanoElapsedTime) call patmo_volc_resetAshState()
    volcanoElapsedTime = newTime
    call patmo_volc_updateAshState(volcanoElapsedTime)

  end subroutine patmo_volc_advanceTime

  !***************
  !Keep volcanic forcing inactive while the background atmosphere spins up.
  !After calling this, start_day in the event file is measured from the later
  !call to patmo_volc_startAfterEquilibrium().
  subroutine patmo_volc_holdUntilEquilibrium()
    implicit none

    volcanoEnabled = .false.
    volcanoClockRunning = .false.
    volcanoElapsedTime = 0d0
    call patmo_volc_resetAshState()

  end subroutine patmo_volc_holdUntilEquilibrium

  !***************
  !Release the volcanic events after the background atmosphere has equilibrated.
  !This resets the volcanic clock, so start_day=0 means immediate post-spinup
  !eruption.
  subroutine patmo_volc_startAfterEquilibrium()
    implicit none

    volcanoEnabled = hasSO2Source .or. hasAshOpacity .or. hasSulfateOpacity
    volcanoClockRunning = .true.
    volcanoElapsedTime = 0d0
    call patmo_volc_resetAshState()
    if(volcanoEnabled) then
       print *,"Volcano forcing starts after background equilibrium."
       print *,"  volcano_time_day=0; event start_day values are relative to this time."
    end if

  end subroutine patmo_volc_startAfterEquilibrium

  !***************
  ! Preferred event file format, one event per line:
  ! event_id=name start_day=... duration_day=... plume_center_km=...
  ! plume_sigma_km=... so2_flux_cm2_s=...
  ! ash_tau_550=... sulfate_tau_550=...
  ! ash_lifetime_day=... sulfate_lifetime_day=...
  ! ash_horizontal_lifetime_day=... sulfate_horizontal_lifetime_day=...
  ! ash_rise_start_km=... ash_rise_time_hour=...
  ! sulfate_formation_day=...
  ! ash_particle_radius_um=... ash_particle_density_g_cm3=...
  ! sulfate_particle_radius_um=... sulfate_particle_density_g_cm3=...
  ! ash_vertical_diffusion_cm2_s=...
  ! sulfate_vertical_diffusion_cm2_s=...
  ! ash_settling_cm_s=... sulfate_settling_cm_s=...
  ! ash_lambda_exponent=... sulfate_lambda_exponent=...
  ! If patmo_volc_holdUntilEquilibrium()/startAfterEquilibrium() are used,
  ! start_day is measured from the post-equilibrium volcano clock.
  ! Convenience inputs such as so2_column_cm2, so2_mass_tg with
  ! plume_radius_km/injection_area_km2, plume_top_km/plume_bottom_km,
  ! plume_fwhm_km, and ash_settling_km_day are also accepted.
  ! Legacy numeric rows are still accepted in the same column order.

  !======================================================================
  ! Event input and forcing parameters (volcano_events.dat)
  !======================================================================
  subroutine patmo_volc_loadEvents(fname)
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname
    character(len=4096)::line
    integer::ios,unitEvent,i,commentPos,lineNumber
    real*8::startDay,durationDay,centerKm,sigmaKm
    real*8::so2ColumnFlux,ashTau550,ashLifetimeDay
    real*8::ashHorizontalLifetimeDay,ashSettling,ashWavelengthExp
    real*8::ashRiseStartKm,ashRiseTimeDay
    real*8::ashVerticalDiffusion,ashParticleRadiusUm
    real*8::ashParticleDensity,ashAirViscosity,ashGravity
    real*8::ashCunninghamFactor
    real*8::sulfateTau550,sulfateLifetimeDay
    real*8::sulfateHorizontalLifetimeDay,sulfateSettling
    real*8::sulfateWavelengthExp,sulfateRiseStartKm
    real*8::sulfateRiseTimeDay,sulfateVerticalDiffusion
    real*8::sulfateParticleRadiusUm,sulfateParticleDensity
    real*8::sulfateAirViscosity,sulfateGravity
    real*8::sulfateCunninghamFactor,sulfateStartDelayDay
    real*8::sulfateFormationTimeDay

    call patmo_volc_reset()

    unitEvent = 91
    open(unitEvent,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problem while opening volcano event file ",trim(fname)
       stop
    end if

    volcanoEventsNumber = 0
    do
       read(unitEvent,'(A)',iostat=ios) line
       if(ios/=0) exit
       commentPos = patmo_volc_commentStart(line)
       if(commentPos>0) line = line(:commentPos-1)
       if(len_trim(line)==0) cycle
       volcanoEventsNumber = volcanoEventsNumber + 1
    end do
    close(unitEvent)

    if(volcanoEventsNumber<=0) then
       print *,"WARNING: volcano event file has no active event: ",trim(fname)
       return
    end if

    allocate(events(volcanoEventsNumber))
    open(unitEvent,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problem while reopening volcano event file ",trim(fname)
       stop
    end if

    i = 0
    lineNumber = 0
    do
       read(unitEvent,'(A)',iostat=ios) line
       if(ios/=0) exit
       lineNumber = lineNumber + 1
       commentPos = patmo_volc_commentStart(line)
       if(commentPos>0) line = line(:commentPos-1)
       if(len_trim(line)==0) cycle

       startDay = 0d0
       durationDay = 0d0
       centerKm = 0d0
       sigmaKm = 1d0
       so2ColumnFlux = 0d0
       ashTau550 = 0d0
       ashLifetimeDay = 0d0
       ashHorizontalLifetimeDay = 0d0
       ashSettling = 0d0
       ashWavelengthExp = 0d0
       ashRiseStartKm = 0d0
       ashRiseTimeDay = -1d0
       ashVerticalDiffusion = 0d0
       ashParticleRadiusUm = 0d0
       ashParticleDensity = 2.5d0
       ashAirViscosity = 1.7d-4
       ashGravity = 980.665d0
       ashCunninghamFactor = 1d0
       sulfateTau550 = 0d0
       sulfateLifetimeDay = 0d0
       sulfateHorizontalLifetimeDay = 0d0
       sulfateSettling = 0d0
       sulfateWavelengthExp = 0d0
       sulfateRiseStartKm = -huge(1d0)
       sulfateRiseTimeDay = -1d0
       sulfateVerticalDiffusion = 0d0
       sulfateParticleRadiusUm = 0d0
       sulfateParticleDensity = 1.6d0
       sulfateAirViscosity = 1.7d-4
       sulfateGravity = 980.665d0
       sulfateCunninghamFactor = 1d0
       sulfateStartDelayDay = 0d0
       sulfateFormationTimeDay = -1d0

       if(index(line,"=")>0) then
          call patmo_volc_parseKeywordEvent(line,startDay,durationDay, &
               centerKm,sigmaKm,so2ColumnFlux,ashTau550,ashLifetimeDay, &
               ashHorizontalLifetimeDay,ashSettling,ashWavelengthExp, &
               ashRiseStartKm,ashRiseTimeDay,ashVerticalDiffusion, &
               ashParticleRadiusUm,ashParticleDensity,ashAirViscosity, &
               ashGravity,ashCunninghamFactor,sulfateTau550, &
               sulfateLifetimeDay,sulfateHorizontalLifetimeDay, &
               sulfateSettling,sulfateWavelengthExp,sulfateRiseStartKm, &
               sulfateRiseTimeDay,sulfateVerticalDiffusion, &
               sulfateParticleRadiusUm,sulfateParticleDensity, &
               sulfateAirViscosity,sulfateGravity,sulfateCunninghamFactor, &
               sulfateStartDelayDay,sulfateFormationTimeDay,ios)
       else
          read(line,*,iostat=ios) startDay,durationDay,centerKm,sigmaKm, &
               so2ColumnFlux,ashTau550,ashLifetimeDay,ashSettling, &
               ashWavelengthExp
       end if
       if(ios/=0) then
          print *,"ERROR: malformed volcano event row at line ",lineNumber
          print *,trim(line)
          stop
       end if

       i = i + 1
       if(ashRiseTimeDay<0d0) ashRiseTimeDay = durationDay
       if(sulfateRiseStartKm<-0.5d0*huge(1d0)) then
          sulfateRiseStartKm = ashRiseStartKm
       end if
       if(sulfateRiseTimeDay<0d0) sulfateRiseTimeDay = ashRiseTimeDay
       if(sulfateFormationTimeDay<0d0) sulfateFormationTimeDay = durationDay
       if(ashSettling<=0d0.and.ashParticleRadiusUm>0d0) then
          ashSettling = patmo_volc_stokesSettling(ashParticleRadiusUm*1d-4, &
               ashParticleDensity,ashAirViscosity,ashGravity, &
               ashCunninghamFactor)
       end if
       if(sulfateSettling<=0d0.and.sulfateParticleRadiusUm>0d0) then
          sulfateSettling = patmo_volc_stokesSettling(sulfateParticleRadiusUm*1d-4, &
               sulfateParticleDensity,sulfateAirViscosity,sulfateGravity, &
               sulfateCunninghamFactor)
       end if
       events(i)%startTime = startDay * secondsPerDay
       events(i)%duration = max(durationDay,0d0) * secondsPerDay
       events(i)%center = centerKm * 1d5
       events(i)%sigma = max(abs(sigmaKm) * 1d5,1d0)
       events(i)%so2ColumnFlux = max(so2ColumnFlux,0d0)
       events(i)%ashTau550 = max(ashTau550,0d0)
       events(i)%ashLifetime = ashLifetimeDay * secondsPerDay
       events(i)%ashHorizontalLifetime = ashHorizontalLifetimeDay * secondsPerDay
       events(i)%ashSettling = max(ashSettling,0d0)
       events(i)%ashWavelengthExp = ashWavelengthExp
       events(i)%ashRiseStart = ashRiseStartKm * 1d5
       events(i)%ashRiseTime = max(ashRiseTimeDay,0d0) * secondsPerDay
       events(i)%ashVerticalDiffusion = ashVerticalDiffusion
       events(i)%ashParticleRadius = max(ashParticleRadiusUm,0d0) * 1d-4
       events(i)%ashParticleDensity = max(ashParticleDensity,0d0)
       events(i)%ashAirViscosity = max(ashAirViscosity,0d0)
       events(i)%ashGravity = max(ashGravity,0d0)
       events(i)%ashCunninghamFactor = max(ashCunninghamFactor,0d0)
       events(i)%sulfateTau550 = max(sulfateTau550,0d0)
       events(i)%sulfateLifetime = sulfateLifetimeDay * secondsPerDay
       events(i)%sulfateHorizontalLifetime = sulfateHorizontalLifetimeDay &
            * secondsPerDay
       events(i)%sulfateSettling = max(sulfateSettling,0d0)
       events(i)%sulfateWavelengthExp = sulfateWavelengthExp
       events(i)%sulfateRiseStart = sulfateRiseStartKm * 1d5
       events(i)%sulfateRiseTime = max(sulfateRiseTimeDay,0d0) * secondsPerDay
       events(i)%sulfateVerticalDiffusion = sulfateVerticalDiffusion
       events(i)%sulfateParticleRadius = max(sulfateParticleRadiusUm,0d0) * 1d-4
       events(i)%sulfateParticleDensity = max(sulfateParticleDensity,0d0)
       events(i)%sulfateAirViscosity = max(sulfateAirViscosity,0d0)
       events(i)%sulfateGravity = max(sulfateGravity,0d0)
       events(i)%sulfateCunninghamFactor = max(sulfateCunninghamFactor,0d0)
       events(i)%sulfateStartDelay = max(sulfateStartDelayDay,0d0) &
            * secondsPerDay
       events(i)%sulfateFormationTime = max(sulfateFormationTimeDay,0d0) &
            * secondsPerDay
       if(events(i)%so2ColumnFlux>0d0) hasSO2Source = .true.
       if(events(i)%ashTau550>0d0) hasAshOpacity = .true.
       if(events(i)%sulfateTau550>0d0) hasSulfateOpacity = .true.
    end do
    close(unitEvent)

    volcanoEnabled = hasSO2Source .or. hasAshOpacity .or. hasSulfateOpacity
    if(hasSO2Source) i = patmo_volc_getSO2Index(required=.true.)
    if(.not.volcanoEnabled) then
       print *,"WARNING: volcano event file has no positive SO2 flux, ash tau, or sulfate tau: ",trim(fname)
    end if

    write(*,'(A)') "Loaded volcano events"
    write(*,'(A,1X,A)') "  file: ",trim(fname)
    write(*,'(A,I0)') "  count: ",volcanoEventsNumber
    do i=1,volcanoEventsNumber
       write(*,'(A)') ""
       write(*,'(A,I0,A)') "  Event ",i,":"
       write(*,'(A,F12.4)') "    start_day:                ", &
            events(i)%startTime/secondsPerDay
       write(*,'(A,F12.4)') "    duration_day:             ", &
            events(i)%duration/secondsPerDay
       write(*,'(A,F12.4)') "    plume_center_km:          ", &
            events(i)%center/1d5
       write(*,'(A,F12.4)') "    plume_sigma_km:           ", &
            events(i)%sigma/1d5
       write(*,'(A,E14.6E3)') "    so2_flux_cm2_s:           ", &
            events(i)%so2ColumnFlux
       write(*,'(A,F12.4)') "    ash_tau_550:              ", &
            events(i)%ashTau550
       write(*,'(A,F12.4)') "    ash_lifetime_day:         ", &
            events(i)%ashLifetime/secondsPerDay
       write(*,'(A,F12.4)') "    ash_horizontal_lifetime_day:", &
            events(i)%ashHorizontalLifetime/secondsPerDay
       write(*,'(A,E14.6E3)') "    ash_settling_cm_s:        ", &
            events(i)%ashSettling
       write(*,'(A,F12.4)') "    ash_settling_km_day:      ", &
            events(i)%ashSettling*secondsPerDay/1d5
       write(*,'(A,F12.4)') "    ash_lambda_exponent:      ", &
            events(i)%ashWavelengthExp
       write(*,'(A,F12.4)') "    ash_rise_start_km:        ", &
            events(i)%ashRiseStart/1d5
       write(*,'(A,F12.4)') "    ash_rise_time_day:        ", &
            events(i)%ashRiseTime/secondsPerDay
       write(*,'(A,E14.6E3)') "    ash_vertical_diffusion_cm2_s:", &
            events(i)%ashVerticalDiffusion
       if(events(i)%ashVerticalDiffusion<0d0) write(*,'(A)') "      negative value = background Kzz from profile.dat"
       write(*,'(A,F12.4)') "    ash_particle_radius_um:   ", &
            events(i)%ashParticleRadius/1d-4
       write(*,'(A,F12.4)') "    ash_particle_density_g_cm3:", &
            events(i)%ashParticleDensity
       write(*,'(A,F12.4)') "    ash_gravity_cm_s2:        ", &
            events(i)%ashGravity
       write(*,'(A,F12.4)') "    sulfate_tau_550:          ", &
            events(i)%sulfateTau550
       write(*,'(A,F12.4)') "    sulfate_formation_efold_day:", &
            events(i)%sulfateFormationTime/secondsPerDay
       write(*,'(A,F12.4)') "    sulfate_start_delay_day:  ", &
            events(i)%sulfateStartDelay/secondsPerDay
       write(*,'(A,F12.4)') "    sulfate_lifetime_day:     ", &
            events(i)%sulfateLifetime/secondsPerDay
       write(*,'(A,F12.4)') "    sulfate_horizontal_lifetime_day:", &
            events(i)%sulfateHorizontalLifetime/secondsPerDay
       write(*,'(A,E14.6E3)') "    sulfate_settling_cm_s:    ", &
            events(i)%sulfateSettling
       write(*,'(A,F12.4)') "    sulfate_settling_km_day:  ", &
            events(i)%sulfateSettling*secondsPerDay/1d5
       write(*,'(A,F12.4)') "    sulfate_lambda_exponent:  ", &
            events(i)%sulfateWavelengthExp
       write(*,'(A,F12.4)') "    sulfate_rise_start_km:    ", &
            events(i)%sulfateRiseStart/1d5
       write(*,'(A,F12.4)') "    sulfate_rise_time_day:    ", &
            events(i)%sulfateRiseTime/secondsPerDay
       write(*,'(A,E14.6E3)') "    sulfate_vertical_diffusion_cm2_s:", &
            events(i)%sulfateVerticalDiffusion
       if(events(i)%sulfateVerticalDiffusion<0d0) write(*,'(A)') "      negative value = background Kzz from profile.dat"
       write(*,'(A,F12.4)') "    sulfate_particle_radius_um:", &
            events(i)%sulfateParticleRadius/1d-4
       write(*,'(A,F12.4)') "    sulfate_particle_density_g_cm3:", &
            events(i)%sulfateParticleDensity
       write(*,'(A,F12.4)') "    sulfate_gravity_cm_s2:    ", &
            events(i)%sulfateGravity
    end do
    write(*,'(A)') ""

  end subroutine patmo_volc_loadEvents

  !***************

  !======================================================================
  ! Chemical SO2 sources
  !======================================================================
  subroutine patmo_volc_injectSO2Column(centerKm,sigmaKm,columnSO2)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::centerKm,sigmaKm,columnSO2
    real*8::w(cellsNumber),norm,center,sigma
    integer::j,idxSO2

    if(columnSO2<=0d0) return

    idxSO2 = patmo_volc_getSO2Index(required=.true.)
    center = centerKm * 1d5
    sigma = max(abs(sigmaKm)*1d5,1d0)
    call patmo_volc_layerWeights(center,sigma,w,norm)
    if(norm<=0d0) return

    do j=1,cellsNumber
       nall(j,idxSO2) = nall(j,idxSO2) + columnSO2 * w(j) / norm
    end do

  end subroutine patmo_volc_injectSO2Column

  !***************
  subroutine patmo_volc_addSources(tlocal,n,dn)
    use patmo_commons
    implicit none
    real*8,intent(in)::tlocal
    real*8,intent(in)::n(cellsNumber,speciesNumber)
    real*8,intent(inout)::dn(cellsNumber,speciesNumber)
    real*8::absoluteTime,w(cellsNumber),norm
    integer::i,j,idxSO2

    if(.not.volcanoEnabled) return
    if(volcanoEventsNumber<=0) return
    if(.not.hasSO2Source) return

    idxSO2 = patmo_volc_getSO2Index(required=.false.)
    if(idxSO2<0) then
       if(.not.warnedMissingSO2) then
          print *,"WARNING: volcano SO2 source requested, but SO2 is absent."
          warnedMissingSO2 = .true.
       end if
       return
    end if

    absoluteTime = volcanoElapsedTime + max(tlocal,0d0)
    do i=1,volcanoEventsNumber
       if(.not.patmo_volc_so2Active(i,absoluteTime)) cycle
       if(events(i)%so2ColumnFlux<=0d0) cycle
       call patmo_volc_layerWeights(events(i)%center,events(i)%sigma,w,norm)
       if(norm<=0d0) cycle
       do j=1,cellsNumber
          dn(j,idxSO2) = dn(j,idxSO2) + events(i)%so2ColumnFlux * w(j) / norm
       end do
    end do

  end subroutine patmo_volc_addSources

  !***************

  !======================================================================
  ! Particle optical shielding and pre-run data writers
  !======================================================================
  subroutine patmo_volc_applyAshOpacity(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(inout)::tau(photoBinsNumber,cellsNumber)

    call patmo_volc_addParticleOpacity(tau,volcComponentTotal)

  end subroutine patmo_volc_applyAshOpacity

  !***************
  subroutine patmo_volc_addParticleOpacity(tau,component)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(inout)::tau(photoBinsNumber,cellsNumber)
    integer,intent(in)::component
    real*8::absoluteTime,columnTau(photoBinsNumber)
    real*8::spectralScale(photoBinsNumber),layerTau550
    integer::i,j,k
    logical::includeAsh,includeSulfate

    if(.not.volcanoEnabled) return
    if(volcanoEventsNumber<=0) return
    if(.not.(hasAshOpacity.or.hasSulfateOpacity)) return
    includeAsh = component==volcComponentAsh .or. component==volcComponentTotal
    includeSulfate = component==volcComponentSulfate &
         .or. component==volcComponentTotal

    absoluteTime = volcanoElapsedTime
    call patmo_volc_updateAshState(absoluteTime)

    if(includeAsh.and.hasAshOpacity.and.allocated(ashLayerTau550)) then
       do k=1,volcanoEventsNumber
          if(events(k)%ashTau550<=0d0) cycle
          do i=1,photoBinsNumber
             spectralScale(i) = patmo_volc_particleSpectralScale(i, &
                  events(k)%ashWavelengthExp)
          end do
          columnTau(:) = 0d0
          do j=cellsNumber,1,-1
             layerTau550 = max(ashLayerTau550(j,k),0d0)
             columnTau(:) = columnTau(:) + layerTau550 * spectralScale(:)
             tau(:,j) = tau(:,j) + columnTau(:)
          end do
       end do
    end if

    if(includeSulfate.and.hasSulfateOpacity.and.allocated(sulfateLayerTau550)) then
       do k=1,volcanoEventsNumber
          if(events(k)%sulfateTau550<=0d0) cycle
          do i=1,photoBinsNumber
             spectralScale(i) = patmo_volc_particleSpectralScale(i, &
                  events(k)%sulfateWavelengthExp)
          end do
          columnTau(:) = 0d0
          do j=cellsNumber,1,-1
             layerTau550 = max(sulfateLayerTau550(j,k),0d0)
             columnTau(:) = columnTau(:) + layerTau550 * spectralScale(:)
             tau(:,j) = tau(:,j) + columnTau(:)
          end do
       end do
    end if

  end subroutine patmo_volc_addParticleOpacity

  !***************
  !Dump only the volcanic ash plus sulfate-aerosol optical depth.
  !The regular gas opacity is not included, and the chemistry solver is not run.
  !endAfterStart<0 selects an automatic end time based on tauFloor.
  !The total tau is the shielding used by the photolysis solver; ash and sulfate
  !components are also written for diagnosis.
  subroutine patmo_volc_dumpOpticalDepth(fname,timeStep,timeUnitSeconds, &
       timeUnitName,wavelengthStepNm,tauFloor,endAfterStart)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname,timeUnitName
    real*8,intent(in)::timeStep,timeUnitSeconds,wavelengthStepNm
    real*8,intent(in)::tauFloor,endAfterStart
    real*8::tauTotal(photoBinsNumber,cellsNumber)
    real*8::tauAsh(photoBinsNumber,cellsNumber)
    real*8::tauSulfate(photoBinsNumber,cellsNumber)
    real*8::originalTime,firstTime,lastTime,nextTime,t,dtOut,tauStop
    real*8::timeScale
    real*8,allocatable::targetWavelength(:),modelWavelength(:)
    integer,allocatable::selectedBins(:)
    integer::unitOut,wavelengthsNumber,iw,j,ibin
    character(len=24)::wavelengthText

    dtOut = max(timeStep,1d0)
    timeScale = max(timeUnitSeconds,1d-99)
    tauStop = max(tauFloor,1d-12)
    firstTime = patmo_volc_firstAshStart()
    if(endAfterStart>0d0) then
       lastTime = firstTime + endAfterStart
    else
       lastTime = patmo_volc_estimateOpticalDepthEnd(tauStop)
    end if
    if(lastTime<firstTime) lastTime = firstTime

    call patmo_volc_selectWavelengths(wavelengthStepNm,selectedBins, &
         targetWavelength,modelWavelength,wavelengthsNumber)

    originalTime = volcanoElapsedTime
    unitOut = 93
    open(unitOut,file=trim(fname),status="replace")
    write(unitOut,'(A)') "# volcanic ash plus sulfate-aerosol optical depth pre-run"
    write(unitOut,'(A)') "# gas opacity and chemistry solver are not included"
    write(unitOut,'(A)') "# ash and sulfate are stored as separate per-layer tau550 budgets"
    write(unitOut,'(A)') "# sulfate forms by first-order conversion; unformed precursor has no opacity or settling"
    write(unitOut,'(A)') "# cumulative tau includes the local layer and every layer above, including the top layer"
    write(unitOut,'(A,2E17.8E3)') "# early_output_step_s until_s ",earlyOutputStep,earlyOutputUntil
    write(unitOut,'(A)') "# total_tau columns are the shielding used by photolysis"
    write(unitOut,'(A)') "# ash_tau and sulfate_tau columns show each particle component"
    write(unitOut,'(A)') "# continuous vent sources; height-limited ascent lasts duration + rise_time"
    write(unitOut,'(A,E17.8E3)') "# output_time_step_s ",dtOut
    write(unitOut,'(A,A)') "# time_unit ",trim(timeUnitName)
    write(unitOut,'(A)') "# time is measured since the first volcanic event"
    write(unitOut,'(A,E17.8E3)') "# requested_wavelength_step_nm ", &
         wavelengthStepNm
    write(unitOut,'(A)') "# wavelengths are nearest PATMO model wavelength bins"
    write(unitOut,'(A,E17.8E3)') "# tau_floor_for_auto_end ",tauStop
    write(unitOut,'(A)',advance="no") "# wavelength_nm:"
    do iw=1,wavelengthsNumber
       write(unitOut,'(1X,E14.6E3)',advance="no") modelWavelength(iw)
    end do
    write(unitOut,*)
    write(unitOut,'(A)',advance="no") "# columns: time layer altitude_km"
    do iw=1,wavelengthsNumber
       write(wavelengthText,'(F12.3)') modelWavelength(iw)
       wavelengthText = adjustl(wavelengthText)
       write(unitOut,'(A)',advance="no") " total_tau_"//trim(wavelengthText)//"nm"
       write(unitOut,'(A)',advance="no") " ash_tau_"//trim(wavelengthText)//"nm"
       write(unitOut,'(A)',advance="no") " sulfate_tau_"//trim(wavelengthText)//"nm"
    end do
    write(unitOut,*)

    t = firstTime
    do
       call patmo_volc_setTime(t)
       tauTotal(:,:) = 0d0
       tauAsh(:,:) = 0d0
       tauSulfate(:,:) = 0d0
       call patmo_volc_addParticleOpacity(tauTotal(:,:),volcComponentTotal)
       call patmo_volc_addParticleOpacity(tauAsh(:,:),volcComponentAsh)
       call patmo_volc_addParticleOpacity(tauSulfate(:,:),volcComponentSulfate)
       do j=1,cellsNumber
          write(unitOut,'(E17.8E3,I6,E14.6E3)',advance="no") &
               (t-firstTime)/timeScale,j,height(j)/1d5
          do iw=1,wavelengthsNumber
             ibin = selectedBins(iw)
             write(unitOut,'(1X,E17.8E3)',advance="no") tauTotal(ibin,j)
             write(unitOut,'(1X,E17.8E3)',advance="no") tauAsh(ibin,j)
             write(unitOut,'(1X,E17.8E3)',advance="no") tauSulfate(ibin,j)
          end do
          write(unitOut,*)
       end do
       write(unitOut,*)

       if(t>=lastTime) exit
       nextTime = patmo_volc_nextOutput(t,firstTime,dtOut)
       if(nextTime>lastTime) nextTime = lastTime
       if(nextTime<=t) exit
       t = nextTime
    end do
    close(unitOut)

    call patmo_volc_setTime(originalTime)
    if(allocated(selectedBins)) deallocate(selectedBins)
    if(allocated(targetWavelength)) deallocate(targetWavelength)
    if(allocated(modelWavelength)) deallocate(modelWavelength)

    write(*,'(A)') "Pre-run output"
    write(*,'(A,1X,A)') "  file:                ",trim(fname)
    write(*,'(A,F12.4)') "  first_eruption_day:  ", &
         firstTime/secondsPerDay
    write(*,'(A,F12.4)') "  final_volcano_day:   ", &
         lastTime/secondsPerDay
    write(*,'(A,F12.4,1X,A)') "  output_time_step:    ", &
         dtOut/timeScale,trim(timeUnitName)
    write(*,'(A,F12.4)') "  output_time_step_day:", &
         dtOut/secondsPerDay
    write(*,'(A,I0)') "  selected_wavelengths: ",wavelengthsNumber
    write(*,'(A)') ""

  end subroutine patmo_volc_dumpOpticalDepth

  !***************
  subroutine patmo_volc_computeLocalState(timeSeconds,ashLocal, &
       ashColumn,sulfateLocal,sulfateColumn,totalLocal,totalColumn)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::timeSeconds
    real*8,intent(out)::ashLocal(cellsNumber)
    real*8,intent(out)::ashColumn(cellsNumber)
    real*8,intent(out)::sulfateLocal(cellsNumber),sulfateColumn(cellsNumber)
    real*8,intent(out)::totalLocal(cellsNumber),totalColumn(cellsNumber)
    real*8::absoluteTime
    integer::j,k

    ashLocal(:) = 0d0
    ashColumn(:) = 0d0
    sulfateLocal(:) = 0d0
    sulfateColumn(:) = 0d0
    totalLocal(:) = 0d0
    totalColumn(:) = 0d0
    absoluteTime = timeSeconds

    if(volcanoEnabled.and.volcanoEventsNumber>0) then
       if(hasAshOpacity.or.hasSulfateOpacity) then
          call patmo_volc_updateAshState(absoluteTime)
          if(hasAshOpacity.and.allocated(ashLayerTau550)) then
             do k=1,volcanoEventsNumber
                do j=1,cellsNumber
                   ashLocal(j) = ashLocal(j) + max(ashLayerTau550(j,k),0d0)
                end do
             end do
          end if
          if(hasSulfateOpacity.and.allocated(sulfateLayerTau550)) then
             do k=1,volcanoEventsNumber
                do j=1,cellsNumber
                   sulfateLocal(j) = sulfateLocal(j) &
                        + max(sulfateLayerTau550(j,k),0d0)
                end do
             end do
          end if
       end if
    end if

    totalLocal(:) = ashLocal(:) + sulfateLocal(:)
    ashColumn(cellsNumber) = ashLocal(cellsNumber)
    sulfateColumn(cellsNumber) = sulfateLocal(cellsNumber)
    totalColumn(cellsNumber) = totalLocal(cellsNumber)
    do j=cellsNumber-1,1,-1
       ashColumn(j) = ashColumn(j+1) + ashLocal(j)
       sulfateColumn(j) = sulfateColumn(j+1) + sulfateLocal(j)
       totalColumn(j) = totalColumn(j+1) + totalLocal(j)
    end do

  end subroutine patmo_volc_computeLocalState

  !***************
  subroutine patmo_volc_computeSO2SourceProfile(timeSeconds,source)
    use patmo_parameters
    implicit none
    real*8,intent(in)::timeSeconds
    real*8,intent(out)::source(cellsNumber)
    real*8::w(cellsNumber),norm
    integer::i,j

    source(:) = 0d0
    if(.not.hasSO2Source) return
    do i=1,volcanoEventsNumber
       if(.not.patmo_volc_so2Active(i,timeSeconds)) cycle
       call patmo_volc_layerWeights(events(i)%center,events(i)%sigma,w,norm)
       if(norm<=0d0) cycle
       do j=1,cellsNumber
          source(j) = source(j) + events(i)%so2ColumnFlux * w(j) / norm
       end do
    end do

  end subroutine patmo_volc_computeSO2SourceProfile

  !***************
  !Dump time-varying local ash and sulfate plume strength.
  !*_local_tau550 is each layer's own 550 nm optical-depth contribution,
  !not the cumulative overlying optical depth used by radiative shielding.
  !Ash and sulfate are tracked as separate particle budgets.
  subroutine patmo_volc_dumpAshProfile(fname,timeStep,timeUnitSeconds, &
       timeUnitName,tauFloor,endAfterStart)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname,timeUnitName
    real*8,intent(in)::timeStep,timeUnitSeconds,tauFloor,endAfterStart
    real*8::ashLocal(cellsNumber),ashColumn(cellsNumber)
    real*8::sulfateLocal(cellsNumber),sulfateColumn(cellsNumber)
    real*8::totalLocal(cellsNumber),totalColumn(cellsNumber),precursorLocal(cellsNumber)
    real*8::originalTime,firstTime,lastTime,nextTime,t,dtOut,tauStop
    real*8::timeScale
    integer::j,unitOut

    dtOut = max(timeStep,1d0)
    timeScale = max(timeUnitSeconds,1d-99)
    tauStop = max(tauFloor,1d-12)
    firstTime = patmo_volc_firstAshStart()
    if(endAfterStart>0d0) then
       lastTime = firstTime + endAfterStart
    else
       lastTime = patmo_volc_estimateOpticalDepthEnd(tauStop)
    end if
    if(lastTime<firstTime) lastTime = firstTime

    originalTime = volcanoElapsedTime
    unitOut = 95
    open(unitOut,file=trim(fname),status="replace")
    write(unitOut,'(A)') "# volcanic ash and sulfate aerosol local plume pre-run"
    write(unitOut,'(A)') "# chemistry solver is not included"
    write(unitOut,'(A)') "# ash and sulfate use separate per-layer tau550 budgets"
    write(unitOut,'(A)') "# continuous vent sources; height-limited ascent lasts duration + rise_time"
    write(unitOut,'(A)') "# rise_time sets nominal transit speed, not the duration of emission"
    write(unitOut,'(A)') "# formed SSA shares precursor airflow, but forms locally and settles independently"
    write(unitOut,'(A)') "# ash_local_tau550 is the per-layer 550 nm optical-depth contribution"
    write(unitOut,'(A)') "# sulfate_local_tau550 is the per-layer 550 nm sulfate-aerosol contribution"
    write(unitOut,'(A)') "# sulfate_precursor_potential is unconverted tau550-equivalent potential, NOT opacity"
    write(unitOut,'(A)') "# precursor mixes and rises with air but does not settle; only formed sulfate settles"
    write(unitOut,'(A)') "# formation_day is a first-order e-folding time, not time to 100 percent conversion"
    write(unitOut,'(A,2E17.8E3)') "# early_output_step_s until_s ",earlyOutputStep,earlyOutputUntil
    write(unitOut,'(A)') "# total_local_tau550 = ash_local_tau550 + sulfate_local_tau550"
    write(unitOut,'(A)') "# *_column_tau550 is cumulative downward optical depth at 550 nm"
    write(unitOut,'(A,E17.8E3)') "# output_time_step_s ",dtOut
    write(unitOut,'(A,A)') "# time_unit ",trim(timeUnitName)
    write(unitOut,'(A)') "# time is measured since the first volcanic event"
    write(unitOut,'(A,E17.8E3)') "# tau_floor_for_auto_end ",tauStop
    write(unitOut,'(A,A,A)') "# columns: time layer altitude_km ash_local_tau550 ", &
         "ash_column_tau550 sulfate_local_tau550 sulfate_column_tau550 ", &
         "total_local_tau550 total_column_tau550 sulfate_precursor_potential"

    t = firstTime
    do
       call patmo_volc_setTime(t)
       call patmo_volc_computeLocalState(t,ashLocal,ashColumn, &
            sulfateLocal,sulfateColumn,totalLocal,totalColumn)
       precursorLocal=0d0
       if(allocated(sulfatePrecursor)) precursorLocal=sum(sulfatePrecursor,2)
       do j=1,cellsNumber
          write(unitOut,'(E17.8E3,I6,8E17.8E3)') &
               (t-firstTime)/timeScale,j,height(j)/1d5,ashLocal(j), &
               ashColumn(j),sulfateLocal(j),sulfateColumn(j), &
               totalLocal(j),totalColumn(j),precursorLocal(j)
       end do
       write(unitOut,*)

       if(t>=lastTime) exit
       nextTime = patmo_volc_nextOutput(t,firstTime,dtOut)
       if(nextTime>lastTime) nextTime = lastTime
       if(nextTime<=t) exit
       t = nextTime
    end do
    close(unitOut)

    call patmo_volc_setTime(originalTime)

    write(*,'(A)') "Ash and sulfate plume profile output"
    write(*,'(A,1X,A)') "  file:                ",trim(fname)
    write(*,'(A,F12.4)') "  first_eruption_day:  ", &
         firstTime/secondsPerDay
    write(*,'(A,F12.4)') "  final_volcano_day:   ", &
         lastTime/secondsPerDay
    write(*,'(A,F12.4,1X,A)') "  output_time_step:    ", &
         dtOut/timeScale,trim(timeUnitName)
    write(*,'(A,F12.4)') "  output_time_step_day:", &
         dtOut/secondsPerDay
    write(*,'(A)') ""

  end subroutine patmo_volc_dumpAshProfile

  !***************
  subroutine patmo_volc_dumpState(fname)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    real*8::source(cellsNumber),ashLocal(cellsNumber),ashColumn(cellsNumber)
    real*8::sulfateLocal(cellsNumber),sulfateColumn(cellsNumber)
    real*8::totalLocal(cellsNumber),totalColumn(cellsNumber)
    integer::j,unitOut

    call patmo_volc_computeSO2SourceProfile(volcanoElapsedTime,source)
    call patmo_volc_computeLocalState(volcanoElapsedTime,ashLocal, &
         ashColumn,sulfateLocal,sulfateColumn,totalLocal,totalColumn)

    unitOut = 92
    open(unitOut,file=trim(fname),status="replace")
    write(unitOut,'(A,A,A)') "# altitude_km so2_source_cm-3_s-1 ", &
         "ash_local_tau550 ash_column_tau550 sulfate_local_tau550 ", &
         "sulfate_column_tau550 total_local_tau550 total_column_tau550"
    do j=1,cellsNumber
       write(unitOut,'(8E17.8E3)') height(j)/1d5,source(j),ashLocal(j), &
            ashColumn(j),sulfateLocal(j),sulfateColumn(j),totalLocal(j), &
            totalColumn(j)
    end do
    close(unitOut)

  end subroutine patmo_volc_dumpState

  !***************
  function patmo_volc_getSO2Index(required)
    use patmo_utils
    implicit none
    logical,intent(in)::required
    integer::patmo_volc_getSO2Index

    if(so2IndexCache==-2) so2IndexCache = getSpeciesIndex("SO2",error=.false.)
    if(required.and.so2IndexCache<0) then
       print *,"ERROR: volcano module requires SO2 in the reaction network."
       stop
    end if
    patmo_volc_getSO2Index = so2IndexCache

  end function patmo_volc_getSO2Index

  !***************
  function patmo_volc_so2Active(eventIndex,timeSeconds)
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(in)::timeSeconds
    logical::patmo_volc_so2Active

    patmo_volc_so2Active = .false.
    if(eventIndex<1.or.eventIndex>volcanoEventsNumber) return
    if(events(eventIndex)%duration<=0d0) return
    if(timeSeconds<events(eventIndex)%startTime) return
    if(timeSeconds>=events(eventIndex)%startTime+events(eventIndex)%duration) return
    patmo_volc_so2Active = .true.

  end function patmo_volc_so2Active


  !======================================================================
  ! Ash, precursor and SSA state evolution
  !======================================================================
  subroutine patmo_volc_initAshState()
    use patmo_parameters
    implicit none

    if(volcanoEventsNumber<=0) return
    if(.not.allocated(ashLayerTau550)) then
       allocate(ashLayerTau550(cellsNumber,volcanoEventsNumber))
    end if
    if(.not.allocated(sulfateLayerTau550)) then
       allocate(sulfateLayerTau550(cellsNumber,volcanoEventsNumber))
    end if
    if(.not.allocated(sulfatePrecursor)) then
       allocate(sulfatePrecursor(cellsNumber,volcanoEventsNumber))
    end if
    ashLayerTau550(:,:) = 0d0
    sulfateLayerTau550(:,:) = 0d0
    sulfatePrecursor(:,:) = 0d0
    ashStateTime = 0d0
    ashStateReady = .true.

  end subroutine patmo_volc_initAshState

  !***************
  subroutine patmo_volc_resetAshState()
    implicit none

    if(allocated(ashLayerTau550)) ashLayerTau550(:,:) = 0d0
    if(allocated(sulfateLayerTau550)) sulfateLayerTau550(:,:) = 0d0
    if(allocated(sulfatePrecursor)) sulfatePrecursor(:,:) = 0d0
    ashStateTime = 0d0
    ashStateReady = .false.

  end subroutine patmo_volc_resetAshState

  !***************
  subroutine patmo_volc_updateAshState(targetTime)
    implicit none
    real*8,intent(in)::targetTime
    real*8::timeNow,timeTarget,dtStep

    if(.not.(hasAshOpacity.or.hasSulfateOpacity)) return
    if(volcanoEventsNumber<=0) return
    timeTarget = max(targetTime,0d0)

    if((.not.ashStateReady).or.timeTarget<ashStateTime) then
       call patmo_volc_initAshState()
    end if
    if(.not.ashStateReady) return
    if(timeTarget<=ashStateTime) return

    timeNow = ashStateTime
    do
       if(timeNow>=timeTarget) exit
       dtStep = min(timeTarget-timeNow,patmo_volc_ashInternalStep(timeNow))
       if(dtStep<=0d0) exit
       call patmo_volc_stepAshState(timeNow,dtStep)
       timeNow = timeNow + dtStep
       ashStateTime = timeNow
    end do

  end subroutine patmo_volc_updateAshState

  !***************
  function patmo_volc_ashInternalStep(atTime) result(dtLimit)
    use patmo_constants, only: secondsPerHour
    implicit none
    real*8,intent(in)::atTime
    real*8::dtLimit,starts(2),rises(2),boundary,age
    integer::i,k,j
    dtLimit = 0.25d0*secondsPerHour
    do i=1,volcanoEventsNumber
       starts = [events(i)%startTime,events(i)%startTime+events(i)%sulfateStartDelay]
       rises = [events(i)%ashRiseTime,events(i)%sulfateRiseTime]
       do k=1,2
          age=atTime-starts(k)
          if(rises(k)>0d0.and.age>=0d0.and.age<events(i)%duration+rises(k)) &
               dtLimit=min(dtLimit,rises(k)/60d0)
          if(age>=0d0.and.age<events(i)%duration) dtLimit=min(dtLimit,120d0)
          do j=1,3
             boundary=starts(k)
             if(j==2) boundary=boundary+events(i)%duration
             if(j==3) boundary=boundary+max(events(i)%duration,0d0)+rises(k)
             if(boundary>atTime) dtLimit=min(dtLimit,boundary-atTime)
          end do
       end do
    end do
  end function patmo_volc_ashInternalStep

  function patmo_volc_limitStep(requested) result(dt)
    real*8,intent(in)::requested
    real*8::dt
    dt=requested
    if(volcanoEnabled.and.volcanoClockRunning) &
         dt=min(dt,patmo_volc_ashInternalStep(volcanoElapsedTime))
  end function patmo_volc_limitStep

  !***************
  subroutine patmo_volc_stepAshState(t0,dt)
    implicit none
    real*8,intent(in)::t0,dt
    real*8::mid,start
    integer::i
    mid=t0+0.5d0*dt
    do i=1,volcanoEventsNumber
       if(events(i)%ashTau550>0d0) then
          call patmo_volc_stepParticleState(i,t0,t0+dt, &
               ashLayerTau550(:,i),events(i)%ashTau550, &
               events(i)%startTime,events(i)%duration, &
               events(i)%ashRiseStart,events(i)%ashRiseTime, &
               events(i)%startTime+events(i)%duration, &
               events(i)%ashLifetime,events(i)%ashHorizontalLifetime, &
               events(i)%ashSettling,events(i)%ashVerticalDiffusion)
       end if
       if(events(i)%sulfateTau550<=0d0) cycle
       start=events(i)%startTime+events(i)%sulfateStartDelay
       if(events(i)%duration<=0d0) call patmo_volc_addParticleSourceToState(i,t0,t0+dt, &
            sulfatePrecursor(:,i),events(i)%sulfateTau550,start,events(i)%duration, &
            events(i)%sulfateRiseStart,events(i)%sulfateRiseTime)
       call sulfate_reaction_loss(t0,mid)
       call sulfate_transport(t0,mid)
       if(events(i)%duration>0d0) call patmo_volc_addParticleSourceToState(i,t0,t0+dt, &
            sulfatePrecursor(:,i),events(i)%sulfateTau550,start,events(i)%duration, &
            events(i)%sulfateRiseStart,events(i)%sulfateRiseTime)
       call sulfate_transport(mid,t0+dt)
       call sulfate_reaction_loss(mid,t0+dt)
    end do
  contains
    subroutine sulfate_reaction_loss(a,b)
      real*8,intent(in)::a,b
      !Common horizontal export commutes with formation. Residual particle loss
      !is split around conversion; a gas-like precursor does not sediment.
      call patmo_volc_applyParticleLossToState(a,b,sulfatePrecursor(:,i), &
           start,start,-1d0,events(i)%sulfateHorizontalLifetime)
      call patmo_volc_applyParticleLossToState(a,b,sulfateLayerTau550(:,i), &
           start,start,-1d0,events(i)%sulfateHorizontalLifetime)
      call patmo_volc_applyParticleLossToState(a,0.5d0*(a+b),sulfateLayerTau550(:,i), &
           start,start,events(i)%sulfateLifetime,-1d0)
      call volc_form_sulfate(sulfatePrecursor(:,i),sulfateLayerTau550(:,i), &
           max(b-max(a,start),0d0),events(i)%sulfateFormationTime)
      call patmo_volc_applyParticleLossToState(0.5d0*(a+b),b,sulfateLayerTau550(:,i), &
           start,start,events(i)%sulfateLifetime,-1d0)
    end subroutine sulfate_reaction_loss
    subroutine sulfate_transport(a,b)
      real*8,intent(in)::a,b
      call patmo_volc_moveParticleState(i,a,b,sulfatePrecursor(:,i),start, &
           events(i)%sulfateRiseStart,events(i)%sulfateRiseTime,0d0, &
           events(i)%sulfateVerticalDiffusion)
      call patmo_volc_moveParticleState(i,a,b,sulfateLayerTau550(:,i),start, &
           events(i)%sulfateRiseStart,events(i)%sulfateRiseTime,events(i)%sulfateSettling, &
           events(i)%sulfateVerticalDiffusion)
    end subroutine sulfate_transport
  end subroutine patmo_volc_stepAshState

  !***************
  subroutine patmo_volc_stepParticleState(eventIndex,t0,t1,state,tau550, &
       sourceStart,sourceDuration,riseStart,riseTime,decayStart,lifetime, &
       horizontalLifetime,settlingVelocity,verticalDiffusion)
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(in)::t0,t1,tau550,sourceStart,sourceDuration
    real*8,intent(in)::riseStart,riseTime,decayStart,lifetime,horizontalLifetime
    real*8,intent(in)::settlingVelocity,verticalDiffusion
    real*8,intent(inout)::state(:)
    real*8::mid
    mid=0.5d0*(t0+t1)
    if(sourceDuration<=0d0) call patmo_volc_addParticleSourceToState(eventIndex,t0,t1,state,tau550, &
         sourceStart,sourceDuration,riseStart,riseTime)
    call patmo_volc_applyParticleLossToState(t0,mid,state,sourceStart, &
         decayStart,lifetime,horizontalLifetime)
    call patmo_volc_moveParticleState(eventIndex,t0,mid,state,sourceStart, &
         riseStart,riseTime,settlingVelocity,verticalDiffusion)
    if(sourceDuration>0d0) call patmo_volc_addParticleSourceToState(eventIndex,t0,t1,state,tau550, &
         sourceStart,sourceDuration,riseStart,riseTime)
    call patmo_volc_moveParticleState(eventIndex,mid,t1,state,sourceStart, &
         riseStart,riseTime,settlingVelocity,verticalDiffusion)
    call patmo_volc_applyParticleLossToState(mid,t1,state,sourceStart, &
         decayStart,lifetime,horizontalLifetime)
  end subroutine patmo_volc_stepParticleState

  subroutine patmo_volc_moveParticleState(eventIndex,t0,t1,state,sourceStart, &
       riseStart,riseTime,settlingVelocity,verticalDiffusion)
    use patmo_parameters
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(inout)::state(:)
    real*8,intent(in)::t0,t1,sourceStart,riseStart,riseTime,settlingVelocity,verticalDiffusion
    real*8::air(cellsNumber),kzz(cellsNumber),v(0:cellsNumber),riseDt
    real*8::riseEnd,speed,lower,upper,face,taper
    integer::j
    if(t1<=t0.or.sum(state)<=0d0) return
    !PATMO's M pseudo-species duplicates the physical total number density.
    air=max(0.5d0*sum(nall(:,1:chemSpeciesNumber),2),1d-99)
    kzz=max(verticalDiffusion,0d0)
    if(verticalDiffusion<0d0) kzz=max(eddyKzz,0d0)
    !Sustained, height-limited plume flow, plus one transit time to clear the
    !last emissions. This is prescribed transport, not a buoyant plume solver.
    riseEnd=sourceStart+max(events(eventIndex)%duration,0d0)+riseTime
    riseDt=max(min(t1,riseEnd)-max(t0,sourceStart),0d0)
    speed=patmo_volc_particleRiseVelocity(events(eventIndex)%center, &
         riseStart,riseTime)*riseDt/(t1-t0)
    lower=max(riseStart,events(eventIndex)%center-2d0*events(eventIndex)%sigma)
    upper=max(lower+1d0,events(eventIndex)%center+2d0*events(eventIndex)%sigma)
    v=-settlingVelocity
    do j=1,cellsNumber-1
       face=0.5d0*(height(j)+height(j+1))
       taper=min(1d0,max((upper-face)/(upper-lower),0d0))
       if(face<riseStart) taper=0d0
       v(j)=v(j)+speed*taper
    end do
    v(0)=min(v(0),0d0)
    v(cellsNumber)=0d0
    call volc_transport(state,gridSpace,air,kzz,v,t1-t0)
  end subroutine patmo_volc_moveParticleState

  !***************
  subroutine patmo_volc_addParticleSourceToState(eventIndex,t0,t1,state, &
       tau550,sourceTime,sourceDuration,riseStart,riseTime)
    use patmo_parameters
    implicit none
    integer,intent(in)::eventIndex
    real*8,intent(in)::t0,t1,tau550,sourceTime,sourceDuration
    real*8,intent(in)::riseStart,riseTime
    real*8,intent(inout)::state(:)
    real*8::sourceStart,sourceEnd
    real*8::dTau,w(cellsNumber),norm,layerFraction
    integer::j

    if(eventIndex<1.or.eventIndex>volcanoEventsNumber) return
    if(tau550<=0d0) return

    if(sourceDuration>0d0) then
       sourceStart = max(t0,sourceTime)
       sourceEnd = min(t1,sourceTime+sourceDuration)
       if(sourceEnd<=sourceStart) return
       dTau = tau550 * (sourceEnd-sourceStart) / sourceDuration
    else
       if(.not.(t0<=sourceTime.and.t1>sourceTime)) return
       dTau = tau550
    end if

    if(riseTime>0d0.and.events(eventIndex)%center>riseStart) then
       !All new material enters the same vent layer, including late emissions.
       !A vent below the grid is represented by the lowest model layer.
       j=minloc(abs(height-riseStart),dim=1)
       state(j)=state(j)+dTau
       return
    end if
    !Zero rise time retains the explicit high-altitude injection option.
    call patmo_volc_layerWeights(events(eventIndex)%center,events(eventIndex)%sigma,w,norm)
    if(norm<=0d0) return

    do j=1,cellsNumber
       layerFraction = w(j) * max(gridSpace(j),1d0) / norm
       state(j) = state(j) + dTau * layerFraction
    end do

  end subroutine patmo_volc_addParticleSourceToState

  !***************


  !***************
  subroutine patmo_volc_applyParticleLossToState(t0,t1,state,sourceStart, &
       decayStart,lifetime,horizontalLifetime)
    implicit none
    real*8,intent(in)::t0,t1,sourceStart,decayStart
    real*8,intent(in)::lifetime,horizontalLifetime
    real*8,intent(inout)::state(:)
    real*8::lossFactor,lossTime

    lossFactor = 1d0

    if(horizontalLifetime>0d0) then
       lossTime = max(t1 - max(t0,sourceStart),0d0)
       lossFactor = lossFactor * exp(-lossTime/horizontalLifetime)
    end if

    if(lifetime>0d0) then
       lossTime = max(t1 - max(t0,decayStart),0d0)
       lossFactor = lossFactor * exp(-lossTime/lifetime)
    elseif(lifetime==0d0.and.t1>decayStart) then
       state(:) = 0d0
       return
    end if

    state(:) = state(:) * lossFactor

  end subroutine patmo_volc_applyParticleLossToState

  !***************


  !***************


  !***************
  function patmo_volc_particleRiseVelocity(center,riseStart,riseTime)
    implicit none
    real*8,intent(in)::center,riseStart,riseTime
    real*8::patmo_volc_particleRiseVelocity

    patmo_volc_particleRiseVelocity = 0d0
    if(riseTime<=0d0) return
    patmo_volc_particleRiseVelocity = max(center-riseStart,0d0) / riseTime

  end function patmo_volc_particleRiseVelocity

  !***************
  !***************


  !***************
  function patmo_volc_stokesSettling(radiusCm,particleDensity,airViscosity, &
       gravityLocal,cunninghamFactor)
    implicit none
    real*8,intent(in)::radiusCm,particleDensity,airViscosity
    real*8,intent(in)::gravityLocal,cunninghamFactor
    real*8::patmo_volc_stokesSettling

    patmo_volc_stokesSettling = 0d0
    if(radiusCm<=0d0) return
    if(particleDensity<=0d0) return
    if(airViscosity<=0d0) return
    if(gravityLocal<=0d0) return

    patmo_volc_stokesSettling = 2d0/9d0 * radiusCm**2 &
         * particleDensity * gravityLocal / airViscosity &
         * max(cunninghamFactor,0d0)

  end function patmo_volc_stokesSettling

  !***************
  subroutine patmo_volc_layerWeights(center,sigma,w,norm)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::center,sigma
    real*8,intent(out)::w(cellsNumber),norm
    real*8::arg,dz,bestDistance,distance,sig
    integer::j,bestCell

    sig = max(abs(sigma),1d0)
    norm = 0d0
    do j=1,cellsNumber
       arg = (height(j)-center)/sig
       if(abs(arg)>40d0) then
          w(j) = 0d0
       else
          w(j) = exp(-0.5d0*arg*arg)
       end if
       dz = max(gridSpace(j),1d0)
       norm = norm + w(j) * dz
    end do

    if(norm>0d0) return

    bestCell = 1
    bestDistance = abs(height(1)-center)
    do j=2,cellsNumber
       distance = abs(height(j)-center)
       if(distance<bestDistance) then
          bestDistance = distance
          bestCell = j
       end if
    end do
    w(:) = 0d0
    w(bestCell) = 1d0
    norm = max(gridSpace(bestCell),1d0)

  end subroutine patmo_volc_layerWeights

  !***************
  function patmo_volc_particleSpectralScale(ibin,wavelengthExp)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer,intent(in)::ibin
    real*8,intent(in)::wavelengthExp
    real*8::patmo_volc_particleSpectralScale,lambdaNm

    patmo_volc_particleSpectralScale = 1d0
    if(abs(wavelengthExp)<=1d-99) return
    if(energyMid(ibin)<=0d0) return

    lambdaNm = 1d7 * planck_eV * clight / energyMid(ibin)
    if(lambdaNm<=0d0) return

    patmo_volc_particleSpectralScale = (lambdaNm/550d0)**(-wavelengthExp)
    patmo_volc_particleSpectralScale = &
         min(max(patmo_volc_particleSpectralScale,1d-6),1d6)

  end function patmo_volc_particleSpectralScale

  !***************
  function patmo_volc_binWavelengthNm(ibin)
    use patmo_commons
    use patmo_parameters
    use patmo_constants
    implicit none
    integer,intent(in)::ibin
    real*8::patmo_volc_binWavelengthNm

    patmo_volc_binWavelengthNm = 0d0
    if(ibin<1.or.ibin>photoBinsNumber) return
    if(energyMid(ibin)<=0d0) return
    patmo_volc_binWavelengthNm = 1d7 * planck_eV * clight / energyMid(ibin)

  end function patmo_volc_binWavelengthNm

  !***************
  function patmo_volc_firstAshStart()
    use patmo_constants
    implicit none
    real*8::patmo_volc_firstAshStart
    integer::i
    logical::found

    patmo_volc_firstAshStart = 0d0
    found = .false.
    do i=1,volcanoEventsNumber
       if(events(i)%so2ColumnFlux<=0d0.and.events(i)%ashTau550<=0d0 &
            .and.events(i)%sulfateTau550<=0d0) cycle
       if((.not.found).or.events(i)%startTime<patmo_volc_firstAshStart) then
          patmo_volc_firstAshStart = events(i)%startTime
          found = .true.
       end if
    end do

  end function patmo_volc_firstAshStart

  !***************
  function patmo_volc_estimateOpticalDepthEnd(tauFloor) result(endTime)
    use patmo_commons
    use patmo_constants, only: secondsPerDay
    implicit none
    real*8,intent(in)::tauFloor
    real*8::endTime,lastSource,bound,ashScale(volcanoEventsNumber)
    real*8::sulfateScale(volcanoEventsNumber),originalTime
    integer::i,j,step
    lastSource=0d0
    do i=1,volcanoEventsNumber
       lastSource=max(lastSource,events(i)%startTime+max(events(i)%duration,0d0) &
            +events(i)%ashRiseTime,events(i)%startTime+events(i)%sulfateStartDelay &
            +max(events(i)%duration,0d0)+events(i)%sulfateRiseTime)
       if(events(i)%duration<=0d0) lastSource=max(lastSource,events(i)%startTime &
            +events(i)%sulfateStartDelay+1d0)
       ashScale(i)=1d0
       sulfateScale(i)=1d0
       do j=1,photoBinsNumber
          ashScale(i)=max(ashScale(i),patmo_volc_particleSpectralScale(j,events(i)%ashWavelengthExp))
          sulfateScale(i)=max(sulfateScale(i),patmo_volc_particleSpectralScale(j,events(i)%sulfateWavelengthExp))
       end do
    end do
    originalTime=volcanoElapsedTime
    call patmo_volc_resetAshState()
    endTime=lastSource
    !The unconverted budget bounds all future SSA opacity, including UV.
    !Never use a ballistic fall time as proof that a diffusing plume is gone.
    do step=0,36500
       call patmo_volc_updateAshState(endTime)
       bound=0d0
       if(allocated(ashLayerTau550)) then
          do i=1,volcanoEventsNumber
             bound=bound+sum(ashLayerTau550(:,i))*ashScale(i) &
                  +sum(sulfateLayerTau550(:,i)+sulfatePrecursor(:,i))*sulfateScale(i)
          end do
       end if
       if(bound<=max(tauFloor,1d-12)) exit
       if(step==36500) error stop 'No optical auto-end within 100 years; set output_end_hour explicitly'
       endTime=endTime+secondsPerDay
    end do
    call patmo_volc_setTime(originalTime)
  end function patmo_volc_estimateOpticalDepthEnd

  !***************


  !***************
  subroutine patmo_volc_selectWavelengths(wavelengthStepNm,selectedBins, &
       targetWavelength,modelWavelength,wavelengthsNumber)
    use patmo_commons
    implicit none
    real*8,intent(in)::wavelengthStepNm
    integer,allocatable,intent(out)::selectedBins(:)
    real*8,allocatable,intent(out)::targetWavelength(:),modelWavelength(:)
    integer,intent(out)::wavelengthsNumber
    real*8::lambdaMin,lambdaMax,lambdaNm,target,step
    real*8::bestDistance,distance
    integer::i,iw,estimateNumber,bestBin,lastBin

    lambdaMin = huge(1d0)
    lambdaMax = 0d0
    do i=1,photoBinsNumber
       lambdaNm = patmo_volc_binWavelengthNm(i)
       if(lambdaNm<=0d0) cycle
       lambdaMin = min(lambdaMin,lambdaNm)
       lambdaMax = max(lambdaMax,lambdaNm)
    end do
    if(lambdaMax<=0d0.or.lambdaMin>=huge(1d0)) then
       print *,"ERROR: wavelength grid is not initialized for volcano pre-run."
       stop
    end if

    if(wavelengthStepNm<=0d0) then
       allocate(selectedBins(photoBinsNumber))
       allocate(targetWavelength(photoBinsNumber))
       allocate(modelWavelength(photoBinsNumber))
       wavelengthsNumber = 0
       do i=photoBinsNumber,1,-1
          lambdaNm = patmo_volc_binWavelengthNm(i)
          if(lambdaNm<=0d0) cycle
          wavelengthsNumber = wavelengthsNumber + 1
          selectedBins(wavelengthsNumber) = i
          targetWavelength(wavelengthsNumber) = lambdaNm
          modelWavelength(wavelengthsNumber) = lambdaNm
       end do
       return
    end if

    step = max(wavelengthStepNm,1d-12)
    estimateNumber = int((lambdaMax-lambdaMin)/step) + 1
    if(lambdaMin + dble(estimateNumber-1)*step < lambdaMax-1d-8) then
       estimateNumber = estimateNumber + 1
    end if
    estimateNumber = max(estimateNumber,1)

    allocate(selectedBins(estimateNumber))
    allocate(targetWavelength(estimateNumber))
    allocate(modelWavelength(estimateNumber))
    wavelengthsNumber = 0
    lastBin = -1
    do iw=1,estimateNumber
       target = lambdaMin + dble(iw-1)*step
       if(iw==estimateNumber.or.target>lambdaMax) target = lambdaMax

       bestBin = 1
       bestDistance = huge(1d0)
       do i=1,photoBinsNumber
          lambdaNm = patmo_volc_binWavelengthNm(i)
          if(lambdaNm<=0d0) cycle
          distance = abs(lambdaNm-target)
          if(distance<bestDistance) then
             bestDistance = distance
             bestBin = i
          end if
       end do

       if(bestBin==lastBin) cycle
       wavelengthsNumber = wavelengthsNumber + 1
       selectedBins(wavelengthsNumber) = bestBin
       targetWavelength(wavelengthsNumber) = target
       modelWavelength(wavelengthsNumber) = patmo_volc_binWavelengthNm(bestBin)
       lastBin = bestBin
    end do

  end subroutine patmo_volc_selectWavelengths

  !***************
  function patmo_volc_commentStart(line)
    implicit none
    character(len=*),intent(in)::line
    integer::patmo_volc_commentStart,hashPos,bangPos

    hashPos = index(line,"#")
    bangPos = index(line,"!")
    if(hashPos>0.and.bangPos>0) then
       patmo_volc_commentStart = min(hashPos,bangPos)
    elseif(hashPos>0) then
       patmo_volc_commentStart = hashPos
    elseif(bangPos>0) then
       patmo_volc_commentStart = bangPos
    else
       patmo_volc_commentStart = 0
    end if

  end function patmo_volc_commentStart

  !***************

  !======================================================================
  ! Event parsing and string helpers
  !======================================================================
  subroutine patmo_volc_parseKeywordEvent(line,startDay,durationDay, &
       centerKm,sigmaKm,so2ColumnFlux,ashTau550,ashLifetimeDay, &
       ashHorizontalLifetimeDay,ashSettling,ashWavelengthExp, &
       ashRiseStartKm,ashRiseTimeDay,ashVerticalDiffusion, &
       ashParticleRadiusUm,ashParticleDensity,ashAirViscosity, &
       ashGravity,ashCunninghamFactor,sulfateTau550,sulfateLifetimeDay, &
       sulfateHorizontalLifetimeDay,sulfateSettling,sulfateWavelengthExp, &
       sulfateRiseStartKm,sulfateRiseTimeDay,sulfateVerticalDiffusion, &
       sulfateParticleRadiusUm,sulfateParticleDensity, &
       sulfateAirViscosity,sulfateGravity,sulfateCunninghamFactor, &
       sulfateStartDelayDay,sulfateFormationTimeDay,ios)
    use patmo_constants, only: secondsPerDay, pi, av
    implicit none
    character(len=*),intent(in)::line
    real*8,intent(inout)::startDay,durationDay,centerKm,sigmaKm
    real*8,intent(inout)::so2ColumnFlux,ashTau550,ashLifetimeDay
    real*8,intent(inout)::ashHorizontalLifetimeDay,ashSettling
    real*8,intent(inout)::ashWavelengthExp
    real*8,intent(inout)::ashRiseStartKm,ashRiseTimeDay
    real*8,intent(inout)::ashVerticalDiffusion,ashParticleRadiusUm
    real*8,intent(inout)::ashParticleDensity,ashAirViscosity,ashGravity
    real*8,intent(inout)::ashCunninghamFactor
    real*8,intent(inout)::sulfateTau550,sulfateLifetimeDay
    real*8,intent(inout)::sulfateHorizontalLifetimeDay,sulfateSettling
    real*8,intent(inout)::sulfateWavelengthExp
    real*8,intent(inout)::sulfateRiseStartKm,sulfateRiseTimeDay
    real*8,intent(inout)::sulfateVerticalDiffusion,sulfateParticleRadiusUm
    real*8,intent(inout)::sulfateParticleDensity,sulfateAirViscosity
    real*8,intent(inout)::sulfateGravity,sulfateCunninghamFactor
    real*8,intent(inout)::sulfateStartDelayDay,sulfateFormationTimeDay
    integer,intent(out)::ios
    real*8,parameter::fwhmToSigma = 1d0/2.3548200450309493d0
    real*8,parameter::so2MolarMass = 64.066d0
    real*8,parameter::tgToG = 1d12
    real*8,parameter::km2ToCm2 = 1d10
    character(len=4096)::work
    character(len=80)::key,value
    integer::eqPos,keyStart,keyEnd,valueStart,valueEnd,nline,readIos
    logical::hasStart,hasDuration,hasCenter,hasSigma,hasSource
    logical::hasSO2Flux,hasSO2Column,hasSO2Mass,hasAsh,hasSulfate
    logical::hasTop,hasBottom,hasFwhm,hasArea,hasRadius
    real*8::tmp,plumeTopKm,plumeBottomKm,plumeFwhmKm
    real*8::so2ColumnAmount,so2MassTg,injectionAreaKm2,plumeRadiusKm
    real*8::durationSeconds,areaKm2

    ios = 0
    hasStart = .false.
    hasDuration = .false.
    hasCenter = .false.
    hasSigma = .false.
    hasSource = .false.
    hasSO2Flux = .false.
    hasSO2Column = .false.
    hasSO2Mass = .false.
    hasAsh = .false.
    hasSulfate = .false.
    hasTop = .false.
    hasBottom = .false.
    hasFwhm = .false.
    hasArea = .false.
    hasRadius = .false.
    plumeTopKm = 0d0
    plumeBottomKm = 0d0
    plumeFwhmKm = 0d0
    so2ColumnAmount = 0d0
    so2MassTg = 0d0
    injectionAreaKm2 = 0d0
    plumeRadiusKm = 0d0

    work = adjustl(line)
    do eqPos=1,len(work)
       if(work(eqPos:eqPos)==",".or.work(eqPos:eqPos)==char(9)) work(eqPos:eqPos) = " "
    end do
    nline = len_trim(work)

    do
       eqPos = index(work,"=")
       if(eqPos<=0) exit

       keyEnd = eqPos - 1
       do while(keyEnd>=1.and.work(keyEnd:keyEnd)==" ")
          keyEnd = keyEnd - 1
       end do

       keyStart = keyEnd
       do while(keyStart>1.and.work(keyStart-1:keyStart-1)/=" ")
          keyStart = keyStart - 1
       end do

       valueStart = eqPos + 1
       do while(valueStart<=nline.and.work(valueStart:valueStart)==" ")
          valueStart = valueStart + 1
       end do

       valueEnd = valueStart
       do while(valueEnd<=nline)
          if(work(valueEnd:valueEnd)==" ".or.work(valueEnd:valueEnd)==",") exit
          valueEnd = valueEnd + 1
       end do
       valueEnd = valueEnd - 1

       key = " "
       value = " "
       if(keyStart<=keyEnd) key = patmo_volc_lower(work(keyStart:keyEnd))
       if(valueStart<=valueEnd) value = adjustl(work(valueStart:valueEnd))

       select case(trim(key))
       case("event_id","event","name","scenario")
          continue
       case("start_day","start_days","start")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          startDay = tmp
          hasStart = .true.
       case("start_hour","start_hours","start_hr")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          startDay = tmp / 24d0
          hasStart = .true.
       case("start_s","start_sec","start_second","start_seconds")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          startDay = tmp / secondsPerDay
          hasStart = .true.
       case("duration_day","duration_days","duration")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          durationDay = tmp
          hasDuration = .true.
       case("duration_hour","duration_hours","duration_hr")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          durationDay = tmp / 24d0
          hasDuration = .true.
       case("duration_s","duration_sec","duration_second","duration_seconds")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          durationDay = tmp / secondsPerDay
          hasDuration = .true.
       case("plume_center_km","center_km","injection_center_km", &
            "injection_altitude_km","altitude_km","z_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          centerKm = tmp
          hasCenter = .true.
       case("plume_sigma_km","sigma_km","injection_sigma_km", &
            "vertical_sigma_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sigmaKm = tmp
          hasSigma = .true.
       case("plume_fwhm_km","fwhm_km","vertical_fwhm_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          plumeFwhmKm = tmp
          hasFwhm = .true.
       case("plume_top_km","top_km","injection_top_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          plumeTopKm = tmp
          hasTop = .true.
       case("plume_bottom_km","bottom_km","injection_bottom_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          plumeBottomKm = tmp
          hasBottom = .true.
       case("ash_rise_start_km","ash_start_altitude_km", &
            "plume_rise_start_km","vent_altitude_km", &
            "eruption_start_altitude_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashRiseStartKm = tmp
       case("ash_rise_time_day","ash_rise_time_days", &
            "plume_rise_time_day","plume_rise_time_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashRiseTimeDay = tmp
       case("ash_rise_time_hour","ash_rise_time_hours","ash_rise_time_hr", &
            "plume_rise_time_hour","plume_rise_time_hours", &
            "plume_rise_time_hr")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashRiseTimeDay = tmp / 24d0
       case("ash_rise_time_s","ash_rise_time_sec", &
            "ash_rise_time_second","ash_rise_time_seconds", &
            "plume_rise_time_s","plume_rise_time_sec", &
            "plume_rise_time_second","plume_rise_time_seconds")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashRiseTimeDay = tmp / secondsPerDay
       case("sulfate_rise_start_km","sulfate_start_altitude_km", &
            "ssa_rise_start_km","ssa_start_altitude_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateRiseStartKm = tmp
       case("sulfate_rise_time_day","sulfate_rise_time_days", &
            "ssa_rise_time_day","ssa_rise_time_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateRiseTimeDay = tmp
       case("sulfate_rise_time_hour","sulfate_rise_time_hours", &
            "sulfate_rise_time_hr","ssa_rise_time_hour", &
            "ssa_rise_time_hours","ssa_rise_time_hr")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateRiseTimeDay = tmp / 24d0
       case("sulfate_rise_time_s","sulfate_rise_time_sec", &
            "sulfate_rise_time_second","sulfate_rise_time_seconds", &
            "ssa_rise_time_s","ssa_rise_time_sec", &
            "ssa_rise_time_second","ssa_rise_time_seconds")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateRiseTimeDay = tmp / secondsPerDay
       case("sulfate_start_delay_day","sulfate_start_delay_days", &
            "sulfate_delay_day","sulfate_delay_days", &
            "ssa_start_delay_day","ssa_delay_day")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateStartDelayDay = tmp
       case("sulfate_start_delay_hour","sulfate_start_delay_hours", &
            "sulfate_delay_hour","sulfate_delay_hours", &
            "ssa_start_delay_hour","ssa_delay_hour")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateStartDelayDay = tmp / 24d0
       case("sulfate_start_delay_s","sulfate_delay_s", &
            "sulfate_delay_sec","ssa_start_delay_s","ssa_delay_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateStartDelayDay = tmp / secondsPerDay
       case("sulfate_formation_day","sulfate_formation_days", &
            "sulfate_growth_day","sulfate_growth_days", &
            "ssa_formation_day","ssa_growth_day")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateFormationTimeDay = tmp
       case("sulfate_formation_hour","sulfate_formation_hours", &
            "sulfate_growth_hour","sulfate_growth_hours", &
            "ssa_formation_hour","ssa_growth_hour")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateFormationTimeDay = tmp / 24d0
       case("sulfate_formation_s","sulfate_formation_sec", &
            "sulfate_growth_s","sulfate_growth_sec", &
            "ssa_formation_s","ssa_growth_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateFormationTimeDay = tmp / secondsPerDay
       case("so2_flux_cm2_s","so2_column_flux_cm2_s","so2_column_flux")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2ColumnFlux = tmp
          hasSO2Flux = .true.
          hasSource = .true.
       case("so2_column_cm2","so2_column_molecules_cm2", &
            "so2_total_column_cm2","so2_burden_cm2")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2ColumnAmount = tmp
          hasSO2Column = .true.
          hasSource = .true.
       case("so2_mass_tg","so2_tg","so2_mass_mt","so2_mt")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2MassTg = tmp
          hasSO2Mass = .true.
          hasSource = .true.
       case("so2_mass_kg","so2_kg")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2MassTg = tmp * 1d-9
          hasSO2Mass = .true.
          hasSource = .true.
       case("so2_mass_g","so2_g")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          so2MassTg = tmp * 1d-12
          hasSO2Mass = .true.
          hasSource = .true.
       case("injection_area_km2","plume_area_km2","area_km2")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          injectionAreaKm2 = tmp
          hasArea = .true.
       case("plume_radius_km","injection_radius_km","radius_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          plumeRadiusKm = tmp
          hasRadius = .true.
       case("plume_diameter_km","injection_diameter_km","diameter_km")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          plumeRadiusKm = 0.5d0 * tmp
          hasRadius = .true.
       case("ash_tau_550","ash_optical_depth_550","ash_aod_550", &
            "aod_550","tau_550","optical_depth_550")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashTau550 = tmp
          hasAsh = .true.
          hasSource = .true.
       case("sulfate_tau_550","sulfate_optical_depth_550", &
            "sulfate_aod_550","sulfate_aerosol_tau_550", &
            "sulfate_aerosol_aod_550","ssa_tau_550","ssa_aod_550")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateTau550 = tmp
          hasSulfate = .true.
          hasSource = .true.
       case("ash_lifetime_day","ash_lifetime_days","ash_decay_day", &
            "ash_decay_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashLifetimeDay = tmp
       case("ash_lifetime_year","ash_lifetime_years","ash_decay_year", &
            "ash_decay_years")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashLifetimeDay = tmp * 365d0
       case("ash_lifetime_s","ash_lifetime_sec","ash_decay_s", &
            "ash_decay_sec")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashLifetimeDay = tmp / secondsPerDay
       case("sulfate_lifetime_day","sulfate_lifetime_days", &
            "sulfate_decay_day","sulfate_decay_days", &
            "sulfate_aerosol_lifetime_day","sulfate_aerosol_lifetime_days", &
            "ssa_lifetime_day","ssa_lifetime_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateLifetimeDay = tmp
       case("sulfate_lifetime_year","sulfate_lifetime_years", &
            "sulfate_decay_year","sulfate_decay_years", &
            "ssa_lifetime_year","ssa_lifetime_years")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateLifetimeDay = tmp * 365d0
       case("sulfate_lifetime_s","sulfate_lifetime_sec", &
            "sulfate_decay_s","sulfate_decay_sec", &
            "ssa_lifetime_s","ssa_lifetime_sec")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateLifetimeDay = tmp / secondsPerDay
       case("ash_horizontal_lifetime_day","ash_horizontal_lifetime_days", &
            "ash_export_lifetime_day","ash_export_lifetime_days", &
            "ash_local_lifetime_day","ash_local_lifetime_days", &
            "aerosol_horizontal_lifetime_day", &
            "aerosol_horizontal_lifetime_days", &
            "aerosol_export_lifetime_day","aerosol_export_lifetime_days", &
            "local_horizontal_lifetime_day","local_horizontal_lifetime_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashHorizontalLifetimeDay = tmp
       case("ash_horizontal_lifetime_year","ash_horizontal_lifetime_years", &
            "ash_export_lifetime_year","ash_export_lifetime_years", &
            "ash_local_lifetime_year","ash_local_lifetime_years", &
            "aerosol_horizontal_lifetime_year", &
            "aerosol_horizontal_lifetime_years", &
            "aerosol_export_lifetime_year", &
            "aerosol_export_lifetime_years", &
            "local_horizontal_lifetime_year", &
            "local_horizontal_lifetime_years")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashHorizontalLifetimeDay = tmp * 365d0
       case("ash_horizontal_lifetime_s","ash_horizontal_lifetime_sec", &
            "ash_export_lifetime_s","ash_export_lifetime_sec", &
            "ash_local_lifetime_s","ash_local_lifetime_sec", &
            "aerosol_horizontal_lifetime_s", &
            "aerosol_horizontal_lifetime_sec", &
            "aerosol_export_lifetime_s","aerosol_export_lifetime_sec", &
            "local_horizontal_lifetime_s","local_horizontal_lifetime_sec")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashHorizontalLifetimeDay = tmp / secondsPerDay
       case("sulfate_horizontal_lifetime_day", &
            "sulfate_horizontal_lifetime_days", &
            "sulfate_export_lifetime_day","sulfate_export_lifetime_days", &
            "sulfate_local_lifetime_day","sulfate_local_lifetime_days", &
            "sulfate_aerosol_horizontal_lifetime_day", &
            "sulfate_aerosol_horizontal_lifetime_days", &
            "ssa_horizontal_lifetime_day","ssa_horizontal_lifetime_days")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateHorizontalLifetimeDay = tmp
       case("sulfate_horizontal_lifetime_year", &
            "sulfate_horizontal_lifetime_years", &
            "sulfate_export_lifetime_year","sulfate_export_lifetime_years", &
            "sulfate_local_lifetime_year","sulfate_local_lifetime_years", &
            "ssa_horizontal_lifetime_year","ssa_horizontal_lifetime_years")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateHorizontalLifetimeDay = tmp * 365d0
       case("sulfate_horizontal_lifetime_s", &
            "sulfate_horizontal_lifetime_sec", &
            "sulfate_export_lifetime_s","sulfate_export_lifetime_sec", &
            "sulfate_local_lifetime_s","sulfate_local_lifetime_sec", &
            "ssa_horizontal_lifetime_s","ssa_horizontal_lifetime_sec")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateHorizontalLifetimeDay = tmp / secondsPerDay
       case("ash_settling_cm_s","ash_fall_speed_cm_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashSettling = tmp
       case("ash_settling_m_s","ash_fall_speed_m_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashSettling = tmp * 1d2
       case("ash_settling_km_day","ash_fall_speed_km_day")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashSettling = tmp * 1d5 / secondsPerDay
       case("sulfate_settling_cm_s","sulfate_fall_speed_cm_s", &
            "sulfate_aerosol_settling_cm_s","ssa_settling_cm_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateSettling = tmp
       case("sulfate_settling_m_s","sulfate_fall_speed_m_s", &
            "sulfate_aerosol_settling_m_s","ssa_settling_m_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateSettling = tmp * 1d2
       case("sulfate_settling_km_day","sulfate_fall_speed_km_day", &
            "sulfate_aerosol_settling_km_day","ssa_settling_km_day")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateSettling = tmp * 1d5 / secondsPerDay
       case("ash_particle_radius_um","ash_radius_um", &
            "ash_effective_radius_um","aerosol_radius_um")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashParticleRadiusUm = tmp
       case("sulfate_particle_radius_um","sulfate_radius_um", &
            "sulfate_effective_radius_um","sulfate_aerosol_radius_um", &
            "ssa_particle_radius_um","ssa_radius_um")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateParticleRadiusUm = tmp
       case("ash_particle_density_g_cm3","ash_density_g_cm3", &
            "particle_density_g_cm3","aerosol_density_g_cm3")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashParticleDensity = tmp
       case("sulfate_particle_density_g_cm3","sulfate_density_g_cm3", &
            "sulfate_aerosol_density_g_cm3","ssa_particle_density_g_cm3", &
            "ssa_density_g_cm3")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateParticleDensity = tmp
       case("ash_air_viscosity_g_cm_s","air_viscosity_g_cm_s", &
            "dynamic_viscosity_g_cm_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashAirViscosity = tmp
       case("sulfate_air_viscosity_g_cm_s", &
            "sulfate_aerosol_air_viscosity_g_cm_s", &
            "ssa_air_viscosity_g_cm_s")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateAirViscosity = tmp
       case("ash_gravity_cm_s2","gravity_cm_s2","g_cm_s2")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashGravity = tmp
       case("sulfate_gravity_cm_s2","sulfate_aerosol_gravity_cm_s2", &
            "ssa_gravity_cm_s2")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateGravity = tmp
       case("ash_cunningham_factor","cunningham_factor", &
            "slip_correction_factor")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashCunninghamFactor = tmp
       case("sulfate_cunningham_factor", &
            "sulfate_aerosol_cunningham_factor", &
            "ssa_cunningham_factor")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateCunninghamFactor = tmp
       case("ash_vertical_diffusion_cm2_s","ash_kzz_cm2_s", &
            "aerosol_vertical_diffusion_cm2_s","aerosol_kzz_cm2_s")
          read(value,*,iostat=readIos) tmp
          if(trim(value)=="background") then
             tmp=-1d0
             readIos=0
          end if
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashVerticalDiffusion = tmp
       case("sulfate_vertical_diffusion_cm2_s","sulfate_kzz_cm2_s", &
            "sulfate_aerosol_vertical_diffusion_cm2_s", &
            "sulfate_aerosol_kzz_cm2_s","ssa_vertical_diffusion_cm2_s", &
            "ssa_kzz_cm2_s")
          read(value,*,iostat=readIos) tmp
          if(trim(value)=="background") then
             tmp=-1d0
             readIos=0
          end if
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateVerticalDiffusion = tmp
       case("ash_lambda_exponent","ash_wavelength_exp","ash_angstrom_exp", &
            "ash_alpha","angstrom_exponent","alpha")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          ashWavelengthExp = tmp
       case("sulfate_lambda_exponent","sulfate_wavelength_exp", &
            "sulfate_angstrom_exp","sulfate_alpha", &
            "sulfate_aerosol_lambda_exponent", &
            "ssa_lambda_exponent","ssa_wavelength_exp","ssa_alpha")
          read(value,*,iostat=readIos) tmp
          if(readIos/=0) then
             ios = 1
             return
          end if
          sulfateWavelengthExp = tmp
       case default
          print *,"WARNING: unknown volcano event key ignored: ",trim(key)
       end select

       if(valueEnd>=nline) exit
       work = adjustl(work(valueEnd+1:))
       nline = len_trim(work)
    end do

    if(.not.hasCenter.and.hasTop.and.hasBottom) then
       centerKm = 0.5d0 * (plumeTopKm + plumeBottomKm)
       hasCenter = .true.
    end if

    if(.not.hasSigma.and.hasFwhm) then
       sigmaKm = abs(plumeFwhmKm) * fwhmToSigma
       hasSigma = .true.
    end if

    if(.not.hasSigma.and.hasTop.and.hasBottom) then
       sigmaKm = max(abs(plumeTopKm - plumeBottomKm) / 4d0,1d-12)
       hasSigma = .true.
    end if

    if((.not.hasSO2Flux).and.hasSO2Column.and.so2ColumnAmount>0d0) then
       durationSeconds = durationDay * secondsPerDay
       if(durationSeconds<=0d0) then
          ios = 1
          return
       end if
       so2ColumnFlux = so2ColumnAmount / durationSeconds
       hasSO2Flux = .true.
    end if

    if((.not.hasSO2Flux).and.hasSO2Mass.and.so2MassTg>0d0) then
       areaKm2 = injectionAreaKm2
       if((.not.hasArea).and.hasRadius) areaKm2 = pi * plumeRadiusKm**2
       durationSeconds = durationDay * secondsPerDay
       if(durationSeconds<=0d0.or.areaKm2<=0d0) then
          ios = 1
          return
       end if
       so2ColumnFlux = (so2MassTg * tgToG / so2MolarMass * av) &
            / (areaKm2 * km2ToCm2) / durationSeconds
       hasSO2Flux = .true.
    end if

    hasSource = hasSO2Flux .or. hasSO2Column .or. hasSO2Mass .or. hasAsh &
         .or. hasSulfate

    if(.not.(hasStart.and.hasDuration.and.hasCenter.and.hasSigma.and.hasSource)) then
       ios = 1
    end if

  end subroutine patmo_volc_parseKeywordEvent

  !***************
  function patmo_volc_lower(text)
    implicit none
    character(len=*),intent(in)::text
    character(len=len(text))::patmo_volc_lower
    integer::i,ich

    do i=1,len(text)
       ich = iachar(text(i:i))
       if(ich>=iachar("A").and.ich<=iachar("Z")) then
          patmo_volc_lower(i:i) = achar(ich + iachar("a") - iachar("A"))
       else
          patmo_volc_lower(i:i) = text(i:i)
       end if
    end do

  end function patmo_volc_lower


  !======================================================================
  ! Finite-volume transport and sulfate formation
  !======================================================================

  !Conserve optical potential, not sulfur mass. Only the product is opacity.
  subroutine volc_form_sulfate(precursor, aerosol, dt, formationTime)
    real*8,intent(inout)::precursor(:),aerosol(:)
    real*8,intent(in)::dt,formationTime
    real*8::remaining,formed(size(precursor))
    if(dt<=0d0) return
    remaining = 0d0
    if(formationTime>0d0) remaining = exp(-dt/formationTime)
    formed = precursor * (1d0-remaining)
    precursor = precursor * remaining
    aerosol = aerosol + formed
  end subroutine volc_form_sulfate

  !Layer-integrated budget; velocity is positive upward at cell faces.
  !TVD reconstruction and SSP-RK2 reduce upwind spreading without clipping mass.
  !The lower boundary allows settling outflow; all other boundary fluxes vanish.
  subroutine volc_transport(q,dz,air,kzz,velocity,dt,bottomLoss)
    real*8,intent(inout)::q(:)
    real*8,intent(in)::dz(:),air(:),kzz(:),velocity(0:),dt
    real*8,intent(out),optional::bottomLoss
    real*8::q1(size(q)),f0(0:size(q)),f1(0:size(q))
    real*8::conductance(0:size(q)),rate,h,left,time,loss
    integer::j,n
    n = size(q)
    if(present(bottomLoss)) bottomLoss=0d0
    if(dt<=0d0) return
    if(any(dz<=0d0).or.any(air<=0d0).or.any(kzz<0d0)) &
         error stop 'Invalid volcanic transport grid or diffusivity'
    conductance=0d0
    do j=1,n-1
       conductance(j) = 0.5d0*(kzz(j)+kzz(j+1)) &
            * 0.5d0*(air(j)+air(j+1)) / (0.5d0*(dz(j)+dz(j+1)))
    end do
    rate=0d0
    do j=1,n
       left=2d0*(max(velocity(j),0d0)+max(-velocity(j-1),0d0))/dz(j) &
            +(conductance(j)+conductance(j-1))/(air(j)*dz(j))
       rate=max(rate,left)
    end do
    time=0d0
    loss=0d0
    do while(time<dt)
       h=dt-time
       if(rate>0d0) h=min(h,0.45d0/rate)
       call flux(q,f0)
       q1=q+h*(f0(0:n-1)-f0(1:n))
       call flux(q1,f1)
       q=0.5d0*(q+q1+h*(f1(0:n-1)-f1(1:n)))
       loss=loss-0.5d0*h*(f0(0)+f1(0))
       time=time+h
    end do
    if(present(bottomLoss)) bottomLoss=loss
  contains
    subroutine flux(state,f)
      real*8,intent(in)::state(:)
      real*8,intent(out)::f(0:)
      real*8::c(n),slope(n),a,b,central,edge
      integer::i
      c=state/dz
      slope=0d0
      do i=2,n-1
         a=(c(i)-c(i-1))/(0.5d0*(dz(i)+dz(i-1)))
         b=(c(i+1)-c(i))/(0.5d0*(dz(i+1)+dz(i)))
         central=(c(i+1)-c(i-1))/(dz(i)+0.5d0*(dz(i-1)+dz(i+1)))
         if(a*b>0d0) slope(i)=sign(min(2d0*abs(a),2d0*abs(b),abs(central), &
              2d0*max(c(i),0d0)/dz(i)),a)
      end do
      f=0d0
      f(0)=min(velocity(0),0d0)*c(1)
      do i=1,n-1
         if(velocity(i)>=0d0) then
            edge=c(i)+0.5d0*dz(i)*slope(i)
         else
            edge=c(i+1)-0.5d0*dz(i+1)*slope(i+1)
         end if
         f(i)=velocity(i)*edge &
              -conductance(i)*(c(i+1)/air(i+1)-c(i)/air(i))
      end do
    end subroutine flux
  end subroutine volc_transport

  !======================================================================
  ! Standalone optical pre-run configuration and execution
  !======================================================================

  subroutine patmo_volc_configurePreRun(config)
    type(volcano_prerun_settings),intent(out)::config
    config=volcano_prerun_settings()
  call read_volcano_prerun_config("volcano_prerun.in",config%eventFile,config%profileFile, &
       config%photoMetricFile,config%outputFile,config%ashProfileFile,config%outputTimeStep,config%outputTimeUnit, &
       config%outputTimeUnitSeconds,config%wavelengthStepNm,config%tauFloor,config%endAfterStart,config%earlyStep,config%earlyUntil)

  write(*,'(A)') ""
  write(*,'(A)') "Volcano pre-run"
  write(*,'(A)') "---------------"
  write(*,'(A)') "Input files"
  write(*,'(A,1X,A)') "  event_file:       ",trim(config%eventFile)
  write(*,'(A,1X,A)') "  profile_file:     ",trim(config%profileFile)
  write(*,'(A,1X,A)') "  photo_metric_file:",trim(config%photoMetricFile)
  write(*,'(A,1X,A)') "  output_file:      ",trim(config%outputFile)
  write(*,'(A,1X,A)') "  ash_profile_file: ",trim(config%ashProfileFile)
  write(*,'(A)') ""

  end subroutine patmo_volc_configurePreRun

  ! The executable loads the atmospheric profile and metric before this call.
  subroutine patmo_volc_runPreRun(config)
    type(volcano_prerun_settings),intent(in)::config
  call patmo_volc_loadEvents(trim(config%eventFile))
  call patmo_volc_setClockRunning(.true.)
  call patmo_volc_setTime(0d0)
  call patmo_volc_setEarlyOutput(config%earlyStep,config%earlyUntil)
  call patmo_volc_dumpOpticalDepth(trim(config%outputFile),config%outputTimeStep, &
       config%outputTimeUnitSeconds,trim(config%outputTimeUnit),config%wavelengthStepNm,config%tauFloor, &
       config%endAfterStart)
  call patmo_volc_dumpAshProfile(trim(config%ashProfileFile),config%outputTimeStep, &
       config%outputTimeUnitSeconds,trim(config%outputTimeUnit),config%tauFloor,config%endAfterStart)
  end subroutine patmo_volc_runPreRun

  !***************
  subroutine read_volcano_prerun_config(fname,eventFile,profileFile, &
       photoMetricFile,outputFile,ashProfileFile,outputTimeStep,outputTimeUnit, &
       outputTimeUnitSeconds,wavelengthStepNm,tauFloor,endAfterStart,earlyStep,earlyUntil)
    use patmo_constants
    implicit none
    character(len=*),intent(in)::fname
    character(len=*),intent(inout)::eventFile,profileFile,photoMetricFile
    character(len=*),intent(inout)::outputFile,ashProfileFile,outputTimeUnit
    real*8,intent(inout)::outputTimeStep,outputTimeUnitSeconds
    real*8,intent(inout)::wavelengthStepNm,tauFloor
    real*8,intent(inout)::endAfterStart,earlyStep,earlyUntil
    character(len=512)::line,value
    character(len=80)::key
    integer::unitIn,ios,eqPos,commentPos,lineNumber
    real*8::tmp
    logical::exists

    inquire(file=trim(fname),exist=exists)
    if(.not.exists) then
       print *,"WARNING: ",trim(fname)," not found; using volcano pre-run defaults."
       return
    end if

    unitIn = 94
    open(unitIn,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problem while opening ",trim(fname)
       stop
    end if

    lineNumber = 0
    do
       read(unitIn,'(A)',iostat=ios) line
       if(ios/=0) exit
       lineNumber = lineNumber + 1
       commentPos = patmo_volc_commentStart(line)
       if(commentPos>0) line = line(:commentPos-1)
       if(len_trim(line)==0) cycle

       eqPos = index(line,"=")
       if(eqPos<=0) then
          print *,"WARNING: volcano_prerun.in line ignored: ",lineNumber
          cycle
       end if

       key = patmo_volc_lower(adjustl(line(:eqPos-1)))
       value = adjustl(line(eqPos+1:))

       select case(trim(key))
       case("event_file","volcano_events_file","volcano_file")
          eventFile = trim(value)
       case("profile_file","atmosphere_profile_file")
          profileFile = trim(value)
       case("photo_metric_file","photometric_file","wavelength_grid_file")
          photoMetricFile = trim(value)
       case("output_file","optical_depth_file","tau_file")
          outputFile = trim(value)
       case("ash_profile_file","output_ash_profile_file", &
            "plume_profile_file","output_plume_profile_file")
          ashProfileFile = trim(value)
       case("output_time_step_day","time_step_day","dt_day")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          outputTimeStep = tmp * secondsPerDay
          outputTimeUnit = "day"
          outputTimeUnitSeconds = secondsPerDay
       case("output_time_step_hour","time_step_hour","dt_hour")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          outputTimeStep = tmp * secondsPerHour
          outputTimeUnit = "hour"
          outputTimeUnitSeconds = secondsPerHour
       case("output_time_step_s","output_time_step_sec", &
            "time_step_s","time_step_sec","dt_s","dt_sec")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          outputTimeStep = tmp
          outputTimeUnit = "s"
          outputTimeUnitSeconds = 1d0
       case("output_wavelength_step_nm","wavelength_step_nm", &
            "lambda_step_nm","dlambda_nm")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          wavelengthStepNm = tmp
       case("output_early_time_step_s")
          read(value,*,iostat=ios) earlyStep
          if(ios/=0.or.earlyStep<0d0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
       case("output_early_until_hour")
          read(value,*,iostat=ios) tmp
          if(ios/=0.or.tmp<0d0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          earlyUntil = tmp * secondsPerHour
       case("tau_floor","tau_threshold","tau_stop")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          tauFloor = tmp
       case("output_end_day","end_day","until_day")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          if(tmp>0d0) then
             endAfterStart = tmp * secondsPerDay
          else
             endAfterStart = -1d0
          end if
       case("output_end_hour","end_hour","until_hour")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          if(tmp>0d0) then
             endAfterStart = tmp * secondsPerHour
          else
             endAfterStart = -1d0
          end if
       case("output_end_s","output_end_sec","end_s","end_sec", &
            "until_s","until_sec")
          read(value,*,iostat=ios) tmp
          if(ios/=0) call patmo_volc_prerunBadValue(lineNumber,trim(key))
          if(tmp>0d0) then
             endAfterStart = tmp
          else
             endAfterStart = -1d0
          end if
       case default
          print *,"WARNING: unknown volcano pre-run key ignored: ",trim(key)
       end select
    end do
    close(unitIn)

    if(outputTimeStep<=0d0) then
       print *,"ERROR: output time step must be positive in ",trim(fname)
       stop
    end if
    if(tauFloor<=0d0) then
       print *,"ERROR: tau_floor must be positive in ",trim(fname)
       stop
    end if

  end subroutine read_volcano_prerun_config

  !***************
  subroutine patmo_volc_prerunBadValue(lineNumber,key)
    implicit none
    integer,intent(in)::lineNumber
    character(len=*),intent(in)::key

    print *,"ERROR: bad value for ",trim(key)," in volcano_prerun.in line ", &
         lineNumber
    stop

  end subroutine patmo_volc_prerunBadValue


  !======================================================================
  ! Formal species, reaction and spectral diagnostics
  !======================================================================
  subroutine patmo_volc_historyConfigure(model_end)
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real*8,intent(in)::model_end
    logical::exists
    integer::u,ios,i,j
    real*8::solar_gb
    character(len=512)::message
    inquire(file='volcano_history.in',exist=exists)
    if(exists) then
       open(newunit=u,file='volcano_history.in',status='old',action='read')
       read(u,nml=history_output,iostat=ios,iomsg=message)
       close(u)
       if(ios/=0) then
          print *,trim(message)
          error stop 'Invalid volcano_history.in namelist'
       end if
    end if
    if(all(species_names=='')) species_names(1)='SO2'
    ns=0
    do i=1,size(species_names)
       if(len_trim(species_names(i))==0) cycle
       ns=ns+1
       species_indices(ns)=getSpeciesIndex(trim(species_names(i)),error=.false.)
       if(species_indices(ns)<1.or.species_indices(ns)>chemSpeciesNumber) then
          print *, 'Unknown chemical species: ',trim(species_names(i))
          error stop 1
       end if
       do j=1,ns-1
          if(species_indices(j)==species_indices(ns)) error stop 'Duplicate output species'
       end do
    end do
    ! Blank reaction selection means all reactions in the current network.
    if(all(reaction_ids==0)) reaction_ids=[(i,i=1,reactionsNumber)]
    nr=count(reaction_ids/=0)
    reaction_ids=pack(reaction_ids,reaction_ids/=0,vector=reaction_ids*0)
    do i=1,nr
       if(reaction_ids(i)<1.or.reaction_ids(i)>reactionsNumber) error stop 'Invalid output reaction ID'
       do j=1,i-1
          if(reaction_ids(j)==reaction_ids(i)) error stop 'Duplicate output reaction ID'
       end do
    end do
    every=[species_every_s,reaction_every_s,solar_every_s]
    finish=[species_end_s,reaction_end_s,solar_end_s]
    if(.not.all(ieee_is_finite(every)).or..not.all(ieee_is_finite(finish))) &
         error stop 'History times must be finite'
    enabled=every>0d0
    if(any(enabled.and.finish<0d0)) error stop 'History end seconds must be nonnegative'
    eruption_start=patmo_volc_historyFirstEventStart('volcano_events.dat')
    if(any(enabled.and.(eruption_start+finish>model_end))) &
         error stop 'History window exceeds volcano_tend: shorten output end or extend model duration'
    if(len_trim(history_prefix)==0) error stop 'Empty history prefix'
    next_sample=huge(1d0)
    where(enabled) next_sample=eruption_start
    print '(A,F14.3)', 'History origin on volcano clock [s]: ',eruption_start
    print '(A,3ES14.5)', 'History intervals, species/reactions/solar [s]: ',every
    print '(A,3ES14.5)', 'History end after eruption [s]: ',finish
    print '(A)', 'History settings: volcano_history.in (optional); output begins after spin-up.'
    if(enabled(3)) then
       solar_gb=(ceiling(finish(3)/every(3),kind=8)+1d0)*cellsNumber &
            *(photoBinsNumber*25d0+60d0)/1d9
       print '(A,F12.3)', 'Estimated solar history text size [GB]: ',solar_gb
       if(solar_gb>10d0) print '(A)', 'WARNING: dense full-spectrum output exceeds 10 GB; consider a longer solar interval.'
    end if
  end subroutine patmo_volc_historyConfigure

  !Read the same start-time aliases as patmo_volc, including legacy numeric rows.
  function patmo_volc_historyFirstEventStart(filename) result(first)
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*),intent(in)::filename
    real*8::first,start,scale
    character(len=8192)::line,work
    character(len=128)::key,value
    integer::u,ios,i,p
    logical::found
    first=huge(1d0)
    open(newunit=u,file=filename,status='old',action='read')
    do
       read(u,'(A)',iostat=ios) line
       if(ios/=0) exit
       p=scan(line,'#!')
       if(p>0) line=line(:p-1)
       if(len_trim(line)==0) cycle
       if(index(line,'=')==0) then
          read(line,*,iostat=ios) start
          if(ios/=0) error stop 'Cannot read legacy event start'
          first=min(first,start*secondsPerDay)
          cycle
       end if
       work=adjustl(line)
       do i=1,len_trim(work)
          if(work(i:i)=='='.or.work(i:i)==','.or.iachar(work(i:i))==9) work(i:i)=' '
          p=iachar(work(i:i))
          if(p>=65.and.p<=90) work(i:i)=achar(p+32)
       end do
       found=.false.
       do while(len_trim(work)>0)
          call patmo_volc_historyToken(work,key)
          call patmo_volc_historyToken(work,value)
          scale=0d0
          select case(trim(key))
          case('start','start_day','start_days')
             scale=secondsPerDay
          case('start_hour','start_hours','start_hr')
             scale=secondsPerHour
          case('start_s','start_sec','start_second','start_seconds')
             scale=1d0
          end select
          if(scale==0d0) cycle
          read(value,*,iostat=ios) start
          if(ios/=0) error stop 'Cannot read event start'
          start=start*scale
          found=.true.
       end do
       if(.not.found) error stop 'Event has no recognized start time'
       first=min(first,start)
    end do
    close(u)
    if(first==huge(1d0).or.first<0d0) error stop 'No nonnegative eruption start found'
  end function patmo_volc_historyFirstEventStart

  subroutine patmo_volc_historyToken(work,value)
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=*),intent(inout)::work
    character(len=*),intent(out)::value
    integer::p
    work=adjustl(work)
    p=index(trim(work),' ')
    if(p==0) then
       value=trim(work)
       work=''
    else
       value=work(:p-1)
       work=adjustl(work(p+1:))
    end if
  end subroutine patmo_volc_historyToken

  subroutine patmo_volc_historyBegin()
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(len=12),parameter::suffix(3)=[character(len=12)::'species','reactions','solar']
    character(len=maxNameLength)::names(speciesNumber)
    integer::g,u,r,m,order,refs(3)
    real*8::width(photoBinsNumber)
    names=getSpeciesNames()
    do g=1,3
       if(.not.enabled(g)) cycle
       open(newunit=units(g),file=trim(history_prefix)//'_'//trim(suffix(g))//'.dat',status='replace')
       write(units(g),'(A)') '# time_s is seconds since first eruption, after background steady state'
       write(units(g),'(A,ES24.15E3)') '# eruption_start_on_volcano_clock_s ',eruption_start
       write(units(g),'(A,2ES24.15E3)') '# interval_s end_s ',every(g),finish(g)
    end do
    if(enabled(1)) write(units(1),'(A)') &
         '# columns: time_s layer altitude_km species number_density_cm-3 mixing_ratio_pptv'
    if(enabled(2)) then
       write(units(2),'(A)') '# instantaneous coefficients and event rates, recomputed at output state'
       write(units(2),'(A)') '# NOT interval-integrated flux or species net loss; see reaction_map for stoichiometry'
       write(units(2),'(A)') '# columns: time_s layer altitude_km reaction_id coefficient_kind k_or_J reaction_rate_cm-3_s-1'
    end if
    if(enabled(3)) then
       write(units(3),'(A)') '# effective direct photon spectral flux used in current J calculation, NOT multiple-scattering actinic flux'
       write(units(3),'(A)') '# units inherit solar_flux.txt; assumed photons cm-2 s-1 nm-1'
       width=photoWavelengthWidths()
       write(units(3),'(A,2ES24.15E3)') '# mu coef ',photo_mu,photo_coef
       write(units(3),'(A)') '# per-bin integration widths are in the companion wavelengths file'
       write(units(3),'(A)',advance='no') '# columns: time_s layer altitude_km'
       do m=1,photoBinsNumber
          write(units(3),'(A,I0)',advance='no') ' flux_bin_',m
       end do
       write(units(3),*)
       open(newunit=u,file=trim(history_prefix)//'_wavelengths.dat',status='replace')
       write(u,'(A)') '# columns: bin wavelength_nm incident_flux model_integration_width_nm'
       do m=1,photoBinsNumber
          write(u,'(I0,3(1X,ES24.15E3))') m,1d7*planck_eV*clight/energyMid(m),photoFlux(m),width(m)
       end do
       close(u)
    end if
    open(newunit=u,file=trim(history_prefix)//'_reaction_map.dat',status='replace')
    write(u,'(A)') '# ID | kind | number of density factors | coefficient units | reaction'
    do r=1,reactionsNumber
       refs=[indexReactants1(r),indexReactants2(r),indexReactants3(r)]
       order=count(refs/=positionDummy)
       if(r>chemReactionsNumber.and.r<=chemReactionsNumber+photoReactionsNumber) then
          write(u,'(I0,A,I0,A,A)') r,' | J | ',order,' | s-1 | ',trim(reactionsVerbatim(r))
       else
          write(u,'(I0,A,I0,A,I0,A,A)') r,' | k | ',order,' | cm^',3*(order-1),' s-1 | ',trim(reactionsVerbatim(r))
       end if
    end do
    close(u)
    opened=.true.
    call patmo_volc_historySample(0d0)
  end subroutine patmo_volc_historyBegin

  function patmo_volc_historyLimit(time,requested) result(step)
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real*8,intent(in)::time,requested
    real*8::step
    step=min(requested,max(minval(next_sample)-time,0d0))
    if(step<=0d0) error stop 'History scheduler has an unconsumed output time'
  end function patmo_volc_historyLimit

  subroutine patmo_volc_historySample(time)
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real*8,intent(in)::time
    real*8::tau(photoBinsNumber,cellsNumber),flux(photoBinsNumber),air(cellsNumber),width(photoBinsNumber)
    real*8::saved_rates(cellsNumber,reactionsNumber),rates(reactionsNumber),jcheck,tolerance
    logical::due(3)
    character(len=maxNameLength)::names(speciesNumber)
    character(len=1)::kind
    integer::g,i,j,r,p
    if(.not.opened) return
    tolerance=1d-8
    due=enabled.and.abs(time-next_sample)<=tolerance
    if(.not.any(due)) return
    if(due(3)) width=photoWavelengthWidths()
    names=getSpeciesNames()
    air=0.5d0*sum(nall(:,1:chemSpeciesNumber),2)
    if(any(air<=0d0)) error stop 'Nonpositive air density at history output'
    if(due(2).or.due(3)) then
       !Reproduce PATMO gas opacity and volcanic opacity at the current state,
       !not the previous outer step. Leave the integrator's stored rates intact.
       tau=0d0
       do j=cellsNumber-1,1,-1
          tau(:,j)=tau(:,j+1)
          do p=1,photoReactionsNumber
             tau(:,j)=tau(:,j)+gridSpace(j)*xsecAll(:,p)*nall(j,photoPartnerIndex(p))
          end do
       end do
       call patmo_volc_applyAshOpacity(tau)
       saved_rates=krate
       call computeRates(TgasAll)
       call computeReverseRates(TgasAll)
       call computePhotoRates(tau)
    end if
    do j=1,cellsNumber
       if(due(1)) then
          do i=1,ns
             r=species_indices(i)
             write(units(1),'(ES24.15E3,1X,I0,1X,ES24.15E3,1X,A,2(1X,ES24.15E3))') &
                  time-eruption_start,j,height(j)/1d5,trim(names(r)),nall(j,r),nall(j,r)/air(j)*1d12
          end do
       end if
       if(due(2)) then
          rates=getFlux(nall,j)
          do i=1,nr
             r=reaction_ids(i)
             kind='k'
             if(r>chemReactionsNumber.and.r<=chemReactionsNumber+photoReactionsNumber) kind='J'
             write(units(2),'(ES24.15E3,1X,I0,1X,ES24.15E3,1X,I0,1X,A,2(1X,ES24.15E3))') &
                  time-eruption_start,j,height(j)/1d5,r,kind,krate(j,r),rates(r)
          end do
       end if
       if(due(3)) then
          flux=photoDirectFlux(tau(:,j))
          do p=1,photoReactionsNumber
             jcheck=sum(xsecAll(:,p)*flux*width)
             if(abs(jcheck-krate(j,chemReactionsNumber+p))> &
                  1d-10*max(abs(jcheck),abs(krate(j,chemReactionsNumber+p)),1d-99)) &
                  error stop 'Solar output formula differs from current photoRates: update diagnostic constants'
          end do
          write(units(3),'(ES24.15E3,1X,I0,1X,ES24.15E3,*(1X,ES24.15E3))') &
               time-eruption_start,j,height(j)/1d5,flux
       end if
    end do
    if(due(2).or.due(3)) krate=saved_rates
    do g=1,3
       if(.not.due(g)) cycle
       flush(units(g))
       if(time>=eruption_start+finish(g)-tolerance) then
          next_sample(g)=huge(1d0)
       else
          sample_index(g)=sample_index(g)+1
          next_sample(g)=eruption_start+min(dble(sample_index(g))*every(g),finish(g))
       end if
    end do
  end subroutine patmo_volc_historySample

  subroutine patmo_volc_historyClose()
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  use patmo_utils
  use patmo_rates
  use patmo_reverseRates
  use patmo_photoRates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer::g
    if(.not.opened) return
    do g=1,3
       if(enabled(g)) close(units(g))
    end do
    opened=.false.
  end subroutine patmo_volc_historyClose

end module patmo_volc
