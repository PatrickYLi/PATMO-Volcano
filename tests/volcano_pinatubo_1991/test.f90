program test
  use patmo
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  use patmo_ode
  use patmo_utils
  use patmo_volc
  implicit none
  real*8::dt,x(speciesNumber),spinup_time,volcano_time
  real*8::spinup_tend,volcano_tend,imass,one_year
  real*8::spinup_change
  real*8::heff(chemSpeciesNumber)
  real*8::dep(chemSpeciesNumber) !cm/s
  real*8::total_elapsed_ms
  integer::icell,i,j
  integer(kind=8)::clock_start,clock_finish,clock_rate
  real*8::convergence = 100.0
  logical::background_equilibrated

  !init photochemistry
  call patmo_volc_readRunSettings()
  call patmo_init()

  !load temperature and density profile
  call patmo_loadInitialProfile("profile.dat",unitH="km",unitX="1/cm3")
  call patmo_volc_loadEvents("volcano_events.dat")
  call patmo_volc_holdUntilEquilibrium()

  !set BB flux (default Sun@1AU)
  call patmo_setFluxBB()
  call patmo_setGravity(9.8d2)
  !call patmo_setEddyKzzAll(1.0d5)
  
  !read initial value 
  call patmo_dumpHydrostaticProfile("hydrostat.out")
  wetdep(:,:) = 0d0 

  !calculate wet deposition
  call computewetdep(patmo_idx_COS,2.0d-2)   !OCS
  call computewetdep(patmo_idx_CS2,5.0d-2)  !CS2
	call computewetdep(patmo_idx_H2S,1.0d-1)  !H2S
	call computewetdep(patmo_idx_SO2,4.0d3)  !SO2
  ! call computewetdep(patmo_idx_H2SO4,5d14)  !H2SO4
  call computewetdep(patmo_idx_SO4,5d14)  !SO4

  !Turco H2SO4
  wetdep(12,patmo_idx_H2SO4) = 1.77E-06
  wetdep(11,patmo_idx_H2SO4) = 3.54E-06
  wetdep(10,patmo_idx_H2SO4) = 5.31E-06
  wetdep(9,patmo_idx_H2SO4) = 7.08E-06
  wetdep(8,patmo_idx_H2SO4) = 8.85E-06
  wetdep(7,patmo_idx_H2SO4) =  1.06E-05
  wetdep(6,patmo_idx_H2SO4) = 1.24E-05
  wetdep(5,patmo_idx_H2SO4) = 1.42E-05
  wetdep(4,patmo_idx_H2SO4) = 1.59E-05
  wetdep(3,patmo_idx_H2SO4) = 1.77E-05
  wetdep(2,patmo_idx_H2SO4) = 1.95E-05
  wetdep(1,patmo_idx_H2SO4) = 2.12E-05

  va(:) = 0d0  
  pa(:) = 0d0
  
  open(60,file="vapor_H2SO4.txt",status="old")  
    do i=13,34
      read(60,*) va(i)
	  end do
  close(60)

  open(61,file="partial_H2SO4.txt",status="old") 
    do i=13,34
      read(61,*) pa(i)
	  end do
  close(61)

  gd(:) = 0d0

  open(62,file="SO4_deposition_rate.txt",status="old")  
    do i=1,60
      read(62,*) gd(i)
	  end do
  close(62)

  !get initial mass, g/cm3
    imass = patmo_getTotalMass()
    print*,"mass:",imass
    print *,"Volcano forcing is held while the background reaches steady state."
 
  !first spin up the no-volcano background atmosphere to steady state
  dt = secondsPerDay 
  call patmo_volc_getRunLimits(spinup_tend,volcano_tend)
  call patmo_volc_historyConfigure(volcano_tend)
  one_year = secondsPerDay*365
  spinup_time = 0d0
  volcano_time = 0d0
  background_equilibrated = .false.
  call patmo_volc_beginEquilibriumCheck()
  spinup_change = 1d99
  call system_clock(clock_start, clock_rate)

  do
    call patmo_run(dt,convergence)
    spinup_time = spinup_time + dt
    background_equilibrated=patmo_volc_checkEquilibrium(spinup_change)
    if(background_equilibrated) exit
    print '(A,F11.2,a,ES12.4)', &
         "spin-up ",spinup_time/spinup_tend*1d2," %, rel_change=",spinup_change
    if(spinup_time>=spinup_tend) exit
  end do

  call patmo_volc_requireEquilibrium(background_equilibrated)

  print *,"Background steady state reached at day ",spinup_time/secondsPerDay
  print *,"Background relative daily change = ",spinup_change
  call patmo_dumpFinalSpeciesDensityCsv("background_equilibrium_species_number_density.csv")

  !now reset the volcano clock: start_day=0 means eruption starts after spin-up
  call patmo_volc_startAfterEquilibrium()
  call patmo_volc_dumpState("volcano_initial_state.dat")
  call patmo_volc_historyBegin()

  do
    dt=patmo_volc_limitStep(min(secondsPerDay,volcano_tend-volcano_time))
    dt=patmo_volc_historyLimit(volcano_time,dt)
    call patmo_run(dt,convergence)
    volcano_time = volcano_time + dt
    call patmo_volc_historySample(volcano_time)
    print '(A,F11.2,a2)',"volcano ",volcano_time/volcano_tend*1d2," %"
    if(volcano_time>=volcano_tend) exit
  end do
  call patmo_volc_historyClose()
  call patmo_dumpFinalSpeciesDensityCsv("final_species_number_density.csv")

  total_elapsed_ms = 0d0
  call system_clock(clock_finish)
  if(clock_rate>0) then
    total_elapsed_ms = 1d3 * dble(clock_finish-clock_start) / dble(clock_rate)
  end if
  call patmo_printElapsedTime("Total model elapsed: ", total_elapsed_ms)
  call patmo_volc_dumpState("volcano_final_state.dat")
  
 !get mass, g/cm3
  imass = patmo_getTotalMass()
  print *,"mass:",imass
  
  !dump final hydrostatic equilibrium
  call patmo_dumpHydrostaticProfile("hydrostatEnd.out")
  call patmo_dumpJValue("jvalue.dat")
  call patmo_dumpOpacity("opacity.dat")
  call patmo_dumpAllRates("rates.dat")
  call patmo_dumpAllMixingRatioToFile("allNDs.dat")

end program test
!**************

subroutine computewetdep(i,heff)
  use patmo_commons
  use patmo_constants
  use patmo_parameters
  implicit none
  real*8::gamma(cellsNumber)  ! Precipation and Nonprecipitation time 
  real*8::wh2o(cellsNumber)   ! Rate of wet removal
  real*8::rkj(cellsNumber,chemSpeciesNumber)    ! Average Remeval Frequency	  
  real*8::y(cellsNumber),fz(cellsNumber),wl,qj(cellsNumber,chemSpeciesNumber)
	real*8::heff !Henry's Constant
  real*8::zkm(cellsNumber),temp(cellsNumber),gam15,gam8
	integer::i,j
	!Gas constant
  real*8,parameter::Rgas_atm = 1.36d-22 !Latm/K/mol

  gam15 = 8.64d5/2d0
  gam8 = 7.0d6/2d0
  wl = 1d0
  zkm(:) = height(:)/1d5
	temp(:) = TgasAll(:)
	 
  do j=1,12

	!FIND APPROPRIATE GAMMA
  if (zkm(j).LE.1.51d0) then
    gamma(j) = gam15
  else if (zkm(j).LT.8d0) then
    gamma(j) = gam15 + (gam8-gam15)*((zkm(j)-1.5d0)/6.5d0)
  else
    gamma(j) = gam8
  end if 

  !FIND WH2O
  if (zkm(j).LE.1d0) then
    y(j) = 11.35d0 + 1d-1*zkm(j)
  else
    y(j) = 11.5444d0 - 0.085333d0*zkm(j) - 9.1111d-03*zkm(j)*zkm(j)
  end if
  wh2o(j) = 10d0**y(j)

	!FIND F(Z)
  if (zkm(j).LE.1.51d0) then
     fz(j) = 1d-1
  else
     fz(j) = 0.16615d0 - 0.04916d0*zkm(j) + 3.37451d-3*zkm(j)*zkm(j)
  end if
		
	   !Raintout rates
        rkj(j,i) = wh2o(j)/55d0 /(av*wl*1.0d-9 + 1d0/(heff*Rgas_atm*temp(j)))
	    qj(j,i) = 1d0 - fz(j) + fz(j)/(gamma(j)*rkj(j,i)) * (1d0 - EXP(-rkj(j,i)*gamma(j)))
        wetdep(j,i) = (1d0 - EXP(-rkj(j,i)*gamma(j)))/(gamma(j)*qj(j,i))
	end do

	!Output 
   if (i==patmo_idx_COS) then
  	do j=1,12
      open  (20,file="Rainout-OCS.txt")
      write (20,*) 'GAMMA', gamma(j)
      write (20,*) 'ZKM', zkm(j)
      write (20,*) 'WH2O', wh2o(j)
	  write (20,*) 'RKJ', rkj(j,i)
 	  write (20,*) 'QJ', qj(j,i)
	  write (20,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if
	 	 
   if (i==patmo_idx_CS2) then
	do j=1,12
      open  (25,file="Rainout-CS2.txt")
      write (25,*) 'ZKM', zkm(j)
	  write (25,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
	 	 
   if (i==patmo_idx_H2S) then
 	do j=1,12
      open  (26,file="Rainout-H2S.txt")
      write (26,*) 'ZKM', zkm(j)
	  write (26,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 	
   
   if (i==patmo_idx_SO2) then
   	do j=1,12
      open  (24,file="Rainout-SO2.txt")
      write (24,*) 'ZKM', zkm(j)
	  write (24,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
   
   if (i==patmo_idx_H2SO4) then
   	do j=1,12
      open  (22,file="Rainout-H2SO4.txt")
      write (22,*) 'ZKM', zkm(j)
	  write (22,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if 
   
   if (i==patmo_idx_SO4) then
	do j=1,12
      open  (21,file="Rainout-SO4.txt")
      write (21,*) 'ZKM', zkm(j)
	  write (21,*) 'K_RAIN',  wetdep(j,i)
    end do
   end if

   end subroutine computewetdep
	!**************

  !***************
  !dump final number densities to CSV in the requested species order.
 subroutine patmo_dumpFinalSpeciesDensityCsv(fname)
    use patmo_commons
    use patmo_parameters
    implicit none
    character(len=*),intent(in)::fname
    integer,parameter::dumpSpeciesNumber = 18
    integer,dimension(dumpSpeciesNumber),parameter::speciesIdx = (/ &
         patmo_idx_COS, patmo_idx_SH, patmo_idx_SO,patmo_idx_CS2, &
         patmo_idx_CS, patmo_idx_S, patmo_idx_S2, patmo_idx_SCSOH, &
         patmo_idx_HSO2, patmo_idx_H2S, patmo_idx_SO2, patmo_idx_SO3, &
         patmo_idx_HSO3, patmo_idx_H2SO4, patmo_idx_SO4, patmo_idx_CH3SCH3, &
         patmo_idx_CH4O3S, patmo_idx_CS2E /)
    character(len=32)::value
    integer::i,j

    open(67,file=trim(fname),status="replace")
    write(67,'(A)') "COS,SH,SO,CS2," // &
         "CS,S,S2,SCSOH," // &
         "HSO2,H2S,SO2,SO3," // &
         "HSO3,H2SO4,SO4,CH3SCH3," // &
         "CH4O3S,CS2E"

    do i=1,cellsNumber
       do j=1,dumpSpeciesNumber
          if(j>1) write(67,'(A)',advance='no') ","
          write(value,'(E17.8E3)') nall(i,speciesIdx(j))
          write(67,'(A)',advance='no') trim(adjustl(value))
       end do
       write(67,*)
    end do
    close(67)

    print *, "Final species number densities dumped in ", trim(fname)

 end subroutine patmo_dumpFinalSpeciesDensityCsv
