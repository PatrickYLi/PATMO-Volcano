program check_equilibrium
  use patmo_commons
  use patmo_parameters
  use patmo_volc
  implicit none
  real*8::change,spinupEnd,volcanoEnd
  real*8::previous(cellsNumber,chemSpeciesNumber),scale(cellsNumber,chemSpeciesNumber),expected
  logical::ready
  integer::i
  character(len=16)::mode

  call get_command_argument(1,mode)
  call patmo_volc_readRunSettings()
  call patmo_volc_getRunLimits(spinupEnd,volcanoEnd)
  if(spinupEnd/=365d0*86400d0.or.volcanoEnd/=2d0*86400d0) error stop 'Wrong run limits'
  call patmo_volc_loadEvents('volcano_events.dat')
  call patmo_volc_holdUntilEquilibrium()
  call patmo_volc_advanceTime(86400d0)
  if(patmo_volc_isActive().or.patmo_volc_getTime()/=0d0) error stop 'Forcing active during spin-up'
  if(trim(mode)=='fail') then
     call patmo_volc_requireEquilibrium(.false.)
     error stop 'Failure gate returned'
  end if

  nall=1d10
  nall(:,patmo_idx_SO2)=1d-30
  call patmo_volc_beginEquilibriumCheck()
  do i=1,2
     ready=patmo_volc_checkEquilibrium(change)
     if(ready.or.change/=0d0) error stop 'Premature equilibrium'
  end do
  previous=nall(:,1:chemSpeciesNumber)
  nall(1,patmo_idx_SO2)=1d-7
  scale=spread(max(0.5d0*sum(nall(:,1:chemSpeciesNumber),2),1d-99), &
       2,chemSpeciesNumber)*1d-15
  expected=maxval(abs(nall(:,1:chemSpeciesNumber)-previous) &
       /max(abs(nall(:,1:chemSpeciesNumber)),abs(previous),scale))
  ready=patmo_volc_checkEquilibrium(change)
  if(abs(change-expected)>1d-14*expected.or.ready) error stop 'Changed trace floor or metric'
  do i=1,3
     ready=patmo_volc_checkEquilibrium(change)
     if(ready.neqv.(i==3)) error stop 'Stable counter did not reset'
  end do
  call patmo_volc_requireEquilibrium(ready)
  call patmo_volc_startAfterEquilibrium()
  if(.not.patmo_volc_isActive().or.patmo_volc_getTime()/=0d0) error stop 'Incorrect eruption start'
  call patmo_volc_beginEquilibriumCheck()
  ready=patmo_volc_checkEquilibrium(change)
  if(ready) error stop 'Monitor did not reset'
  print *, 'Equilibrium gate checks passed.'
end program check_equilibrium
