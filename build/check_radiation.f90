program check_radiation
  use patmo
  use patmo_photoRates
  use patmo_photo
  use patmo_commons
  use patmo_parameters
  use patmo_constants
  implicit none
  real*8::width(photoBinsNumber),tau(photoBinsNumber,cellsNumber)
  real*8::flux(photoBinsNumber),expected,actual(cellsNumber)
  integer::j,p
  call patmo_init()
  call patmo_setFluxBB()
  width=photoWavelengthWidths()
  if(any(width<=0d0)) error stop 'Nonpositive wavelength width'
  expected=planck_eV*clight*1d7*(1d0/energyLeft(1)-1d0/energyRight(photoBinsNumber))
  call close(sum(width),expected,'Wavelength widths telescope')
  do j=1,cellsNumber
     tau(:,j)=0.01d0*j*energyMid/energyMid(1)
  end do
  call computePhotoRates(tau)
  do j=1,cellsNumber
     flux=photoDirectFlux(tau(:,j))
     if(photo_mu>0d0) then
        if(any(flux>photo_coef*photoFlux*(1d0+1d-12))) error stop 'Attenuated flux exceeds incident flux'
     else
        if(any(flux/=0d0)) error stop 'Nonzero nighttime direct flux'
     end if
     do p=1,photoReactionsNumber
        expected=sum(xsecAll(:,p)*flux*width)
        call close(krate(j,chemReactionsNumber+p),expected,'J from exported spectrum')
     end do
  end do
  !Independent constant-spectrum analytic integral, not just closure of two helpers.
  photoFlux=2d0
  xsecAll(:,1)=3d0
  tau=0d0
  expected=0d0
  if(photo_mu>0d0) expected=6d0*photo_coef*sum(width)
  actual=integrateXsec(1,tau)
  call close(actual(1),expected,'Constant spectrum integral')
  print *, 'All radiation checks passed.'
contains
  subroutine close(a,b,label)
    real*8,intent(in)::a,b
    character(len=*),intent(in)::label
    if(abs(a-b)>1d-10*max(abs(a),abs(b),1d-99)) then
       print *,label,a,b
       error stop 1
    end if
  end subroutine close
end program check_radiation
