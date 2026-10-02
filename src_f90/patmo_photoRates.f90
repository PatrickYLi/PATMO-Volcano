module patmo_photoRates
  implicit none
  real*8,parameter :: photo_mu = #PATMO_zenith_mu
  real*8,parameter :: photo_coef = #PATMO_TOA_coef
contains

  function photoWavelengthWidths() result(width)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    real*8::width(photoBinsNumber)
    if(any(energyLeft<=0d0).or.any(energyRight<=energyLeft)) &
         error stop 'Invalid photochemistry bin edges'
    width=planck_eV*clight*1d7*(1d0/energyLeft-1d0/energyRight)
  end function photoWavelengthWidths

  !Effective direct photon spectrum [photons cm-2 s-1 nm-1], shared with diagnostics.
  function photoDirectFlux(tauLayer) result(flux)
    use patmo_parameters
    use patmo_commons
    implicit none
    real*8,intent(in)::tauLayer(photoBinsNumber)
    real*8::flux(photoBinsNumber)
    flux=0d0
    if(photo_mu>0d0) flux=photo_coef*photoFlux*exp(-tauLayer/photo_mu)
  end function photoDirectFlux

  !**************
  subroutine computePhotoRates(tau)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)

#PATMO_photoRates

  end subroutine computePhotoRates

  !*************
  function integrateXsec(index,tau)
    use patmo_parameters
    use patmo_commons
    use patmo_constants
    implicit none
    integer,intent(in)::index
    real*8,intent(in)::tau(photoBinsNumber,cellsNumber)
    real*8::integrateXsec(cellsNumber), width(photoBinsNumber)
    integer::j

    width=photoWavelengthWidths()
    !loop on cells (stride photobins)
    do j=1,cellsNumber
          integrateXsec(j) = sum(xsecAll(:,index)*photoDirectFlux(tau(:,j))*width)
    end do


  end function integrateXsec

end module patmo_photoRates
