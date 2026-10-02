module patmo_photoRates
  implicit none
  real*8,parameter :: photo_mu = 5.0000000000000011d-01
  real*8,parameter :: photo_coef = 5.0000000000000000d-01
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

    !COS -> CO + S
    krate(:,44) = integrateXsec(1, tau(:,:))

    !OH -> O + H
    krate(:,45) = integrateXsec(2, tau(:,:))

    !CO2 -> CO + O
    krate(:,46) = integrateXsec(3, tau(:,:))

    !SO -> S + O
    krate(:,47) = integrateXsec(4, tau(:,:))

    !CS2 -> CS + S
    krate(:,48) = integrateXsec(5, tau(:,:))

    !O2 -> O + O
    krate(:,49) = integrateXsec(6, tau(:,:))

    !O3 -> O2 + O
    krate(:,50) = integrateXsec(7, tau(:,:))

    !H2S -> SH + H
    krate(:,51) = integrateXsec(8, tau(:,:))

    !H2O -> OH + H
    krate(:,52) = integrateXsec(9, tau(:,:))

    !H2O -> H2 + O
    krate(:,53) = integrateXsec(10, tau(:,:))

    !H2 -> H + H
    krate(:,54) = integrateXsec(11, tau(:,:))

    !HO2 -> OH + O
    krate(:,55) = integrateXsec(12, tau(:,:))

    !SO2 -> SO + O
    krate(:,56) = integrateXsec(13, tau(:,:))

    !SO3 -> SO2 + O
    krate(:,57) = integrateXsec(14, tau(:,:))

    !H2SO4 -> SO2 + OH + OH
    krate(:,58) = integrateXsec(15, tau(:,:))

    !N2 -> N + N
    krate(:,59) = integrateXsec(16, tau(:,:))

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
