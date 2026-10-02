module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    use patmo_volc, only: patmo_volc_addSources
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0
    dn(:,patmo_idx_COS) = &
        - krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,8)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        + krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,44)*n(:,patmo_idx_COS) &
        + krate(:,60)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        + krate(:,61)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        - krate(:,64)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,67)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        - krate(:,69)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,71)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_SH) = &
        + krate(:,1)*n(:,patmo_idx_COS)*n(:,patmo_idx_OH) &
        + krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        + krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        + krate(:,15)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        - krate(:,18)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,19)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        + krate(:,30)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        + krate(:,51)*n(:,patmo_idx_H2S) &
        - krate(:,60)*n(:,patmo_idx_CO2)*n(:,patmo_idx_SH) &
        - krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        - krate(:,72)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        - krate(:,73)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        - krate(:,74)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,76)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        + krate(:,77)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,78)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,89)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO) = &
        + krate(:,2)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_SH)*n(:,patmo_idx_O) &
        + krate(:,18)*n(:,patmo_idx_SH)*n(:,patmo_idx_O2) &
        - krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,22)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        + krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        + krate(:,25)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,41)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,47)*n(:,patmo_idx_SO) &
        + krate(:,56)*n(:,patmo_idx_SO2) &
        - krate(:,61)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        - krate(:,63)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        - krate(:,70)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO) &
        - krate(:,76)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        - krate(:,77)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,79)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,81)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        - krate(:,82)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        - krate(:,83)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,84)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        - krate(:,100)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CS2) = &
        - krate(:,3)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,6)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,7)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_M) &
        - krate(:,48)*n(:,patmo_idx_CS2) &
        + krate(:,62)*n(:,patmo_idx_COS)*n(:,patmo_idx_SH) &
        + krate(:,63)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,64)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        + krate(:,65)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO) &
        + krate(:,66)*n(:,patmo_idx_SCSOH) &
        - krate(:,101)*n(:,patmo_idx_CS2)*n(:,patmo_idx_M)

    dn(:,patmo_idx_CS) = &
        + krate(:,4)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,10)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,11)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,12)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        + krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,48)*n(:,patmo_idx_CS2) &
        - krate(:,63)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO) &
        + krate(:,68)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
        + krate(:,69)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,70)*n(:,patmo_idx_SO)*n(:,patmo_idx_CO) &
        + krate(:,71)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        - krate(:,102)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_S) = &
        + krate(:,5)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        + krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O) &
        - krate(:,23)*n(:,patmo_idx_S)*n(:,patmo_idx_O2) &
        - krate(:,24)*n(:,patmo_idx_S)*n(:,patmo_idx_O3) &
        - krate(:,25)*n(:,patmo_idx_S)*n(:,patmo_idx_OH) &
        + krate(:,41)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        + krate(:,44)*n(:,patmo_idx_COS) &
        + krate(:,47)*n(:,patmo_idx_SO) &
        + krate(:,48)*n(:,patmo_idx_CS2) &
        - krate(:,64)*n(:,patmo_idx_COS)*n(:,patmo_idx_S) &
        - krate(:,68)*n(:,patmo_idx_S)*n(:,patmo_idx_CO) &
        + krate(:,82)*n(:,patmo_idx_SO)*n(:,patmo_idx_O) &
        + krate(:,83)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,84)*n(:,patmo_idx_SO)*n(:,patmo_idx_H) &
        - krate(:,100)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_S2) = &
        + krate(:,6)*n(:,patmo_idx_CS2)*n(:,patmo_idx_O) &
        - krate(:,41)*n(:,patmo_idx_S2)*n(:,patmo_idx_O) &
        - krate(:,65)*n(:,patmo_idx_S2)*n(:,patmo_idx_CO) &
        + krate(:,100)*n(:,patmo_idx_S)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_SCSOH) = &
        + krate(:,7)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,8)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,66)*n(:,patmo_idx_SCSOH) &
        + krate(:,67)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2)

    dn(:,patmo_idx_HSO2) = &
        + krate(:,8)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,31)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,67)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        + krate(:,90)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2)

    dn(:,patmo_idx_H2S) = &
        - krate(:,13)*n(:,patmo_idx_H2S)*n(:,patmo_idx_OH) &
        - krate(:,14)*n(:,patmo_idx_H2S)*n(:,patmo_idx_O) &
        - krate(:,15)*n(:,patmo_idx_H2S)*n(:,patmo_idx_H) &
        - krate(:,16)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        - krate(:,51)*n(:,patmo_idx_H2S) &
        + krate(:,72)*n(:,patmo_idx_H2O)*n(:,patmo_idx_SH) &
        + krate(:,73)*n(:,patmo_idx_OH)*n(:,patmo_idx_SH) &
        + krate(:,74)*n(:,patmo_idx_H2)*n(:,patmo_idx_SH) &
        + krate(:,75)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO)

    dn(:,patmo_idx_HSO) = &
        + krate(:,16)*n(:,patmo_idx_H2S)*n(:,patmo_idx_HO2) &
        + krate(:,19)*n(:,patmo_idx_SH)*n(:,patmo_idx_O3) &
        - krate(:,29)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        - krate(:,30)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O3) &
        - krate(:,75)*n(:,patmo_idx_H2O)*n(:,patmo_idx_HSO) &
        - krate(:,78)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,88)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        + krate(:,89)*n(:,patmo_idx_O2)*n(:,patmo_idx_O2)*n(:,patmo_idx_SH)

    dn(:,patmo_idx_SO2) = &
        + krate(:,20)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,21)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,22)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,27)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,29)*n(:,patmo_idx_HSO)*n(:,patmo_idx_O2) &
        + krate(:,31)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,33)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,37)*n(:,patmo_idx_SO2) &
        + krate(:,38)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        + krate(:,39)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,56)*n(:,patmo_idx_SO2) &
        + krate(:,57)*n(:,patmo_idx_SO3) &
        + krate(:,58)*n(:,patmo_idx_H2SO4) &
        - krate(:,79)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,80)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,81)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H) &
        + krate(:,85)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        + krate(:,86)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        + krate(:,87)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO3) &
        - krate(:,88)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH) &
        - krate(:,90)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,92)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M) &
        + krate(:,96)*n(:,patmo_idx_SO4) &
        - krate(:,97)*n(:,patmo_idx_SO2) &
        - krate(:,98)*n(:,patmo_idx_SO2) &
        - krate(:,99)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S) &
        - krate(:,102)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_M) = &
        - krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        + krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        - krate(:,33)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,33)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        - krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_M) &
        + krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_M) &
        - krate(:,85)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        + krate(:,85)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        - krate(:,92)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M) &
        + krate(:,92)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M) &
        - krate(:,101)*n(:,patmo_idx_CS2)*n(:,patmo_idx_M) &
        + krate(:,101)*n(:,patmo_idx_CS2)*n(:,patmo_idx_M)

    dn(:,patmo_idx_SO3) = &
        + krate(:,26)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)*n(:,patmo_idx_M) &
        + krate(:,27)*n(:,patmo_idx_SO2)*n(:,patmo_idx_HO2) &
        + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O3) &
        + krate(:,32)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        - krate(:,34)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,57)*n(:,patmo_idx_SO3) &
        - krate(:,85)*n(:,patmo_idx_SO3)*n(:,patmo_idx_M) &
        - krate(:,86)*n(:,patmo_idx_OH)*n(:,patmo_idx_SO3) &
        - krate(:,87)*n(:,patmo_idx_O2)*n(:,patmo_idx_SO3) &
        - krate(:,91)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        + krate(:,93)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_HSO3) = &
        - krate(:,32)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_O2) &
        + krate(:,33)*n(:,patmo_idx_SO2)*n(:,patmo_idx_OH)*n(:,patmo_idx_M) &
        + krate(:,91)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO3) &
        - krate(:,92)*n(:,patmo_idx_HSO3)*n(:,patmo_idx_M)

    dn(:,patmo_idx_H2SO4) = &
        + krate(:,34)*n(:,patmo_idx_SO3)*n(:,patmo_idx_H2O) &
        - krate(:,58)*n(:,patmo_idx_H2SO4) &
        - krate(:,93)*n(:,patmo_idx_H2SO4)

    dn(:,patmo_idx_SO4) = &
        + krate(:,37)*n(:,patmo_idx_SO2) &
        - krate(:,96)*n(:,patmo_idx_SO4)

    dn(:,patmo_idx_CH3SCH3) = &
        - krate(:,38)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_O) &
        - krate(:,39)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        + krate(:,97)*n(:,patmo_idx_SO2) &
        + krate(:,98)*n(:,patmo_idx_SO2) &
        + krate(:,99)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_CH4O3S) = &
        + krate(:,40)*n(:,patmo_idx_CH3SCH3)*n(:,patmo_idx_OH) &
        - krate(:,99)*n(:,patmo_idx_SO2)*n(:,patmo_idx_CH4O3S)

    dn(:,patmo_idx_CS2E) = &
        - krate(:,42)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_M) &
        - krate(:,43)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,101)*n(:,patmo_idx_CS2)*n(:,patmo_idx_M) &
        + krate(:,102)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

    do i=1,chemSpeciesNumber
      dn(:,i) = dn(:,i) &
          + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
          - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
          + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
          + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
    end do

    !Chemical Species with constant concentration
    dn(:,patmo_idx_HO2) = 0d0
    dn(:,patmo_idx_N) = 0d0
    dn(:,patmo_idx_CO2) = 0d0
    dn(:,patmo_idx_H2O) = 0d0
    dn(:,patmo_idx_CO) = 0d0
    dn(:,patmo_idx_O2) = 0d0
    dn(:,patmo_idx_N2) = 0d0
    dn(:,patmo_idx_OH) = 0d0
    dn(:,patmo_idx_O) = 0d0
    dn(:,patmo_idx_H2) = 0d0
    dn(:,patmo_idx_H) = 0d0
    dn(:,patmo_idx_O3) = 0d0

    ! Gravity Settling
    do j = cellsNumber, 2, -1
      dn(j    , patmo_idx_SO4) = dn(j    , patmo_idx_SO4) - gd(j) * n(j, patmo_idx_SO4)
      dn(j - 1, patmo_idx_SO4) = dn(j - 1, patmo_idx_SO4) + gd(j) * n(j, patmo_idx_SO4)
    end do
    SO4SurFall = gd(j) * n(1, patmo_idx_SO4)
    dn(1, patmo_idx_SO4) = dn(1, patmo_idx_SO4) - SO4SurFall

    ! Dry Deposition: assumed a deposition rate of 0.1 cm/s
    !dn(1,patmo_idx_A)=dn(1,patmo_idx_A) - 0.1/(layer_thickness(in cm))*n(1,patmo_idx_A)
    if (n(1,patmo_idx_COS) > 9.5d-3/(1000*1d2)) then
      dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - (9.5d-3/(1000*1d2)) * n(1,patmo_idx_COS)
    end if
    if (n(1,patmo_idx_CS2) > 4.48d-2/(1000*1d2)) then
      dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) - (4.48d-2/(1000*1d2)) * n(1,patmo_idx_CS2)
    end if
    if (n(1,patmo_idx_SO2) > 1/(1000*1d2)) then
      dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - (1/(1000*1d2)) * n(1,patmo_idx_SO2)
    end if
    if (n(1,patmo_idx_H2S) > 1.7d-1/(1000*1d2)) then
      dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) - (1.7d-1/(1000*1d2)) * n(1,patmo_idx_H2S)
    end if
    if (n(1,patmo_idx_CH3SCH3) > 1.48d-1/(1000*1d2)) then
      dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) - (1.48d-1/(1000*1d2)) * n(1,patmo_idx_CH3SCH3)
    end if

    ! Emission
    dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 8.1001d2
    dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 5.9886d2
    dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 84.7910d2
    dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 615.84d2
    dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 39.50274d3

    ! Volcanic emission
    call patmo_volc_addSources(tt,n(:,:),dn(:,:))

    ! Wet Deposition
    do j=12, 2, -1
      do i = 1, chemSpeciesNumber
        dn(j,     i) = dn(j,     i) - wetdep(j, i) * n(j, i)
        dn(j - 1, i) = dn(j - 1, i) + wetdep(j, i) * n(j, i)
      end do
    end do
    do i = 1, chemSpeciesNumber
      dn(1, i) = dn(1, i) - wetdep(1, i) * n(1, i)
    end do
    !aerosol formation
    do i=13,34
      if (va(i) <= n(i, patmo_idx_H2SO4) .AND. pa(i) >= n(i, patmo_idx_H2SO4)) then
        dn(i, patmo_idx_H2SO4) = dn(i, patmo_idx_H2SO4) - (n(i, patmo_idx_H2SO4) - va(i))
        dn(i, patmo_idx_SO4)   = dn(i, patmo_idx_SO4)   + (n(i, patmo_idx_H2SO4) - va(i))
      end if
    end do

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
