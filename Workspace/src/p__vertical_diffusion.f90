module p__vertical_diffusion
  
  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private

  ! Module-level variables ------------------------------------------------------------------------
  real(dp), private :: z0, zu, zl, zp, zm
  real(dp), private :: m0, mu, ml, mp, mm, qi
  real(dp), private :: dzp, dz0, dzm, mi, niu, ni0, nil, nip, nim, neu, ne0, nel, nep, nem
  real(dp), private :: Diu, Di0, Dil, Dip, Dim, Kiu, Ki0, Kil, Kip, Kim
  real(dp), private :: Tiu, Ti0, Til, Teu, Te0, Tel, Tep, Tem, Tip, Tim
  real(dp), private :: Peu, Pe0, Pel, Pep, Pem
  real(dp), private :: dni_dzp, dni_dzm, dne_dzp, dne_dzm
  real(dp), private :: dTi_dzp, dTi_dzm, dTe_dzp, dTe_dzm, dPe_dzp, dPe_dzm
  real(dp), private :: Hiu, Hi0, Hil, Hip, Him, Hau, Ha0, Hal, Hap, Ham
  real(dp), private :: gradPep, gradPem
  real(dp), private :: Thermp, Thermm
  real(dp), private :: zetap, zetam
  real(dp), parameter, private :: k_B = 1.38064852e-16_dp ! erg K^-1 = g cm^2 s^-2 K^-1
  real(dp), parameter, private :: G_grav = 6.67408e-8_dp  ! cm^3 s^-2 g^-1

  ! Public interface ------------------------------------------------------------------------------
  public :: p__vertical_diffusion__exe

contains


  subroutine p__vertical_diffusion__exe(species,                                      & ! in:  name of variable species
    &                                   nsp, nz,                                      & ! in:  number of species and grids in z direction
    &                                   alt,                                          & ! in:  altitude
    &                                   mass, mass_mean, charge, alpha,               & ! in:  mass and charge of species, and thermal diffusion coefficient
    &                                   D_binary, K_eddy,                             & ! in:  binary diffusion coefficient and eddy diffusion coefficient
    &                                   ni, ne,                                       & ! in:  number density of variable species and electron density
    &                                   Tn, Ti, Te,                                   & ! in:  temperature of neutrals, ions and electrons
    &                                   phi_lower_boundary, phi_upper_boundary,       & ! in:  lower and upper boundary of fluxes
    &                                   v_lower_boundary, v_upper_boundary,           & ! in:  lower and upper boundary of velocities
    &                                   M_planet, R_planet,                           & ! in:  mass and radius of planet
    &                                   phi, dphi_dz,                                 & ! out: vertical flux, vertical gradient of vertical flux
    &                                   d_dniu_dphi_dz, d_dni0_dphi_dz, d_dnil_dphi_dz) ! out: Jacobian term of vertical flux
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - species(nsp)
    !   * Strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of the species string array.
    !
    ! - nsp
    !   * The number of chemical species to be treated as variables. 
    !
    ! - nz
    !   * The number of grids in vertical direction. 
    !
    ! - alt(nz)
    !   * Altitude in [cm] in the simulation vertical grids.
    !
    ! - mass(nsp), mass_mean(nz), charge(nsp), alpha(nsp)
    !   * Mass of species in [g], mean mass of species in [g] in vertical grids, 
    !     charge of species devided by elementary charge (no unit), and thermal diffusion coefficient.
    !
    ! - D_binary(nz,nsp), K_eddy(nz)
    !   * Binary diffusion coefficient in [cm^2 s^-1] and eddy diffusion coefficient in [cm^2 s^-1] 
    !     in the simulation vertical grids.
    !
    ! - ni(nz,nsp), ne(nz)
    !   * Number density of variable species in [cm^-3] and electron density in [cm^-3] 
    !     in the simulation vertical grids. 
    !
    ! - Tn(nz), Ti(nz,nsp), Te(nz)
    !   * Tn, Ti, and Te are the neutral, ion and electron temperatures in [K] in the simulation vertical grids. 
    !   * Ion temperature Ti has the number of elements equal to the number of variable species and simulation grids,  
    !     and should be given for each variable species ALSO FOR NEUTRALS AND ELECTRONS IF THEY ARE VARIABLES  
    !     because PROTEUS does not give a dedicated arrays only for ions to avoid the complexity of the code interface.
    !   * You can simply define "0.0 [K]" for neutral and electrons in Ti array. 
    !
    ! - phi_lower_boundary(nsp), phi_upper_boundary(nsp)
    !   * Lower and upper boundary of fluxes in [cm^-2 s^-1] for variable species.
    !
    ! - v_lower_boundary(nsp), v_upper_boundary(nsp)
    !   * Lower and upper boundary of velocities in [cm s^-1] for variable species.
    !
    ! - M_planet, R_planet
    !   * Mass of planet in [g] and radius of planet in [cm].
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - phi(nz+1,nsp)
    !   * Vertical flux of variable species in [cm^-2 s^-1] in the simulation vertical grids.
    !
    ! - dphi_dz(nz,nsp)
    !   * Vertical gradient of vertical flux in [cm^-3 s^-1] in the simulation vertical grids.
    !
    ! - d_dniu_dphi_dz(nz,nsp), d_dni0_dphi_dz(nz,nsp), d_dnil_dphi_dz(nz,nsp)
    !   * Jacobian term of vertical flux in [s^-1] in the simulation vertical grids.
    !     
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer,          intent(in)  :: nsp, nz 
    character(len=*), intent(in)  :: species(nsp)
    real(dp),         intent(in)  :: alt(nz)
    real(dp),         intent(in)  :: mass(nsp), mass_mean(nz), charge(nsp), alpha(nsp)
    real(dp),         intent(in)  :: Tn(nz), Ti(nz,nsp), Te(nz)
    real(dp),         intent(in)  :: ni(nz,nsp), ne(nz)
    real(dp),         intent(in)  :: D_binary(nz,nsp), K_eddy(nz)
    real(dp),         intent(in)  :: phi_lower_boundary(nsp), phi_upper_boundary(nsp)
    real(dp),         intent(in)  :: v_lower_boundary(nsp), v_upper_boundary(nsp)
    real(dp),         intent(in)  :: M_planet, R_planet
    real(dp),         intent(out) :: phi(nz+1,nsp)
    real(dp),         intent(out) :: dphi_dz(nz,nsp)
    real(dp),         intent(out) :: d_dniu_dphi_dz(nz,nsp), d_dni0_dphi_dz(nz,nsp), d_dnil_dphi_dz(nz,nsp)
    integer isp, iz
    real(dp) dalt(nz)

    dphi_dz        = 0.0_dp ! dphi/dz
    d_dniu_dphi_dz = 0.0_dp ! d/dniu(dphi/dz)
    d_dni0_dphi_dz = 0.0_dp ! d/dni0(dphi/dz)
    d_dnil_dphi_dz = 0.0_dp ! d/dnil(dphi/dz)

    ! altitude difference
    do iz = 2, nz
      dalt(iz) = alt(iz) - alt(iz-1)
    end do
    dalt(1) = dalt(2)

    !
    ! Space derivative of vertical flux
    !
    do isp = 1, nsp

      mi  = mass(isp)
      qi  = charge(isp)

      if (species(isp) /= 'e-' .and. species(isp) /= 'M' .and. species(isp) /= 'products' .and. species(isp) /= 'hv') then

        do iz = 2, nz-1

          ! grid definition

          ! iz+1    : u 
          ! iz+1/2  : p 
          ! iz      : 0 
          ! iz-1/2  : m  
          ! iz-1    : l

          zu  = alt(iz+1)
          z0  = alt(iz)
          zl  = alt(iz-1)
          zp  = ( zu + z0 ) / 2.0_dp
          zm  = ( z0 + zl ) / 2.0_dp

          dzp = dalt(iz)
          dzm = dalt(iz-1)
          dz0 = ( dzp + dzm ) / 2.0_dp

          mu  = mass_mean(iz+1)
          m0  = mass_mean(iz)
          ml  = mass_mean(iz-1)
          mp  = ( mu + m0 ) / 2.0_dp
          mm  = ( m0 + ml ) / 2.0_dp

          niu = ni(iz+1,isp)
          ni0 = ni(iz,isp)
          nil = ni(iz-1,isp)
          nip = ( niu + ni0 ) / 2.0_dp
          nim = ( ni0 + nil ) / 2.0_dp
          neu = ne(iz+1)
          ne0 = ne(iz)
          nel = ne(iz-1)
          nep = ( neu + ne0 ) / 2.0_dp
          nem = ( ne0 + nel ) / 2.0_dp

          Diu = D_binary(iz+1,isp)
          Di0 = D_binary(iz,isp)
          Dil = D_binary(iz-1,isp)
          Dip = ( Diu + Di0 ) / 2.0_dp
          Dim = ( Di0 + Dil ) / 2.0_dp
          Kiu = K_eddy(iz+1)
          Ki0 = K_eddy(iz)
          Kil = K_eddy(iz-1)
          Kip = ( Kiu + Ki0 ) / 2.0_dp
          Kim = ( Ki0 + Kil ) / 2.0_dp

          if ( charge(isp) /= 0.0_dp ) then
            Tiu = Ti(iz+1,isp)
            Ti0 = Ti(iz,isp)
            Til = Ti(iz-1,isp)
          else if ( charge(isp) == 0.0_dp ) then
            Tiu = Tn(iz+1)
            Ti0 = Tn(iz)
            Til = Tn(iz-1)
          end if
          Tip = ( Tiu + Ti0 ) / 2.0_dp
          Tim = ( Ti0 + Til ) / 2.0_dp

          Teu = Te(iz+1)
          Te0 = Te(iz)
          Tel = Te(iz-1)
          Tep = ( Teu + Te0 ) / 2.0_dp
          Tem = ( Te0 + Tel ) / 2.0_dp

          Peu = neu * k_B * Teu
          Pe0 = ne0 * k_B * Te0
          Pel = nel * k_B * Tel
          Pep = ( Peu + Pe0 ) / 2.0_dp
          Pem = ( Pe0 + Pel ) / 2.0_dp

          dni_dzp = ( niu - ni0 ) / dzp
          dni_dzm = ( ni0 - nil ) / dzm
          dne_dzp = ( neu - ne0 ) / dzp
          dne_dzm = ( ne0 - nel ) / dzm
          dTi_dzp = ( Tiu - Ti0 ) / dzp
          dTi_dzm = ( Ti0 - Til ) / dzm
          dTe_dzp = ( Teu - Te0 ) / dzp
          dTe_dzm = ( Te0 - Tel ) / dzm
          dPe_dzp = ( Peu - Pe0 ) / dzp
          dPe_dzm = ( Pe0 - Pel ) / dzm

          if ( neu > 0.0_dp .and. ne0 > 0.0_dp .and. nel > 0.0_dp &
              & .and. charge(isp) /= 0.0_dp ) then
            gradPep = (Tep/Tip) / Pep * dPe_dzp
            gradPem = (Tem/Tim) / Pem * dPe_dzm
          else
            gradPep = 0.0_dp
            gradPem = 0.0_dp
          end if
          
          Thermp = 1.0_dp / Tip * dTi_dzp
          Thermm = 1.0_dp / Tim * dTi_dzm
          Hiu = k_B * Tiu / mi / (G_grav*M_planet/(R_planet+zu)**2.0_dp)
          Hi0 = k_B * Ti0 / mi / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
          Hil = k_B * Til / mi / (G_grav*M_planet/(R_planet+zl)**2.0_dp)
          Hau = k_B * Tiu / mu / (G_grav*M_planet/(R_planet+zu)**2.0_dp)
          Ha0 = k_B * Ti0 / m0 / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
          Hal = k_B * Til / ml / (G_grav*M_planet/(R_planet+zl)**2.0_dp)
          Hip = ( Hiu + Hi0 ) / 2.0_dp
          Him = ( Hi0 + Hil ) / 2.0_dp
          Hap = ( Hau + Ha0 ) / 2.0_dp
          Ham = ( Ha0 + Hal ) / 2.0_dp

          if ( charge(isp) == 0.0_dp ) then
            zetap = Dip * ( 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
            zetam = Dim * ( 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
          else if ( charge(isp) > 0.0_dp ) then 
            zetap = Dip * (qi*gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
            zetam = Dim * (qi*gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
          else if ( charge(isp) < 0.0_dp ) then 
            zetap = Dip * (qi*gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
            zetam = Dim * (qi*gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
          end if

          dphi_dz(iz,isp) = niu * ( - (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0 ) & ! iz+1 -> iz
            &             + ni0 * (   (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0   & ! iz   -> iz+1
            &                       + (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0 ) & ! iz   -> iz-1
            &             + nil * ( - (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0 )   ! iz-1 -> iz
          d_dniu_dphi_dz(iz,isp) =  - (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0
          d_dni0_dphi_dz(iz,isp) =    (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0   &
            &                       + (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0
          d_dnil_dphi_dz(iz,isp) =  - (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0

          phi(iz,isp)   = - (Dim+Kim)*dni_dzm - zetam*nim
          phi(iz+1,isp) = - (Dip+Kip)*dni_dzp - zetap*nip

        end do

        ! Lower Boundary Condition
        zu  = alt(2)
        z0  = alt(1)
        zp  = ( zu + z0 ) / 2.0_dp

        dzp = dalt(1)
        dz0 = dzp

        mu  = mass_mean(2)
        m0  = mass_mean(1)
        mp  = ( mu + m0 ) / 2.0_dp

        niu = ni(2,isp)
        ni0 = ni(1,isp)
        nip = ( niu + ni0 ) / 2.0_dp
        neu = ne(2)
        ne0 = ne(1)
        nep = ( neu + ne0 ) / 2.0_dp

        Diu = D_binary(2,isp)
        Di0 = D_binary(1,isp)
        Dip = ( Diu + Di0 ) / 2.0_dp
        Kiu = K_eddy(2)
        Ki0 = K_eddy(1)
        Kip = ( Kiu + Ki0 ) / 2.0_dp

        if ( charge(isp) /= 0.0_dp ) then
          Tiu = Ti(2,isp)
          Ti0 = Ti(1,isp)
        else if ( charge(isp) == 0.0_dp ) then
          Tiu = Tn(2)
          Ti0 = Tn(1)
        end if
        Tip = ( Tiu + Ti0 ) / 2.0_dp

        Teu = Te(2)
        Te0 = Te(1)
        Tep = ( Teu + Te0 ) / 2.0_dp

        Peu = neu * k_B * Teu
        Pe0 = ne0 * k_B * Te0
        Pep = nep * k_B * Tep

        dni_dzp = ( niu - ni0 ) / dzp
        dne_dzp = ( neu - ne0 ) / dzp
        dTi_dzp = ( Tiu - Ti0 ) / dzp
        dTe_dzp = ( Teu - Te0 ) / dzp
        dPe_dzp = ( Peu - Pe0 ) / dzp

        if ( neu > 0.0_dp .and. ne0 > 0.0_dp &
            & .and. charge(isp) /= 0.0_dp ) then
          gradPep = (Tep/Tip) / Pep * dPe_dzp
        else
          gradPep = 0.0_dp
        end if
        Thermp = 1.0_dp / Tip * dTi_dzp
        Hiu = k_B * Tiu / mi / (G_grav*M_planet/(R_planet+zu)**2.0_dp)
        Hi0 = k_B * Ti0 / mi / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
        Hau = k_B * Tiu / mu / (G_grav*M_planet/(R_planet+zu)**2.0_dp)
        Ha0 = k_B * Ti0 / m0 / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
        Hip = ( Hiu + Hi0 ) / 2.0_dp
        Hap = ( Hau + Ha0 ) / 2.0_dp

        if ( charge(isp) == 0.0_dp ) then
          zetap = Dip * ( 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
        else if ( charge(isp) > 0.0_dp ) then 
          zetap = Dip * (qi*gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
        else if ( charge(isp) < 0.0_dp ) then 
          zetap = Dip * (qi*gradPep + 1.0_dp/Hip + (1.0_dp+alpha(isp))*Thermp) + Kip * (1.0_dp/Hap + Thermp)
        end if

        dphi_dz(1,isp) = niu * ( - (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0 ) & ! iz+1 -> iz
          &            + ni0 * (   (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0 ) & ! iz   -> iz+1
          &            - phi_lower_boundary(isp)/dz0 &                      ! Lower boundary: flux 
          &            - ni0 * v_lower_boundary(isp)/dz0                    ! Lower boundary: velocity
        d_dniu_dphi_dz(1,isp) =  - (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0
        d_dni0_dphi_dz(1,isp) =    (Dip+Kip)/dz0/dzp - zetap/2.0_dp/dz0 &
          &                      -  v_lower_boundary(isp)/dz0                ! Lower boundary: velocity
        d_dnil_dphi_dz(1,isp) = 0.0_dp

        phi(1,isp) = phi_lower_boundary(isp) + ni0 * v_lower_boundary(isp)

        ! Upper Boundary Condition
        z0  = alt(nz)
        zl  = alt(nz-1)
        zm  = ( z0 + zl ) / 2.0_dp

        dzm = dalt(nz)
        dz0 = dzm

        m0  = mass_mean(nz)
        ml  = mass_mean(nz-1)
        mm  = ( m0 + ml ) / 2.0_dp

        mi  = mass(isp)

        ni0 = ni(nz,isp)
        nil = ni(nz-1,isp)
        nim = ( ni0 + nil ) / 2.0_dp
        ne0 = ne(nz)
        nel = ne(nz-1)
        nem = ( ne0 + nel ) / 2.0_dp

        Di0 = D_binary(nz,isp)
        Dil = D_binary(nz-1,isp)
        Dim = ( Di0 + Dil ) / 2.0_dp
        Ki0 = K_eddy(nz)
        Kil = K_eddy(nz-1)
        Kim = ( Ki0 + Kil ) / 2.0_dp

        if ( charge(isp) /= 0.0_dp ) then
          Ti0 = Ti(nz,isp)
          Til = Ti(nz-1,isp)
        else if ( charge(isp) == 0.0_dp ) then
          Ti0 = Tn(nz)
          Til = Tn(nz-1)
        end if
        Tim = ( Ti0 + Til ) / 2.0_dp

        Te0 = Te(nz)
        Tel = Te(nz-1)
        Tem = ( Te0 + Tel ) / 2.0_dp

        Pe0 = ne0 * k_B * Te0
        Pel = nel * k_B * Tel
        Pem = nem * k_B * Tem

        dni_dzm = ( ni0 - nil ) / dzm
        dne_dzm = ( ne0 - nel ) / dzm
        dTi_dzm = ( Ti0 - Til ) / dzm
        dTe_dzm = ( Te0 - Tel ) / dzm
        dPe_dzm = ( Pe0 - Pel ) / dzm

        if ( ne0 > 0.0_dp .and. nel > 0.0_dp &
            & .and. charge(isp) /= 0.0_dp ) then
          gradPem = (Tem/Tim) / Pem * dPe_dzm
        else
          gradPem = 0.0_dp
        end if
        Thermm = 1.0_dp / Tim * dTi_dzm
        Hi0 = k_B * Ti0 / mi / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
        Hil = k_B * Til / mi / (G_grav*M_planet/(R_planet+zl)**2.0_dp)
        Ha0 = k_B * Ti0 / m0 / (G_grav*M_planet/(R_planet+z0)**2.0_dp)
        Hal = k_B * Til / ml / (G_grav*M_planet/(R_planet+zl)**2.0_dp)
        Him = ( Hi0 + Hil ) / 2.0_dp
        Ham = ( Ha0 + Hal ) / 2.0_dp

        if ( charge(isp) == 0.0_dp ) then
          zetam = Dim * ( 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
        else if ( charge(isp) > 0.0_dp ) then 
          zetam = Dim * (qi*gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
        else if ( charge(isp) < 0.0_dp ) then 
          zetam = Dim * (qi*gradPem + 1.0_dp/Him + (1.0_dp+alpha(isp))*Thermm) + Kim * (1.0_dp/Ham + Thermm)
        end if

        dphi_dz(nz,isp) = ni0 * (   (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0 ) & ! iz   -> iz-1
          &             + nil * ( - (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0 ) & ! iz-1 -> iz
          &             + phi_upper_boundary(isp)/dz0 &                      ! Upper boundary: flux 
          &             + ni0 * v_upper_boundary(isp)/dz0                    ! Upper boundary: velocity
        d_dniu_dphi_dz(nz,isp) = 0.0_dp
        d_dni0_dphi_dz(nz,isp) =    (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0 &
          &                       +  v_upper_boundary(isp)/dz0                ! Upper boundary: velocity
        d_dnil_dphi_dz(nz,isp) =  - (Dim+Kim)/dz0/dzm + zetam/2.0_dp/dz0

        phi(nz+1,isp) = phi_upper_boundary(isp) + ni0 * v_upper_boundary(isp)

      end if
      
    end do


  end subroutine p__vertical_diffusion__exe

end module p__vertical_diffusion