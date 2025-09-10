module p__photochem_transport
  
  use v__tdec,                only : spl_, var_, grd_, cst_, set_
  use p__search,              only : p__search_reactant, p__search_product, sp_index
  use p__eddy_diffusion,      only : p__eddy_diffusion__exe
  use p__molecular_diffusion, only : p__molecular_diffusion__exe
  use p__vertical_diffusion,  only : p__vertical_diffusion__exe

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private

  ! Module-level variables ----------------------------------------------------------------------------------
  character(len=256), allocatable, private :: species(:)
  real(dp),           allocatable, private :: alt(:)
  real(dp),           allocatable, private :: mass(:), mass_mean(:), charge(:), alpha(:)
  real(dp),           allocatable, private :: Tn(:), Ti(:,:), Te(:)
  real(dp),           allocatable, private :: ni(:,:), ne(:)
  real(dp),           allocatable, private :: D_binary(:,:), K_eddy(:)
  real(dp),           allocatable, private :: phi_lower_boundary(:), phi_upper_boundary(:)
  real(dp),           allocatable, private :: v_lower_boundary(:), v_upper_boundary(:)
  real(dp),           allocatable, private :: phi(:,:), dphi_dz(:,:)
  real(dp),           allocatable, private :: d_dniu_dphi_dz(:,:), d_dni0_dphi_dz(:,:), d_dnil_dphi_dz(:,:)

  ! Public interface ----------------------------------------------------------------------------------------
  public :: p__photochem_transport__ini, p__photochem_transport__exe

contains

  subroutine p__photochem_transport__ini(spl, grd)
    implicit none
    type(spl_), intent(in) :: spl
    type(grd_), intent(in) :: grd

    ! Allocate arrays for diffusion transport
    allocate(species(spl%nsp_i), mass(spl%nsp_i), charge(spl%nsp_i), alpha(spl%nsp_i))
    allocate(alt(grd%nz), mass_mean(grd%nz), Tn(grd%nz), Ti(grd%nz,spl%nsp_i), Te(grd%nz))
    allocate(ni(grd%nz,spl%nsp_i), ne(grd%nz))
    allocate(D_binary(grd%nz,spl%nsp_i), K_eddy(grd%nz))
    allocate(phi_lower_boundary(spl%nsp_i), phi_upper_boundary(spl%nsp_i))
    allocate(v_lower_boundary(spl%nsp_i), v_upper_boundary(spl%nsp_i))
    allocate(phi(grd%nz+1,spl%nsp_i), dphi_dz(grd%nz,spl%nsp_i))
    allocate(d_dniu_dphi_dz(grd%nz,spl%nsp_i), d_dni0_dphi_dz(grd%nz,spl%nsp_i), d_dnil_dphi_dz(grd%nz,spl%nsp_i))

  end subroutine p__photochem_transport__ini


  subroutine p__photochem_transport__exe(spl, cst, grd, set, & ! in
    &                                    var                 ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var

    integer isp, jsp, iz
    real(dp) Texo, lambda, vth
    real(dp) M_planet, R_planet

    ! Vertical diffusion ---------------------------------------------------------------------------

      !--------------------------------------
      ! Eddy diffusion coefficient
      !--------------------------------------
      var%K_eddy = 0.0_dp
      call p__eddy_diffusion__exe(spl, cst, grd, & ! in
        &                         var            ) ! inout

      !--------------------------------------
      ! Molecular diffusion coefficient
      !--------------------------------------
      var%D_mol = 0.0_dp
      call p__molecular_diffusion__exe(spl, cst, grd, & ! in
        &                              var            ) ! inout

      !--------------------------------------
      ! Boundary conditions
      !--------------------------------------
      v_lower_boundary = 0.0_dp
      v_upper_boundary = 0.0_dp
      phi_lower_boundary = 0.0_dp
      phi_upper_boundary = 0.0_dp

      do isp = 1, spl%nsp
        jsp = spl%all_to_var(isp)
        if (nint(var%Upper_v(isp,1)) == 10) then 
          Texo = var%Tn(grd%nz)
          lambda = cst%BigG * cst%Mplanet * var%m(isp) / ( cst%k_B * Texo * (cst%R+grd%alt(grd%nz)) )
          vth = dsqrt(2.0_dp*cst%k_B*Texo/var%m(isp))
          v_upper_boundary(jsp) = dexp(-lambda)*vth*(lambda+1.0_dp)/(2.0_dp*cst%pi**0.5_dp) * 1.0e2_dp
        end if
        if (nint(var%Upper_v(isp,1)) == 20) then
          v_upper_boundary(jsp) = var%D_mol(grd%nz,isp)*&
          &                     (1.0_dp/(cst%k_B * var%Tn(grd%nz) / var%m_mean(grd%nz) &
          &                           / (cst%BigG*cst%Mplanet/(cst%R+grd%alt(grd%nz))**2.0_dp))&
          &                    - 1.0_dp/(cst%k_B * var%Tn(grd%nz) / var%m(isp) &
          &                           / (cst%BigG*cst%Mplanet/(cst%R+grd%alt(grd%nz))**2.0_dp)))
        end if
        if (nint(var%Upper_f(isp,1)) == 1) then 
          phi_upper_boundary(jsp) = var%Upper_f(isp,2) * 1.0e-4_dp
        end if
        if (nint(var%Upper_v(isp,1)) == 1) then 
          v_upper_boundary(jsp) = var%Upper_v(isp,2) * 1.0e2_dp
        end if
        if (nint(var%Lower_f(isp,1)) == 1) then 
          phi_lower_boundary(jsp) = var%Lower_f(isp,2) * 1.0e-4_dp
        end if
        if (nint(var%Lower_v(isp,1)) == 1) then 
          v_lower_boundary(jsp) = var%Lower_v(isp,2) * 1.0e2_dp
        end if
      end do

      ! Thermal diffusion coefficients
      alpha = 0.0_dp
      isp = sp_index(spl, 'H')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        isp = spl%all_to_var(isp)
        if (isp >= 1) alpha(isp) = -0.25_dp
      end if
      isp = sp_index(spl, 'H2')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        isp = spl%all_to_var(isp)
        if (isp >= 1) alpha(isp) = -0.25_dp
      end if

      ! Temperature
      do iz = 1, grd%nz
        ne(iz) = 1.0e-20_dp
        do isp = 1, spl%nsp
          if (spl%species(isp) == 'e-') then
            ne(iz) = var%ni(iz,isp) * 1.0e-6_dp
          end if
        end do
        Tn(iz) = var%Tn(iz)
        Ti(iz,:) = var%Ti(iz)
        Te(iz) = var%Te(iz)
      end do

      do isp = 1, spl%nsp_i
        jsp = spl%var_to_all(isp)
        species(isp) = spl%species(jsp)
        mass(isp) = var%m(jsp) * 1.0e3_dp
        charge(isp) = var%q(jsp)/cst%q_e
        ni(1:grd%nz,isp) = var%ni(1:grd%nz,jsp) * 1.0e-6_dp
        D_binary(1:grd%nz,isp) = var%D_mol(1:grd%nz,jsp) * 1.0e4_dp
      end do 
      alt(1:grd%nz) = grd%alt(1:grd%nz) * 1.0e2_dp
      mass_mean(1:grd%nz) = var%m_mean(1:grd%nz) * 1.0e3_dp
      K_eddy(1:grd%nz) = var%K_eddy(1:grd%nz) * 1.0e4_dp

      M_planet = cst%Mplanet * 1.0e3_dp ! kg -> g
      R_planet = cst%R * 1.0e2_dp ! m -> cm

      call p__vertical_diffusion__exe(species,                                      & ! in:  name of variable species
        &                             spl%nsp_i, grd%nz,                            & ! in:  number of species and grids in z direction
        &                             alt,                                          & ! in:  altitude
        &                             mass, mass_mean, charge, alpha,               & ! in:  mass and charge of species, and thermal diffusion coefficient
        &                             D_binary, K_eddy,                             & ! in:  binary diffusion coefficient and eddy diffusion coefficient
        &                             ni, ne,                                       & ! in:  number density of variable
        &                             Tn, Ti, Te,                                   & ! in:  temperature of neutrals, ions and electrons
        &                             phi_lower_boundary, phi_upper_boundary,       & ! in:  lower and upper boundary of fluxes
        &                             v_lower_boundary, v_upper_boundary,           & ! in:  lower and upper boundary of velocities
        &                             M_planet, R_planet,                           & ! in:  mass and radius of planet
        &                             phi, dphi_dz,                                 & ! out:
        &                             d_dniu_dphi_dz, d_dni0_dphi_dz, d_dnil_dphi_dz) ! out: 

        var%dphi_dz = dphi_dz * 1.0e6_dp
        var%d_dniu_dphi_dz = d_dniu_dphi_dz
        var%d_dni0_dphi_dz = d_dni0_dphi_dz
        var%d_dnil_dphi_dz = d_dnil_dphi_dz

        var%Fluxup  = 0.0_dp
        var%Fluxdwn = 0.0_dp
        do iz = 1, grd%nz+1
          do isp = 1, spl%nsp
            jsp = spl%all_to_var(isp)
            if (jsp >= 1 .and. jsp <= spl%nsp_i) then 
              if (Phi(iz,jsp) >= 0.0_dp) then
                var%Fluxup(iz,jsp) = Phi(iz,jsp) * 1.0e4_dp
              else if (Phi(iz,jsp) < 0.0_dp) then
                var%Fluxdwn(iz,jsp) = - Phi(iz,jsp) * 1.0e4_dp
              end if
            end if

          end do
        end do

  end subroutine p__photochem_transport__exe


end module p__photochem_transport
