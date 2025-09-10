module p__PROTEUS
  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: get_number_of_reaction, get_number_of_all_species, get_number_of_var_species, get_number_of_fix_species, &
            get_all_species_name, get_var_species_name, get_fix_species_name, get_mass_of_species, &
            get_reaction_type_list, get_reactant_product_list, &
            get_number_of_wavelength_bins, get_wavelength_bin, get_solar_flux, &
            p__PROTEUS_source, p__PROTEUS_Jacobian 

contains


  integer function get_number_of_reaction()
    get_number_of_reaction = 64
    ! R1: CO2 + hv -> CO + O 
    ! R2: CO2 + hv -> CO + O(1D) 
    ! R3: H2O + hv -> H + OH 
    ! R4: H2O + hv -> H2 + O(1D) 
    ! R5: O3 + hv -> O2 + O 
    ! R6: O3 + hv -> O2 + O(1D) 
    ! R7: O2 + hv -> O + O 
    ! R8: O2 + hv -> O + O(1D) 
    ! R9: H2 + hv -> H + H 
    ! R10: OH + hv -> O + H 
    ! R11: OH + hv -> O(1D) + H 
    ! R12: HO2 + hv -> OH + O 
    ! R13: H2O2 + hv -> OH + OH 
    ! R14: H2O2 + hv -> HO2 + H 
    ! R15: H2CO + hv -> HCO + H 
    ! R16: H2CO + hv -> CO + H2 
    ! R17: O + O + M -> O2 + M 
    ! R18: O + O2 + CO2 -> O3 + CO2 
    ! R19: O + O3 -> O2 + O2 
    ! R20: O + CO + M -> CO2 + M 
    ! R21: O(1D) + O2 -> O + O2 
    ! R22: O(1D) + O3 -> O2 + O2 
    ! R23: O(1D) + O3 -> O + O + O2 
    ! R24: O(1D) + H2 -> H + OH 
    ! R25: O(1D) + CO2 -> O + CO2 
    ! R26: O(1D) + H2O -> OH + OH 
    ! R27: H2 + O -> OH + H 
    ! R28: OH + H2 -> H2O + H 
    ! R29: H + H + CO2 -> H2 + CO2 
    ! R30: H + OH + CO2 -> H2O + CO2 
    ! R31: H + HO2 -> OH + OH 
    ! R32: H + HO2 -> H2O + O(1D) 
    ! R33: H + HO2 -> H2 + O2 
    ! R34: H + H2O2 -> HO2 + H2 
    ! R35: H + H2O2 -> H2O + OH 
    ! R36: H + O2 + M -> HO2 + M 
    ! R37: H + O3 -> OH + O2 
    ! R38: O + OH -> O2 + H 
    ! R39: O + HO2 -> OH + O2 
    ! R40: O + H2O2 -> OH + HO2 
    ! R41: OH + OH -> H2O + O 
    ! R42: OH + OH + M -> H2O2 + M 
    ! R43: OH + O3 -> HO2 + O2 
    ! R44: OH + HO2 -> H2O + O2 
    ! R45: OH + H2O2 -> H2O + HO2 
    ! R46: HO2 + O3 -> OH + O2 + O2 
    ! R47: HO2 + HO2 -> H2O2 + O2 
    ! R48: HO2 + HO2 + M -> H2O2 + O2 + M 
    ! R49: CO + OH + M -> CO2 + H + M 
    ! R50: CO + OH + M -> HOCO + M 
    ! R51: HOCO + O2 -> HO2 + CO2 
    ! R52: CO2+ + H2 -> CO2 + H + H 
    ! R53: H + CO + M -> HCO + M 
    ! R54: H + HCO -> H2 + CO 
    ! R55: HCO + HCO -> H2CO + CO 
    ! R56: OH + HCO -> H2O + CO 
    ! R57: O + HCO -> H + CO2 
    ! R58: O + HCO -> OH + CO 
    ! R59: H2CO + H -> H2 + HCO 
    ! R60: HCO + O2 -> HO2 + CO 
    ! R61: H2CO + OH -> H2O + HCO 
    ! R62: H2CO + O -> OH + HCO 
    ! R63: HCO + HCO -> CO + CO + H2 
    ! R64: HO2 + HCO -> H2O2 + CO 
  end function


  integer function get_number_of_all_species()
    get_number_of_all_species = 17
  end function


  integer function get_number_of_var_species()
    get_number_of_var_species = 13
  end function


  integer function get_number_of_fix_species()
    get_number_of_fix_species = 4
  end function


  function get_all_species_name()
    implicit none
    character(len=256), dimension(17) :: get_all_species_name
    get_all_species_name(1) = 'CO2'
    get_all_species_name(2) = 'CO'
    get_all_species_name(3) = 'O'
    get_all_species_name(4) = 'O(1D)'
    get_all_species_name(5) = 'H2O'
    get_all_species_name(6) = 'H'
    get_all_species_name(7) = 'OH'
    get_all_species_name(8) = 'H2'
    get_all_species_name(9) = 'O3'
    get_all_species_name(10) = 'O2'
    get_all_species_name(11) = 'HO2'
    get_all_species_name(12) = 'H2O2'
    get_all_species_name(13) = 'H2CO'
    get_all_species_name(14) = 'HCO'
    get_all_species_name(15) = 'M'
    get_all_species_name(16) = 'HOCO'
    get_all_species_name(17) = 'CO2+'
  end function


  function get_var_species_name()
    implicit none
    character(len=256), dimension(13) :: get_var_species_name
    get_var_species_name(1) = 'CO'
    get_var_species_name(2) = 'O'
    get_var_species_name(3) = 'O(1D)'
    get_var_species_name(4) = 'H'
    get_var_species_name(5) = 'OH'
    get_var_species_name(6) = 'H2'
    get_var_species_name(7) = 'O3'
    get_var_species_name(8) = 'O2'
    get_var_species_name(9) = 'HO2'
    get_var_species_name(10) = 'H2O2'
    get_var_species_name(11) = 'H2CO'
    get_var_species_name(12) = 'HCO'
    get_var_species_name(13) = 'HOCO'
  end function


  function get_fix_species_name()
    implicit none
    character(len=256), dimension(4) :: get_fix_species_name
    get_fix_species_name(1) = 'CO2'
    get_fix_species_name(2) = 'H2O'
    get_fix_species_name(3) = 'M'
    get_fix_species_name(4) = 'CO2+'
  end function


  function get_mass_of_species()
    implicit none
    integer, parameter :: nsp = 17
    real(dp), dimension(1:nsp) :: get_mass_of_species
    real(dp), parameter :: m_u = 1.660538921e-27_dp ! kg
    get_mass_of_species = 0.0_dp
    get_mass_of_species(1) = 44.0_dp * m_u ! CO2
    get_mass_of_species(2) = 28.0_dp * m_u ! CO
    get_mass_of_species(3) = 16.0_dp * m_u ! O
    get_mass_of_species(4) = 16.0_dp * m_u ! O(1D)
    get_mass_of_species(5) = 18.0_dp * m_u ! H2O
    get_mass_of_species(6) = 1.0_dp * m_u ! H
    get_mass_of_species(7) = 17.0_dp * m_u ! OH
    get_mass_of_species(8) = 2.0_dp * m_u ! H2
    get_mass_of_species(9) = 48.0_dp * m_u ! O3
    get_mass_of_species(10) = 32.0_dp * m_u ! O2
    get_mass_of_species(11) = 33.0_dp * m_u ! HO2
    get_mass_of_species(12) = 34.0_dp * m_u ! H2O2
    get_mass_of_species(13) = 30.0_dp * m_u ! H2CO
    get_mass_of_species(14) = 29.0_dp * m_u ! HCO
    get_mass_of_species(15) = 1000000.0_dp * m_u ! M
    get_mass_of_species(16) = 45.0_dp * m_u ! HOCO
    get_mass_of_species(17) = 44.0_dp * m_u ! CO2+
  end function get_mass_of_species


  subroutine get_reactant_product_list(reactant_list, product_list)
    implicit none
    integer, parameter :: nch = 64
    integer, intent(out) :: reactant_list(1:nch,0:20)
    integer, intent(out) :: product_list(1:nch,0:20)

    ! reactant list ----------------------------
    reactant_list(1,0:1) = (/1,1/)
    reactant_list(2,0:1) = (/1,1/)
    reactant_list(3,0:1) = (/1,5/)
    reactant_list(4,0:1) = (/1,5/)
    reactant_list(5,0:1) = (/1,9/)
    reactant_list(6,0:1) = (/1,9/)
    reactant_list(7,0:1) = (/1,10/)
    reactant_list(8,0:1) = (/1,10/)
    reactant_list(9,0:1) = (/1,8/)
    reactant_list(10,0:1) = (/1,7/)
    reactant_list(11,0:1) = (/1,7/)
    reactant_list(12,0:1) = (/1,11/)
    reactant_list(13,0:1) = (/1,12/)
    reactant_list(14,0:1) = (/1,12/)
    reactant_list(15,0:1) = (/1,13/)
    reactant_list(16,0:1) = (/1,13/)
    reactant_list(17,0:3) = (/3,3,3,15/)
    reactant_list(18,0:3) = (/3,3,10,1/)
    reactant_list(19,0:2) = (/2,3,9/)
    reactant_list(20,0:3) = (/3,3,2,15/)
    reactant_list(21,0:2) = (/2,4,10/)
    reactant_list(22,0:2) = (/2,4,9/)
    reactant_list(23,0:2) = (/2,4,9/)
    reactant_list(24,0:2) = (/2,4,8/)
    reactant_list(25,0:2) = (/2,4,1/)
    reactant_list(26,0:2) = (/2,4,5/)
    reactant_list(27,0:2) = (/2,8,3/)
    reactant_list(28,0:2) = (/2,7,8/)
    reactant_list(29,0:3) = (/3,6,6,1/)
    reactant_list(30,0:3) = (/3,6,7,1/)
    reactant_list(31,0:2) = (/2,6,11/)
    reactant_list(32,0:2) = (/2,6,11/)
    reactant_list(33,0:2) = (/2,6,11/)
    reactant_list(34,0:2) = (/2,6,12/)
    reactant_list(35,0:2) = (/2,6,12/)
    reactant_list(36,0:3) = (/3,6,10,15/)
    reactant_list(37,0:2) = (/2,6,9/)
    reactant_list(38,0:2) = (/2,3,7/)
    reactant_list(39,0:2) = (/2,3,11/)
    reactant_list(40,0:2) = (/2,3,12/)
    reactant_list(41,0:2) = (/2,7,7/)
    reactant_list(42,0:3) = (/3,7,7,15/)
    reactant_list(43,0:2) = (/2,7,9/)
    reactant_list(44,0:2) = (/2,7,11/)
    reactant_list(45,0:2) = (/2,7,12/)
    reactant_list(46,0:2) = (/2,11,9/)
    reactant_list(47,0:2) = (/2,11,11/)
    reactant_list(48,0:3) = (/3,11,11,15/)
    reactant_list(49,0:3) = (/3,2,7,15/)
    reactant_list(50,0:3) = (/3,2,7,15/)
    reactant_list(51,0:2) = (/2,16,10/)
    reactant_list(52,0:2) = (/2,17,8/)
    reactant_list(53,0:3) = (/3,6,2,15/)
    reactant_list(54,0:2) = (/2,6,14/)
    reactant_list(55,0:2) = (/2,14,14/)
    reactant_list(56,0:2) = (/2,7,14/)
    reactant_list(57,0:2) = (/2,3,14/)
    reactant_list(58,0:2) = (/2,3,14/)
    reactant_list(59,0:2) = (/2,13,6/)
    reactant_list(60,0:2) = (/2,14,10/)
    reactant_list(61,0:2) = (/2,13,7/)
    reactant_list(62,0:2) = (/2,13,3/)
    reactant_list(63,0:2) = (/2,14,14/)
    reactant_list(64,0:2) = (/2,11,14/)

    ! product list ----------------------------
    product_list(1,0:2) = (/2,2,3/)
    product_list(2,0:2) = (/2,2,4/)
    product_list(3,0:2) = (/2,6,7/)
    product_list(4,0:2) = (/2,8,4/)
    product_list(5,0:2) = (/2,10,3/)
    product_list(6,0:2) = (/2,10,4/)
    product_list(7,0:2) = (/2,3,3/)
    product_list(8,0:2) = (/2,3,4/)
    product_list(9,0:2) = (/2,6,6/)
    product_list(10,0:2) = (/2,3,6/)
    product_list(11,0:2) = (/2,4,6/)
    product_list(12,0:2) = (/2,7,3/)
    product_list(13,0:2) = (/2,7,7/)
    product_list(14,0:2) = (/2,11,6/)
    product_list(15,0:2) = (/2,14,6/)
    product_list(16,0:2) = (/2,2,8/)
    product_list(17,0:2) = (/2,10,15/)
    product_list(18,0:2) = (/2,9,1/)
    product_list(19,0:2) = (/2,10,10/)
    product_list(20,0:2) = (/2,1,15/)
    product_list(21,0:2) = (/2,3,10/)
    product_list(22,0:2) = (/2,10,10/)
    product_list(23,0:3) = (/3,3,3,10/)
    product_list(24,0:2) = (/2,6,7/)
    product_list(25,0:2) = (/2,3,1/)
    product_list(26,0:2) = (/2,7,7/)
    product_list(27,0:2) = (/2,7,6/)
    product_list(28,0:2) = (/2,5,6/)
    product_list(29,0:2) = (/2,8,1/)
    product_list(30,0:2) = (/2,5,1/)
    product_list(31,0:2) = (/2,7,7/)
    product_list(32,0:2) = (/2,5,4/)
    product_list(33,0:2) = (/2,8,10/)
    product_list(34,0:2) = (/2,11,8/)
    product_list(35,0:2) = (/2,5,7/)
    product_list(36,0:2) = (/2,11,15/)
    product_list(37,0:2) = (/2,7,10/)
    product_list(38,0:2) = (/2,10,6/)
    product_list(39,0:2) = (/2,7,10/)
    product_list(40,0:2) = (/2,7,11/)
    product_list(41,0:2) = (/2,5,3/)
    product_list(42,0:2) = (/2,12,15/)
    product_list(43,0:2) = (/2,11,10/)
    product_list(44,0:2) = (/2,5,10/)
    product_list(45,0:2) = (/2,5,11/)
    product_list(46,0:3) = (/3,7,10,10/)
    product_list(47,0:2) = (/2,12,10/)
    product_list(48,0:3) = (/3,12,10,15/)
    product_list(49,0:3) = (/3,1,6,15/)
    product_list(50,0:2) = (/2,16,15/)
    product_list(51,0:2) = (/2,11,1/)
    product_list(52,0:3) = (/3,1,6,6/)
    product_list(53,0:2) = (/2,14,15/)
    product_list(54,0:2) = (/2,8,2/)
    product_list(55,0:2) = (/2,13,2/)
    product_list(56,0:2) = (/2,5,2/)
    product_list(57,0:2) = (/2,6,1/)
    product_list(58,0:2) = (/2,7,2/)
    product_list(59,0:2) = (/2,8,14/)
    product_list(60,0:2) = (/2,11,2/)
    product_list(61,0:2) = (/2,5,14/)
    product_list(62,0:2) = (/2,7,14/)
    product_list(63,0:3) = (/3,2,2,8/)
    product_list(64,0:2) = (/2,12,2/)

  end subroutine get_reactant_product_list


  subroutine get_reaction_type_list(reaction_type)
    implicit none
    integer, parameter :: nch = 64
    integer, intent(out) :: reaction_type(1:nch)

    ! reaction type
    reaction_type = 0
    reaction_type(1) = 2
    reaction_type(2) = 2
    reaction_type(3) = 2
    reaction_type(4) = 2
    reaction_type(5) = 2
    reaction_type(6) = 2
    reaction_type(7) = 2
    reaction_type(8) = 2
    reaction_type(9) = 2
    reaction_type(10) = 2
    reaction_type(11) = 2
    reaction_type(12) = 2
    reaction_type(13) = 2
    reaction_type(14) = 2
    reaction_type(15) = 2
    reaction_type(16) = 2

  end subroutine get_reaction_type_list


  integer function get_number_of_wavelength_bins()
    get_number_of_wavelength_bins = 1000
  end function


  subroutine get_wavelength_bin(lambda, dlambda)
    implicit none
    real(dp), intent(inout) :: lambda(1:), dlambda(1:)
    integer i
    do i = 1, 1000
      dlambda(i) = 1.0_dp
      lambda(i) = 0.0_dp + 0.5_dp * dlambda(i) + dble(i-1) * 1.0_dp
    end do
  end subroutine get_wavelength_bin


  subroutine get_solar_flux(lambda, dlambda, solar_flux)
    implicit none
    real(dp), intent(in) :: lambda(1:), dlambda(1:)
    real(dp), intent(out) :: solar_flux(1:)
    integer, parameter :: nwl = 1000
    character(len=256) fname, unit1, unit2

    fname = './UV/ref_solar_irradiance_whi-2008_ver2_1.dat'
    unit1 = 'nm'
    unit2 = 'W/m^2/nm'

    call get_solar_flux_data(lambda, dlambda, nwl, fname, unit1, unit2, & ! in
      &                      solar_flux                                 ) ! out

  end subroutine get_solar_flux


  subroutine p__PROTEUS_source(var_species_list, fix_species_list, & ! in:  name of variable species and fixed species
    &                          nsp_var_out, nsp_fix_out,           & ! in:  number of variable species and fixed species
    &                          nx, ny, nz,                         & ! in:  number of grids in x, y, z direction
    &                          T_n, T_i, T_e,                      & ! in:  temperature of neutrals, ions and electrons
    &                          v_var, v_fix,                       & ! in:  three dimensional velocity vector of variable and fixed species
    &                          n_var, n_fix,                       & ! in:  number density of variable and fixed species
    &                          J_rate_in,                          & ! in:  Reaction rate input, e.g., photolysis rate, rate calculated by other model
    &                          prod, loss, k_coef                  ) ! out: production and loss rates for variable species, and reaction rate coefficients
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    !   * Both arrays are strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables and fixed, respectively. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of species_list string arrays.
    !
    ! - nsp_var_out, nsp_fix_out
    !   * The number of chemical species to be treated as variables and fixed, respectively. 
    !
    ! - nx, ny, nz
    !   * The number of grids in x, y, z direction. 
    !   * Please input "1, 1, nz" if the simulation is vertical 1D, and "1, 1, 1" if it is 0D.
    !
    ! - T_n(nx,ny,nz), T_i(nx,ny,nz,nsp_fix_out+nsp_var_out), T_e(nx,ny,nz)
    !   * T_n, T_i, and T_e are the neutral, ion and electron temperatures in [K] in the simulation grids. 
    !   * Ion temperature Ti has the number of elements equal to the number of variable species and simulation grids,  
    !     and should be given for each variable species ALSO FOR NEUTRALS AND ELECTRONS IF THEY ARE VARIABLES  
    !     because PROTEUS does not give a dedicated arrays only for ions to avoid the complexity of the code interface.
    !   * You can simply define "0.0 [K]" for neutral and electrons in T_i array. 
    !
    ! - v_var(3,nx,ny,nz,nsp_var_out), v_fix(3,nx,ny,nz,nsp_fix_out)
    !   * Three dimensional velocity vector of variable and fixed species in [cm s^-1] in simulation grids. 
    !     v_*(1,ix,iy,iz,i) = x component of the velocity of ith species at (ix,iy,iz) grid. 
    !     v_*(2,ix,iy,iz,i) = y component of the velocity of ith species at (ix,iy,iz) grid. 
    !     v_*(3,ix,iy,iz,i) = z component of the velocity of ith species at (ix,iy,iz) grid. 
    !
    ! - n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    !   * Number density of variable and fixed species in [cm^-3] in simulation grids. 
    !
    ! - J_rate_in(nx,ny,nz,nch)
    !   * Users can input reaction rate coefficient calculated by other modules, such as photolysis rate. Units are arbitrary.
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - prod(nx,ny,nz,nsp_var_out), loss(nx,ny,nz,nsp_var_out)
    !   * Production and loss rates for variable species in [cm^-3 s^-1]. 
    !
    ! - k_coef(nx,ny,nz,nch)
    !   * Reaction rate coefficients for each reaction in the simulation grids in the units as follows:
    !     [s^-1] for reactions with only one reactant 
    !     [cm^3 s^-1] for two-body reactions
    !     [cm^6 s^-1] for three-body reactions
    !
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer, parameter :: nsp = 17
    integer, parameter :: nsp_fix_in = 4
    integer, parameter :: nsp_var_in = 13
    integer, parameter :: nch = 64
    integer,          intent(in)  :: nsp_var_out, nsp_fix_out, nx, ny, nz 
    character(len=*), intent(in)  :: var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    real(dp),         intent(in)  :: T_n(nx,ny,nz), T_i(nx,ny,nz,nsp_fix_out+nsp_var_out), T_e(nx,ny,nz)
    real(dp),         intent(in)  :: v_var(3,nx,ny,nz,nsp_var_out), v_fix(3,nx,ny,nz,nsp_fix_out)
    real(dp),         intent(in)  :: n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    real(dp),         intent(in)  :: J_rate_in(nx,ny,nz,nch)
    real(dp),         intent(out) :: prod(nx,ny,nz,nsp_var_out), loss(nx,ny,nz,nsp_var_out), k_coef(nx,ny,nz,nch)
    integer isp, jsp, ksp 
    integer ix, iy, iz 
    integer  all_in_var_out(nsp), var_out_all_in(nsp_var_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_fix_out(nsp), fix_out_all_in(nsp_fix_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_var_in(nsp) ! 
    real(dp) P(nsp_var_in), L(nsp_var_in), k(nch), k0, kinf, k2, k3, n(nsp), n_tot, Tn, Te, Ti(nsp), T_tmp, v(3,nsp)
    character(len = 256) species(nsp)

    ! species list ----------------------------
    ! 
    ! CO2, CO, O, O(1D), H2O, H, OH, H2, O3, O2, HO2, H2O2, H2CO, HCO, M, HOCO, CO2+
    species(1) = "CO2"
    species(2) = "CO"
    species(3) = "O"
    species(4) = "O(1D)"
    species(5) = "H2O"
    species(6) = "H"
    species(7) = "OH"
    species(8) = "H2"
    species(9) = "O3"
    species(10) = "O2"
    species(11) = "HO2"
    species(12) = "H2O2"
    species(13) = "H2CO"
    species(14) = "HCO"
    species(15) = "M"
    species(16) = "HOCO"
    species(17) = "CO2+"

    ! Converting index --------------------------
    fix_out_all_in = 0
    all_in_fix_out = 0
    var_out_all_in = 0
    all_in_var_out = 0
    all_in_var_in  = 0

    all_in_var_in(2) = 1 ! CO: variable
    all_in_var_in(3) = 2 ! O: variable
    all_in_var_in(4) = 3 ! O(1D): variable
    all_in_var_in(6) = 4 ! H: variable
    all_in_var_in(7) = 5 ! OH: variable
    all_in_var_in(8) = 6 ! H2: variable
    all_in_var_in(9) = 7 ! O3: variable
    all_in_var_in(10) = 8 ! O2: variable
    all_in_var_in(11) = 9 ! HO2: variable
    all_in_var_in(12) = 10 ! H2O2: variable
    all_in_var_in(13) = 11 ! H2CO: variable
    all_in_var_in(14) = 12 ! HCO: variable
    all_in_var_in(16) = 13 ! HOCO: variable

    do isp = 1, nsp_fix_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(fix_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          fix_out_all_in(isp) = jsp 
          all_in_fix_out(jsp) = isp 
        end if 
      end do 
    end do 

    do isp = 1, nsp_var_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(var_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          var_out_all_in(isp) = jsp 
          all_in_var_out(jsp) = isp 
        end if 
      end do 
    end do 

    ! Calculating production and loss rates ---

    prod = 0.0_dp
    loss = 0.0_dp

    do iz = 1, nz
    do iy = 1, ny
    do ix = 1, nx

      Tn = T_n(ix,iy,iz)
      Te = T_e(ix,iy,iz)

      n = 0.0_dp
      n_tot = 0.0_dp
      do isp = 1, nsp_fix_out
        jsp = fix_out_all_in(isp)
        if (jsp/=0) then 
          n(jsp) = n_fix(ix,iy,iz,isp) 
          v(1:3,jsp) = v_fix(1:3,ix,iy,iz,isp) 
          Ti(jsp) = T_i(ix,iy,iz,isp) 
          if (trim(ADJUSTL(species(jsp)))/='M') n_tot = n_tot + n(jsp)
        end if 
      end do 
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        if (jsp/=0) then
          n(jsp) = n_var(ix,iy,iz,isp) 
          v(1:3,jsp) = v_var(1:3,ix,iy,iz,isp) 
          Ti(jsp) = T_i(ix,iy,iz,nsp_fix_out+isp)
          if (trim(ADJUSTL(species(jsp)))/='M') n_tot = n_tot + n(jsp)
        end if 
      end do 

      ! Reaction rate coefficient -----------------
      k = 0.0_dp
      k(1) = J_rate_in(ix,iy,iz,1)
      k(2) = J_rate_in(ix,iy,iz,2)
      k(3) = J_rate_in(ix,iy,iz,3)
      k(4) = J_rate_in(ix,iy,iz,4)
      k(5) = J_rate_in(ix,iy,iz,5)
      k(6) = J_rate_in(ix,iy,iz,6)
      k(7) = J_rate_in(ix,iy,iz,7)
      k(8) = J_rate_in(ix,iy,iz,8)
      k(9) = J_rate_in(ix,iy,iz,9)
      k(10) = J_rate_in(ix,iy,iz,10)
      k(11) = J_rate_in(ix,iy,iz,11)
      k(12) = J_rate_in(ix,iy,iz,12)
      k(13) = J_rate_in(ix,iy,iz,13)
      k(14) = J_rate_in(ix,iy,iz,14)
      k(15) = J_rate_in(ix,iy,iz,15)
      k(16) = J_rate_in(ix,iy,iz,16)
      k(17) = 5.4e-33_dp*(300.0_dp/Tn)**3.25_dp
      k(18) = 1.5e-33_dp*(300.0_dp/Tn)**2.4_dp
      k(19) = 8.0e-12_dp*dexp(-2060.0_dp/Tn)
      k(20) = 2.2e-33_dp*dexp(-1780.0_dp/Tn)
      k(21) = 3.2e-11_dp*dexp(70.0_dp/Tn)
      k(22) = 1.2e-10_dp
      k(23) = 1.2e-10_dp
      k(24) = 1.2e-10_dp
      k(25) = 7.5e-11_dp*dexp(115.0_dp/Tn)
      k(26) = 1.63e-10_dp*dexp(60.0_dp/Tn)
      k(27) = 6.34e-12_dp*dexp(-4000.0_dp/Tn)
      k(28) = 9.01e-13_dp*dexp(-1526.0_dp/Tn)
      k(29) = 1.6e-32_dp*(298.0_dp/Tn)**2.27_dp
      k(30) = 1.292e-30_dp*(300.0_dp/Tn)**2.0_dp
      k(31) = 7.2e-11_dp
      k(32) = 1.6e-12_dp
      k(33) = 3.45e-12_dp
      k(34) = 2.8e-12_dp*dexp(-1890.0_dp/Tn)
      k(35) = 1.7e-11_dp*dexp(-1800.0_dp/Tn)
      k0   = 8.8e-32_dp*(300.0_dp/Tn)**1.3_dp
      kinf = 7.5e-11_dp*(300.0_dp/Tn)**(-0.2_dp)
      k(36) = k0 / (1.0_dp + k0*n_tot/kinf) &
        & * 0.6_dp ** (1.0_dp / (1.0_dp + (dlog10(k0*n_tot/kinf))**2.0_dp) )
      k(37) = 1.4e-10_dp*dexp(-470.0_dp/Tn)
      k(38) = 1.8e-11_dp*dexp(180.0_dp/Tn)
      k(39) = 3.0e-11_dp*dexp(200.0_dp/Tn)
      k(40) = 1.4e-12_dp*dexp(-2000.0_dp/Tn)
      k(41) = 1.8e-12_dp
      k0   = 8.97e-31_dp*(300.0_dp/Tn)**1.0_dp
      kinf = 2.6e-11_dp
      k(42) = k0 / (1.0_dp + k0*n_tot/kinf) &
        & * 0.6_dp ** (1.0_dp / (1.0_dp + (dlog10(k0*n_tot/kinf))**2.0_dp) )
      k(43) = 1.7e-12_dp*dexp(-940.0_dp/Tn)
      k(44) = 4.8e-11_dp*dexp(250.0_dp/Tn)
      k(45) = 1.8e-12_dp
      k(46) = 1.0e-14_dp*dexp(-490.0_dp/Tn)
      k(47) = 3.0e-13_dp*dexp(460.0_dp/Tn)
      k(48) = 4.2e-33_dp*dexp(920.0_dp/Tn)
      k0   = 1.5e-13_dp*(300.0_dp/Tn)**(-0.6_dp)
      kinf = 2.1e9_dp*(300.0_dp/Tn)**(-6.1_dp)
      k(49) = k0 / n_tot / (1.0_dp + k0*n_tot/kinf) &
        & * 0.6_dp ** (1.0_dp / (1.0_dp + (dlog10(k0*n_tot/kinf))**2.0_dp) )
      k0   = 5.9e-33_dp*(300.0_dp/Tn)**1.4_dp
      kinf = 1.1e-12_dp*(300.0_dp/Tn)**(-1.3_dp)
      k(50) = k0 / (1.0_dp + k0*n_tot/kinf) &
        & * 0.6_dp ** (1.0_dp / (1.0_dp + (dlog10(k0*n_tot/kinf))**2.0_dp) )
      k(51) = 2.0e-12_dp
      k(52) = 8.7e-10_dp
      k0   = 1.4e-34_dp*dexp(-100.0_dp/Tn)
      kinf = 1.96e-13_dp*dexp(-1366.0_dp/Tn)
      k(53) = k0 / (1.0_dp + k0*n_tot/kinf) &
        & * 0.6_dp ** (1.0_dp / (1.0_dp + (dlog10(k0*n_tot/kinf))**2.0_dp) )
      k(54) = 1.8e-10_dp
      k(55) = 4.5e-11_dp
      k(56) = 1.7e-10_dp
      k(57) = 5.0e-11_dp
      k(58) = 5.0e-11_dp
      k(59) = 2.1e-16_dp*(Tn)**(1.62_dp)*dexp(-1090.0_dp/Tn)
      k(60) = 5.6e-12_dp*(Tn/298.0_dp)**(-0.4_dp)
      k(61) = 8.2e-12_dp*dexp(40.0_dp/Tn)
      k(62) = 3.4e-11_dp*dexp(-1600.0_dp/Tn)
      k(63) = 3.63e-11_dp
      k(64) = 1.0e-11_dp

      k_coef(ix,iy,iz,:) = k(:)

      ! Production rate ---------------------------
      P = 0.0_dp
      P(1) = k(1)*n(1) &
             & + k(2)*n(1) &
             & + k(16)*n(13) &
             & + k(54)*n(6)*n(14) &
             & + k(55)*n(14)*n(14) &
             & + k(56)*n(7)*n(14) &
             & + k(58)*n(3)*n(14) &
             & + k(60)*n(14)*n(10) &
             & + 2.0_dp*k(63)*n(14)*n(14) &
             & + k(64)*n(11)*n(14)
      P(2) = k(1)*n(1) &
             & + k(5)*n(9) &
             & + 2.0_dp*k(7)*n(10) &
             & + k(8)*n(10) &
             & + k(10)*n(7) &
             & + k(12)*n(11) &
             & + k(21)*n(4)*n(10) &
             & + 2.0_dp*k(23)*n(4)*n(9) &
             & + k(25)*n(4)*n(1) &
             & + k(41)*n(7)*n(7)
      P(3) = k(2)*n(1) &
             & + k(4)*n(5) &
             & + k(6)*n(9) &
             & + k(8)*n(10) &
             & + k(11)*n(7) &
             & + k(32)*n(6)*n(11)
      P(4) = k(3)*n(5) &
             & + 2.0_dp*k(9)*n(8) &
             & + k(10)*n(7) &
             & + k(11)*n(7) &
             & + k(14)*n(12) &
             & + k(15)*n(13) &
             & + k(24)*n(4)*n(8) &
             & + k(27)*n(8)*n(3) &
             & + k(28)*n(7)*n(8) &
             & + k(38)*n(3)*n(7) &
             & + k(49)*n(2)*n(7)*n(15) &
             & + 2.0_dp*k(52)*n(17)*n(8) &
             & + k(57)*n(3)*n(14)
      P(5) = k(3)*n(5) &
             & + k(12)*n(11) &
             & + 2.0_dp*k(13)*n(12) &
             & + k(24)*n(4)*n(8) &
             & + 2.0_dp*k(26)*n(4)*n(5) &
             & + k(27)*n(8)*n(3) &
             & + 2.0_dp*k(31)*n(6)*n(11) &
             & + k(35)*n(6)*n(12) &
             & + k(37)*n(6)*n(9) &
             & + k(39)*n(3)*n(11) &
             & + k(40)*n(3)*n(12) &
             & + k(46)*n(11)*n(9) &
             & + k(58)*n(3)*n(14) &
             & + k(62)*n(13)*n(3)
      P(6) = k(4)*n(5) &
             & + k(16)*n(13) &
             & + k(29)*n(6)*n(6)*n(1) &
             & + k(33)*n(6)*n(11) &
             & + k(34)*n(6)*n(12) &
             & + k(54)*n(6)*n(14) &
             & + k(59)*n(13)*n(6) &
             & + k(63)*n(14)*n(14)
      P(7) = k(18)*n(3)*n(10)*n(1)
      P(8) = k(5)*n(9) &
             & + k(6)*n(9) &
             & + k(17)*n(3)*n(3)*n(15) &
             & + 2.0_dp*k(19)*n(3)*n(9) &
             & + k(21)*n(4)*n(10) &
             & + 2.0_dp*k(22)*n(4)*n(9) &
             & + k(23)*n(4)*n(9) &
             & + k(33)*n(6)*n(11) &
             & + k(37)*n(6)*n(9) &
             & + k(38)*n(3)*n(7) &
             & + k(39)*n(3)*n(11) &
             & + k(43)*n(7)*n(9) &
             & + k(44)*n(7)*n(11) &
             & + 2.0_dp*k(46)*n(11)*n(9) &
             & + k(47)*n(11)*n(11) &
             & + k(48)*n(11)*n(11)*n(15)
      P(9) = k(14)*n(12) &
             & + k(34)*n(6)*n(12) &
             & + k(36)*n(6)*n(10)*n(15) &
             & + k(40)*n(3)*n(12) &
             & + k(43)*n(7)*n(9) &
             & + k(45)*n(7)*n(12) &
             & + k(51)*n(16)*n(10) &
             & + k(60)*n(14)*n(10)
      P(10) = k(42)*n(7)*n(7)*n(15) &
             & + k(47)*n(11)*n(11) &
             & + k(48)*n(11)*n(11)*n(15) &
             & + k(64)*n(11)*n(14)
      P(11) = k(55)*n(14)*n(14)
      P(12) = k(15)*n(13) &
             & + k(53)*n(6)*n(2)*n(15) &
             & + k(59)*n(13)*n(6) &
             & + k(61)*n(13)*n(7) &
             & + k(62)*n(13)*n(3)
      P(13) = k(50)*n(2)*n(7)*n(15)

      ! Loss rate ---------------------------
      L = 0.0_dp
      L(1) = k(20)*n(3)*n(2)*n(15) &
             & + k(49)*n(2)*n(7)*n(15) &
             & + k(50)*n(2)*n(7)*n(15) &
             & + k(53)*n(6)*n(2)*n(15)
      L(2) = 2.0_dp*k(17)*n(3)*n(3)*n(15) &
             & + k(18)*n(3)*n(10)*n(1) &
             & + k(19)*n(3)*n(9) &
             & + k(20)*n(3)*n(2)*n(15) &
             & + k(27)*n(8)*n(3) &
             & + k(38)*n(3)*n(7) &
             & + k(39)*n(3)*n(11) &
             & + k(40)*n(3)*n(12) &
             & + k(57)*n(3)*n(14) &
             & + k(58)*n(3)*n(14) &
             & + k(62)*n(13)*n(3)
      L(3) = k(21)*n(4)*n(10) &
             & + k(22)*n(4)*n(9) &
             & + k(23)*n(4)*n(9) &
             & + k(24)*n(4)*n(8) &
             & + k(25)*n(4)*n(1) &
             & + k(26)*n(4)*n(5)
      L(4) = 2.0_dp*k(29)*n(6)*n(6)*n(1) &
             & + k(30)*n(6)*n(7)*n(1) &
             & + k(31)*n(6)*n(11) &
             & + k(32)*n(6)*n(11) &
             & + k(33)*n(6)*n(11) &
             & + k(34)*n(6)*n(12) &
             & + k(35)*n(6)*n(12) &
             & + k(36)*n(6)*n(10)*n(15) &
             & + k(37)*n(6)*n(9) &
             & + k(53)*n(6)*n(2)*n(15) &
             & + k(54)*n(6)*n(14) &
             & + k(59)*n(13)*n(6)
      L(5) = k(10)*n(7) &
             & + k(11)*n(7) &
             & + k(28)*n(7)*n(8) &
             & + k(30)*n(6)*n(7)*n(1) &
             & + k(38)*n(3)*n(7) &
             & + 2.0_dp*k(41)*n(7)*n(7) &
             & + 2.0_dp*k(42)*n(7)*n(7)*n(15) &
             & + k(43)*n(7)*n(9) &
             & + k(44)*n(7)*n(11) &
             & + k(45)*n(7)*n(12) &
             & + k(49)*n(2)*n(7)*n(15) &
             & + k(50)*n(2)*n(7)*n(15) &
             & + k(56)*n(7)*n(14) &
             & + k(61)*n(13)*n(7)
      L(6) = k(9)*n(8) &
             & + k(24)*n(4)*n(8) &
             & + k(27)*n(8)*n(3) &
             & + k(28)*n(7)*n(8) &
             & + k(52)*n(17)*n(8)
      L(7) = k(5)*n(9) &
             & + k(6)*n(9) &
             & + k(19)*n(3)*n(9) &
             & + k(22)*n(4)*n(9) &
             & + k(23)*n(4)*n(9) &
             & + k(37)*n(6)*n(9) &
             & + k(43)*n(7)*n(9) &
             & + k(46)*n(11)*n(9)
      L(8) = k(7)*n(10) &
             & + k(8)*n(10) &
             & + k(18)*n(3)*n(10)*n(1) &
             & + k(21)*n(4)*n(10) &
             & + k(36)*n(6)*n(10)*n(15) &
             & + k(51)*n(16)*n(10) &
             & + k(60)*n(14)*n(10)
      L(9) = k(12)*n(11) &
             & + k(31)*n(6)*n(11) &
             & + k(32)*n(6)*n(11) &
             & + k(33)*n(6)*n(11) &
             & + k(39)*n(3)*n(11) &
             & + k(44)*n(7)*n(11) &
             & + k(46)*n(11)*n(9) &
             & + 2.0_dp*k(47)*n(11)*n(11) &
             & + 2.0_dp*k(48)*n(11)*n(11)*n(15) &
             & + k(64)*n(11)*n(14)
      L(10) = k(13)*n(12) &
             & + k(14)*n(12) &
             & + k(34)*n(6)*n(12) &
             & + k(35)*n(6)*n(12) &
             & + k(40)*n(3)*n(12) &
             & + k(45)*n(7)*n(12)
      L(11) = k(15)*n(13) &
             & + k(16)*n(13) &
             & + k(59)*n(13)*n(6) &
             & + k(61)*n(13)*n(7) &
             & + k(62)*n(13)*n(3)
      L(12) = k(54)*n(6)*n(14) &
             & + 2.0_dp*k(55)*n(14)*n(14) &
             & + k(56)*n(7)*n(14) &
             & + k(57)*n(3)*n(14) &
             & + k(58)*n(3)*n(14) &
             & + k(60)*n(14)*n(10) &
             & + 2.0_dp*k(63)*n(14)*n(14) &
             & + k(64)*n(11)*n(14)
      L(13) = k(51)*n(16)*n(10)


      ! Converting index --------------------------
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        ksp = all_in_var_in(jsp)
        prod(ix,iy,iz,isp) = P(ksp)
        loss(ix,iy,iz,isp) = L(ksp)
      end do 

    end do 
    end do 
    end do 

  end subroutine p__PROTEUS_source



  subroutine p__PROTEUS_Jacobian(var_species_list, fix_species_list, & ! in:  name of variable species and fixed species
    &                            nsp_var_out, nsp_fix_out,           & ! in:  number of variable species and fixed species
    &                            nz,                                 & ! in:  number of vertical grids
    &                            n_var, n_fix,                       & ! in:  number density of variable and fixed species
    &                            k_coef,                             & ! in:  reaction rate coefficients
    &                            Jmtx                                ) ! out: chemical Jacobian matrix
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    !   * Both arrays are strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables and fixed, respectively. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of species_list string arrays.
    !
    ! - nsp_var_out, nsp_fix_out
    !   * The number of chemical species to be treated as variables and fixed, respectively. 
    !
    ! - nz
    !   * The number of vertical grids. 
    !
    ! - n_var(nx,ny,nz,nsp_var_out), n_fix(nx,ny,nz,nsp_fix_out)
    !   * Number density of variable and fixed species in [cm^-3] in simulation grids. 
    !
    ! - k_coef(nx,ny,nz,nch)
    !   * Reaction rate coefficients for each reaction in the simulation grids in the units as follows:
    !     [s^-1] for reactions with only one reactant 
    !     [cm^3 s^-1] for two-body reactions
    !     [cm^6 s^-1] for three-body reactions
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - Jmtx(1:,1:)
    !   * If nz >= 2, Jmtx is the transposed chemical Jacobian matrix with a dimension (2*nsp_var_out+1,nsp_var_out*nz).
    !   * If nz == 1, Jmtx is the chemical Jacobian matrix with a dimension (nsp_var_out,nsp_var_out).
    !
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer, parameter :: nsp = 17
    integer, parameter :: nsp_fix_in = 4
    integer, parameter :: nsp_var_in = 13
    integer, parameter :: nch = 64
    integer,          intent(in)  :: nsp_var_out, nsp_fix_out, nz
    character(len=*), intent(in)  :: var_species_list(nsp_var_out), fix_species_list(nsp_fix_out)
    real(dp),         intent(in)  :: n_var(nz,nsp_var_out), n_fix(nz,nsp_fix_out)
    real(dp),         intent(in)  :: k_coef(nz,nch)
    real(dp),         intent(out) :: Jmtx(1:,1:)
    integer isp, jsp, i0, i1, j0, j1, i, j, iz
    
    integer  all_in_var_out(nsp), var_out_all_in(nsp_var_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_fix_out(nsp), fix_out_all_in(nsp_fix_out) ! in: PROTEUS index, out: outside model index
    integer  all_in_var_in(nsp) ! 
    real(dp) k(nch), n(nsp), Jmtx_tmp(nsp_var_in,nsp_var_in)
    character(len = 256) species(nsp)
    
    ! species list ----------------------------
    ! 
    ! CO2, CO, O, O(1D), H2O, H, OH, H2, O3, O2, HO2, H2O2, H2CO, HCO, M, HOCO, CO2+
    species(1) = "CO2"
    species(2) = "CO"
    species(3) = "O"
    species(4) = "O(1D)"
    species(5) = "H2O"
    species(6) = "H"
    species(7) = "OH"
    species(8) = "H2"
    species(9) = "O3"
    species(10) = "O2"
    species(11) = "HO2"
    species(12) = "H2O2"
    species(13) = "H2CO"
    species(14) = "HCO"
    species(15) = "M"
    species(16) = "HOCO"
    species(17) = "CO2+"

    ! Converting index --------------------------
    fix_out_all_in = 0
    all_in_fix_out = 0
    var_out_all_in = 0
    all_in_var_out = 0
    all_in_var_in  = 0

    all_in_var_in(2) = 1 ! CO: variable
    all_in_var_in(3) = 2 ! O: variable
    all_in_var_in(4) = 3 ! O(1D): variable
    all_in_var_in(6) = 4 ! H: variable
    all_in_var_in(7) = 5 ! OH: variable
    all_in_var_in(8) = 6 ! H2: variable
    all_in_var_in(9) = 7 ! O3: variable
    all_in_var_in(10) = 8 ! O2: variable
    all_in_var_in(11) = 9 ! HO2: variable
    all_in_var_in(12) = 10 ! H2O2: variable
    all_in_var_in(13) = 11 ! H2CO: variable
    all_in_var_in(14) = 12 ! HCO: variable
    all_in_var_in(16) = 13 ! HOCO: variable

    do isp = 1, nsp_fix_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(fix_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          fix_out_all_in(isp) = jsp 
          all_in_fix_out(jsp) = isp 
        end if 
      end do 
    end do 

    do isp = 1, nsp_var_out
      do jsp = 1, nsp
        if (trim(ADJUSTL(var_species_list(isp)))==trim(ADJUSTL(species(jsp)))) then 
          var_out_all_in(isp) = jsp 
          all_in_var_out(jsp) = isp 
        end if 
      end do 
    end do 

    ! Calculating chemical Jacobian matrix ---
    Jmtx = 0.0_dp

    do iz = 1, nz

      n = 0.0_dp
      do isp = 1, nsp_fix_out
        jsp = fix_out_all_in(isp)
        if (jsp/=0) then 
          n(jsp) = n_fix(iz,isp) 
        end if 
      end do 
      do isp = 1, nsp_var_out
        jsp = var_out_all_in(isp)
        if (jsp/=0) then
          n(jsp) = n_var(iz,isp) 
        end if 
      end do 

      k(1:nch) = k_coef(iz,1:nch)

      Jmtx_tmp = 0.0_dp

      Jmtx_tmp(1,1) = - k(20)*n(3)*n(15) &
                  & - k(49)*n(7)*n(15) &
                  & - k(50)*n(7)*n(15) &
                  & - k(53)*n(6)*n(15)

      Jmtx_tmp(1,2) = k(58)*n(14) &
                  & - k(20)*n(2)*n(15)

      Jmtx_tmp(1,4) = k(54)*n(14) &
                  & - k(53)*n(2)*n(15)

      Jmtx_tmp(1,5) = k(56)*n(14) &
                  & - k(49)*n(2)*n(15) &
                  & - k(50)*n(2)*n(15)

      Jmtx_tmp(1,8) = k(60)*n(14)

      Jmtx_tmp(1,9) = k(64)*n(14)

      Jmtx_tmp(1,11) = k(16)

      Jmtx_tmp(1,12) = k(54)*n(6) &
                  & + k(55)*2.0_dp*n(14) &
                  & + k(56)*n(7) &
                  & + k(58)*n(3) &
                  & + k(60)*n(10) &
                  & + 2.0_dp*k(63)*2.0_dp*n(14) &
                  & + k(64)*n(11)

      Jmtx_tmp(2,1) = - k(20)*n(3)*n(15)

      Jmtx_tmp(2,2) = - 2.0_dp*k(17)*n(15)*2.0_dp*n(3) &
                  & - k(18)*n(10)*n(1) &
                  & - k(19)*n(9) &
                  & - k(20)*n(2)*n(15) &
                  & - k(27)*n(8) &
                  & - k(38)*n(7) &
                  & - k(39)*n(11) &
                  & - k(40)*n(12) &
                  & - k(57)*n(14) &
                  & - k(58)*n(14) &
                  & - k(62)*n(13)

      Jmtx_tmp(2,3) = k(21)*n(10) &
                  & + 2.0_dp*k(23)*n(9) &
                  & + k(25)*n(1)

      Jmtx_tmp(2,5) = k(10) &
                  & + k(41)*2.0_dp*n(7) &
                  & - k(38)*n(3)

      Jmtx_tmp(2,6) = - k(27)*n(3)

      Jmtx_tmp(2,7) = k(5) &
                  & + 2.0_dp*k(23)*n(4) &
                  & - k(19)*n(3)

      Jmtx_tmp(2,8) = 2.0_dp*k(7) &
                  & + k(8) &
                  & + k(21)*n(4) &
                  & - k(18)*n(3)*n(1)

      Jmtx_tmp(2,9) = k(12) &
                  & - k(39)*n(3)

      Jmtx_tmp(2,10) = - k(40)*n(3)

      Jmtx_tmp(2,11) = - k(62)*n(3)

      Jmtx_tmp(2,12) = - k(57)*n(3) &
                  & - k(58)*n(3)

      Jmtx_tmp(3,3) = - k(21)*n(10) &
                  & - k(22)*n(9) &
                  & - k(23)*n(9) &
                  & - k(24)*n(8) &
                  & - k(25)*n(1) &
                  & - k(26)*n(5)

      Jmtx_tmp(3,4) = k(32)*n(11)

      Jmtx_tmp(3,5) = k(11)

      Jmtx_tmp(3,6) = - k(24)*n(4)

      Jmtx_tmp(3,7) = k(6) &
                  & - k(22)*n(4) &
                  & - k(23)*n(4)

      Jmtx_tmp(3,8) = k(8) &
                  & - k(21)*n(4)

      Jmtx_tmp(3,9) = k(32)*n(6)

      Jmtx_tmp(4,1) = k(49)*n(7)*n(15) &
                  & - k(53)*n(6)*n(15)

      Jmtx_tmp(4,2) = k(27)*n(8) &
                  & + k(38)*n(7) &
                  & + k(57)*n(14)

      Jmtx_tmp(4,3) = k(24)*n(8)

      Jmtx_tmp(4,4) = - 2.0_dp*k(29)*n(1)*2.0_dp*n(6) &
                  & - k(30)*n(7)*n(1) &
                  & - k(31)*n(11) &
                  & - k(32)*n(11) &
                  & - k(33)*n(11) &
                  & - k(34)*n(12) &
                  & - k(35)*n(12) &
                  & - k(36)*n(10)*n(15) &
                  & - k(37)*n(9) &
                  & - k(53)*n(2)*n(15) &
                  & - k(54)*n(14) &
                  & - k(59)*n(13)

      Jmtx_tmp(4,5) = k(10) &
                  & + k(11) &
                  & + k(28)*n(8) &
                  & + k(38)*n(3) &
                  & + k(49)*n(2)*n(15) &
                  & - k(30)*n(6)*n(1)

      Jmtx_tmp(4,6) = 2.0_dp*k(9) &
                  & + k(24)*n(4) &
                  & + k(27)*n(3) &
                  & + k(28)*n(7) &
                  & + 2.0_dp*k(52)*n(17)

      Jmtx_tmp(4,7) = - k(37)*n(6)

      Jmtx_tmp(4,8) = - k(36)*n(6)*n(15)

      Jmtx_tmp(4,9) = - k(31)*n(6) &
                  & - k(32)*n(6) &
                  & - k(33)*n(6)

      Jmtx_tmp(4,10) = k(14) &
                  & - k(34)*n(6) &
                  & - k(35)*n(6)

      Jmtx_tmp(4,11) = k(15) &
                  & - k(59)*n(6)

      Jmtx_tmp(4,12) = k(57)*n(3) &
                  & - k(54)*n(6)

      Jmtx_tmp(5,1) = - k(49)*n(7)*n(15) &
                  & - k(50)*n(7)*n(15)

      Jmtx_tmp(5,2) = k(27)*n(8) &
                  & + k(39)*n(11) &
                  & + k(40)*n(12) &
                  & + k(58)*n(14) &
                  & + k(62)*n(13) &
                  & - k(38)*n(7)

      Jmtx_tmp(5,3) = k(24)*n(8) &
                  & + 2.0_dp*k(26)*n(5)

      Jmtx_tmp(5,4) = 2.0_dp*k(31)*n(11) &
                  & + k(35)*n(12) &
                  & + k(37)*n(9) &
                  & - k(30)*n(7)*n(1)

      Jmtx_tmp(5,5) = - k(10) &
                  & - k(11) &
                  & - k(28)*n(8) &
                  & - k(30)*n(6)*n(1) &
                  & - k(38)*n(3) &
                  & - 2.0_dp*k(41)*2.0_dp*n(7) &
                  & - 2.0_dp*k(42)*n(15)*2.0_dp*n(7) &
                  & - k(43)*n(9) &
                  & - k(44)*n(11) &
                  & - k(45)*n(12) &
                  & - k(49)*n(2)*n(15) &
                  & - k(50)*n(2)*n(15) &
                  & - k(56)*n(14) &
                  & - k(61)*n(13)

      Jmtx_tmp(5,6) = k(24)*n(4) &
                  & + k(27)*n(3) &
                  & - k(28)*n(7)

      Jmtx_tmp(5,7) = k(37)*n(6) &
                  & + k(46)*n(11) &
                  & - k(43)*n(7)

      Jmtx_tmp(5,9) = k(12) &
                  & + 2.0_dp*k(31)*n(6) &
                  & + k(39)*n(3) &
                  & + k(46)*n(9) &
                  & - k(44)*n(7)

      Jmtx_tmp(5,10) = 2.0_dp*k(13) &
                  & + k(35)*n(6) &
                  & + k(40)*n(3) &
                  & - k(45)*n(7)

      Jmtx_tmp(5,11) = k(62)*n(3) &
                  & - k(61)*n(7)

      Jmtx_tmp(5,12) = k(58)*n(3) &
                  & - k(56)*n(7)

      Jmtx_tmp(6,2) = - k(27)*n(8)

      Jmtx_tmp(6,3) = - k(24)*n(8)

      Jmtx_tmp(6,4) = k(29)*n(1)*2.0_dp*n(6) &
                  & + k(33)*n(11) &
                  & + k(34)*n(12) &
                  & + k(54)*n(14) &
                  & + k(59)*n(13)

      Jmtx_tmp(6,5) = - k(28)*n(8)

      Jmtx_tmp(6,6) = - k(9) &
                  & - k(24)*n(4) &
                  & - k(27)*n(3) &
                  & - k(28)*n(7) &
                  & - k(52)*n(17)

      Jmtx_tmp(6,9) = k(33)*n(6)

      Jmtx_tmp(6,10) = k(34)*n(6)

      Jmtx_tmp(6,11) = k(16) &
                  & + k(59)*n(6)

      Jmtx_tmp(6,12) = k(54)*n(6) &
                  & + k(63)*2.0_dp*n(14)

      Jmtx_tmp(7,2) = k(18)*n(10)*n(1) &
                  & - k(19)*n(9)

      Jmtx_tmp(7,3) = - k(22)*n(9) &
                  & - k(23)*n(9)

      Jmtx_tmp(7,4) = - k(37)*n(9)

      Jmtx_tmp(7,5) = - k(43)*n(9)

      Jmtx_tmp(7,7) = - k(5) &
                  & - k(6) &
                  & - k(19)*n(3) &
                  & - k(22)*n(4) &
                  & - k(23)*n(4) &
                  & - k(37)*n(6) &
                  & - k(43)*n(7) &
                  & - k(46)*n(11)

      Jmtx_tmp(7,8) = k(18)*n(3)*n(1)

      Jmtx_tmp(7,9) = - k(46)*n(9)

      Jmtx_tmp(8,2) = k(17)*n(15)*2.0_dp*n(3) &
                  & + 2.0_dp*k(19)*n(9) &
                  & + k(38)*n(7) &
                  & + k(39)*n(11) &
                  & - k(18)*n(10)*n(1)

      Jmtx_tmp(8,3) = k(21)*n(10) &
                  & + 2.0_dp*k(22)*n(9) &
                  & + k(23)*n(9) &
                  & - k(21)*n(10)

      Jmtx_tmp(8,4) = k(33)*n(11) &
                  & + k(37)*n(9) &
                  & - k(36)*n(10)*n(15)

      Jmtx_tmp(8,5) = k(38)*n(3) &
                  & + k(43)*n(9) &
                  & + k(44)*n(11)

      Jmtx_tmp(8,7) = k(5) &
                  & + k(6) &
                  & + 2.0_dp*k(19)*n(3) &
                  & + 2.0_dp*k(22)*n(4) &
                  & + k(23)*n(4) &
                  & + k(37)*n(6) &
                  & + k(43)*n(7) &
                  & + 2.0_dp*k(46)*n(11)

      Jmtx_tmp(8,8) = k(21)*n(4) &
                  & - k(7) &
                  & - k(8) &
                  & - k(18)*n(3)*n(1) &
                  & - k(21)*n(4) &
                  & - k(36)*n(6)*n(15) &
                  & - k(51)*n(16) &
                  & - k(60)*n(14)

      Jmtx_tmp(8,9) = k(33)*n(6) &
                  & + k(39)*n(3) &
                  & + k(44)*n(7) &
                  & + 2.0_dp*k(46)*n(9) &
                  & + k(47)*2.0_dp*n(11) &
                  & + k(48)*n(15)*2.0_dp*n(11)

      Jmtx_tmp(8,12) = - k(60)*n(10)

      Jmtx_tmp(8,13) = - k(51)*n(10)

      Jmtx_tmp(9,2) = k(40)*n(12) &
                  & - k(39)*n(11)

      Jmtx_tmp(9,4) = k(34)*n(12) &
                  & + k(36)*n(10)*n(15) &
                  & - k(31)*n(11) &
                  & - k(32)*n(11) &
                  & - k(33)*n(11)

      Jmtx_tmp(9,5) = k(43)*n(9) &
                  & + k(45)*n(12) &
                  & - k(44)*n(11)

      Jmtx_tmp(9,7) = k(43)*n(7) &
                  & - k(46)*n(11)

      Jmtx_tmp(9,8) = k(36)*n(6)*n(15) &
                  & + k(51)*n(16) &
                  & + k(60)*n(14)

      Jmtx_tmp(9,9) = - k(12) &
                  & - k(31)*n(6) &
                  & - k(32)*n(6) &
                  & - k(33)*n(6) &
                  & - k(39)*n(3) &
                  & - k(44)*n(7) &
                  & - k(46)*n(9) &
                  & - 2.0_dp*k(47)*2.0_dp*n(11) &
                  & - 2.0_dp*k(48)*n(15)*2.0_dp*n(11) &
                  & - k(64)*n(14)

      Jmtx_tmp(9,10) = k(14) &
                  & + k(34)*n(6) &
                  & + k(40)*n(3) &
                  & + k(45)*n(7)

      Jmtx_tmp(9,12) = k(60)*n(10) &
                  & - k(64)*n(11)

      Jmtx_tmp(9,13) = k(51)*n(10)

      Jmtx_tmp(10,2) = - k(40)*n(12)

      Jmtx_tmp(10,4) = - k(34)*n(12) &
                  & - k(35)*n(12)

      Jmtx_tmp(10,5) = k(42)*n(15)*2.0_dp*n(7) &
                  & - k(45)*n(12)

      Jmtx_tmp(10,9) = k(47)*2.0_dp*n(11) &
                  & + k(48)*n(15)*2.0_dp*n(11) &
                  & + k(64)*n(14)

      Jmtx_tmp(10,10) = - k(13) &
                  & - k(14) &
                  & - k(34)*n(6) &
                  & - k(35)*n(6) &
                  & - k(40)*n(3) &
                  & - k(45)*n(7)

      Jmtx_tmp(10,12) = k(64)*n(11)

      Jmtx_tmp(11,2) = - k(62)*n(13)

      Jmtx_tmp(11,4) = - k(59)*n(13)

      Jmtx_tmp(11,5) = - k(61)*n(13)

      Jmtx_tmp(11,11) = - k(15) &
                  & - k(16) &
                  & - k(59)*n(6) &
                  & - k(61)*n(7) &
                  & - k(62)*n(3)

      Jmtx_tmp(11,12) = k(55)*2.0_dp*n(14)

      Jmtx_tmp(12,1) = k(53)*n(6)*n(15)

      Jmtx_tmp(12,2) = k(62)*n(13) &
                  & - k(57)*n(14) &
                  & - k(58)*n(14)

      Jmtx_tmp(12,4) = k(53)*n(2)*n(15) &
                  & + k(59)*n(13) &
                  & - k(54)*n(14)

      Jmtx_tmp(12,5) = k(61)*n(13) &
                  & - k(56)*n(14)

      Jmtx_tmp(12,8) = - k(60)*n(14)

      Jmtx_tmp(12,9) = - k(64)*n(14)

      Jmtx_tmp(12,11) = k(15) &
                  & + k(59)*n(6) &
                  & + k(61)*n(7) &
                  & + k(62)*n(3)

      Jmtx_tmp(12,12) = - k(54)*n(6) &
                  & - 2.0_dp*k(55)*2.0_dp*n(14) &
                  & - k(56)*n(7) &
                  & - k(57)*n(3) &
                  & - k(58)*n(3) &
                  & - k(60)*n(10) &
                  & - 2.0_dp*k(63)*2.0_dp*n(14) &
                  & - k(64)*n(11)

      Jmtx_tmp(13,1) = k(50)*n(7)*n(15)

      Jmtx_tmp(13,5) = k(50)*n(2)*n(15)

      Jmtx_tmp(13,8) = - k(51)*n(16)

      Jmtx_tmp(13,13) = - k(51)*n(10)


      ! Converting index --------------------------
      if (nz >= 2) then
        do i1 = 1, nsp_var_out
          jsp = var_out_all_in(i1)
          i0 = all_in_var_in(jsp)
          do j1 = 1, nsp_var_out
            jsp = var_out_all_in(j1)
            j0 = all_in_var_in(jsp)
            i = (iz-1)*nsp_var_out+i1
            j = (iz-1)*nsp_var_out+j1 + nsp_var_out + 1 - i
            Jmtx(j,i) = - Jmtx_tmp(i0,j0)
          end do 
        end do 
      else if (nz == 1) then
        do i1 = 1, nsp_var_out
          jsp = var_out_all_in(i1)
          i0 = all_in_var_in(jsp)
          do j1 = 1, nsp_var_out
            jsp = var_out_all_in(j1)
            j0 = all_in_var_in(jsp)
            Jmtx(j0,i0) = - Jmtx_tmp(i0,j0)
          end do 
        end do 
      end if

    end do 

  end subroutine p__PROTEUS_Jacobian


  subroutine get_solar_flux_data(lambda, dlambda, nwl, fname, unit1, unit2, & ! in
    &                            solar_flux                                 ) ! out
    implicit none
    real(dp), intent(in)  :: lambda(1:), dlambda(1:)
    integer,  intent(in)  :: nwl
    character(len=*), intent(in) :: unit1, unit2, fname
    real(dp), intent(out) :: solar_flux(1:)
    integer i, nh, il, nl
    real(dp), allocatable :: idata(:,:), odata(:)

    write(*,'(a)',advance='no')  '  Reading datafile: '//trim(ADJUSTL(fname))//'...'

    allocate(odata(nwl))

    call get_header_line_number(nh, fname)
    call get_data_line_number(nl, fname)
    allocate(idata(nl,2))

    open(11, file = fname, status = 'old')
      do i = 1, nh; read(11,*); end do
      do il = 1, nl
        read(11,*) idata(il,1), idata(il,2) 

        ! Unit conversion: column 1
        if (unit1/='nm' .and. unit1/='A' .and. unit1/='/cm' .and. unit1/='cm^-1') then 
          print *, ''
          print *, 'error!'
          print *, '  The unit of the datafile "'//trim(adjustl(fname))//'" column 1 is not recognized.'
          print *, '  The unit of column 1 should be "nm", "A", "/cm" or "cm^-1".'
          stop
        end if
        if (unit1 == 'A') then 
          idata(il,1) = idata(il,1) * 1.0e-1_dp ! [A -> nm]
        else if (unit1 == '/cm' .or. unit1 == 'cm^-1') then 
          idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
        end if

        ! Unit conversion: column 2
        if (unit2 == '/cm^2/s/nm') then 
          idata(il,2) = idata(il,2) * 1.0e4_dp ! [/cm^2/s/nm -> /m^2/s/nm]
        else if (unit2 == 'W/m^2/nm') then 
          idata(il,2) = idata(il,2) * idata(il,1) * 1.0e-9_dp / 6.626e-34_dp / 2.99792458e8_dp ! [W/m^2/nm  ->  /m^2/s/nm]
        end if

      end do
    close(11)

    write(*,*) 'done.'

    call binning_flux(lambda, dlambda, idata, nl, nwl, odata)
    solar_flux(1:nwl) = odata(1:nwl)
    deallocate(idata)
    deallocate(odata)

  end subroutine get_solar_flux_data


  subroutine binning_flux(lambda, dlambda, idata, nl, nwl, & ! in
    &                     odata                            ) ! out
    implicit none
    real(dp),   intent(in)  :: lambda(nwl), dlambda(nwl)
    integer,    intent(in)  :: nl, nwl
    real(dp),   intent(in)  :: idata(nl,2)
    real(dp),   intent(out) :: odata(nwl)
    integer iwl, il, label, il0, ndata
    real(dp) idata_tmp(nl,2)
    real(dp) i0, ip, im, dip, dim, o0, dop, dom, sum_data, sum_wl, tmp
  
    if (idata(1,1) > idata(nl,1)) then 
      do il = 1, nl
        idata_tmp(il,1) = idata(nl-il+1,1)
        idata_tmp(il,2) = idata(nl-il+1,2)
      end do 
    else 
      idata_tmp(:,:) = idata(:,:)
    end if
  
    odata = 0.0_dp
    il0 = 1
  
    do iwl = 1, nwl 
  
      label = -1
      ndata = 0
      sum_data = 0.0_dp
      sum_wl   = 0.0_dp
      tmp = 0.0_dp
  
      loop: do il = il0, nl 
        
        ! input bin
        if (il == 1) then 
          i0 = idata_tmp(il,1)
          ip = idata_tmp(il+1,1)
          dip = ip - i0
          dim = dip
        else if (il > 1 .and. il < nl) then 
          im = idata_tmp(il-1,1)
          i0 = idata_tmp(il,1)
          ip = idata_tmp(il+1,1)
          dim = i0 - im
          dip = ip - i0
        else if (il == nl) then 
          im = idata_tmp(il-1,1)
          i0 = idata_tmp(il,1)
          dim = i0 - im
          dip = dim
        end if
  
        ! model bin
        if (iwl == 1) then 
          o0 = lambda(1)
          dop = dlambda(1) * 0.5d0
          dom = dlambda(1) * 0.5d0
        else if (iwl > 1 .and. iwl < nwl) then 
          o0 = lambda(iwl)
          dom = dlambda(iwl) * 0.5d0
          dop = dlambda(iwl) * 0.5d0
        else if (iwl == nwl) then 
          o0 = lambda(iwl)
          dom = dlambda(iwl) * 0.5d0
          dop = 1.0e10_dp ! to take large enough value
        end if 
        
        ! for binning
        if (o0-dom <= i0 .and. i0 <= o0+dop) then 
          sum_data = sum_data + idata_tmp(il,2)*(dim+dip)
          sum_wl   = sum_wl + (dim+dip)
          ndata = ndata+1
          if (ndata >= 2) label = 1
          if (iwl==nwl .or. il==nl) label = 1
          if (iwl==1 .or. il==1) label = 1
        end if
        
        ! for interpolate
        if (i0 <= o0 .and. o0 < ip) then 
          tmp = (idata_tmp(il,2)*(ip-o0) + idata_tmp(il+1,2)*(o0-i0))/(ip-i0)
          if (label /= 1) label = 0
        end if
  
        ! to escape from this loop
        if (label == 1 .and. i0 > o0+dop) then 
          if (il >= 2) il0 = il-1
          exit loop
        else if (label == 0 .and. ip < o0) then 
          if (il >= 2) il0 = il-1
          exit loop
        end if
        
      end do loop ! il
  
      ! binning if an input bin is smaller than a model bin
      if (label == 1) then 
        odata(iwl) = sum_data / sum_wl
      end if
  
      ! interpolate if an input bin is larger than a model bin
      if (label == 0) then 
        odata(iwl) = tmp
      end if
  
    end do ! iwl
  
  end subroutine binning_flux
  
  
  subroutine get_header_line_number(nh, fname)
    implicit none
    integer nh
    character(len=256) fname
    character(len=256) strm
    open(11, file = fname, status = 'old' )
      nh = 0
      do
        read(11,'(A)') strm
        if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
          nh = nh + 1
        else
          exit
        end if
      end do
    close(11)
  end subroutine get_header_line_number
  
  
  subroutine get_data_line_number(nl, fname)
    implicit none
    integer nl, ios
    character(len=256) fname
    character(len=256) strm
    open(11, file = fname, status = 'old' )
      nl = 0
      do
        read(11,'(A)',iostat = ios) strm
        if (ios < 0) exit
        if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
          cycle
        else
          nl = nl + 1
        endif
      end do
    close(11)
  end subroutine get_data_line_number
  
  
end module p__PROTEUS
