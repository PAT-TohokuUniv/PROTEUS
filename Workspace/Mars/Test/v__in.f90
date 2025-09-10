module v__in
  use v__tdec,   only : set_, var_, cst_, spl_, grd_
  use p__search, only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__in__ini, v__in__exe

contains


  subroutine v__in__ini(spl, set) ! out
    type(spl_),           intent(out) :: spl
    type(set_),           intent(out) :: set

    ! Planet and project
    spl%planet = 'Mars'
    set%dir_name = './Mars/Test'
    set%dir_name_win = '.\Mars\Test'

    ! Calculation settings
    set%mode = '1D'
    set%use_1d = 1
    set%use_2d = 0
    set%F107 = 140.0_dp
    set%nstep = 10000
    set%fin_sec = 3.0e15_dp
    set%dtime_limit = 1.0e14_dp
    set%latitude = 0.0_dp
    set%sza = 0.0_dp
    set%lat_res = 1.0_dp
    set%lt_res = 0.1_dp
    set%nx = 1
    set%ny = 1
    set%fix_sza = 1
    set%diurnal_ave = 1
    set%Ls = 0.0_dp
    set%nday = 2.0_dp
    set%scheme = 'implicit'
    set%inversion = 'default'
    set%dt_rate = 10.0_dp
    set%dt_inc_eps = 0.1_dp
    set%max_eps = 1.0e4_dp
    set%solar_flux = './UV/ref_solar_irradiance_whi-2008_ver2_1.dat'
    set%solar_flux_unit1 = 'nm'
    set%solar_flux_unit2 = 'W/m^2/nm'
    set%euv_input = 'EUV-DATA'
    set%euv_factor = 1.0_dp
    set%n_wl_bin = 1
    allocate(set%wl_bin(set%n_wl_bin,3))
    set%wl_bin(1,1) = 0.0_dp
    set%wl_bin(1,2) = 1000.0_dp
    set%wl_bin(1,3) = 1.0_dp
    set%nsp_tout = 3
    allocate(set%species_tout(set%nsp_tout))
    set%species_tout(1) = "O"
    set%species_tout(2) = "H"
    set%species_tout(3) = "H2"
    set%dt_out = 100.0_dp

  end subroutine v__in__ini


  subroutine v__in__exe(cst,      & ! in
    &                   set, spl, & ! inout
    &                   var, grd  ) ! out
    type(cst_),           intent(in)    :: cst
    type(set_),           intent(inout) :: set
    type(spl_),           intent(inout) :: spl
    type(var_),           intent(out)   :: var
    type(grd_),           intent(out)   :: grd
    integer i, j, ip, isp, jsp, ich, jch, iz, nh
    real(dp) tmp
    character(len = 256) strm, fname


    ! grid setting
    grd%nx    = set%nx
    grd%ny    = set%ny
    grd%nz    = 101
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:grd%nz) = 2.0e3_dp ! [m]
    grd%alt(1)      = 0.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 17
    spl%nsp_i   = 13
    spl%nch     = 64

    ! allocate
    allocate(set%rate_from_datafile(grd%nz,0:spl%nch))
    allocate(var%ni(grd%nz,0:spl%nsp))
    allocate(var%ni_0(grd%nz,0:spl%nsp))
    allocate(var%ni_new(grd%nz,spl%nsp))
    allocate(var%ni_stable(grd%ny,grd%nz,spl%nsp))
    allocate(var%ni_3d(grd%nx,grd%ny,grd%nz,spl%nsp))
    allocate(var%clm_ni(grd%nz,spl%nsp))
    allocate(var%Ti(grd%nz),var%Te(grd%nz),var%Tn(grd%nz))
    allocate(var%Ti_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Te_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Tn_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%m(spl%nsp), var%q(spl%nsp))
    allocate(spl%reactant_list(spl%nch, 0:20))
    allocate(spl%product_list(spl%nch, 0:20))
    allocate(spl%species(0:spl%nsp))
    allocate(spl%label_fix(0:spl%nsp))
    allocate(spl%all_to_var(0:spl%nsp))
    allocate(spl%var_to_all(0:spl%nsp_i))
    allocate(spl%reaction_type_list(spl%nch))
    allocate(spl%reaction_type_char(spl%nch))
    allocate(var%ki(grd%nz,spl%nch))
    allocate(var%rate(grd%nz,spl%nch))
    allocate(var%Pi(grd%nz,spl%nsp_i))
    allocate(var%Pij(grd%nz,spl%nsp_i,spl%nch))
    allocate(var%Li(grd%nz,spl%nsp_i))
    allocate(var%Lij(grd%nz,spl%nsp_i,spl%nch))
    allocate(var%K_eddy(grd%nz),var%D_mol(grd%nz,spl%nsp),var%nu(grd%nz,spl%nsp))
    allocate(var%Fluxup(0:grd%nz+1,0:spl%nsp_i),var%Fluxdwn(0:grd%nz+1,0:spl%nsp_i))
    allocate(var%vFluxup(grd%nz+1,0:spl%nsp_i),var%vFluxdwn(grd%nz+1,0:spl%nsp_i))
    allocate(var%Upper_n(0:spl%nsp,2),var%Upper_f(0:spl%nsp,2),var%Upper_v(0:spl%nsp,2))
    allocate(var%Lower_n(0:spl%nsp,2),var%Lower_f(0:spl%nsp,2),var%Lower_v(0:spl%nsp,2))
    allocate(var%Jmtx(spl%nsp_i, spl%nsp_i))
    allocate(spl%rate_cases(spl%nch))
    allocate(spl%T_range(spl%nch,3,3))
    allocate(spl%major_species(grd%nz))
    allocate(var%n_tot(grd%nz),var%m_mean(grd%nz))
    allocate(var%Phip(grd%nz,spl%nsp_i))
    allocate(var%Phim(grd%nz,spl%nsp_i))
    allocate(var%dPhi_dz(grd%nz,spl%nsp_i))
    allocate(var%d_dniu_dPhi_dz(grd%nz,spl%nsp_i))
    allocate(var%d_dni0_dPhi_dz(grd%nz,spl%nsp_i))
    allocate(var%d_dnil_dPhi_dz(grd%nz,spl%nsp_i))
    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))
    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))
    allocate(var%tAmtx(2*spl%nsp_i+1,spl%nsp_i*grd%nz))
    allocate(var%tLmtx(spl%nsp_i+1,spl%nsp_i*grd%nz))
    allocate(var%Umtx(spl%nsp_i+1,spl%nsp_i*grd%nz))

    ! species list ----------------------------
    ! 
    ! CO2, CO, O, O(1D), H2O, H, OH, H2, O3, O2, HO2, H2O2, H2CO, HCO, M, HOCO, CO2+
    spl%species(1) = "CO2"
    spl%species(2) = "CO"
    spl%species(3) = "O"
    spl%species(4) = "O(1D)"
    spl%species(5) = "H2O"
    spl%species(6) = "H"
    spl%species(7) = "OH"
    spl%species(8) = "H2"
    spl%species(9) = "O3"
    spl%species(10) = "O2"
    spl%species(11) = "HO2"
    spl%species(12) = "H2O2"
    spl%species(13) = "H2CO"
    spl%species(14) = "HCO"
    spl%species(15) = "M"
    spl%species(16) = "HOCO"
    spl%species(17) = "CO2+"

    ! label_fix
    spl%label_fix(1) = 1
    spl%label_fix(2) = 0
    spl%label_fix(3) = 0
    spl%label_fix(4) = 0
    spl%label_fix(5) = 1
    spl%label_fix(6) = 0
    spl%label_fix(7) = 0
    spl%label_fix(8) = 0
    spl%label_fix(9) = 0
    spl%label_fix(10) = 0
    spl%label_fix(11) = 0
    spl%label_fix(12) = 0
    spl%label_fix(13) = 0
    spl%label_fix(14) = 0
    spl%label_fix(15) = 1
    spl%label_fix(16) = 0
    spl%label_fix(17) = 1

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(2) = 1
    spl%all_to_var(3) = 2
    spl%all_to_var(4) = 3
    spl%all_to_var(6) = 4
    spl%all_to_var(7) = 5
    spl%all_to_var(8) = 6
    spl%all_to_var(9) = 7
    spl%all_to_var(10) = 8
    spl%all_to_var(11) = 9
    spl%all_to_var(12) = 10
    spl%all_to_var(13) = 11
    spl%all_to_var(14) = 12
    spl%all_to_var(16) = 13

    ! var_to_all
    spl%var_to_all = 0
    spl%var_to_all(1) = 2
    spl%var_to_all(2) = 3
    spl%var_to_all(3) = 4
    spl%var_to_all(4) = 6
    spl%var_to_all(5) = 7
    spl%var_to_all(6) = 8
    spl%var_to_all(7) = 9
    spl%var_to_all(8) = 10
    spl%var_to_all(9) = 11
    spl%var_to_all(10) = 12
    spl%var_to_all(11) = 13
    spl%var_to_all(12) = 14
    spl%var_to_all(13) = 16

    ! Mass
    var%m(1) = 44.0_dp * cst%m_u ! CO2
    var%m(2) = 28.0_dp * cst%m_u ! CO
    var%m(3) = 16.0_dp * cst%m_u ! O
    var%m(4) = 16.0_dp * cst%m_u ! O(1D)
    var%m(5) = 18.0_dp * cst%m_u ! H2O
    var%m(6) = 1.0_dp * cst%m_u ! H
    var%m(7) = 17.0_dp * cst%m_u ! OH
    var%m(8) = 2.0_dp * cst%m_u ! H2
    var%m(9) = 48.0_dp * cst%m_u ! O3
    var%m(10) = 32.0_dp * cst%m_u ! O2
    var%m(11) = 33.0_dp * cst%m_u ! HO2
    var%m(12) = 34.0_dp * cst%m_u ! H2O2
    var%m(13) = 30.0_dp * cst%m_u ! H2CO
    var%m(14) = 29.0_dp * cst%m_u ! HCO
    var%m(15) = 1000000.0_dp * cst%m_u ! M
    var%m(16) = 45.0_dp * cst%m_u ! HOCO
    var%m(17) = 44.0_dp * cst%m_u ! CO2+

    ! Charge
    var%q(1) = 0.0_dp * cst%q_e ! CO2
    var%q(2) = 0.0_dp * cst%q_e ! CO
    var%q(3) = 0.0_dp * cst%q_e ! O
    var%q(4) = 0.0_dp * cst%q_e ! O(1D)
    var%q(5) = 0.0_dp * cst%q_e ! H2O
    var%q(6) = 0.0_dp * cst%q_e ! H
    var%q(7) = 0.0_dp * cst%q_e ! OH
    var%q(8) = 0.0_dp * cst%q_e ! H2
    var%q(9) = 0.0_dp * cst%q_e ! O3
    var%q(10) = 0.0_dp * cst%q_e ! O2
    var%q(11) = 0.0_dp * cst%q_e ! HO2
    var%q(12) = 0.0_dp * cst%q_e ! H2O2
    var%q(13) = 0.0_dp * cst%q_e ! H2CO
    var%q(14) = 0.0_dp * cst%q_e ! HCO
    var%q(15) = 0.0_dp * cst%q_e ! M
    var%q(16) = 0.0_dp * cst%q_e ! HOCO
    var%q(17) = 1.0_dp * cst%q_e ! CO2+

    ! reaction type characters
    var%nspecial = 0
    spl%reaction_type_list(1) = 2
    spl%reaction_type_char(1) = 'photodissociation'
    spl%reaction_type_list(2) = 2
    spl%reaction_type_char(2) = 'photodissociation'
    spl%reaction_type_list(3) = 2
    spl%reaction_type_char(3) = 'photodissociation'
    spl%reaction_type_list(4) = 2
    spl%reaction_type_char(4) = 'photodissociation'
    spl%reaction_type_list(5) = 2
    spl%reaction_type_char(5) = 'photodissociation'
    spl%reaction_type_list(6) = 2
    spl%reaction_type_char(6) = 'photodissociation'
    spl%reaction_type_list(7) = 2
    spl%reaction_type_char(7) = 'photodissociation'
    spl%reaction_type_list(8) = 2
    spl%reaction_type_char(8) = 'photodissociation'
    spl%reaction_type_list(9) = 2
    spl%reaction_type_char(9) = 'photodissociation'
    spl%reaction_type_list(10) = 2
    spl%reaction_type_char(10) = 'photodissociation'
    spl%reaction_type_list(11) = 2
    spl%reaction_type_char(11) = 'photodissociation'
    spl%reaction_type_list(12) = 2
    spl%reaction_type_char(12) = 'photodissociation'
    spl%reaction_type_list(13) = 2
    spl%reaction_type_char(13) = 'photodissociation'
    spl%reaction_type_list(14) = 2
    spl%reaction_type_char(14) = 'photodissociation'
    spl%reaction_type_list(15) = 2
    spl%reaction_type_char(15) = 'photodissociation'
    spl%reaction_type_list(16) = 2
    spl%reaction_type_char(16) = 'photodissociation'
    spl%reaction_type_list(36) = 30
    spl%reaction_type_char(36) = 'pressure_dependent_3body'
    spl%reaction_type_list(42) = 30
    spl%reaction_type_char(42) = 'pressure_dependent_3body'
    spl%reaction_type_list(49) = 31
    spl%reaction_type_char(49) = 'pressure_dependent_3bodyM'
    spl%reaction_type_list(50) = 30
    spl%reaction_type_char(50) = 'pressure_dependent_3body'
    spl%reaction_type_list(53) = 30
    spl%reaction_type_char(53) = 'pressure_dependent_3body'

    ! reactant list ----------------------------
    spl%reactant_list(1,0:1) = (/1,1/)
    spl%reactant_list(2,0:1) = (/1,1/)
    spl%reactant_list(3,0:1) = (/1,5/)
    spl%reactant_list(4,0:1) = (/1,5/)
    spl%reactant_list(5,0:1) = (/1,9/)
    spl%reactant_list(6,0:1) = (/1,9/)
    spl%reactant_list(7,0:1) = (/1,10/)
    spl%reactant_list(8,0:1) = (/1,10/)
    spl%reactant_list(9,0:1) = (/1,8/)
    spl%reactant_list(10,0:1) = (/1,7/)
    spl%reactant_list(11,0:1) = (/1,7/)
    spl%reactant_list(12,0:1) = (/1,11/)
    spl%reactant_list(13,0:1) = (/1,12/)
    spl%reactant_list(14,0:1) = (/1,12/)
    spl%reactant_list(15,0:1) = (/1,13/)
    spl%reactant_list(16,0:1) = (/1,13/)
    spl%reactant_list(17,0:3) = (/3,3,3,15/)
    spl%reactant_list(18,0:3) = (/3,3,10,1/)
    spl%reactant_list(19,0:2) = (/2,3,9/)
    spl%reactant_list(20,0:3) = (/3,3,2,15/)
    spl%reactant_list(21,0:2) = (/2,4,10/)
    spl%reactant_list(22,0:2) = (/2,4,9/)
    spl%reactant_list(23,0:2) = (/2,4,9/)
    spl%reactant_list(24,0:2) = (/2,4,8/)
    spl%reactant_list(25,0:2) = (/2,4,1/)
    spl%reactant_list(26,0:2) = (/2,4,5/)
    spl%reactant_list(27,0:2) = (/2,8,3/)
    spl%reactant_list(28,0:2) = (/2,7,8/)
    spl%reactant_list(29,0:3) = (/3,6,6,1/)
    spl%reactant_list(30,0:3) = (/3,6,7,1/)
    spl%reactant_list(31,0:2) = (/2,6,11/)
    spl%reactant_list(32,0:2) = (/2,6,11/)
    spl%reactant_list(33,0:2) = (/2,6,11/)
    spl%reactant_list(34,0:2) = (/2,6,12/)
    spl%reactant_list(35,0:2) = (/2,6,12/)
    spl%reactant_list(36,0:3) = (/3,6,10,15/)
    spl%reactant_list(37,0:2) = (/2,6,9/)
    spl%reactant_list(38,0:2) = (/2,3,7/)
    spl%reactant_list(39,0:2) = (/2,3,11/)
    spl%reactant_list(40,0:2) = (/2,3,12/)
    spl%reactant_list(41,0:2) = (/2,7,7/)
    spl%reactant_list(42,0:3) = (/3,7,7,15/)
    spl%reactant_list(43,0:2) = (/2,7,9/)
    spl%reactant_list(44,0:2) = (/2,7,11/)
    spl%reactant_list(45,0:2) = (/2,7,12/)
    spl%reactant_list(46,0:2) = (/2,11,9/)
    spl%reactant_list(47,0:2) = (/2,11,11/)
    spl%reactant_list(48,0:3) = (/3,11,11,15/)
    spl%reactant_list(49,0:3) = (/3,2,7,15/)
    spl%reactant_list(50,0:3) = (/3,2,7,15/)
    spl%reactant_list(51,0:2) = (/2,16,10/)
    spl%reactant_list(52,0:2) = (/2,17,8/)
    spl%reactant_list(53,0:3) = (/3,6,2,15/)
    spl%reactant_list(54,0:2) = (/2,6,14/)
    spl%reactant_list(55,0:2) = (/2,14,14/)
    spl%reactant_list(56,0:2) = (/2,7,14/)
    spl%reactant_list(57,0:2) = (/2,3,14/)
    spl%reactant_list(58,0:2) = (/2,3,14/)
    spl%reactant_list(59,0:2) = (/2,13,6/)
    spl%reactant_list(60,0:2) = (/2,14,10/)
    spl%reactant_list(61,0:2) = (/2,13,7/)
    spl%reactant_list(62,0:2) = (/2,13,3/)
    spl%reactant_list(63,0:2) = (/2,14,14/)
    spl%reactant_list(64,0:2) = (/2,11,14/)

    ! product list ----------------------------
    spl%product_list(1,0:2) = (/2,2,3/)
    spl%product_list(2,0:2) = (/2,2,4/)
    spl%product_list(3,0:2) = (/2,6,7/)
    spl%product_list(4,0:2) = (/2,8,4/)
    spl%product_list(5,0:2) = (/2,10,3/)
    spl%product_list(6,0:2) = (/2,10,4/)
    spl%product_list(7,0:2) = (/2,3,3/)
    spl%product_list(8,0:2) = (/2,3,4/)
    spl%product_list(9,0:2) = (/2,6,6/)
    spl%product_list(10,0:2) = (/2,3,6/)
    spl%product_list(11,0:2) = (/2,4,6/)
    spl%product_list(12,0:2) = (/2,7,3/)
    spl%product_list(13,0:2) = (/2,7,7/)
    spl%product_list(14,0:2) = (/2,11,6/)
    spl%product_list(15,0:2) = (/2,14,6/)
    spl%product_list(16,0:2) = (/2,2,8/)
    spl%product_list(17,0:2) = (/2,10,15/)
    spl%product_list(18,0:2) = (/2,9,1/)
    spl%product_list(19,0:2) = (/2,10,10/)
    spl%product_list(20,0:2) = (/2,1,15/)
    spl%product_list(21,0:2) = (/2,3,10/)
    spl%product_list(22,0:2) = (/2,10,10/)
    spl%product_list(23,0:3) = (/3,3,3,10/)
    spl%product_list(24,0:2) = (/2,6,7/)
    spl%product_list(25,0:2) = (/2,3,1/)
    spl%product_list(26,0:2) = (/2,7,7/)
    spl%product_list(27,0:2) = (/2,7,6/)
    spl%product_list(28,0:2) = (/2,5,6/)
    spl%product_list(29,0:2) = (/2,8,1/)
    spl%product_list(30,0:2) = (/2,5,1/)
    spl%product_list(31,0:2) = (/2,7,7/)
    spl%product_list(32,0:2) = (/2,5,4/)
    spl%product_list(33,0:2) = (/2,8,10/)
    spl%product_list(34,0:2) = (/2,11,8/)
    spl%product_list(35,0:2) = (/2,5,7/)
    spl%product_list(36,0:2) = (/2,11,15/)
    spl%product_list(37,0:2) = (/2,7,10/)
    spl%product_list(38,0:2) = (/2,10,6/)
    spl%product_list(39,0:2) = (/2,7,10/)
    spl%product_list(40,0:2) = (/2,7,11/)
    spl%product_list(41,0:2) = (/2,5,3/)
    spl%product_list(42,0:2) = (/2,12,15/)
    spl%product_list(43,0:2) = (/2,11,10/)
    spl%product_list(44,0:2) = (/2,5,10/)
    spl%product_list(45,0:2) = (/2,5,11/)
    spl%product_list(46,0:3) = (/3,7,10,10/)
    spl%product_list(47,0:2) = (/2,12,10/)
    spl%product_list(48,0:3) = (/3,12,10,15/)
    spl%product_list(49,0:3) = (/3,1,6,15/)
    spl%product_list(50,0:2) = (/2,16,15/)
    spl%product_list(51,0:2) = (/2,11,1/)
    spl%product_list(52,0:3) = (/3,1,6,6/)
    spl%product_list(53,0:2) = (/2,14,15/)
    spl%product_list(54,0:2) = (/2,8,2/)
    spl%product_list(55,0:2) = (/2,13,2/)
    spl%product_list(56,0:2) = (/2,5,2/)
    spl%product_list(57,0:2) = (/2,6,1/)
    spl%product_list(58,0:2) = (/2,7,2/)
    spl%product_list(59,0:2) = (/2,8,14/)
    spl%product_list(60,0:2) = (/2,11,2/)
    spl%product_list(61,0:2) = (/2,5,14/)
    spl%product_list(62,0:2) = (/2,7,14/)
    spl%product_list(63,0:3) = (/3,2,2,8/)
    spl%product_list(64,0:2) = (/2,12,2/)

    ! input Temperature profiles

    open(11, file = './Mars/Test/input/temperature/T_e.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Mars/Test/input/temperature/T_i.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Mars/Test/input/temperature/T_n.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-50_dp

    isp = sp_index(spl, 'CO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Test/input/density/CO2.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Test/input/density/H2O.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CO2+')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Test/input/density/CO2+.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    var%ni_0 = var%ni


    ! Reaction rate from external datafile in cgs unit



    ! Lower boundary condition
    var%Lower_n = 0.0_dp
    var%Lower_f = 0.0_dp
    var%Lower_v = 0.0_dp


    ! Upper boundary condition
    var%Upper_n = 0.0_dp
    var%Upper_f = 0.0_dp
    var%Upper_v = 0.0_dp

    isp = sp_index(spl, 'H')
    var%Upper_v(isp,1) = 10.0_dp ! Jeans escape
    var%Upper_v(isp,2) = 0.0_dp

    isp = sp_index(spl, 'H2')
    var%Upper_v(isp,1) = 10.0_dp ! Jeans escape
    var%Upper_v(isp,2) = 0.0_dp

    isp = sp_index(spl, 'O')
    var%Upper_f(isp,1) = 1.0_dp
    var%Upper_f(isp,2) = 1.2e12_dp


  end subroutine v__in__exe

end module v__in


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
