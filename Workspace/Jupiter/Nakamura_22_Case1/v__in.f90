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
    spl%planet = 'Jupiter'
    set%dir_name = './Jupiter/Nakamura_22_Case1'
    set%dir_name_win = '.\Jupiter\Nakamura_22_Case1'

    ! Calculation settings
    set%mode = '1D'
    set%use_1d = 0
    set%use_2d = 1
    set%F107 = 140.0_dp
    set%nstep = 10000
    set%fin_sec = 30000.0_dp
    set%dtime_limit = 10000.0_dp
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
    set%euv_input = 'HEUVAC'
    set%euv_factor = 1.0_dp
    set%n_wl_bin = 1
    allocate(set%wl_bin(set%n_wl_bin,3))
    set%wl_bin(1,1) = 0.0_dp
    set%wl_bin(1,2) = 1000.0_dp
    set%wl_bin(1,3) = 1.0_dp
    set%nsp_tout = 1
    allocate(set%species_tout(set%nsp_tout))
    set%species_tout(1) = ""
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
    grd%nz    = 141
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:grd%nz) = 20.0e3_dp ! [m]
    grd%alt(1)      = 200.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 52
    spl%nsp_i   = 35
    spl%nch     = 218

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
    ! H, H+, e-, H2, H2+, He, He+, CH4, CH4+, CH3+, C2H2, C2H2+, C2H4, C2H4+, C2H3+, C2H+, C2H6, C2H6+, C2H5+, HeH+, H3+, C+, C, CH+, CH2+, CH, CH2, CH3, CH5+, C2+, C2, C2H, C2H3, C2H5, C2H7+, C3H+, C3H2+, C3H3+, C3H4+, C3H5+, C3H6+, C3H7+, C3H8+, C3H9+, C4H+, C4H2+, C4H3+, C4H5+, C4H7+, C4H9+, H2(v2), H2(v4)
    spl%species(1) = "H"
    spl%species(2) = "H+"
    spl%species(3) = "e-"
    spl%species(4) = "H2"
    spl%species(5) = "H2+"
    spl%species(6) = "He"
    spl%species(7) = "He+"
    spl%species(8) = "CH4"
    spl%species(9) = "CH4+"
    spl%species(10) = "CH3+"
    spl%species(11) = "C2H2"
    spl%species(12) = "C2H2+"
    spl%species(13) = "C2H4"
    spl%species(14) = "C2H4+"
    spl%species(15) = "C2H3+"
    spl%species(16) = "C2H+"
    spl%species(17) = "C2H6"
    spl%species(18) = "C2H6+"
    spl%species(19) = "C2H5+"
    spl%species(20) = "HeH+"
    spl%species(21) = "H3+"
    spl%species(22) = "C+"
    spl%species(23) = "C"
    spl%species(24) = "CH+"
    spl%species(25) = "CH2+"
    spl%species(26) = "CH"
    spl%species(27) = "CH2"
    spl%species(28) = "CH3"
    spl%species(29) = "CH5+"
    spl%species(30) = "C2+"
    spl%species(31) = "C2"
    spl%species(32) = "C2H"
    spl%species(33) = "C2H3"
    spl%species(34) = "C2H5"
    spl%species(35) = "C2H7+"
    spl%species(36) = "C3H+"
    spl%species(37) = "C3H2+"
    spl%species(38) = "C3H3+"
    spl%species(39) = "C3H4+"
    spl%species(40) = "C3H5+"
    spl%species(41) = "C3H6+"
    spl%species(42) = "C3H7+"
    spl%species(43) = "C3H8+"
    spl%species(44) = "C3H9+"
    spl%species(45) = "C4H+"
    spl%species(46) = "C4H2+"
    spl%species(47) = "C4H3+"
    spl%species(48) = "C4H5+"
    spl%species(49) = "C4H7+"
    spl%species(50) = "C4H9+"
    spl%species(51) = "H2(v2)"
    spl%species(52) = "H2(v4)"

    ! label_fix
    spl%label_fix(1) = 0
    spl%label_fix(2) = 0
    spl%label_fix(3) = 1
    spl%label_fix(4) = 1
    spl%label_fix(5) = 0
    spl%label_fix(6) = 1
    spl%label_fix(7) = 0
    spl%label_fix(8) = 1
    spl%label_fix(9) = 0
    spl%label_fix(10) = 0
    spl%label_fix(11) = 1
    spl%label_fix(12) = 0
    spl%label_fix(13) = 1
    spl%label_fix(14) = 0
    spl%label_fix(15) = 0
    spl%label_fix(16) = 0
    spl%label_fix(17) = 1
    spl%label_fix(18) = 0
    spl%label_fix(19) = 0
    spl%label_fix(20) = 0
    spl%label_fix(21) = 0
    spl%label_fix(22) = 0
    spl%label_fix(23) = 1
    spl%label_fix(24) = 0
    spl%label_fix(25) = 0
    spl%label_fix(26) = 1
    spl%label_fix(27) = 1
    spl%label_fix(28) = 1
    spl%label_fix(29) = 0
    spl%label_fix(30) = 0
    spl%label_fix(31) = 1
    spl%label_fix(32) = 1
    spl%label_fix(33) = 1
    spl%label_fix(34) = 1
    spl%label_fix(35) = 0
    spl%label_fix(36) = 0
    spl%label_fix(37) = 0
    spl%label_fix(38) = 0
    spl%label_fix(39) = 0
    spl%label_fix(40) = 0
    spl%label_fix(41) = 0
    spl%label_fix(42) = 0
    spl%label_fix(43) = 0
    spl%label_fix(44) = 0
    spl%label_fix(45) = 0
    spl%label_fix(46) = 0
    spl%label_fix(47) = 0
    spl%label_fix(48) = 0
    spl%label_fix(49) = 0
    spl%label_fix(50) = 0
    spl%label_fix(51) = 1
    spl%label_fix(52) = 1

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1
    spl%all_to_var(2) = 2
    spl%all_to_var(5) = 3
    spl%all_to_var(7) = 4
    spl%all_to_var(9) = 5
    spl%all_to_var(10) = 6
    spl%all_to_var(12) = 7
    spl%all_to_var(14) = 8
    spl%all_to_var(15) = 9
    spl%all_to_var(16) = 10
    spl%all_to_var(18) = 11
    spl%all_to_var(19) = 12
    spl%all_to_var(20) = 13
    spl%all_to_var(21) = 14
    spl%all_to_var(22) = 15
    spl%all_to_var(24) = 16
    spl%all_to_var(25) = 17
    spl%all_to_var(29) = 18
    spl%all_to_var(30) = 19
    spl%all_to_var(35) = 20
    spl%all_to_var(36) = 21
    spl%all_to_var(37) = 22
    spl%all_to_var(38) = 23
    spl%all_to_var(39) = 24
    spl%all_to_var(40) = 25
    spl%all_to_var(41) = 26
    spl%all_to_var(42) = 27
    spl%all_to_var(43) = 28
    spl%all_to_var(44) = 29
    spl%all_to_var(45) = 30
    spl%all_to_var(46) = 31
    spl%all_to_var(47) = 32
    spl%all_to_var(48) = 33
    spl%all_to_var(49) = 34
    spl%all_to_var(50) = 35

    ! var_to_all
    spl%var_to_all = 0
    spl%var_to_all(1) = 1
    spl%var_to_all(2) = 2
    spl%var_to_all(3) = 5
    spl%var_to_all(4) = 7
    spl%var_to_all(5) = 9
    spl%var_to_all(6) = 10
    spl%var_to_all(7) = 12
    spl%var_to_all(8) = 14
    spl%var_to_all(9) = 15
    spl%var_to_all(10) = 16
    spl%var_to_all(11) = 18
    spl%var_to_all(12) = 19
    spl%var_to_all(13) = 20
    spl%var_to_all(14) = 21
    spl%var_to_all(15) = 22
    spl%var_to_all(16) = 24
    spl%var_to_all(17) = 25
    spl%var_to_all(18) = 29
    spl%var_to_all(19) = 30
    spl%var_to_all(20) = 35
    spl%var_to_all(21) = 36
    spl%var_to_all(22) = 37
    spl%var_to_all(23) = 38
    spl%var_to_all(24) = 39
    spl%var_to_all(25) = 40
    spl%var_to_all(26) = 41
    spl%var_to_all(27) = 42
    spl%var_to_all(28) = 43
    spl%var_to_all(29) = 44
    spl%var_to_all(30) = 45
    spl%var_to_all(31) = 46
    spl%var_to_all(32) = 47
    spl%var_to_all(33) = 48
    spl%var_to_all(34) = 49
    spl%var_to_all(35) = 50

    ! Mass
    var%m(1) = 1.0_dp * cst%m_u ! H
    var%m(2) = 1.0_dp * cst%m_u ! H+
    var%m(3) = 0.00054858_dp * cst%m_u ! e-
    var%m(4) = 2.0_dp * cst%m_u ! H2
    var%m(5) = 2.0_dp * cst%m_u ! H2+
    var%m(6) = 4.0_dp * cst%m_u ! He
    var%m(7) = 4.0_dp * cst%m_u ! He+
    var%m(8) = 16.0_dp * cst%m_u ! CH4
    var%m(9) = 16.0_dp * cst%m_u ! CH4+
    var%m(10) = 15.0_dp * cst%m_u ! CH3+
    var%m(11) = 26.0_dp * cst%m_u ! C2H2
    var%m(12) = 26.0_dp * cst%m_u ! C2H2+
    var%m(13) = 28.0_dp * cst%m_u ! C2H4
    var%m(14) = 28.0_dp * cst%m_u ! C2H4+
    var%m(15) = 27.0_dp * cst%m_u ! C2H3+
    var%m(16) = 25.0_dp * cst%m_u ! C2H+
    var%m(17) = 30.0_dp * cst%m_u ! C2H6
    var%m(18) = 30.0_dp * cst%m_u ! C2H6+
    var%m(19) = 29.0_dp * cst%m_u ! C2H5+
    var%m(20) = 5.0_dp * cst%m_u ! HeH+
    var%m(21) = 3.0_dp * cst%m_u ! H3+
    var%m(22) = 12.0_dp * cst%m_u ! C+
    var%m(23) = 12.0_dp * cst%m_u ! C
    var%m(24) = 13.0_dp * cst%m_u ! CH+
    var%m(25) = 14.0_dp * cst%m_u ! CH2+
    var%m(26) = 13.0_dp * cst%m_u ! CH
    var%m(27) = 14.0_dp * cst%m_u ! CH2
    var%m(28) = 15.0_dp * cst%m_u ! CH3
    var%m(29) = 17.0_dp * cst%m_u ! CH5+
    var%m(30) = 24.0_dp * cst%m_u ! C2+
    var%m(31) = 24.0_dp * cst%m_u ! C2
    var%m(32) = 25.0_dp * cst%m_u ! C2H
    var%m(33) = 27.0_dp * cst%m_u ! C2H3
    var%m(34) = 29.0_dp * cst%m_u ! C2H5
    var%m(35) = 31.0_dp * cst%m_u ! C2H7+
    var%m(36) = 37.0_dp * cst%m_u ! C3H+
    var%m(37) = 38.0_dp * cst%m_u ! C3H2+
    var%m(38) = 39.0_dp * cst%m_u ! C3H3+
    var%m(39) = 40.0_dp * cst%m_u ! C3H4+
    var%m(40) = 41.0_dp * cst%m_u ! C3H5+
    var%m(41) = 42.0_dp * cst%m_u ! C3H6+
    var%m(42) = 43.0_dp * cst%m_u ! C3H7+
    var%m(43) = 44.0_dp * cst%m_u ! C3H8+
    var%m(44) = 45.0_dp * cst%m_u ! C3H9+
    var%m(45) = 49.0_dp * cst%m_u ! C4H+
    var%m(46) = 50.0_dp * cst%m_u ! C4H2+
    var%m(47) = 51.0_dp * cst%m_u ! C4H3+
    var%m(48) = 53.0_dp * cst%m_u ! C4H5+
    var%m(49) = 55.0_dp * cst%m_u ! C4H7+
    var%m(50) = 57.0_dp * cst%m_u ! C4H9+
    var%m(51) = 2.0_dp * cst%m_u ! H2(v2)
    var%m(52) = 2.0_dp * cst%m_u ! H2(v4)

    ! Charge
    var%q(1) = 0.0_dp * cst%q_e ! H
    var%q(2) = 1.0_dp * cst%q_e ! H+
    var%q(3) = -1.0_dp * cst%q_e ! e-
    var%q(4) = 0.0_dp * cst%q_e ! H2
    var%q(5) = 1.0_dp * cst%q_e ! H2+
    var%q(6) = 0.0_dp * cst%q_e ! He
    var%q(7) = 1.0_dp * cst%q_e ! He+
    var%q(8) = 0.0_dp * cst%q_e ! CH4
    var%q(9) = 1.0_dp * cst%q_e ! CH4+
    var%q(10) = 1.0_dp * cst%q_e ! CH3+
    var%q(11) = 0.0_dp * cst%q_e ! C2H2
    var%q(12) = 1.0_dp * cst%q_e ! C2H2+
    var%q(13) = 0.0_dp * cst%q_e ! C2H4
    var%q(14) = 1.0_dp * cst%q_e ! C2H4+
    var%q(15) = 1.0_dp * cst%q_e ! C2H3+
    var%q(16) = 1.0_dp * cst%q_e ! C2H+
    var%q(17) = 0.0_dp * cst%q_e ! C2H6
    var%q(18) = 1.0_dp * cst%q_e ! C2H6+
    var%q(19) = 1.0_dp * cst%q_e ! C2H5+
    var%q(20) = 1.0_dp * cst%q_e ! HeH+
    var%q(21) = 1.0_dp * cst%q_e ! H3+
    var%q(22) = 1.0_dp * cst%q_e ! C+
    var%q(23) = 0.0_dp * cst%q_e ! C
    var%q(24) = 1.0_dp * cst%q_e ! CH+
    var%q(25) = 1.0_dp * cst%q_e ! CH2+
    var%q(26) = 0.0_dp * cst%q_e ! CH
    var%q(27) = 0.0_dp * cst%q_e ! CH2
    var%q(28) = 0.0_dp * cst%q_e ! CH3
    var%q(29) = 1.0_dp * cst%q_e ! CH5+
    var%q(30) = 1.0_dp * cst%q_e ! C2+
    var%q(31) = 0.0_dp * cst%q_e ! C2
    var%q(32) = 0.0_dp * cst%q_e ! C2H
    var%q(33) = 0.0_dp * cst%q_e ! C2H3
    var%q(34) = 0.0_dp * cst%q_e ! C2H5
    var%q(35) = 1.0_dp * cst%q_e ! C2H7+
    var%q(36) = 1.0_dp * cst%q_e ! C3H+
    var%q(37) = 1.0_dp * cst%q_e ! C3H2+
    var%q(38) = 1.0_dp * cst%q_e ! C3H3+
    var%q(39) = 1.0_dp * cst%q_e ! C3H4+
    var%q(40) = 1.0_dp * cst%q_e ! C3H5+
    var%q(41) = 1.0_dp * cst%q_e ! C3H6+
    var%q(42) = 1.0_dp * cst%q_e ! C3H7+
    var%q(43) = 1.0_dp * cst%q_e ! C3H8+
    var%q(44) = 1.0_dp * cst%q_e ! C3H9+
    var%q(45) = 1.0_dp * cst%q_e ! C4H+
    var%q(46) = 1.0_dp * cst%q_e ! C4H2+
    var%q(47) = 1.0_dp * cst%q_e ! C4H3+
    var%q(48) = 1.0_dp * cst%q_e ! C4H5+
    var%q(49) = 1.0_dp * cst%q_e ! C4H7+
    var%q(50) = 1.0_dp * cst%q_e ! C4H9+
    var%q(51) = 0.0_dp * cst%q_e ! H2(v2)
    var%q(52) = 0.0_dp * cst%q_e ! H2(v4)

    ! reaction type characters
    var%nspecial = 0
    spl%reaction_type_list(1) = 1
    spl%reaction_type_char(1) = 'photoionization'
    spl%reaction_type_list(2) = 1
    spl%reaction_type_char(2) = 'photoionization'
    spl%reaction_type_list(3) = 1
    spl%reaction_type_char(3) = 'photoionization'
    spl%reaction_type_list(4) = 1
    spl%reaction_type_char(4) = 'photoionization'
    spl%reaction_type_list(5) = 1
    spl%reaction_type_char(5) = 'photoionization'
    spl%reaction_type_list(6) = 1
    spl%reaction_type_char(6) = 'photoionization'
    spl%reaction_type_list(7) = 1
    spl%reaction_type_char(7) = 'photoionization'
    spl%reaction_type_list(8) = 1
    spl%reaction_type_char(8) = 'photoionization'
    spl%reaction_type_list(9) = 1
    spl%reaction_type_char(9) = 'photoionization'
    spl%reaction_type_list(10) = 1
    spl%reaction_type_char(10) = 'photoionization'
    spl%reaction_type_list(11) = 1
    spl%reaction_type_char(11) = 'photoionization'
    spl%reaction_type_list(12) = 1
    spl%reaction_type_char(12) = 'photoionization'
    spl%reaction_type_list(13) = 1
    spl%reaction_type_char(13) = 'photoionization'
    spl%reaction_type_list(14) = 1
    spl%reaction_type_char(14) = 'photoionization'
    spl%reaction_type_list(15) = 1
    spl%reaction_type_char(15) = 'photoionization'
    spl%reaction_type_list(16) = 1
    spl%reaction_type_char(16) = 'photoionization'
    spl%reaction_type_list(17) = 1
    spl%reaction_type_char(17) = 'photoionization'

    ! reactant list ----------------------------
    spl%reactant_list(1,0:1) = (/1,1/)
    spl%reactant_list(2,0:1) = (/1,4/)
    spl%reactant_list(3,0:1) = (/1,4/)
    spl%reactant_list(4,0:1) = (/1,6/)
    spl%reactant_list(5,0:1) = (/1,8/)
    spl%reactant_list(6,0:1) = (/1,8/)
    spl%reactant_list(7,0:1) = (/1,8/)
    spl%reactant_list(8,0:1) = (/1,11/)
    spl%reactant_list(9,0:1) = (/1,13/)
    spl%reactant_list(10,0:1) = (/1,13/)
    spl%reactant_list(11,0:1) = (/1,13/)
    spl%reactant_list(12,0:1) = (/1,13/)
    spl%reactant_list(13,0:1) = (/1,17/)
    spl%reactant_list(14,0:1) = (/1,17/)
    spl%reactant_list(15,0:1) = (/1,17/)
    spl%reactant_list(16,0:1) = (/1,17/)
    spl%reactant_list(17,0:1) = (/1,17/)
    spl%reactant_list(18,0:2) = (/2,2,3/)
    spl%reactant_list(19,0:2) = (/2,7,3/)
    spl%reactant_list(20,0:2) = (/2,20,3/)
    spl%reactant_list(21,0:2) = (/2,5,3/)
    spl%reactant_list(22,0:2) = (/2,21,3/)
    spl%reactant_list(23,0:2) = (/2,21,3/)
    spl%reactant_list(24,0:2) = (/2,22,3/)
    spl%reactant_list(25,0:2) = (/2,24,3/)
    spl%reactant_list(26,0:2) = (/2,25,3/)
    spl%reactant_list(27,0:2) = (/2,10,3/)
    spl%reactant_list(28,0:2) = (/2,9,3/)
    spl%reactant_list(29,0:2) = (/2,9,3/)
    spl%reactant_list(30,0:2) = (/2,29,3/)
    spl%reactant_list(31,0:2) = (/2,29,3/)
    spl%reactant_list(32,0:2) = (/2,30,3/)
    spl%reactant_list(33,0:2) = (/2,16,3/)
    spl%reactant_list(34,0:2) = (/2,16,3/)
    spl%reactant_list(35,0:2) = (/2,12,3/)
    spl%reactant_list(36,0:2) = (/2,12,3/)
    spl%reactant_list(37,0:2) = (/2,15,3/)
    spl%reactant_list(38,0:2) = (/2,15,3/)
    spl%reactant_list(39,0:2) = (/2,14,3/)
    spl%reactant_list(40,0:2) = (/2,14,3/)
    spl%reactant_list(41,0:2) = (/2,19,3/)
    spl%reactant_list(42,0:2) = (/2,19,3/)
    spl%reactant_list(43,0:2) = (/2,18,3/)
    spl%reactant_list(44,0:2) = (/2,18,3/)
    spl%reactant_list(45,0:2) = (/2,35,3/)
    spl%reactant_list(46,0:2) = (/2,36,3/)
    spl%reactant_list(47,0:2) = (/2,37,3/)
    spl%reactant_list(48,0:2) = (/2,38,3/)
    spl%reactant_list(49,0:2) = (/2,39,3/)
    spl%reactant_list(50,0:2) = (/2,40,3/)
    spl%reactant_list(51,0:2) = (/2,41,3/)
    spl%reactant_list(52,0:2) = (/2,42,3/)
    spl%reactant_list(53,0:2) = (/2,43,3/)
    spl%reactant_list(54,0:2) = (/2,44,3/)
    spl%reactant_list(55,0:2) = (/2,45,3/)
    spl%reactant_list(56,0:2) = (/2,46,3/)
    spl%reactant_list(57,0:2) = (/2,47,3/)
    spl%reactant_list(58,0:2) = (/2,48,3/)
    spl%reactant_list(59,0:2) = (/2,49,3/)
    spl%reactant_list(60,0:2) = (/2,50,3/)
    spl%reactant_list(61,0:2) = (/2,5,4/)
    spl%reactant_list(62,0:2) = (/2,5,1/)
    spl%reactant_list(63,0:2) = (/2,5,6/)
    spl%reactant_list(64,0:2) = (/2,5,8/)
    spl%reactant_list(65,0:2) = (/2,5,8/)
    spl%reactant_list(66,0:2) = (/2,5,8/)
    spl%reactant_list(67,0:2) = (/2,5,11/)
    spl%reactant_list(68,0:2) = (/2,5,11/)
    spl%reactant_list(69,0:2) = (/2,5,17/)
    spl%reactant_list(70,0:2) = (/2,5,17/)
    spl%reactant_list(71,0:2) = (/2,5,17/)
    spl%reactant_list(72,0:2) = (/2,5,17/)
    spl%reactant_list(73,0:2) = (/2,5,17/)
    spl%reactant_list(74,0:2) = (/2,7,4/)
    spl%reactant_list(75,0:2) = (/2,7,4/)
    spl%reactant_list(76,0:2) = (/2,7,51/)
    spl%reactant_list(77,0:2) = (/2,7,8/)
    spl%reactant_list(78,0:2) = (/2,7,8/)
    spl%reactant_list(79,0:2) = (/2,7,8/)
    spl%reactant_list(80,0:2) = (/2,7,8/)
    spl%reactant_list(81,0:2) = (/2,7,8/)
    spl%reactant_list(82,0:2) = (/2,7,11/)
    spl%reactant_list(83,0:2) = (/2,7,11/)
    spl%reactant_list(84,0:2) = (/2,7,11/)
    spl%reactant_list(85,0:2) = (/2,7,11/)
    spl%reactant_list(86,0:2) = (/2,7,13/)
    spl%reactant_list(87,0:2) = (/2,7,13/)
    spl%reactant_list(88,0:2) = (/2,7,13/)
    spl%reactant_list(89,0:2) = (/2,7,13/)
    spl%reactant_list(90,0:2) = (/2,7,13/)
    spl%reactant_list(91,0:2) = (/2,7,17/)
    spl%reactant_list(92,0:2) = (/2,7,17/)
    spl%reactant_list(93,0:2) = (/2,7,17/)
    spl%reactant_list(94,0:2) = (/2,2,8/)
    spl%reactant_list(95,0:2) = (/2,2,8/)
    spl%reactant_list(96,0:2) = (/2,2,17/)
    spl%reactant_list(97,0:2) = (/2,2,17/)
    spl%reactant_list(98,0:2) = (/2,2,17/)
    spl%reactant_list(99,0:2) = (/2,2,52/)
    spl%reactant_list(100,0:2) = (/2,20,4/)
    spl%reactant_list(101,0:2) = (/2,20,1/)
    spl%reactant_list(102,0:2) = (/2,20,13/)
    spl%reactant_list(103,0:2) = (/2,20,13/)
    spl%reactant_list(104,0:2) = (/2,20,17/)
    spl%reactant_list(105,0:2) = (/2,20,17/)
    spl%reactant_list(106,0:2) = (/2,21,8/)
    spl%reactant_list(107,0:2) = (/2,21,11/)
    spl%reactant_list(108,0:2) = (/2,21,13/)
    spl%reactant_list(109,0:2) = (/2,21,13/)
    spl%reactant_list(110,0:2) = (/2,21,17/)
    spl%reactant_list(111,0:2) = (/2,22,8/)
    spl%reactant_list(112,0:2) = (/2,22,8/)
    spl%reactant_list(113,0:2) = (/2,22,11/)
    spl%reactant_list(114,0:2) = (/2,22,13/)
    spl%reactant_list(115,0:2) = (/2,22,13/)
    spl%reactant_list(116,0:2) = (/2,22,13/)
    spl%reactant_list(117,0:2) = (/2,22,13/)
    spl%reactant_list(118,0:2) = (/2,22,13/)
    spl%reactant_list(119,0:2) = (/2,22,17/)
    spl%reactant_list(120,0:2) = (/2,22,17/)
    spl%reactant_list(121,0:2) = (/2,22,17/)
    spl%reactant_list(122,0:2) = (/2,22,17/)
    spl%reactant_list(123,0:2) = (/2,24,1/)
    spl%reactant_list(124,0:2) = (/2,24,4/)
    spl%reactant_list(125,0:2) = (/2,24,8/)
    spl%reactant_list(126,0:2) = (/2,24,8/)
    spl%reactant_list(127,0:2) = (/2,24,8/)
    spl%reactant_list(128,0:2) = (/2,24,11/)
    spl%reactant_list(129,0:2) = (/2,25,4/)
    spl%reactant_list(130,0:2) = (/2,25,8/)
    spl%reactant_list(131,0:2) = (/2,25,8/)
    spl%reactant_list(132,0:2) = (/2,25,11/)
    spl%reactant_list(133,0:2) = (/2,10,8/)
    spl%reactant_list(134,0:2) = (/2,10,11/)
    spl%reactant_list(135,0:2) = (/2,10,13/)
    spl%reactant_list(136,0:2) = (/2,10,13/)
    spl%reactant_list(137,0:2) = (/2,10,13/)
    spl%reactant_list(138,0:2) = (/2,10,17/)
    spl%reactant_list(139,0:2) = (/2,10,17/)
    spl%reactant_list(140,0:2) = (/2,10,17/)
    spl%reactant_list(141,0:2) = (/2,9,4/)
    spl%reactant_list(142,0:2) = (/2,9,8/)
    spl%reactant_list(143,0:2) = (/2,9,11/)
    spl%reactant_list(144,0:2) = (/2,9,11/)
    spl%reactant_list(145,0:2) = (/2,9,11/)
    spl%reactant_list(146,0:2) = (/2,9,13/)
    spl%reactant_list(147,0:2) = (/2,9,13/)
    spl%reactant_list(148,0:2) = (/2,9,13/)
    spl%reactant_list(149,0:2) = (/2,9,17/)
    spl%reactant_list(150,0:2) = (/2,29,1/)
    spl%reactant_list(151,0:2) = (/2,29,11/)
    spl%reactant_list(152,0:2) = (/2,29,13/)
    spl%reactant_list(153,0:2) = (/2,29,17/)
    spl%reactant_list(154,0:2) = (/2,29,17/)
    spl%reactant_list(155,0:2) = (/2,30,4/)
    spl%reactant_list(156,0:2) = (/2,30,8/)
    spl%reactant_list(157,0:2) = (/2,30,8/)
    spl%reactant_list(158,0:2) = (/2,30,8/)
    spl%reactant_list(159,0:2) = (/2,30,8/)
    spl%reactant_list(160,0:2) = (/2,30,8/)
    spl%reactant_list(161,0:2) = (/2,30,11/)
    spl%reactant_list(162,0:2) = (/2,30,13/)
    spl%reactant_list(163,0:2) = (/2,16,4/)
    spl%reactant_list(164,0:2) = (/2,16,8/)
    spl%reactant_list(165,0:2) = (/2,16,8/)
    spl%reactant_list(166,0:2) = (/2,16,8/)
    spl%reactant_list(167,0:2) = (/2,16,8/)
    spl%reactant_list(168,0:2) = (/2,16,11/)
    spl%reactant_list(169,0:2) = (/2,16,13/)
    spl%reactant_list(170,0:2) = (/2,12,4/)
    spl%reactant_list(171,0:2) = (/2,12,8/)
    spl%reactant_list(172,0:2) = (/2,12,8/)
    spl%reactant_list(173,0:2) = (/2,12,11/)
    spl%reactant_list(174,0:2) = (/2,12,11/)
    spl%reactant_list(175,0:2) = (/2,12,13/)
    spl%reactant_list(176,0:2) = (/2,12,13/)
    spl%reactant_list(177,0:2) = (/2,12,13/)
    spl%reactant_list(178,0:2) = (/2,12,17/)
    spl%reactant_list(179,0:2) = (/2,12,17/)
    spl%reactant_list(180,0:2) = (/2,12,17/)
    spl%reactant_list(181,0:2) = (/2,12,17/)
    spl%reactant_list(182,0:2) = (/2,12,17/)
    spl%reactant_list(183,0:2) = (/2,12,17/)
    spl%reactant_list(184,0:2) = (/2,12,17/)
    spl%reactant_list(185,0:2) = (/2,15,1/)
    spl%reactant_list(186,0:2) = (/2,15,8/)
    spl%reactant_list(187,0:2) = (/2,15,11/)
    spl%reactant_list(188,0:2) = (/2,15,13/)
    spl%reactant_list(189,0:2) = (/2,15,17/)
    spl%reactant_list(190,0:2) = (/2,15,17/)
    spl%reactant_list(191,0:2) = (/2,15,17/)
    spl%reactant_list(192,0:2) = (/2,14,1/)
    spl%reactant_list(193,0:2) = (/2,14,11/)
    spl%reactant_list(194,0:2) = (/2,14,11/)
    spl%reactant_list(195,0:2) = (/2,14,13/)
    spl%reactant_list(196,0:2) = (/2,14,13/)
    spl%reactant_list(197,0:2) = (/2,14,17/)
    spl%reactant_list(198,0:2) = (/2,14,17/)
    spl%reactant_list(199,0:2) = (/2,19,1/)
    spl%reactant_list(200,0:2) = (/2,19,8/)
    spl%reactant_list(201,0:2) = (/2,19,11/)
    spl%reactant_list(202,0:2) = (/2,19,11/)
    spl%reactant_list(203,0:2) = (/2,19,13/)
    spl%reactant_list(204,0:2) = (/2,19,17/)
    spl%reactant_list(205,0:2) = (/2,18,1/)
    spl%reactant_list(206,0:2) = (/2,18,11/)
    spl%reactant_list(207,0:2) = (/2,18,11/)
    spl%reactant_list(208,0:2) = (/2,18,11/)
    spl%reactant_list(209,0:2) = (/2,18,13/)
    spl%reactant_list(210,0:2) = (/2,18,17/)
    spl%reactant_list(211,0:2) = (/2,18,17/)
    spl%reactant_list(212,0:2) = (/2,35,11/)
    spl%reactant_list(213,0:2) = (/2,35,13/)
    spl%reactant_list(214,0:3) = (/3,2,4,4/)
    spl%reactant_list(215,0:3) = (/3,10,4,4/)
    spl%reactant_list(216,0:3) = (/3,10,4,6/)
    spl%reactant_list(217,0:3) = (/3,12,4,6/)
    spl%reactant_list(218,0:3) = (/3,1,1,4/)

    ! product list ----------------------------
    spl%product_list(1,0:2) = (/2,2,3/)
    spl%product_list(2,0:2) = (/2,5,3/)
    spl%product_list(3,0:3) = (/3,2,3,1/)
    spl%product_list(4,0:2) = (/2,7,3/)
    spl%product_list(5,0:2) = (/2,9,3/)
    spl%product_list(6,0:3) = (/3,10,3,1/)
    spl%product_list(7,0:2) = (/2,5,3/)
    spl%product_list(8,0:2) = (/2,12,3/)
    spl%product_list(9,0:2) = (/2,14,3/)
    spl%product_list(10,0:3) = (/3,15,3,1/)
    spl%product_list(11,0:2) = (/2,12,3/)
    spl%product_list(12,0:2) = (/2,16,3/)
    spl%product_list(13,0:2) = (/2,18,3/)
    spl%product_list(14,0:3) = (/3,19,3,1/)
    spl%product_list(15,0:2) = (/2,14,3/)
    spl%product_list(16,0:2) = (/2,15,3/)
    spl%product_list(17,0:2) = (/2,12,3/)
    spl%product_list(18,0:1) = (/1,1/)
    spl%product_list(19,0:1) = (/1,6/)
    spl%product_list(20,0:2) = (/2,1,6/)
    spl%product_list(21,0:2) = (/2,1,1/)
    spl%product_list(22,0:2) = (/2,4,1/)
    spl%product_list(23,0:3) = (/3,1,1,1/)
    spl%product_list(24,0:1) = (/1,23/)
    spl%product_list(25,0:2) = (/2,23,1/)
    spl%product_list(26,0:2) = (/2,26,1/)
    spl%product_list(27,0:2) = (/2,27,1/)
    spl%product_list(28,0:2) = (/2,28,1/)
    spl%product_list(29,0:3) = (/3,27,1,1/)
    spl%product_list(30,0:3) = (/3,27,1,4/)
    spl%product_list(31,0:3) = (/3,28,1,1/)
    spl%product_list(32,0:2) = (/2,23,23/)
    spl%product_list(33,0:2) = (/2,31,1/)
    spl%product_list(34,0:2) = (/2,26,23/)
    spl%product_list(35,0:2) = (/2,32,1/)
    spl%product_list(36,0:2) = (/2,26,26/)
    spl%product_list(37,0:2) = (/2,11,1/)
    spl%product_list(38,0:2) = (/2,27,26/)
    spl%product_list(39,0:2) = (/2,33,1/)
    spl%product_list(40,0:2) = (/2,27,27/)
    spl%product_list(41,0:2) = (/2,13,1/)
    spl%product_list(42,0:2) = (/2,28,27/)
    spl%product_list(43,0:2) = (/2,34,1/)
    spl%product_list(44,0:2) = (/2,28,28/)
    spl%product_list(45,0:2) = (/2,17,1/)
    spl%product_list(46,0:1) = (/1,0/)
    spl%product_list(47,0:1) = (/1,0/)
    spl%product_list(48,0:1) = (/1,0/)
    spl%product_list(49,0:1) = (/1,0/)
    spl%product_list(50,0:1) = (/1,0/)
    spl%product_list(51,0:1) = (/1,0/)
    spl%product_list(52,0:1) = (/1,0/)
    spl%product_list(53,0:1) = (/1,0/)
    spl%product_list(54,0:1) = (/1,0/)
    spl%product_list(55,0:1) = (/1,0/)
    spl%product_list(56,0:1) = (/1,0/)
    spl%product_list(57,0:1) = (/1,0/)
    spl%product_list(58,0:1) = (/1,0/)
    spl%product_list(59,0:1) = (/1,0/)
    spl%product_list(60,0:1) = (/1,0/)
    spl%product_list(61,0:2) = (/2,21,1/)
    spl%product_list(62,0:2) = (/2,2,4/)
    spl%product_list(63,0:2) = (/2,20,1/)
    spl%product_list(64,0:2) = (/2,29,1/)
    spl%product_list(65,0:2) = (/2,9,4/)
    spl%product_list(66,0:3) = (/3,10,1,4/)
    spl%product_list(67,0:2) = (/2,15,1/)
    spl%product_list(68,0:2) = (/2,12,4/)
    spl%product_list(69,0:2) = (/2,18,4/)
    spl%product_list(70,0:3) = (/3,19,1,4/)
    spl%product_list(71,0:3) = (/3,14,4,4/)
    spl%product_list(72,0:4) = (/4,15,1,4,4/)
    spl%product_list(73,0:4) = (/4,12,4,4,4/)
    spl%product_list(74,0:2) = (/2,5,6/)
    spl%product_list(75,0:3) = (/3,2,1,6/)
    spl%product_list(76,0:3) = (/3,2,1,6/)
    spl%product_list(77,0:3) = (/3,2,28,6/)
    spl%product_list(78,0:4) = (/4,24,4,1,6/)
    spl%product_list(79,0:3) = (/3,25,4,6/)
    spl%product_list(80,0:3) = (/3,10,1,6/)
    spl%product_list(81,0:2) = (/2,9,6/)
    spl%product_list(82,0:2) = (/2,12,6/)
    spl%product_list(83,0:3) = (/3,16,1,6/)
    spl%product_list(84,0:3) = (/3,30,4,6/)
    spl%product_list(85,0:3) = (/3,24,26,6/)
    spl%product_list(86,0:2) = (/2,14,6/)
    spl%product_list(87,0:3) = (/3,15,1,6/)
    spl%product_list(88,0:3) = (/3,12,4,6/)
    spl%product_list(89,0:4) = (/4,16,1,4,6/)
    spl%product_list(90,0:3) = (/3,25,27,6/)
    spl%product_list(91,0:3) = (/3,14,4,6/)
    spl%product_list(92,0:4) = (/4,15,1,4,6/)
    spl%product_list(93,0:4) = (/4,12,4,4,6/)
    spl%product_list(94,0:2) = (/2,10,4/)
    spl%product_list(95,0:2) = (/2,9,1/)
    spl%product_list(96,0:3) = (/3,15,4,4/)
    spl%product_list(97,0:3) = (/3,14,1,4/)
    spl%product_list(98,0:2) = (/2,19,4/)
    spl%product_list(99,0:2) = (/2,5,1/)
    spl%product_list(100,0:2) = (/2,21,6/)
    spl%product_list(101,0:2) = (/2,5,6/)
    spl%product_list(102,0:3) = (/3,14,1,6/)
    spl%product_list(103,0:3) = (/3,15,4,6/)
    spl%product_list(104,0:4) = (/4,15,4,4,6/)
    spl%product_list(105,0:3) = (/3,19,4,6/)
    spl%product_list(106,0:2) = (/2,29,4/)
    spl%product_list(107,0:2) = (/2,15,4/)
    spl%product_list(108,0:2) = (/2,19,4/)
    spl%product_list(109,0:3) = (/3,15,4,4/)
    spl%product_list(110,0:3) = (/3,19,4,4/)
    spl%product_list(111,0:2) = (/2,12,4/)
    spl%product_list(112,0:2) = (/2,15,1/)
    spl%product_list(113,0:2) = (/2,36,1/)
    spl%product_list(114,0:2) = (/2,15,26/)
    spl%product_list(115,0:2) = (/2,14,23/)
    spl%product_list(116,0:3) = (/3,36,1,4/)
    spl%product_list(117,0:2) = (/2,37,4/)
    spl%product_list(118,0:2) = (/2,38,1/)
    spl%product_list(119,0:2) = (/2,12,8/)
    spl%product_list(120,0:2) = (/2,15,28/)
    spl%product_list(121,0:2) = (/2,19,26/)
    spl%product_list(122,0:3) = (/3,38,1,4/)
    spl%product_list(123,0:2) = (/2,22,4/)
    spl%product_list(124,0:2) = (/2,25,1/)
    spl%product_list(125,0:2) = (/2,14,1/)
    spl%product_list(126,0:2) = (/2,15,4/)
    spl%product_list(127,0:3) = (/3,12,1,4/)
    spl%product_list(128,0:2) = (/2,37,1/)
    spl%product_list(129,0:2) = (/2,10,1/)
    spl%product_list(130,0:2) = (/2,19,1/)
    spl%product_list(131,0:2) = (/2,14,4/)
    spl%product_list(132,0:2) = (/2,38,1/)
    spl%product_list(133,0:2) = (/2,19,4/)
    spl%product_list(134,0:2) = (/2,38,4/)
    spl%product_list(135,0:2) = (/2,15,8/)
    spl%product_list(136,0:3) = (/3,38,4,4/)
    spl%product_list(137,0:2) = (/2,40,4/)
    spl%product_list(138,0:2) = (/2,19,8/)
    spl%product_list(139,0:3) = (/3,40,4,4/)
    spl%product_list(140,0:2) = (/2,42,4/)
    spl%product_list(141,0:2) = (/2,29,1/)
    spl%product_list(142,0:2) = (/2,29,28/)
    spl%product_list(143,0:2) = (/2,12,8/)
    spl%product_list(144,0:2) = (/2,15,28/)
    spl%product_list(145,0:3) = (/3,38,1,4/)
    spl%product_list(146,0:2) = (/2,14,8/)
    spl%product_list(147,0:2) = (/2,19,28/)
    spl%product_list(148,0:3) = (/3,40,1,4/)
    spl%product_list(149,0:3) = (/3,14,8,4/)
    spl%product_list(150,0:2) = (/2,9,4/)
    spl%product_list(151,0:2) = (/2,15,8/)
    spl%product_list(152,0:2) = (/2,19,8/)
    spl%product_list(153,0:3) = (/3,19,8,4/)
    spl%product_list(154,0:2) = (/2,35,8/)
    spl%product_list(155,0:2) = (/2,16,1/)
    spl%product_list(156,0:2) = (/2,16,28/)
    spl%product_list(157,0:2) = (/2,12,27/)
    spl%product_list(158,0:3) = (/3,36,1,4/)
    spl%product_list(159,0:2) = (/2,37,4/)
    spl%product_list(160,0:2) = (/2,38,1/)
    spl%product_list(161,0:2) = (/2,45,1/)
    spl%product_list(162,0:1) = (/1,0/)
    spl%product_list(163,0:2) = (/2,12,1/)
    spl%product_list(164,0:2) = (/2,12,28/)
    spl%product_list(165,0:2) = (/2,38,4/)
    spl%product_list(166,0:2) = (/2,39,1/)
    spl%product_list(167,0:1) = (/1,40/)
    spl%product_list(168,0:2) = (/2,46,1/)
    spl%product_list(169,0:1) = (/1,0/)
    spl%product_list(170,0:2) = (/2,15,1/)
    spl%product_list(171,0:2) = (/2,39,4/)
    spl%product_list(172,0:2) = (/2,40,1/)
    spl%product_list(173,0:2) = (/2,46,4/)
    spl%product_list(174,0:2) = (/2,47,1/)
    spl%product_list(175,0:2) = (/2,14,11/)
    spl%product_list(176,0:2) = (/2,38,28/)
    spl%product_list(177,0:2) = (/2,48,1/)
    spl%product_list(178,0:2) = (/2,14,13/)
    spl%product_list(179,0:2) = (/2,19,33/)
    spl%product_list(180,0:3) = (/3,38,28,4/)
    spl%product_list(181,0:2) = (/2,39,8/)
    spl%product_list(182,0:2) = (/2,40,28/)
    spl%product_list(183,0:3) = (/3,48,1,4/)
    spl%product_list(184,0:2) = (/2,49,1/)
    spl%product_list(185,0:2) = (/2,12,4/)
    spl%product_list(186,0:2) = (/2,40,4/)
    spl%product_list(187,0:2) = (/2,47,4/)
    spl%product_list(188,0:2) = (/2,19,11/)
    spl%product_list(189,0:2) = (/2,19,13/)
    spl%product_list(190,0:2) = (/2,40,8/)
    spl%product_list(191,0:2) = (/2,49,4/)
    spl%product_list(192,0:2) = (/2,15,4/)
    spl%product_list(193,0:2) = (/2,38,28/)
    spl%product_list(194,0:2) = (/2,48,1/)
    spl%product_list(195,0:2) = (/2,40,28/)
    spl%product_list(196,0:2) = (/2,49,1/)
    spl%product_list(197,0:2) = (/2,41,8/)
    spl%product_list(198,0:2) = (/2,42,28/)
    spl%product_list(199,0:2) = (/2,14,4/)
    spl%product_list(200,0:2) = (/2,42,4/)
    spl%product_list(201,0:2) = (/2,38,8/)
    spl%product_list(202,0:2) = (/2,48,4/)
    spl%product_list(203,0:2) = (/2,40,8/)
    spl%product_list(204,0:2) = (/2,50,4/)
    spl%product_list(205,0:2) = (/2,19,4/)
    spl%product_list(206,0:2) = (/2,19,33/)
    spl%product_list(207,0:2) = (/2,40,28/)
    spl%product_list(208,0:2) = (/2,49,1/)
    spl%product_list(209,0:2) = (/2,14,17/)
    spl%product_list(210,0:2) = (/2,43,8/)
    spl%product_list(211,0:2) = (/2,44,28/)
    spl%product_list(212,0:2) = (/2,15,17/)
    spl%product_list(213,0:2) = (/2,19,17/)
    spl%product_list(214,0:2) = (/2,21,4/)
    spl%product_list(215,0:2) = (/2,29,4/)
    spl%product_list(216,0:2) = (/2,29,6/)
    spl%product_list(217,0:2) = (/2,14,6/)
    spl%product_list(218,0:2) = (/2,4,4/)

    ! input Temperature profiles

    open(11, file = './Jupiter/Nakamura_22_Case1/input/temperature/T_e.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/temperature/T_i.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/temperature/T_n.dat', status = 'old' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-50_dp

    isp = sp_index(spl, 'H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'He')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/He.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CH4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/CH4.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H2.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H4.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H6')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H6.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v2)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2_v2.dat', status = 'old' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(iz,isp)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v4)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2_v4.dat', status = 'old' )
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


  end subroutine v__in__exe

end module v__in


! R1: H + hv -> H+ + e- 
! R2: H2 + hv -> H2+ + e- 
! R3: H2 + hv -> H+ + e- + H 
! R4: He + hv -> He+ + e- 
! R5: CH4 + hv -> CH4+ + e- 
! R6: CH4 + hv -> CH3+ + e- + H 
! R7: CH4 + hv -> H2+ + e- + products 
! R8: C2H2 + hv -> C2H2+ + e- 
! R9: C2H4 + hv -> C2H4+ + e- 
! R10: C2H4 + hv -> C2H3+ + e- + H 
! R11: C2H4 + hv -> C2H2+ + e- + products 
! R12: C2H4 + hv -> C2H+ + e- + products 
! R13: C2H6 + hv -> C2H6+ + e- 
! R14: C2H6 + hv -> C2H5+ + e- + H 
! R15: C2H6 + hv -> C2H4+ + e- + products 
! R16: C2H6 + hv -> C2H3+ + e- + products 
! R17: C2H6 + hv -> C2H2+ + e- + products 
! R18: H+ + e- -> H 
! R19: He+ + e- -> He 
! R20: HeH+ + e- -> H + He 
! R21: H2+ + e- -> H + H 
! R22: H3+ + e- -> H2 + H 
! R23: H3+ + e- -> H + H + H 
! R24: C+ + e- -> C 
! R25: CH+ + e- -> C + H 
! R26: CH2+ + e- -> CH + H 
! R27: CH3+ + e- -> CH2 + H 
! R28: CH4+ + e- -> CH3 + H 
! R29: CH4+ + e- -> CH2 + H + H 
! R30: CH5+ + e- -> CH2 + H + H2 
! R31: CH5+ + e- -> CH3 + H + H 
! R32: C2+ + e- -> C + C 
! R33: C2H+ + e- -> C2 + H 
! R34: C2H+ + e- -> CH + C 
! R35: C2H2+ + e- -> C2H + H 
! R36: C2H2+ + e- -> CH + CH 
! R37: C2H3+ + e- -> C2H2 + H 
! R38: C2H3+ + e- -> CH2 + CH 
! R39: C2H4+ + e- -> C2H3 + H 
! R40: C2H4+ + e- -> CH2 + CH2 
! R41: C2H5+ + e- -> C2H4 + H 
! R42: C2H5+ + e- -> CH3 + CH2 
! R43: C2H6+ + e- -> C2H5 + H 
! R44: C2H6+ + e- -> CH3 + CH3 
! R45: C2H7+ + e- -> C2H6 + H 
! R46: C3H+ + e- -> products 
! R47: C3H2+ + e- -> products 
! R48: C3H3+ + e- -> products 
! R49: C3H4+ + e- -> products 
! R50: C3H5+ + e- -> products 
! R51: C3H6+ + e- -> products 
! R52: C3H7+ + e- -> products 
! R53: C3H8+ + e- -> products 
! R54: C3H9+ + e- -> products 
! R55: C4H+ + e- -> products 
! R56: C4H2+ + e- -> products 
! R57: C4H3+ + e- -> products 
! R58: C4H5+ + e- -> products 
! R59: C4H7+ + e- -> products 
! R60: C4H9+ + e- -> products 
! R61: H2+ + H2 -> H3+ + H 
! R62: H2+ + H -> H+ + H2 
! R63: H2+ + He -> HeH+ + H 
! R64: H2+ + CH4 -> CH5+ + H 
! R65: H2+ + CH4 -> CH4+ + H2 
! R66: H2+ + CH4 -> CH3+ + H + H2 
! R67: H2+ + C2H2 -> C2H3+ + H 
! R68: H2+ + C2H2 -> C2H2+ + H2 
! R69: H2+ + C2H6 -> C2H6+ + H2 
! R70: H2+ + C2H6 -> C2H5+ + H + H2 
! R71: H2+ + C2H6 -> C2H4+ + H2 + H2 
! R72: H2+ + C2H6 -> C2H3+ + H + H2 + H2 
! R73: H2+ + C2H6 -> C2H2+ + H2 + H2 + H2 
! R74: He+ + H2 -> H2+ + He 
! R75: He+ + H2 -> H+ + H + He 
! R76: He+ + H2(v2) -> H+ + H + He 
! R77: He+ + CH4 -> H+ + CH3 + He 
! R78: He+ + CH4 -> CH+ + H2 + H + He 
! R79: He+ + CH4 -> CH2+ + H2 + He 
! R80: He+ + CH4 -> CH3+ + H + He 
! R81: He+ + CH4 -> CH4+ + He 
! R82: He+ + C2H2 -> C2H2+ + He 
! R83: He+ + C2H2 -> C2H+ + H + He 
! R84: He+ + C2H2 -> C2+ + H2 + He 
! R85: He+ + C2H2 -> CH+ + CH + He 
! R86: He+ + C2H4 -> C2H4+ + He 
! R87: He+ + C2H4 -> C2H3+ + H + He 
! R88: He+ + C2H4 -> C2H2+ + H2 + He 
! R89: He+ + C2H4 -> C2H+ + H + H2 + He 
! R90: He+ + C2H4 -> CH2+ + CH2 + He 
! R91: He+ + C2H6 -> C2H4+ + H2 + He 
! R92: He+ + C2H6 -> C2H3+ + H + H2 + He 
! R93: He+ + C2H6 -> C2H2+ + H2 + H2 + He 
! R94: H+ + CH4 -> CH3+ + H2 
! R95: H+ + CH4 -> CH4+ + H 
! R96: H+ + C2H6 -> C2H3+ + H2 + H2 
! R97: H+ + C2H6 -> C2H4+ + H + H2 
! R98: H+ + C2H6 -> C2H5+ + H2 
! R99: H+ + H2(v4) -> H2+ + H 
! R100: HeH+ + H2 -> H3+ + He 
! R101: HeH+ + H -> H2+ + He 
! R102: HeH+ + C2H4 -> C2H4+ + H + He 
! R103: HeH+ + C2H4 -> C2H3+ + H2 + He 
! R104: HeH+ + C2H6 -> C2H3+ + H2 + H2 + He 
! R105: HeH+ + C2H6 -> C2H5+ + H2 + He 
! R106: H3+ + CH4 -> CH5+ + H2 
! R107: H3+ + C2H2 -> C2H3+ + H2 
! R108: H3+ + C2H4 -> C2H5+ + H2 
! R109: H3+ + C2H4 -> C2H3+ + H2 + H2 
! R110: H3+ + C2H6 -> C2H5+ + H2 + H2 
! R111: C+ + CH4 -> C2H2+ + H2 
! R112: C+ + CH4 -> C2H3+ + H 
! R113: C+ + C2H2 -> C3H+ + H 
! R114: C+ + C2H4 -> C2H3+ + CH 
! R115: C+ + C2H4 -> C2H4+ + C 
! R116: C+ + C2H4 -> C3H+ + H + H2 
! R117: C+ + C2H4 -> C3H2+ + H2 
! R118: C+ + C2H4 -> C3H3+ + H 
! R119: C+ + C2H6 -> C2H2+ + CH4 
! R120: C+ + C2H6 -> C2H3+ + CH3 
! R121: C+ + C2H6 -> C2H5+ + CH 
! R122: C+ + C2H6 -> C3H3+ + H + H2 
! R123: CH+ + H -> C+ + H2 
! R124: CH+ + H2 -> CH2+ + H 
! R125: CH+ + CH4 -> C2H4+ + H 
! R126: CH+ + CH4 -> C2H3+ + H2 
! R127: CH+ + CH4 -> C2H2+ + H + H2 
! R128: CH+ + C2H2 -> C3H2+ + H 
! R129: CH2+ + H2 -> CH3+ + H 
! R130: CH2+ + CH4 -> C2H5+ + H 
! R131: CH2+ + CH4 -> C2H4+ + H2 
! R132: CH2+ + C2H2 -> C3H3+ + H 
! R133: CH3+ + CH4 -> C2H5+ + H2 
! R134: CH3+ + C2H2 -> C3H3+ + H2 
! R135: CH3+ + C2H4 -> C2H3+ + CH4 
! R136: CH3+ + C2H4 -> C3H3+ + H2 + H2 
! R137: CH3+ + C2H4 -> C3H5+ + H2 
! R138: CH3+ + C2H6 -> C2H5+ + CH4 
! R139: CH3+ + C2H6 -> C3H5+ + H2 + H2 
! R140: CH3+ + C2H6 -> C3H7+ + H2 
! R141: CH4+ + H2 -> CH5+ + H 
! R142: CH4+ + CH4 -> CH5+ + CH3 
! R143: CH4+ + C2H2 -> C2H2+ + CH4 
! R144: CH4+ + C2H2 -> C2H3+ + CH3 
! R145: CH4+ + C2H2 -> C3H3+ + H + H2 
! R146: CH4+ + C2H4 -> C2H4+ + CH4 
! R147: CH4+ + C2H4 -> C2H5+ + CH3 
! R148: CH4+ + C2H4 -> C3H5+ + H + H2 
! R149: CH4+ + C2H6 -> C2H4+ + CH4 + H2 
! R150: CH5+ + H -> CH4+ + H2 
! R151: CH5+ + C2H2 -> C2H3+ + CH4 
! R152: CH5+ + C2H4 -> C2H5+ + CH4 
! R153: CH5+ + C2H6 -> C2H5+ + CH4 + H2 
! R154: CH5+ + C2H6 -> C2H7+ + CH4 
! R155: C2+ + H2 -> C2H+ + H 
! R156: C2+ + CH4 -> C2H+ + CH3 
! R157: C2+ + CH4 -> C2H2+ + CH2 
! R158: C2+ + CH4 -> C3H+ + H + H2 
! R159: C2+ + CH4 -> C3H2+ + H2 
! R160: C2+ + CH4 -> C3H3+ + H 
! R161: C2+ + C2H2 -> C4H+ + H 
! R162: C2+ + C2H4 -> products 
! R163: C2H+ + H2 -> C2H2+ + H 
! R164: C2H+ + CH4 -> C2H2+ + CH3 
! R165: C2H+ + CH4 -> C3H3+ + H2 
! R166: C2H+ + CH4 -> C3H4+ + H 
! R167: C2H+ + CH4 -> C3H5+ 
! R168: C2H+ + C2H2 -> C4H2+ + H 
! R169: C2H+ + C2H4 -> products 
! R170: C2H2+ + H2 -> C2H3+ + H 
! R171: C2H2+ + CH4 -> C3H4+ + H2 
! R172: C2H2+ + CH4 -> C3H5+ + H 
! R173: C2H2+ + C2H2 -> C4H2+ + H2 
! R174: C2H2+ + C2H2 -> C4H3+ + H 
! R175: C2H2+ + C2H4 -> C2H4+ + C2H2 
! R176: C2H2+ + C2H4 -> C3H3+ + CH3 
! R177: C2H2+ + C2H4 -> C4H5+ + H 
! R178: C2H2+ + C2H6 -> C2H4+ + C2H4 
! R179: C2H2+ + C2H6 -> C2H5+ + C2H3 
! R180: C2H2+ + C2H6 -> C3H3+ + CH3 + H2 
! R181: C2H2+ + C2H6 -> C3H4+ + CH4 
! R182: C2H2+ + C2H6 -> C3H5+ + CH3 
! R183: C2H2+ + C2H6 -> C4H5+ + H + H2 
! R184: C2H2+ + C2H6 -> C4H7+ + H 
! R185: C2H3+ + H -> C2H2+ + H2 
! R186: C2H3+ + CH4 -> C3H5+ + H2 
! R187: C2H3+ + C2H2 -> C4H3+ + H2 
! R188: C2H3+ + C2H4 -> C2H5+ + C2H2 
! R189: C2H3+ + C2H6 -> C2H5+ + C2H4 
! R190: C2H3+ + C2H6 -> C3H5+ + CH4 
! R191: C2H3+ + C2H6 -> C4H7+ + H2 
! R192: C2H4+ + H -> C2H3+ + H2 
! R193: C2H4+ + C2H2 -> C3H3+ + CH3 
! R194: C2H4+ + C2H2 -> C4H5+ + H 
! R195: C2H4+ + C2H4 -> C3H5+ + CH3 
! R196: C2H4+ + C2H4 -> C4H7+ + H 
! R197: C2H4+ + C2H6 -> C3H6+ + CH4 
! R198: C2H4+ + C2H6 -> C3H7+ + CH3 
! R199: C2H5+ + H -> C2H4+ + H2 
! R200: C2H5+ + CH4 -> C3H7+ + H2 
! R201: C2H5+ + C2H2 -> C3H3+ + CH4 
! R202: C2H5+ + C2H2 -> C4H5+ + H2 
! R203: C2H5+ + C2H4 -> C3H5+ + CH4 
! R204: C2H5+ + C2H6 -> C4H9+ + H2 
! R205: C2H6+ + H -> C2H5+ + H2 
! R206: C2H6+ + C2H2 -> C2H5+ + C2H3 
! R207: C2H6+ + C2H2 -> C3H5+ + CH3 
! R208: C2H6+ + C2H2 -> C4H7+ + H 
! R209: C2H6+ + C2H4 -> C2H4+ + C2H6 
! R210: C2H6+ + C2H6 -> C3H8+ + CH4 
! R211: C2H6+ + C2H6 -> C3H9+ + CH3 
! R212: C2H7+ + C2H2 -> C2H3+ + C2H6 
! R213: C2H7+ + C2H4 -> C2H5+ + C2H6 
! R214: H+ + H2 + H2 -> H3+ + H2 
! R215: CH3+ + H2 + H2 -> CH5+ + H2 
! R216: CH3+ + H2 + He -> CH5+ + He 
! R217: C2H2+ + H2 + He -> C2H4+ + He 
! R218: H + H + H2 -> H2 + H2 
