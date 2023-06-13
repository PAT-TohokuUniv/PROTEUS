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

    ! Planet type
    spl%planet = 'Jupiter'

    ! Calculation settings
    set%mode = '1D'
    set%F107 = 80.0_dp
    set%nstep = 30000
    set%fin_sec = 35729.685e2_dp
    set%dtime_limit = 1.0e5_dp
    set%latitude = 0.0_dp
    set%Ls = 0.0_dp
    set%nday = 3.0_dp
    set%scheme = 'implicit'
    set%inversion = 'default'
    ! directory setting
    set%dir_name = './Jupiter/Nakamura_22_Case1'

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
    grd%dalt(1:141) = 20.0e3_dp ! [m]
    grd%alt(1)      = 200.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 52
    spl%nsp_i   = 36
    spl%nch     = 218
    spl%nch_P   = 147
    spl%nch_L   = 87
    spl%n_Jlist = 255
    spl%nch_J   = 89
    spl%nrpn    = 8

    ! allocate
    allocate(var%ni(0:spl%nsp,grd%nz))
    allocate(var%ni_0(0:spl%nsp,grd%nz))
    allocate(var%ni_new(spl%nsp,grd%nz))
    allocate(var%ni_stable(spl%nsp,grd%ny,grd%nz))
    allocate(var%ni_3d(spl%nsp,grd%nx,grd%ny,grd%nz))
    allocate(var%clm_ni(spl%nsp,grd%nz))
    allocate(var%Ti(grd%nz),var%Te(grd%nz),var%Tn(grd%nz))
    allocate(var%Ti_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Te_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%Tn_3d(grd%nx,grd%ny,grd%nz))
    allocate(var%m(spl%nsp), var%q(spl%nsp))
    allocate(spl%reactant_list(spl%nch, 10))
    allocate(spl%product_list(spl%nch, 10))
    allocate(spl%species(0:spl%nsp))
    allocate(spl%label_fix(spl%nsp))
    allocate(spl%all_to_var(spl%nsp))
    allocate(spl%var_to_all(spl%nsp_i))
    allocate(spl%Prod_list(spl%nsp_i, spl%nch_P))
    allocate(spl%Loss_list(spl%nsp_i, spl%nch_L))
    allocate(spl%Jmtx_list(spl%n_Jlist, spl%nch_J))
    allocate(spl%reaction_type_list(spl%nch))
    allocate(spl%reaction_type_char(spl%nch))
    allocate(var%ki(spl%nch,grd%nz))
    allocate(var%rate(spl%nch,grd%nz))
    allocate(var%Pi(spl%nsp_i,grd%nz))
    allocate(var%Pij(spl%nsp_i,grd%nz,spl%nch))
    allocate(var%Li(spl%nsp_i,grd%nz))
    allocate(var%K_eddy(grd%nz),var%D_mol(spl%nsp,grd%nz))
    allocate(var%Fluxup(spl%nsp_i,0:grd%nz),var%Fluxdwn(spl%nsp_i,0:grd%nz))
    allocate(var%vFluxup(spl%nsp_i,grd%nz),var%vFluxdwn(spl%nsp_i,grd%nz))
    allocate(var%UpperBC(0:spl%nsp,3),var%LowerBC(0:spl%nsp,3))
    allocate(var%Jmtx(spl%nsp_i, spl%nsp_i))
    allocate(spl%rate_rpn_token(spl%nch,3,spl%nrpn))
    allocate(spl%rate_rpn_label(spl%nch,3,spl%nrpn))
    allocate(spl%rate_cases(spl%nch))
    allocate(spl%T_range(spl%nch,3,3))
    allocate(spl%major_species(grd%nz))
    allocate(var%n_tot(grd%nz),var%m_mean(grd%nz))
    allocate(var%Phip(spl%nsp_i,grd%nz))
    allocate(var%Phim(spl%nsp_i,grd%nz))
    allocate(var%dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dniu_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dni0_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dnil_dPhi_dz(spl%nsp_i,grd%nz))
    allocate(var%d_dneu_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dne0_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%d_dnel_dPhi_dz_add(spl%nsp_i,grd%nz))
    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))
    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))
    allocate(var%tAmtx(2*spl%nsp_i+1,spl%nsp_i*grd%nz))
    allocate(var%tLmtx(spl%nsp_i+1,spl%nsp_i*grd%nz))
    allocate(var%Umtx(spl%nsp_i+1,spl%nsp_i*grd%nz))

    ! species
    ! 
    ! H, H+, e-, H2, H2+, He, He+, CH4, CH4+, CH3+, C2H2, C2H2+, C2H4, C2H4+, C2H3+, C2H+, C2H6, C2H6+, C2H5+, HeH+, H3+, C+, C, CH+, CH2+, CH, CH2, CH3, CH5+, C2+, C2, C2H, C2H3, C2H5, C2H7+, C3H+, C3H2+, C3H3+, C3H4+, C3H5+, C3H6+, C3H7+, C3H8+, C3H9+, C4H+, C4H2+, C4H3+, C4H5+, C4H7+, C4H9+, H2(v>=2), H2(v>=4)
    ! 
    spl%species(1) = 'H'
    spl%species(2) = 'H+'
    spl%species(3) = 'e-'
    spl%species(4) = 'H2'
    spl%species(5) = 'H2+'
    spl%species(6) = 'He'
    spl%species(7) = 'He+'
    spl%species(8) = 'CH4'
    spl%species(9) = 'CH4+'
    spl%species(10) = 'CH3+'
    spl%species(11) = 'C2H2'
    spl%species(12) = 'C2H2+'
    spl%species(13) = 'C2H4'
    spl%species(14) = 'C2H4+'
    spl%species(15) = 'C2H3+'
    spl%species(16) = 'C2H+'
    spl%species(17) = 'C2H6'
    spl%species(18) = 'C2H6+'
    spl%species(19) = 'C2H5+'
    spl%species(20) = 'HeH+'
    spl%species(21) = 'H3+'
    spl%species(22) = 'C+'
    spl%species(23) = 'C'
    spl%species(24) = 'CH+'
    spl%species(25) = 'CH2+'
    spl%species(26) = 'CH'
    spl%species(27) = 'CH2'
    spl%species(28) = 'CH3'
    spl%species(29) = 'CH5+'
    spl%species(30) = 'C2+'
    spl%species(31) = 'C2'
    spl%species(32) = 'C2H'
    spl%species(33) = 'C2H3'
    spl%species(34) = 'C2H5'
    spl%species(35) = 'C2H7+'
    spl%species(36) = 'C3H+'
    spl%species(37) = 'C3H2+'
    spl%species(38) = 'C3H3+'
    spl%species(39) = 'C3H4+'
    spl%species(40) = 'C3H5+'
    spl%species(41) = 'C3H6+'
    spl%species(42) = 'C3H7+'
    spl%species(43) = 'C3H8+'
    spl%species(44) = 'C3H9+'
    spl%species(45) = 'C4H+'
    spl%species(46) = 'C4H2+'
    spl%species(47) = 'C4H3+'
    spl%species(48) = 'C4H5+'
    spl%species(49) = 'C4H7+'
    spl%species(50) = 'C4H9+'
    spl%species(51) = 'H2(v>=2)'
    spl%species(52) = 'H2(v>=4)'

    ! label_fix
    spl%label_fix(1) = 0 ! H: variable
    spl%label_fix(2) = 0 ! H+: variable
    spl%label_fix(3) = 0 ! e-: variable
    spl%label_fix(4) = 1 ! H2: fixed
    spl%label_fix(5) = 0 ! H2+: variable
    spl%label_fix(6) = 1 ! He: fixed
    spl%label_fix(7) = 0 ! He+: variable
    spl%label_fix(8) = 1 ! CH4: fixed
    spl%label_fix(9) = 0 ! CH4+: variable
    spl%label_fix(10) = 0 ! CH3+: variable
    spl%label_fix(11) = 1 ! C2H2: fixed
    spl%label_fix(12) = 0 ! C2H2+: variable
    spl%label_fix(13) = 1 ! C2H4: fixed
    spl%label_fix(14) = 0 ! C2H4+: variable
    spl%label_fix(15) = 0 ! C2H3+: variable
    spl%label_fix(16) = 0 ! C2H+: variable
    spl%label_fix(17) = 1 ! C2H6: fixed
    spl%label_fix(18) = 0 ! C2H6+: variable
    spl%label_fix(19) = 0 ! C2H5+: variable
    spl%label_fix(20) = 0 ! HeH+: variable
    spl%label_fix(21) = 0 ! H3+: variable
    spl%label_fix(22) = 0 ! C+: variable
    spl%label_fix(23) = 1 ! C: fixed
    spl%label_fix(24) = 0 ! CH+: variable
    spl%label_fix(25) = 0 ! CH2+: variable
    spl%label_fix(26) = 1 ! CH: fixed
    spl%label_fix(27) = 1 ! CH2: fixed
    spl%label_fix(28) = 1 ! CH3: fixed
    spl%label_fix(29) = 0 ! CH5+: variable
    spl%label_fix(30) = 0 ! C2+: variable
    spl%label_fix(31) = 1 ! C2: fixed
    spl%label_fix(32) = 1 ! C2H: fixed
    spl%label_fix(33) = 1 ! C2H3: fixed
    spl%label_fix(34) = 1 ! C2H5: fixed
    spl%label_fix(35) = 0 ! C2H7+: variable
    spl%label_fix(36) = 0 ! C3H+: variable
    spl%label_fix(37) = 0 ! C3H2+: variable
    spl%label_fix(38) = 0 ! C3H3+: variable
    spl%label_fix(39) = 0 ! C3H4+: variable
    spl%label_fix(40) = 0 ! C3H5+: variable
    spl%label_fix(41) = 0 ! C3H6+: variable
    spl%label_fix(42) = 0 ! C3H7+: variable
    spl%label_fix(43) = 0 ! C3H8+: variable
    spl%label_fix(44) = 0 ! C3H9+: variable
    spl%label_fix(45) = 0 ! C4H+: variable
    spl%label_fix(46) = 0 ! C4H2+: variable
    spl%label_fix(47) = 0 ! C4H3+: variable
    spl%label_fix(48) = 0 ! C4H5+: variable
    spl%label_fix(49) = 0 ! C4H7+: variable
    spl%label_fix(50) = 0 ! C4H9+: variable
    spl%label_fix(51) = 1 ! H2(v>=2): fixed
    spl%label_fix(52) = 1 ! H2(v>=4): fixed

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1 ! H: variable
    spl%all_to_var(2) = 2 ! H+: variable
    spl%all_to_var(3) = 3 ! e-: variable
    spl%all_to_var(5) = 4 ! H2+: variable
    spl%all_to_var(7) = 5 ! He+: variable
    spl%all_to_var(9) = 6 ! CH4+: variable
    spl%all_to_var(10) = 7 ! CH3+: variable
    spl%all_to_var(12) = 8 ! C2H2+: variable
    spl%all_to_var(14) = 9 ! C2H4+: variable
    spl%all_to_var(15) = 10 ! C2H3+: variable
    spl%all_to_var(16) = 11 ! C2H+: variable
    spl%all_to_var(18) = 12 ! C2H6+: variable
    spl%all_to_var(19) = 13 ! C2H5+: variable
    spl%all_to_var(20) = 14 ! HeH+: variable
    spl%all_to_var(21) = 15 ! H3+: variable
    spl%all_to_var(22) = 16 ! C+: variable
    spl%all_to_var(24) = 17 ! CH+: variable
    spl%all_to_var(25) = 18 ! CH2+: variable
    spl%all_to_var(29) = 19 ! CH5+: variable
    spl%all_to_var(30) = 20 ! C2+: variable
    spl%all_to_var(35) = 21 ! C2H7+: variable
    spl%all_to_var(36) = 22 ! C3H+: variable
    spl%all_to_var(37) = 23 ! C3H2+: variable
    spl%all_to_var(38) = 24 ! C3H3+: variable
    spl%all_to_var(39) = 25 ! C3H4+: variable
    spl%all_to_var(40) = 26 ! C3H5+: variable
    spl%all_to_var(41) = 27 ! C3H6+: variable
    spl%all_to_var(42) = 28 ! C3H7+: variable
    spl%all_to_var(43) = 29 ! C3H8+: variable
    spl%all_to_var(44) = 30 ! C3H9+: variable
    spl%all_to_var(45) = 31 ! C4H+: variable
    spl%all_to_var(46) = 32 ! C4H2+: variable
    spl%all_to_var(47) = 33 ! C4H3+: variable
    spl%all_to_var(48) = 34 ! C4H5+: variable
    spl%all_to_var(49) = 35 ! C4H7+: variable
    spl%all_to_var(50) = 36 ! C4H9+: variable

    ! var_to_all
    spl%var_to_all(1) = 1 ! H: variable
    spl%var_to_all(2) = 2 ! H+: variable
    spl%var_to_all(3) = 3 ! e-: variable
    spl%var_to_all(4) = 5 ! H2+: variable
    spl%var_to_all(5) = 7 ! He+: variable
    spl%var_to_all(6) = 9 ! CH4+: variable
    spl%var_to_all(7) = 10 ! CH3+: variable
    spl%var_to_all(8) = 12 ! C2H2+: variable
    spl%var_to_all(9) = 14 ! C2H4+: variable
    spl%var_to_all(10) = 15 ! C2H3+: variable
    spl%var_to_all(11) = 16 ! C2H+: variable
    spl%var_to_all(12) = 18 ! C2H6+: variable
    spl%var_to_all(13) = 19 ! C2H5+: variable
    spl%var_to_all(14) = 20 ! HeH+: variable
    spl%var_to_all(15) = 21 ! H3+: variable
    spl%var_to_all(16) = 22 ! C+: variable
    spl%var_to_all(17) = 24 ! CH+: variable
    spl%var_to_all(18) = 25 ! CH2+: variable
    spl%var_to_all(19) = 29 ! CH5+: variable
    spl%var_to_all(20) = 30 ! C2+: variable
    spl%var_to_all(21) = 35 ! C2H7+: variable
    spl%var_to_all(22) = 36 ! C3H+: variable
    spl%var_to_all(23) = 37 ! C3H2+: variable
    spl%var_to_all(24) = 38 ! C3H3+: variable
    spl%var_to_all(25) = 39 ! C3H4+: variable
    spl%var_to_all(26) = 40 ! C3H5+: variable
    spl%var_to_all(27) = 41 ! C3H6+: variable
    spl%var_to_all(28) = 42 ! C3H7+: variable
    spl%var_to_all(29) = 43 ! C3H8+: variable
    spl%var_to_all(30) = 44 ! C3H9+: variable
    spl%var_to_all(31) = 45 ! C4H+: variable
    spl%var_to_all(32) = 46 ! C4H2+: variable
    spl%var_to_all(33) = 47 ! C4H3+: variable
    spl%var_to_all(34) = 48 ! C4H5+: variable
    spl%var_to_all(35) = 49 ! C4H7+: variable
    spl%var_to_all(36) = 50 ! C4H9+: variable

    ! mass
    var%m(1) = 1.00000000_dp * cst%m_u !H
    var%m(2) = 0.99945142_dp * cst%m_u !H+
    var%m(3) = 0.00054858_dp * cst%m_u !e-
    var%m(4) = 2.00000000_dp * cst%m_u !H2
    var%m(5) = 1.99945142_dp * cst%m_u !H2+
    var%m(6) = 4.00000000_dp * cst%m_u !He
    var%m(7) = 3.99945142_dp * cst%m_u !He+
    var%m(8) = 16.00000000_dp * cst%m_u !CH4
    var%m(9) = 15.99945142_dp * cst%m_u !CH4+
    var%m(10) = 14.99945142_dp * cst%m_u !CH3+
    var%m(11) = 26.00000000_dp * cst%m_u !C2H2
    var%m(12) = 25.99945142_dp * cst%m_u !C2H2+
    var%m(13) = 28.00000000_dp * cst%m_u !C2H4
    var%m(14) = 27.99945142_dp * cst%m_u !C2H4+
    var%m(15) = 26.99945142_dp * cst%m_u !C2H3+
    var%m(16) = 24.99945142_dp * cst%m_u !C2H+
    var%m(17) = 30.00000000_dp * cst%m_u !C2H6
    var%m(18) = 29.99945142_dp * cst%m_u !C2H6+
    var%m(19) = 28.99945142_dp * cst%m_u !C2H5+
    var%m(20) = 4.99945142_dp * cst%m_u !HeH+
    var%m(21) = 2.99945142_dp * cst%m_u !H3+
    var%m(22) = 11.99945142_dp * cst%m_u !C+
    var%m(23) = 12.00000000_dp * cst%m_u !C
    var%m(24) = 12.99945142_dp * cst%m_u !CH+
    var%m(25) = 13.99945142_dp * cst%m_u !CH2+
    var%m(26) = 13.00000000_dp * cst%m_u !CH
    var%m(27) = 14.00000000_dp * cst%m_u !CH2
    var%m(28) = 15.00000000_dp * cst%m_u !CH3
    var%m(29) = 16.99945142_dp * cst%m_u !CH5+
    var%m(30) = 23.99945142_dp * cst%m_u !C2+
    var%m(31) = 24.00000000_dp * cst%m_u !C2
    var%m(32) = 25.00000000_dp * cst%m_u !C2H
    var%m(33) = 27.00000000_dp * cst%m_u !C2H3
    var%m(34) = 29.00000000_dp * cst%m_u !C2H5
    var%m(35) = 30.99945142_dp * cst%m_u !C2H7+
    var%m(36) = 36.99945142_dp * cst%m_u !C3H+
    var%m(37) = 37.99945142_dp * cst%m_u !C3H2+
    var%m(38) = 38.99945142_dp * cst%m_u !C3H3+
    var%m(39) = 39.99945142_dp * cst%m_u !C3H4+
    var%m(40) = 40.99945142_dp * cst%m_u !C3H5+
    var%m(41) = 41.99945142_dp * cst%m_u !C3H6+
    var%m(42) = 42.99945142_dp * cst%m_u !C3H7+
    var%m(43) = 43.99945142_dp * cst%m_u !C3H8+
    var%m(44) = 44.99945142_dp * cst%m_u !C3H9+
    var%m(45) = 48.99945142_dp * cst%m_u !C4H+
    var%m(46) = 49.99945142_dp * cst%m_u !C4H2+
    var%m(47) = 50.99945142_dp * cst%m_u !C4H3+
    var%m(48) = 52.99945142_dp * cst%m_u !C4H5+
    var%m(49) = 54.99945142_dp * cst%m_u !C4H7+
    var%m(50) = 56.99945142_dp * cst%m_u !C4H9+
    var%m(51) = 2.00000000_dp * cst%m_u !H2(v>=2)
    var%m(52) = 2.00000000_dp * cst%m_u !H2(v>=4)

    ! mass zero error
    do isp = 1, spl%nsp
      if ( var%m(isp) == 0.0_dp ) then 
        write(*,*) 'mass zero error!'
        write(*,*) 'mass of ',trim(ADJUSTL(spl%species(isp))),' is zero.'
        write(*,*) 'Calculation stopped.'
        stop
      end  if
    end do

    ! charge
    var%q(1) = 0.0_dp * cst%q_e !H
    var%q(2) = 1.0_dp * cst%q_e !H+
    var%q(3) = -1.0_dp * cst%q_e !e-
    var%q(4) = 0.0_dp * cst%q_e !H2
    var%q(5) = 1.0_dp * cst%q_e !H2+
    var%q(6) = 0.0_dp * cst%q_e !He
    var%q(7) = 1.0_dp * cst%q_e !He+
    var%q(8) = 0.0_dp * cst%q_e !CH4
    var%q(9) = 1.0_dp * cst%q_e !CH4+
    var%q(10) = 1.0_dp * cst%q_e !CH3+
    var%q(11) = 0.0_dp * cst%q_e !C2H2
    var%q(12) = 1.0_dp * cst%q_e !C2H2+
    var%q(13) = 0.0_dp * cst%q_e !C2H4
    var%q(14) = 1.0_dp * cst%q_e !C2H4+
    var%q(15) = 1.0_dp * cst%q_e !C2H3+
    var%q(16) = 1.0_dp * cst%q_e !C2H+
    var%q(17) = 0.0_dp * cst%q_e !C2H6
    var%q(18) = 1.0_dp * cst%q_e !C2H6+
    var%q(19) = 1.0_dp * cst%q_e !C2H5+
    var%q(20) = 1.0_dp * cst%q_e !HeH+
    var%q(21) = 1.0_dp * cst%q_e !H3+
    var%q(22) = 1.0_dp * cst%q_e !C+
    var%q(23) = 0.0_dp * cst%q_e !C
    var%q(24) = 1.0_dp * cst%q_e !CH+
    var%q(25) = 1.0_dp * cst%q_e !CH2+
    var%q(26) = 0.0_dp * cst%q_e !CH
    var%q(27) = 0.0_dp * cst%q_e !CH2
    var%q(28) = 0.0_dp * cst%q_e !CH3
    var%q(29) = 1.0_dp * cst%q_e !CH5+
    var%q(30) = 1.0_dp * cst%q_e !C2+
    var%q(31) = 0.0_dp * cst%q_e !C2
    var%q(32) = 0.0_dp * cst%q_e !C2H
    var%q(33) = 0.0_dp * cst%q_e !C2H3
    var%q(34) = 0.0_dp * cst%q_e !C2H5
    var%q(35) = 1.0_dp * cst%q_e !C2H7+
    var%q(36) = 1.0_dp * cst%q_e !C3H+
    var%q(37) = 1.0_dp * cst%q_e !C3H2+
    var%q(38) = 1.0_dp * cst%q_e !C3H3+
    var%q(39) = 1.0_dp * cst%q_e !C3H4+
    var%q(40) = 1.0_dp * cst%q_e !C3H5+
    var%q(41) = 1.0_dp * cst%q_e !C3H6+
    var%q(42) = 1.0_dp * cst%q_e !C3H7+
    var%q(43) = 1.0_dp * cst%q_e !C3H8+
    var%q(44) = 1.0_dp * cst%q_e !C3H9+
    var%q(45) = 1.0_dp * cst%q_e !C4H+
    var%q(46) = 1.0_dp * cst%q_e !C4H2+
    var%q(47) = 1.0_dp * cst%q_e !C4H3+
    var%q(48) = 1.0_dp * cst%q_e !C4H5+
    var%q(49) = 1.0_dp * cst%q_e !C4H7+
    var%q(50) = 1.0_dp * cst%q_e !C4H9+
    var%q(51) = 0.0_dp * cst%q_e !H2(v>=2)
    var%q(52) = 0.0_dp * cst%q_e !H2(v>=4)

    ! read P, L, J list
    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/PLJ_list/T_range.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%T_range(ich,i,j), j = 1, 3)
        end do
      end do
    close(11)

    ! reaction type characters
    var%nspecial = 0
    do ich = 1, spl%nch
      if (      spl%reaction_type_list(ich) == 1 ) then
        spl%reaction_type_char(ich) = 'photoionization'
      else if ( spl%reaction_type_list(ich) == 2 ) then
        spl%reaction_type_char(ich) = 'photodissociation'
      else if ( spl%reaction_type_list(ich) == 11 ) then
        spl%reaction_type_char(ich) = 'electron impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 12 ) then
        spl%reaction_type_char(ich) = 'proton impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 13 ) then
        spl%reaction_type_char(ich) = 'H impact'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 30 ) then
        spl%reaction_type_char(ich) = 'pressure_dependent_3body'
      else if ( spl%reaction_type_list(ich) == 31 ) then
        spl%reaction_type_char(ich) = 'pressure_dependent_3bodyM'
      else if ( spl%reaction_type_list(ich) == 32 ) then
        spl%reaction_type_char(ich) = 'Lindemann-Hinshelwood'
      else if ( spl%reaction_type_list(ich) == 41 ) then
        spl%reaction_type_char(ich) = 'Meteoroid ablation'
        var%nspecial = var%nspecial + 1
      else if ( spl%reaction_type_list(ich) == 47 ) then
        spl%reaction_type_char(ich) = 'Rainout'
        var%nspecial = var%nspecial + 1
      end if
    end do

    ! input Temperature profiles

    open(11, file = './Jupiter/Nakamura_22_Case1/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Jupiter/Nakamura_22_Case1/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'He')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/He.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CH4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/CH4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H4')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'C2H6')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/C2H6.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v>=2)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2_v2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2(v>=4)')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Jupiter/Nakamura_22_Case1/input/density/H2_v4.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    var%ni_0 = var%ni

    ! Lower boundary condition
    var%LowerBC = 0.0_dp


    ! Upper boundary condition
    var%UpperBC = 0.0_dp


  end subroutine v__in__exe

end module v__in
