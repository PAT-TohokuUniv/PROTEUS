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
    spl%planet = 'Venus'

    ! Calculation settings
    set%mode = '1D'
    set%F107 = 80.0_dp
    set%nstep = 30000
    set%fin_sec = 3.0e15_dp
    set%dtime_limit = 1.0e14_dp
    set%latitude = 0.0_dp
    set%Ls = 0.0_dp
    set%nday = 3.0_dp
    set%scheme = 'implicit'
    set%inversion = 'Catling'
    ! directory setting
    set%dir_name = './Venus/Krasnopolsky_2007_no-S'

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
    grd%nz    = 48
    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))
    grd%dalt(1:48) = 1.0e3_dp ! [m]
    grd%alt(1)      = 0.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 24
    spl%nsp_i   = 20
    spl%nch     = 42
    spl%nch_P   = 17
    spl%nch_L   = 25
    spl%n_Jlist = 171
    spl%nch_J   = 27
    spl%nrpn    = 9

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
    allocate(spl%species(spl%nsp))
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
    allocate(var%Li(spl%nsp_i,grd%nz))
    allocate(var%K_eddy(grd%nz),var%D_mol(spl%nsp,grd%nz))
    allocate(var%Fluxup(spl%nsp_i,0:grd%nz),var%Fluxdwn(spl%nsp_i,0:grd%nz))
    allocate(var%vFluxup(spl%nsp_i,grd%nz),var%vFluxdwn(spl%nsp_i,grd%nz))
    allocate(var%UpperBC(spl%nsp,3),var%LowerBC(spl%nsp,3))
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
    allocate(var%Amtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Lmtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))
    allocate(var%Umtx(spl%nsp_i*grd%nz,2*spl%nsp_i+1))

    ! species
    ! 
    ! SO3, H2O, H2SO4, CO, CO2, SO2, OCS, (SO)2, SO, H, SH, HCl, H2, Cl, H2S, SNO, NO, M, ClSO2, SO2Cl2, Cl2, HSCl, SCl, OH
    ! 
    spl%species(1) = 'SO3'
    spl%species(2) = 'H2O'
    spl%species(3) = 'H2SO4'
    spl%species(4) = 'CO'
    spl%species(5) = 'CO2'
    spl%species(6) = 'SO2'
    spl%species(7) = 'OCS'
    spl%species(8) = '(SO)2'
    spl%species(9) = 'SO'
    spl%species(10) = 'H'
    spl%species(11) = 'SH'
    spl%species(12) = 'HCl'
    spl%species(13) = 'H2'
    spl%species(14) = 'Cl'
    spl%species(15) = 'H2S'
    spl%species(16) = 'SNO'
    spl%species(17) = 'NO'
    spl%species(18) = 'M'
    spl%species(19) = 'ClSO2'
    spl%species(20) = 'SO2Cl2'
    spl%species(21) = 'Cl2'
    spl%species(22) = 'HSCl'
    spl%species(23) = 'SCl'
    spl%species(24) = 'OH'

    ! label_fix
    spl%label_fix(1) = 0 ! SO3: variable
    spl%label_fix(2) = 0 ! H2O: variable
    spl%label_fix(3) = 0 ! H2SO4: variable
    spl%label_fix(4) = 0 ! CO: variable
    spl%label_fix(5) = 0 ! CO2: variable
    spl%label_fix(6) = 0 ! SO2: variable
    spl%label_fix(7) = 0 ! OCS: variable
    spl%label_fix(8) = 1 ! (SO)2: fixed
    spl%label_fix(9) = 0 ! SO: variable
    spl%label_fix(10) = 0 ! H: variable
    spl%label_fix(11) = 0 ! SH: variable
    spl%label_fix(12) = 0 ! HCl: variable
    spl%label_fix(13) = 0 ! H2: variable
    spl%label_fix(14) = 0 ! Cl: variable
    spl%label_fix(15) = 0 ! H2S: variable
    spl%label_fix(16) = 0 ! SNO: variable
    spl%label_fix(17) = 1 ! NO: fixed
    spl%label_fix(18) = 1 ! M: fixed
    spl%label_fix(19) = 0 ! ClSO2: variable
    spl%label_fix(20) = 0 ! SO2Cl2: variable
    spl%label_fix(21) = 0 ! Cl2: variable
    spl%label_fix(22) = 0 ! HSCl: variable
    spl%label_fix(23) = 1 ! SCl: fixed
    spl%label_fix(24) = 0 ! OH: variable

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1 ! SO3: variable
    spl%all_to_var(2) = 2 ! H2O: variable
    spl%all_to_var(3) = 3 ! H2SO4: variable
    spl%all_to_var(4) = 4 ! CO: variable
    spl%all_to_var(5) = 5 ! CO2: variable
    spl%all_to_var(6) = 6 ! SO2: variable
    spl%all_to_var(7) = 7 ! OCS: variable
    spl%all_to_var(9) = 8 ! SO: variable
    spl%all_to_var(10) = 9 ! H: variable
    spl%all_to_var(11) = 10 ! SH: variable
    spl%all_to_var(12) = 11 ! HCl: variable
    spl%all_to_var(13) = 12 ! H2: variable
    spl%all_to_var(14) = 13 ! Cl: variable
    spl%all_to_var(15) = 14 ! H2S: variable
    spl%all_to_var(16) = 15 ! SNO: variable
    spl%all_to_var(19) = 16 ! ClSO2: variable
    spl%all_to_var(20) = 17 ! SO2Cl2: variable
    spl%all_to_var(21) = 18 ! Cl2: variable
    spl%all_to_var(22) = 19 ! HSCl: variable
    spl%all_to_var(24) = 20 ! OH: variable

    ! var_to_all
    spl%var_to_all(1) = 1 ! SO3: variable
    spl%var_to_all(2) = 2 ! H2O: variable
    spl%var_to_all(3) = 3 ! H2SO4: variable
    spl%var_to_all(4) = 4 ! CO: variable
    spl%var_to_all(5) = 5 ! CO2: variable
    spl%var_to_all(6) = 6 ! SO2: variable
    spl%var_to_all(7) = 7 ! OCS: variable
    spl%var_to_all(8) = 9 ! SO: variable
    spl%var_to_all(9) = 10 ! H: variable
    spl%var_to_all(10) = 11 ! SH: variable
    spl%var_to_all(11) = 12 ! HCl: variable
    spl%var_to_all(12) = 13 ! H2: variable
    spl%var_to_all(13) = 14 ! Cl: variable
    spl%var_to_all(14) = 15 ! H2S: variable
    spl%var_to_all(15) = 16 ! SNO: variable
    spl%var_to_all(16) = 19 ! ClSO2: variable
    spl%var_to_all(17) = 20 ! SO2Cl2: variable
    spl%var_to_all(18) = 21 ! Cl2: variable
    spl%var_to_all(19) = 22 ! HSCl: variable
    spl%var_to_all(20) = 24 ! OH: variable

    ! mass
    var%m(1) = 80.00000000_dp * cst%m_u !SO3
    var%m(2) = 18.00000000_dp * cst%m_u !H2O
    var%m(3) = 98.00000000_dp * cst%m_u !H2SO4
    var%m(4) = 28.00000000_dp * cst%m_u !CO
    var%m(5) = 44.00000000_dp * cst%m_u !CO2
    var%m(6) = 64.00000000_dp * cst%m_u !SO2
    var%m(7) = 60.00000000_dp * cst%m_u !OCS
    var%m(8) = 96.00000000_dp * cst%m_u !(SO)2
    var%m(9) = 48.00000000_dp * cst%m_u !SO
    var%m(10) = 1.00000000_dp * cst%m_u !H
    var%m(11) = 33.00000000_dp * cst%m_u !SH
    var%m(12) = 36.00000000_dp * cst%m_u !HCl
    var%m(13) = 2.00000000_dp * cst%m_u !H2
    var%m(14) = 35.00000000_dp * cst%m_u !Cl
    var%m(15) = 34.00000000_dp * cst%m_u !H2S
    var%m(16) = 62.00000000_dp * cst%m_u !SNO
    var%m(17) = 30.00000000_dp * cst%m_u !NO
    var%m(18) = 10.00000000_dp * cst%m_u !M
    var%m(19) = 99.00000000_dp * cst%m_u !ClSO2
    var%m(20) = 134.00000000_dp * cst%m_u !SO2Cl2
    var%m(21) = 70.00000000_dp * cst%m_u !Cl2
    var%m(22) = 68.00000000_dp * cst%m_u !HSCl
    var%m(23) = 67.00000000_dp * cst%m_u !SCl
    var%m(24) = 17.00000000_dp * cst%m_u !OH

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
    var%q(1) = 0.0_dp * cst%q_e !SO3
    var%q(2) = 0.0_dp * cst%q_e !H2O
    var%q(3) = 0.0_dp * cst%q_e !H2SO4
    var%q(4) = 0.0_dp * cst%q_e !CO
    var%q(5) = 0.0_dp * cst%q_e !CO2
    var%q(6) = 0.0_dp * cst%q_e !SO2
    var%q(7) = 0.0_dp * cst%q_e !OCS
    var%q(8) = 0.0_dp * cst%q_e !(SO)2
    var%q(9) = 0.0_dp * cst%q_e !SO
    var%q(10) = 0.0_dp * cst%q_e !H
    var%q(11) = 0.0_dp * cst%q_e !SH
    var%q(12) = 0.0_dp * cst%q_e !HCl
    var%q(13) = 0.0_dp * cst%q_e !H2
    var%q(14) = 0.0_dp * cst%q_e !Cl
    var%q(15) = 0.0_dp * cst%q_e !H2S
    var%q(16) = 0.0_dp * cst%q_e !SNO
    var%q(17) = 0.0_dp * cst%q_e !NO
    var%q(18) = 0.0_dp * cst%q_e !M
    var%q(19) = 0.0_dp * cst%q_e !ClSO2
    var%q(20) = 0.0_dp * cst%q_e !SO2Cl2
    var%q(21) = 0.0_dp * cst%q_e !Cl2
    var%q(22) = 0.0_dp * cst%q_e !HSCl
    var%q(23) = 0.0_dp * cst%q_e !SCl
    var%q(24) = 0.0_dp * cst%q_e !OH

    ! read P, L, J list
    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/PLJ_list/T_range.dat', status = 'unknown' )
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
      else if ( spl%reaction_type_list(ich) == 333 ) then
        spl%reaction_type_char(ich) = 'S3 photolysis'
        var%nspecial = var%nspecial + 1
      end if
    end do

    ! input Temperature profiles

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Venus/Krasnopolsky_2007_no-S/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'H2O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Venus/Krasnopolsky_2007_no-S/input/density/H2O.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Venus/Krasnopolsky_2007_no-S/input/density/CO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'SO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Venus/Krasnopolsky_2007_no-S/input/density/SO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'NO')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Venus/Krasnopolsky_2007_no-S/input/density/NO.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    var%ni_0 = var%ni

    ! Lower boundary condition
    var%LowerBC = 0.0_dp

    isp = sp_index(spl, 'CO2')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 9.07e26_dp

    isp = sp_index(spl, 'SO2')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 1.18e23_dp

    isp = sp_index(spl, 'H2O')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 2.72e22_dp

    isp = sp_index(spl, 'HCl')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 4.54e20_dp

    isp = sp_index(spl, 'OCS')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 2.46e22_dp

    isp = sp_index(spl, 'NO')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 4.99e18_dp


    ! Upper boundary condition
    var%UpperBC = 0.0_dp

    isp = sp_index(spl, 'H2SO4')
    var%UpperBC(isp,1) = 2.0_dp
    var%UpperBC(isp,2) = -2.2e16_dp


  end subroutine v__in__exe

end module v__in
