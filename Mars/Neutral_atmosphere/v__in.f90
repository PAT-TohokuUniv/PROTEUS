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
    spl%planet = 'Mars'

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
    set%inversion = 'default'
    ! directory setting
    set%dir_name = './Mars/Neutral_atmosphere'

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
    grd%dalt(1:101) = 2.0e3_dp ! [m]
    grd%alt(1)      = 0.0e3_dp ! [m]
    do iz = 2, grd%nz
      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)
    end do

    ! reactions, chemical species
    spl%nsp     = 18
    spl%nsp_i   = 15
    spl%nch     = 62
    spl%nch_P   = 33
    spl%nch_L   = 29
    spl%n_Jlist = 128
    spl%nch_J   = 31
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
    ! CO2, CO, O, O(1D), H2O, H, OH, H2, O3, O2, HO2, H2O2, M, N2, HOCO, CO2+, HCO, H2CO
    ! 
    spl%species(1) = 'CO2'
    spl%species(2) = 'CO'
    spl%species(3) = 'O'
    spl%species(4) = 'O(1D)'
    spl%species(5) = 'H2O'
    spl%species(6) = 'H'
    spl%species(7) = 'OH'
    spl%species(8) = 'H2'
    spl%species(9) = 'O3'
    spl%species(10) = 'O2'
    spl%species(11) = 'HO2'
    spl%species(12) = 'H2O2'
    spl%species(13) = 'M'
    spl%species(14) = 'N2'
    spl%species(15) = 'HOCO'
    spl%species(16) = 'CO2+'
    spl%species(17) = 'HCO'
    spl%species(18) = 'H2CO'

    ! label_fix
    spl%label_fix(1) = 0 ! CO2: variable
    spl%label_fix(2) = 0 ! CO: variable
    spl%label_fix(3) = 0 ! O: variable
    spl%label_fix(4) = 0 ! O(1D): variable
    spl%label_fix(5) = 1 ! H2O: fixed
    spl%label_fix(6) = 0 ! H: variable
    spl%label_fix(7) = 0 ! OH: variable
    spl%label_fix(8) = 0 ! H2: variable
    spl%label_fix(9) = 0 ! O3: variable
    spl%label_fix(10) = 0 ! O2: variable
    spl%label_fix(11) = 0 ! HO2: variable
    spl%label_fix(12) = 0 ! H2O2: variable
    spl%label_fix(13) = 1 ! M: fixed
    spl%label_fix(14) = 0 ! N2: variable
    spl%label_fix(15) = 0 ! HOCO: variable
    spl%label_fix(16) = 1 ! CO2+: fixed
    spl%label_fix(17) = 0 ! HCO: variable
    spl%label_fix(18) = 0 ! H2CO: variable

    ! all_to_var
    spl%all_to_var = 0
    spl%all_to_var(1) = 1 ! CO2: variable
    spl%all_to_var(2) = 2 ! CO: variable
    spl%all_to_var(3) = 3 ! O: variable
    spl%all_to_var(4) = 4 ! O(1D): variable
    spl%all_to_var(6) = 5 ! H: variable
    spl%all_to_var(7) = 6 ! OH: variable
    spl%all_to_var(8) = 7 ! H2: variable
    spl%all_to_var(9) = 8 ! O3: variable
    spl%all_to_var(10) = 9 ! O2: variable
    spl%all_to_var(11) = 10 ! HO2: variable
    spl%all_to_var(12) = 11 ! H2O2: variable
    spl%all_to_var(14) = 12 ! N2: variable
    spl%all_to_var(15) = 13 ! HOCO: variable
    spl%all_to_var(17) = 14 ! HCO: variable
    spl%all_to_var(18) = 15 ! H2CO: variable

    ! var_to_all
    spl%var_to_all(1) = 1 ! CO2: variable
    spl%var_to_all(2) = 2 ! CO: variable
    spl%var_to_all(3) = 3 ! O: variable
    spl%var_to_all(4) = 4 ! O(1D): variable
    spl%var_to_all(5) = 6 ! H: variable
    spl%var_to_all(6) = 7 ! OH: variable
    spl%var_to_all(7) = 8 ! H2: variable
    spl%var_to_all(8) = 9 ! O3: variable
    spl%var_to_all(9) = 10 ! O2: variable
    spl%var_to_all(10) = 11 ! HO2: variable
    spl%var_to_all(11) = 12 ! H2O2: variable
    spl%var_to_all(12) = 14 ! N2: variable
    spl%var_to_all(13) = 15 ! HOCO: variable
    spl%var_to_all(14) = 17 ! HCO: variable
    spl%var_to_all(15) = 18 ! H2CO: variable

    ! mass
    var%m(1) = 44.00000000_dp * cst%m_u !CO2
    var%m(2) = 28.00000000_dp * cst%m_u !CO
    var%m(3) = 16.00000000_dp * cst%m_u !O
    var%m(4) = 16.00000000_dp * cst%m_u !O(1D)
    var%m(5) = 18.00000000_dp * cst%m_u !H2O
    var%m(6) = 1.00000000_dp * cst%m_u !H
    var%m(7) = 17.00000000_dp * cst%m_u !OH
    var%m(8) = 2.00000000_dp * cst%m_u !H2
    var%m(9) = 48.00000000_dp * cst%m_u !O3
    var%m(10) = 32.00000000_dp * cst%m_u !O2
    var%m(11) = 33.00000000_dp * cst%m_u !HO2
    var%m(12) = 34.00000000_dp * cst%m_u !H2O2
    var%m(13) = 10.00000000_dp * cst%m_u !M
    var%m(14) = 28.00000000_dp * cst%m_u !N2
    var%m(15) = 45.00000000_dp * cst%m_u !HOCO
    var%m(16) = 43.99945142_dp * cst%m_u !CO2+
    var%m(17) = 29.00000000_dp * cst%m_u !HCO
    var%m(18) = 30.00000000_dp * cst%m_u !H2CO

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
    var%q(1) = 0.0_dp * cst%q_e !CO2
    var%q(2) = 0.0_dp * cst%q_e !CO
    var%q(3) = 0.0_dp * cst%q_e !O
    var%q(4) = 0.0_dp * cst%q_e !O(1D)
    var%q(5) = 0.0_dp * cst%q_e !H2O
    var%q(6) = 0.0_dp * cst%q_e !H
    var%q(7) = 0.0_dp * cst%q_e !OH
    var%q(8) = 0.0_dp * cst%q_e !H2
    var%q(9) = 0.0_dp * cst%q_e !O3
    var%q(10) = 0.0_dp * cst%q_e !O2
    var%q(11) = 0.0_dp * cst%q_e !HO2
    var%q(12) = 0.0_dp * cst%q_e !H2O2
    var%q(13) = 0.0_dp * cst%q_e !M
    var%q(14) = 0.0_dp * cst%q_e !N2
    var%q(15) = 0.0_dp * cst%q_e !HOCO
    var%q(16) = 1.0_dp * cst%q_e !CO2+
    var%q(17) = 0.0_dp * cst%q_e !HCO
    var%q(18) = 0.0_dp * cst%q_e !H2CO

    ! read P, L, J list
    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/Production_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/Loss_list.dat', status = 'unknown' )
      do isp = 1, spl%nsp_i
        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/Jacobian_list.dat', status = 'unknown' )
      do i = 1, spl%n_Jlist
        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/reactant_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/product_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/reaction_type_list.dat', status = 'unknown' )
      do ich = 1, spl%nch
        read(11,*) spl%reaction_type_list(ich)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/rate_rpn_token.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/rate_rpn_label.dat', status = 'unknown' )
      do ich = 1, spl%nch
        do i = 1, 3
          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)
        end do
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/rate_cases.dat', status = 'unknown' )
      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/PLJ_list/T_range.dat', status = 'unknown' )
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

    open(11, file = './Mars/Neutral_atmosphere/input/Temperature/T_e.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Te(iz)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/Temperature/T_i.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Ti(iz)
      end do
    close(11)

    open(11, file = './Mars/Neutral_atmosphere/input/Temperature/T_n.dat', status = 'unknown' )
      do iz = 1, grd%nz
        read(11,*) tmp,var%Tn(iz)
      end do
    close(11)

    ! input density profiles
    var%ni   = 1.0e-20_dp

    isp = sp_index(spl, 'CO2')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Neutral_atmosphere/input/density/CO2.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'H2O')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Neutral_atmosphere/input/density/H2O.dat', status = 'unknown' )
        do iz = 1, grd%nz
          read(11,*) tmp, var%ni(isp,iz)
        end do
      close(11)
    end if

    isp = sp_index(spl, 'CO2+')
    if (isp >= 1 .and. isp <= spl%nsp) then
      open(11, file = './Mars/Neutral_atmosphere/input/density/CO2+.dat', status = 'unknown' )
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
    var%LowerBC(isp,2) = 2.1e23_dp

    isp = sp_index(spl, 'N2')
    var%LowerBC(isp,1) = 1.0_dp
    var%LowerBC(isp,2) = 4.0e21_dp


    ! Upper boundary condition
    var%UpperBC = 0.0_dp

    isp = sp_index(spl, 'O')
    var%UpperBC(isp,1) = 2.0_dp
    var%UpperBC(isp,2) = 1.2e12_dp

    isp = sp_index(spl, 'H')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape

    isp = sp_index(spl, 'H2')
    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape


  end subroutine v__in__exe

end module v__in
