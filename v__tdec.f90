module v__tdec
  !---------------------------------------------------------------------
  ! Type declaration statements
  !---------------------------------------------------------------------
  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: grd_, var_, cst_, xct_, spl_, flx_, set_

  ! settings
  type set_
    integer(sp)            :: nx, ny, resx, resy            ! number of horizontal grids, resolution
    real(dp)                  latitude, sza, Ls, DOY, F107  ! latitude, solar zenith angle, solar longitude, day of year, F10.7
    real(dp)                  fin_sec, dtime_limit          ! Calculation time, upper limit of delta t
    integer                   nday, nlat, nstep             
    integer                   calc_stable, start_rot, test_loc, read_stable
    character(len=256)        scheme, inversion, mode
    character(len=256)        fnamestable, dir_name
  end type set_

  ! grid
  type grd_
    integer(sp)            :: nx, ny, nz               ! number of grid index
    integer(sp)            :: xs, xe, ys, ye, zs, ze   ! grid index range
    integer(sp)            :: ix, iy, iday
    real(dp), allocatable  :: dalt(:), alt(:), lt(:), lat(:), sza(:,:), sza_xact(:,:)   ! altitude, longitude, latitude
    real(dp)               :: tilt, season, latitude ! solar zenith angle, tilt angle, season
  end type grd_

  ! variables
  type var_
    ! 種類ごとに分けて欲しい
    !-------------------------------------------------------------------
    ! nstep     : total number of time steps
    ! istep     : index of time step
    ! species_i : species for searching cross sections
    ! reactants : reactant list
    ! products  : product list
    ! m         : mass of species
    ! m_major   : mass of major species
    ! q         : charge of species
    !
    !-------------------------------------------------------------------
    integer                :: nsteps, istep, nspecial, iter ! number and index of time steps
    real(dp)               :: t1, t2, t3, t4
    character(len=256)     :: species_i, reactants(10), products(10)
    real(dp), allocatable  :: m(:), m_mean(:), q(:)
    real(dp), allocatable  :: ni(:,:), ni_stable(:,:,:), ni_3d(:,:,:,:), ni_new(:,:), n_tot(:)   ! density
    real(dp), allocatable  :: ni_0(:,:)
    real(dp), allocatable  :: clm_ni(:,:) ! column density
    real(dp), allocatable  :: Ti(:), Te(:), Tn(:)
    real(dp), allocatable  :: Ti_3d(:,:,:), Te_3d(:,:,:), Tn_3d(:,:,:)
    real(dp), allocatable  :: tau_EUV(:,:), tau_UV(:,:), tau_EUV_subsolar(:,:), tau_UV_subsolar(:,:), tau_RS(:,:)
    real(dp), allocatable  :: E_fld(:,:,:), B_fld(:,:,:) ! electric and magnetic field
    real(dp)               :: dtime, sum_time
    real(dp), allocatable  :: ki(:,:), ki_special(:,:,:,:), ich_special(:), Pi(:,:), Li(:,:), Jmtx(:,:), rate(:,:), Pij(:,:,:)
    real(dp), allocatable  :: Fluxup(:,:), Fluxdwn(:,:), vFluxup(:,:), vFluxdwn(:,:)
    real(dp), allocatable  :: UpperBC(:,:), LowerBC(:,:) ! upper and lower boundary condition: label, (1:density, 2:flux, 3:velocity)
    real(dp), allocatable  :: K_eddy(:), D_mol(:,:)
    real(dp), allocatable  :: I_EUV(:,:), I_UV(:,:)
    real(dp), allocatable  :: Phip(:,:), Phim(:,:) ! 
    real(dp), allocatable  :: dPhi_dz(:,:)         ! dPhi/dz
    real(dp), allocatable  :: d_dniu_dPhi_dz(:,:)  ! d/dniu * dPhi/dz
    real(dp), allocatable  :: d_dni0_dPhi_dz(:,:)  ! d/dni0 * dPhi/dz
    real(dp), allocatable  :: d_dnil_dPhi_dz(:,:)  ! d/dnil * dPhi/dz
    real(dp), allocatable  :: d_dneu_dPhi_dz_add(:,:)  ! d/dneu * dPhi/dz
    real(dp), allocatable  :: d_dne0_dPhi_dz_add(:,:)  ! d/dne0 * dPhi/dz
    real(dp), allocatable  :: d_dnel_dPhi_dz_add(:,:)  ! d/dnel * dPhi/dz
    real(dp), allocatable  :: tAmtx(:,:), tLmtx(:,:), Umtx(:,:)
    real(dp), allocatable  :: barr(:), xarr(:), yarr(:), dxarr(:), rarr(:)
    real(dp)               :: max_dn_n(3), i_dn_n
  end type var_

  type LU_
    real(dp), allocatable :: tLmtx(:,:), Umtx(:,:), LUmtx
  end type LU_

  ! species list and variables associated with chemical reactions
  type spl_
    !-------------------------------------------------------------------
    ! planet        : selected planet e.g. Mars, Jupiter, ...
    ! nsp           : total number of species
    ! nsp_i         : total number of variable species
    ! nsp_f         : total number of fixed species
    ! nch           : total number of chemical reactions
    ! nch_P         : total number of columns for reading Production list
    ! nch_L         : total number of columns for reading Loss list
    ! nch_J         : total number of columns for reading Jacobian list
    ! n_Jlist       : total number of lines for reading Jacobian list
    !
    ! species       : species list
    ! reactant_list :
    !
    !-------------------------------------------------------------------
    character(len=256)     :: planet
    integer                :: nsp, nsp_i, nch, nch_P, nch_L, nch_J, n_Jlist, nrpn
    character(len=256), allocatable :: species(:), reaction_type_char(:)
    integer,  allocatable  :: reactant_list(:,:), product_list(:,:), Prod_list(:,:), Loss_list(:,:)
    integer,  allocatable  :: reaction_type_list(:), label_fix(:), all_to_var(:), var_to_all(:), major_species(:)
    integer,  allocatable  :: Jmtx_list(:,:)
    integer,  allocatable  :: rate_rpn_label(:,:,:)
    real(dp), allocatable  :: rate_rpn_token(:,:,:), T_range(:,:,:)
    integer,  allocatable  :: rate_cases(:)
  end type spl_

  ! physical constant
  type cst_
    real(dp)               :: pi, c, h
    real(dp)               :: NA  ! Avogadro constant
    real(dp)               :: k_B ! Boltzmann constant
    real(dp)               :: g   ! gravitational acceleration
    real(dp)               :: BigG ! gravitational constant
    real(dp)               :: Mplanet ! gravitational acceleration
    real(dp)               :: R   ! planetary radius
    real(dp)               :: m_u ! atomic mass unit
    real(dp)               :: q_e
    real(dp)               :: eV  ! electron volt
    real(dp)               :: planck 
    real(dp)               :: daysec
    real(dp)               :: R_gas ! Gas constant
  end type cst_

  ! cross sections
  type xct_
    character(len=256)    :: type          ! type of cross section
    integer,  allocatable :: label_sigma_a_UV(:), label_sigma_a_EUV(:), label_sigma_a_UV_EUV_conv(:)
    integer,  allocatable :: label_sigma_a_RS(:)
    real(dp), allocatable :: sigma_a_EUV(:,:)    ! EUV absorption cross section
    real(dp), allocatable :: sigma_i_EUV(:,:)    ! EUV ionization cross section
    real(dp), allocatable :: sigma_a_UV(:,:,:)   ! UV absorption cross section
    real(dp), allocatable :: sigma_d_UV(:,:,:)   ! UV dissociation cross section
    integer, allocatable  :: wlrange_a_UV(:,:), wlrange_d_UV(:,:) ! wavelength index range of absorption & dissociation cross sections
    real(dp), allocatable :: sigma_e_UV(:) ! UV excitation cross section
    real(dp), allocatable :: sigma_a_UV_EUV(:,:)    ! EUV absorption cross section
    real(dp), allocatable :: sigma_a_RS(:,:)     ! Rayleigh scattering cross section

    ! CO2
    real(dp), allocatable :: sigma_a_CO2(:,:), sigma_a_u13CO2(:,:)
    real(dp), allocatable :: sigma_d_CO2_O1D(:,:)
    real(dp), allocatable :: T_x_CO2(:)
    real(dp), allocatable :: qy_CO2_O(:,:)
    ! H2O
    real(dp), allocatable :: sigma_a_H2O(:,:)
    real(dp), allocatable :: sigma_d_H2O_OH(:,:), sigma_d_H2O_O1D(:,:)
    real(dp), allocatable :: qy_H2O_OH(:,:), qy_H2O_O1D(:,:)
    ! H2O2 
    real(dp), allocatable :: sigma_a_H2O2(:,:)
    real(dp), allocatable :: qy_H2O2_OH(:,:), qy_H2O2_HO2(:,:)
    real(dp), allocatable :: sigma_d_H2O2_OH(:,:), sigma_d_H2O2_HO2(:,:)
    ! H2
    real(dp), allocatable :: sigma_a_H2(:,:)
    real(dp), allocatable :: qy_H2_H(:,:)
    real(dp), allocatable :: sigma_d_H2_H(:,:)
    ! OH
    real(dp), allocatable :: sigma_a_OH(:,:)
    real(dp), allocatable :: sigma_d_OH_O1D(:,:), sigma_d_OH_O(:,:)
    real(dp), allocatable :: qy_OH_O(:,:)
    ! HO2
    real(dp), allocatable :: sigma_a_HO2(:,:)
    real(dp), allocatable :: qy_HO2_OH(:,:)
    real(dp), allocatable :: sigma_d_HO2_OH(:,:)
    ! O3
    real(dp), allocatable :: sigma_a_O3(:,:)
    real(dp), allocatable :: T_x_O3(:)
    real(dp), allocatable :: qy_O3_O1D(:,:), T_qy_O3_O1D(:), qy_O3_O(:,:), T_qy_O3_O(:)
    ! O2
    real(dp), allocatable :: sigma_a_O2(:,:)
    real(dp), allocatable :: T_x_O2(:)
    real(dp), allocatable :: qy_O2_O(:,:), qy_O2_O1D(:,:)
    real(dp), allocatable :: sigma_a_O2_SR_bands_130K(:,:), sigma_a_O2_SR_bands_190K(:,:), sigma_a_O2_SR_bands_280K(:,:)

    ! H2CO
    real(dp), allocatable :: sigma_a_H2CO(:,:), sigma_gamma_H2CO(:,:)
    real(dp), allocatable :: qy_H2CO_CO_300K(:,:)

    ! N2
    real(dp), allocatable :: sigma_a_N2(:,:)
    real(dp), allocatable :: sigma_d_N2_N(:,:)

    ! NO
    real(dp), allocatable :: sigma_a_NO(:,:)
    real(dp), allocatable :: sigma_d_NO_N(:,:)

    ! NO2
    real(dp), allocatable :: sigma_a_NO2(:,:)
    real(dp), allocatable :: T_x_NO2(:)
    real(dp), allocatable :: qy_NO2_O(:,:)
    real(dp), allocatable :: T_qy_NO2(:)
    real(dp), allocatable :: sigma_d_NO2_O1D(:,:)

    ! NO3
    real(dp), allocatable :: sigma_a_NO3(:,:)
    real(dp), allocatable :: T_x_NO3(:)
    real(dp), allocatable :: qy_NO3_NO2(:,:), qy_NO3_NO(:,:)
    real(dp), allocatable :: T_qy_NO3(:)

    ! N2O
    real(dp), allocatable :: sigma_a_N2O(:,:)
    real(dp), allocatable :: T_x_N2O(:)
    real(dp), allocatable :: qy_N2O_O1D(:,:)

    ! N2O5
    real(dp), allocatable :: sigma_a_N2O5(:,:), sigma_a_N2O5_above260nm(:,:)
    real(dp), allocatable :: qy_N2O5_NO2(:,:), qy_N2O5_O(:,:)

    ! HONO
    real(dp), allocatable :: sigma_a_HONO(:,:)
    real(dp), allocatable :: qy_HONO_OH(:,:)

    ! HNO3
    real(dp), allocatable :: sigma_a_HNO3(:,:),sigma_a_HNO3_Tdep(:,:)
    real(dp), allocatable :: qy_HNO3_OH(:,:), qy_HNO3_O(:,:), qy_HNO3_O1D(:,:)

    ! HONO
    real(dp), allocatable :: sigma_a_HO2NO2(:,:), sigma_a_HO2NO2_Tdep(:,:)
    real(dp), allocatable :: qy_HO2NO2_OH(:,:), qy_HO2NO2_HO2(:,:)

    ! HCO
    real(dp), allocatable :: sigma_a_HCO(:,:)
    real(dp), allocatable :: qy_HCO_H(:,:)

    ! CH4
    real(dp), allocatable :: sigma_a_CH4(:,:)
    real(dp), allocatable :: qy_CH4_1CH2(:,:), qy_CH4_3CH2(:,:), qy_CH4_CH3(:,:)
    real(dp), allocatable :: sigma_d_CH4_1CH2(:,:), sigma_d_CH4_3CH2(:,:), sigma_d_CH4_CH3(:,:)

    ! C2H6
    real(dp), allocatable :: sigma_a_C2H6(:,:)
    real(dp), allocatable :: qy_C2H6_3CH2(:,:), qy_C2H6_1CH2(:,:), qy_C2H6_C2H2(:,:), qy_C2H6_C2H4_H(:,:), qy_C2H6_C2H4_H2(:,:)
    real(dp), allocatable :: qy_C2H6_CH3(:,:)

    ! HNO
    real(dp), allocatable :: sigma_a_HNO(:,:)
    real(dp), allocatable :: qy_HNO_NO(:,:)

    ! HNO2
    real(dp), allocatable :: sigma_a_HNO2(:,:)
    real(dp), allocatable :: qy_HNO2_NO(:,:)

    ! CH3
    real(dp), allocatable :: sigma_a_CH3(:,:)
    real(dp), allocatable :: qy_CH3_1CH2(:,:)

    ! NH3
    real(dp), allocatable :: sigma_a_NH3(:,:)
    real(dp), allocatable :: qy_NH3_NH2(:,:)

    ! N2H4
    real(dp), allocatable :: sigma_a_N2H4(:,:)
    real(dp), allocatable :: qy_N2H4_N2H3(:,:)

    ! NH
    real(dp), allocatable :: sigma_a_NH(:,:)
    real(dp), allocatable :: qy_NH_N(:,:)

    ! NH2
    real(dp), allocatable :: sigma_a_NH2(:,:)
    real(dp), allocatable :: qy_NH2_NH(:,:)

    ! C2H2
    real(dp), allocatable :: sigma_a_C2H2(:,:)
    real(dp), allocatable :: qy_C2H2_C2H(:,:), qy_C2H2_C2(:,:)
    real(dp), allocatable :: sigma_d_C2H2_C2H(:,:), sigma_d_C2H2_C2(:,:)

    ! C2H4
    real(dp), allocatable :: sigma_a_C2H4(:,:)
    real(dp), allocatable :: qy_C2H4_C2H2_H2(:,:), qy_C2H4_C2H2_H(:,:)
    real(dp), allocatable :: sigma_d_C2H4_C2H2_H2(:,:), sigma_d_C2H4_C2H2_H(:,:)

    ! C3H8
    real(dp), allocatable :: sigma_a_C3H8(:,:)
    real(dp), allocatable :: qy_C3H8_C3H6(:,:), qy_C3H8_C2H6(:,:), qy_C3H8_C2H4(:,:), qy_C3H8_C2H5(:,:)

    ! C3H6
    real(dp), allocatable :: sigma_a_C3H6(:,:)
    real(dp), allocatable :: qy_C3H6_C2H2(:,:), qy_C3H6_CH2CCH2(:,:), qy_C3H6_C2H4(:,:), qy_C3H6_C2H(:,:)
    
    ! CH
    real(dp), allocatable :: sigma_a_CH(:,:)
    real(dp), allocatable :: qy_CH_C(:,:)

    ! CH2CO
    real(dp), allocatable :: sigma_a_CH2CO(:,:)
    real(dp), allocatable :: qy_CH2CO_CH2(:,:)

    ! CH3CHO
    real(dp), allocatable :: sigma_a_CH3CHO(:,:)
    real(dp), allocatable :: qy_CH3CHO_CH3(:,:), qy_CH3CHO_CH4(:,:)

    ! C2H5CHO
    real(dp), allocatable :: sigma_a_C2H5CHO(:,:)
    real(dp), allocatable :: qy_C2H5CHO_C2H5(:,:)

    ! C3H3
    real(dp), allocatable :: sigma_a_C3H3(:,:)
    real(dp), allocatable :: qy_C3H3_C3H2(:,:)

    ! CH3C2H
    real(dp), allocatable :: sigma_a_CH3C2H(:,:)
    real(dp), allocatable :: qy_CH3C2H_C3H3(:,:), qy_CH3C2H_C3H2(:,:), qy_CH3C2H_CH3(:,:)

    ! CH2CCH2
    real(dp), allocatable :: sigma_a_CH2CCH2(:,:)
    real(dp), allocatable :: qy_CH2CCH2_C3H3(:,:), qy_CH2CCH2_C3H2(:,:), qy_CH2CCH2_C2H2(:,:)
       
    ! HCN
    real(dp), allocatable :: sigma_a_HCN(:,:)
    real(dp), allocatable :: qy_HCN_H(:,:)

    ! HNCO
    real(dp), allocatable :: sigma_a_HNCO(:,:)
    real(dp), allocatable :: qy_HNCO_H(:,:)

    ! NCO
    real(dp), allocatable :: sigma_a_NCO(:,:)
    real(dp), allocatable :: qy_NCO_N(:,:)     

  end type xct_

  ! solar flux
  type flx_
    real(dp)              :: orbit ! orbit radius [AU]
    real(dp)              :: mode_factor

    integer               :: nwl_EUV
    real(dp), allocatable :: lambda_EUV(:), F74113(:), Ai(:), F107, solar_EUV(:)
    real(dp)              :: multiplying_factor_EUV

    integer               :: nwl_UV
    real(dp), allocatable :: lambda_UV(:), dlambda_UV(:), solar_UV(:)
    real(dp)              :: multiplying_factor_UV

  end type flx_


end module v__tdec
