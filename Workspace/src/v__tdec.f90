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
    real(dp)                  fin_sec, dtime_limit, dt_out, dt_rate, dt_inc_eps, max_eps  ! Calculation time, upper limit of delta t
    real(dp)                  lat_res, lt_res
    integer                   use_1d, use_2d
    integer                   nday, nlat, nstep             
    integer                   calc_stable, start_rot, fix_sza, read_stable, diurnal_ave
    real(dp)                  euv_factor
    character(len=256)        scheme, inversion, mode, euv_input
    character(len=256)        fnamestable, dir_name, dir_name_win, solar_flux, solar_flux_unit1, solar_flux_unit2
    integer                   n_wl_bin, nsp_tout, rate_tout
    real(dp), allocatable  :: wl_bin(:,:), rate_from_datafile(:,:)
    character(len=256), allocatable :: species_tout(:)
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
    integer                :: nsteps, istep, nspecial, iter ! number and index of time steps
    real(dp)               :: t1, t2, t3, t4
    character(len=256)     :: species_i, reactants(10), products(10)
    real(dp), allocatable  :: m(:), m_mean(:), q(:)
    real(dp), allocatable  :: ni(:,:), ni_stable(:,:,:), ni_3d(:,:,:,:), ni_new(:,:), n_tot(:)   ! density
    real(dp), allocatable  :: ni_0(:,:)
    real(dp), allocatable  :: clm_ni(:,:) ! column density
    real(dp), allocatable  :: Ti(:), Te(:), Tn(:)
    real(dp), allocatable  :: Ti_3d(:,:,:), Te_3d(:,:,:), Tn_3d(:,:,:)
    real(dp), allocatable  :: tau_EUV(:,:), tau_RS(:,:)
    real(dp), allocatable  :: E_fld(:,:,:), B_fld(:,:,:) ! electric and magnetic field
    real(dp)               :: dtime, sum_time
    real(dp), allocatable  :: ki(:,:), ki_special(:,:,:,:), ich_special(:)
    real(dp), allocatable  :: Pi(:,:), Li(:,:), Jmtx(:,:), rate(:,:), Pij(:,:,:), Lij(:,:,:)
    real(dp), allocatable  :: Fluxup(:,:), Fluxdwn(:,:), vFluxup(:,:), vFluxdwn(:,:)
    real(dp), allocatable  :: Upper_n(:,:), Upper_f(:,:), Upper_v(:,:), Lower_n(:,:), Lower_f(:,:), Lower_v(:,:) ! upper and lower boundary condition: label, (density, flux, velocity)
    real(dp), allocatable  :: K_eddy(:), D_mol(:,:), nu(:,:)
    real(dp), allocatable  :: I_EUV(:,:)
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
    character(len=256)     :: planet
    integer                :: nsp, nsp_i, nch
    character(len=256), allocatable :: species(:), reaction_type_char(:)
    integer,  allocatable  :: reactant_list(:,:), product_list(:,:)
    integer,  allocatable  :: reaction_type_list(:), label_fix(:), all_to_var(:), var_to_all(:), major_species(:)
    real(dp), allocatable  :: T_range(:,:,:)
    integer,  allocatable  :: rate_cases(:)
  end type spl_

  ! physical constant
  type cst_
    real(dp)               :: pi, c, h ! pi, speed of light, Planck constant
    real(dp)               :: NA  ! Avogadro constant
    real(dp)               :: k_B ! Boltzmann constant
    real(dp)               :: g   ! gravitational acceleration
    real(dp)               :: BigG ! gravitational constant
    real(dp)               :: Mplanet ! gravitational acceleration
    real(dp)               :: R   ! planetary radius
    real(dp)               :: m_u ! atomic mass unit
    real(dp)               :: q_e ! elementary charge
    real(dp)               :: eV  ! electron volt
    real(dp)               :: daysec ! seconds in a day
    real(dp)               :: R_gas ! Gas constant
  end type cst_

  ! cross sections
  type xct_
    character(len=256)    :: type          ! type of cross section
    integer,  allocatable :: label_sigma_a_EUV(:)
    real(dp), allocatable :: sigma_a_EUV(:,:)    ! EUVAC absorption cross section
    real(dp), allocatable :: sigma_i_EUV(:,:)    ! EUVAC ionization cross section
  end type xct_

  ! solar flux
  type flx_
    real(dp)              :: orbit ! orbit radius [AU]
    real(dp)              :: irradiance_factor

    integer               :: nwl_EUV
    real(dp), allocatable :: lambda_EUV(:), dlambda_EUV(:), F107, solar_EUV(:)
    real(dp), allocatable :: F74113(:), Ai(:)
    real(dp)              :: multiplying_factor_EUV

    integer               :: nwl_UV
    real(dp), allocatable :: lambda_UV(:), dlambda_UV(:), solar_UV(:)
    real(dp)              :: multiplying_factor_UV

  end type flx_


end module v__tdec
