module wrapper_PROTEUS

  ! This module serves as a wrapper for the PROTEUS photochemistry modules, 
  ! providing an interface to be called from an external model. 
  !
  ! Required fortran modules to use this wrapper:
  ! - p__photolysis_rate : This module calculates photodissociation rates.
  ! - p__photoionization_rate : This module calculates photoionization rates.
  ! - p__PROTEUS : This module is output from the PROTEUS GUI app. 
  !               Please select reactions that you want to include in your model 
  !               and generate this module using the GUI app.
  ! 
  ! Yuki Nakamura, 10 February 2026

  use p__photolysis_rate
  use p__photoionization_rate
  use p__PROTEUS

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private

  ! Module-level variables --------------------------------------------------------------
  integer,  private :: nx, ny, nz, nwl, nwl_EUV
  integer,  private :: nsp_v, nsp_f, nsp, nch
  integer,  private :: nsp_v_ex_sav, nsp_f_ex_sav
  integer,  allocatable, private :: reaction_type(:), product_list(:,:), reactant_list(:,:)
  real(dp), allocatable, private :: lambda(:), dlambda(:), solar_flux(:)
  real(dp), allocatable, private :: mass(:), J_rate_in(:,:,:,:), density(:,:)
  real(dp)          :: rplanet, mplanet
  character(len=256), allocatable, private :: species_list(:), species_var(:), species_fix(:)
  character(len=256), allocatable, private :: species_list_ex_sav(:), species_var_ex_sav(:), species_fix_ex_sav(:)
  
  ! Public interface --------------------------------------------------------------------
  public :: wrapper_PROTEUS__ini, wrapper_PROTEUS__exe

contains
  
  subroutine wrapper_PROTEUS__ini(flag,                           & ! in:  0: load all, 1: load UV only, 2: load EUV only
    &                             rplanet_ex, mplanet_ex,         & ! in:  Radius [m] and mass [kg] of the planet
    &                             nx_ex, ny_ex, nz_ex,            & ! in:  Number of grid points in x, y, z direction. Note that +z directs vertically upward.
    &                             nsp_v_ex, nsp_f_ex,             & ! in:  Number of variable and fixed species in the external model.
    &                             species_var_ex, species_fix_ex, & ! in:  List of variable and fixed species in the external model. The name of species should be exactly the same as the name used in PROTEUS (e.g., "CO2", "O", "CO2+", "O2+", etc.).
    &                             dir_xsect_dis, dir_xsect_ion    ) ! in:  Directory path to the cross section data files for photolysis and photoionization
    implicit none
    integer,  intent(in) :: flag 
    real(dp), intent(in) :: rplanet_ex, mplanet_ex 
    integer,  intent(in) :: nx_ex, ny_ex, nz_ex 
    integer,  intent(in) :: nsp_v_ex, nsp_f_ex
    character(len=256), intent(in) :: species_var_ex(:), species_fix_ex(:)
    character(len=256), intent(in) :: dir_xsect_dis, dir_xsect_ion
    integer iwl, i

    !
    ! set planetary parameters
    !
    rplanet = rplanet_ex 
    mplanet = mplanet_ex

    !
    ! set grid size
    !
    nx = nx_ex
    ny = ny_ex
    nz = nz_ex

    !
    ! set species list from external model
    !
    nsp_v_ex_sav = nsp_v_ex
    nsp_f_ex_sav = nsp_f_ex
    species_var_ex_sav = species_var_ex
    species_fix_ex_sav = species_fix_ex
    
    !
    ! get number of species, reactions, and wavelength bins
    !
    nsp = get_number_of_all_species()
    nsp_v = get_number_of_var_species()
    nsp_f = nsp - nsp_v
    nch = get_number_of_reaction()
    nwl = get_number_of_wavelength_bins()

    allocate(reaction_type(nch))
    allocate(product_list(nch,0:20), reactant_list(nch,0:20))
    allocate(species_list(nsp), species_var(nsp_v), species_fix(nsp-nsp_v))
    allocate(mass(nsp))
    allocate(lambda(nwl), dlambda(nwl), solar_flux(nwl))
    allocate(J_rate_in(nx,ny,nz,nch))
    allocate(density(nz,nsp))

    ! 
    ! Get all necessary data
    !
    species_list(1:nsp) = get_all_species_name()
    species_var(1:nsp_v) = get_var_species_name()
    species_fix(1:nsp_f) = get_fix_species_name()
    mass(1:nsp) = get_mass_of_species()
    call get_reaction_type_list(reaction_type)
    call get_reactant_product_list(reactant_list, product_list)
    call get_wavelength_bin(lambda, dlambda)
    call get_solar_flux(lambda, dlambda, solar_flux)

    ! for EUV wagelength bin limit
    do iwl = 1, nwl
      if (lambda(iwl) > 105.0_dp) then
        nwl_EUV = iwl - 1
        exit
      end if
    end do

    !
    ! Initialize variables & Load cross section data
    !
    if (flag==0 .or. flag==1) then
      call p__photolysis_rate__ini(nwl, nz, nsp, nch,                                              & ! in
        &                          species_list, mass, reaction_type, reactant_list, product_list, & ! in
        &                          lambda, dlambda, solar_flux                                     ) ! in
      call load_cross_section_dat(dir_xsect_dis) ! in
    end if
    if (flag==0 .or. flag==2) then
      call p__photoionization_rate__ini(nwl_EUV, nz, nsp, nch,                                          & ! in
        &                               species_list, mass, reaction_type, reactant_list, product_list, & ! in
        &                               lambda(1:nwl_EUV), dlambda(1:nwl_EUV), solar_flux(1:nwl_EUV)    ) ! in
      call load_EUV_cross_section_dat(dir_xsect_ion) ! in
    end if

  end subroutine wrapper_PROTEUS__ini


  subroutine wrapper_PROTEUS__exe(flag,            & ! in:  0: calculate all, 1: photolysis only, 2: photoionization only
    &                             alt_grid,        & ! in:  grid of altitude [m] with dimension (nx,ny,nz)
    &                             T_n, T_i, T_e,   & ! in:  temperature of neutrals and electrons (nx,ny,nz), and 4D array of ion temperature (nx,ny,nz,nsp) [K]
    &                             v_var, v_fix,    & ! in:  three-dimensional velocity vector of variable and fixed species with dimension (3,nx,ny,nz,nsp_v) and (3,nx,ny,nz,nsp_f) [m/s]. Note that the last dimension corresponds to x, y, z component of velocity vector.
    &                             n_var, n_fix,    & ! in:  number density of variable and fixed species with dimension (nx,ny,nz,nsp_v) and (nx,ny,nz,nsp_f) [cm^-3]
    &                             sza_array,       & ! in:  solar zenith angle with dimension (nx,ny) [rad]
    &                             prod, loss, kout ) ! out: production and loss rates for variable species with dimension (nx,ny,nz,nsp_v) [cm^-3 s^-1] and rate coefficient for each reaction with dimension (nx,ny,nz,nch)
    implicit none
    integer, intent(in)     :: flag 
    real(dp), intent(in)    :: alt_grid(:,:,:)
    real(dp), intent(in)    :: T_n(:,:,:), T_i(:,:,:,:), T_e(:,:,:)
    real(dp), intent(in)    :: v_var(:,:,:,:,:), v_fix(:,:,:,:,:)
    real(dp), intent(in)    :: n_var(:,:,:,:), n_fix(:,:,:,:)
    real(dp), intent(in)    :: sza_array(:,:)
    real(dp), intent(out)   :: prod(:,:,:,:), loss(:,:,:,:), kout(:,:,:,:)

    integer ix, iy, iz, isp, jsp
    real(dp) T_z(nz), alt(nz), dalt(nz), ntot(nz), sza
    
    do ix = 1, nx
    do iy = 1, ny

      sza = sza_array(ix,iy)

      T_z(:) = T_n(ix,iy,:)
      alt(:) = alt_grid(ix,iy,:)
      do iz = 2, nz
        dalt(iz) = alt(iz) - alt(iz-1)
      end do
      dalt(1) = dalt(2)

      density(:,:) = 0.0_dp
      ntot(:) = 0.0_dp
      do isp = 1, nsp_v_ex_sav
        do jsp = 1, nsp
          if (trim(ADJUSTL(species_var_ex_sav(isp)))==trim(ADJUSTL(species_list(jsp)))) then 
            density(:,jsp) = n_var(ix,iy,1:,isp) * 1.0e6_dp
            if (trim(ADJUSTL(species_list(jsp))) /= 'M') ntot(:) = ntot(:) + density(:,jsp)
          end if
        end do 
      end do 
      do isp = 1, nsp_f_ex_sav
        do jsp = 1, nsp
          if (trim(ADJUSTL(species_fix_ex_sav(isp)))==trim(ADJUSTL(species_list(jsp)))) then 
            density(:,jsp) = n_fix(ix,iy,1:,isp) * 1.0e6_dp
            if (trim(ADJUSTL(species_list(jsp))) /= 'M') ntot(:) = ntot(:) + density(:,jsp)
          end if
        end do 
      end do


      !
      ! Calculate cross sections
      !
      if (flag==0 .or. flag==1) then
        call get_cross_section(T_z, ntot, 'absorption', 0) ! in
        call get_cross_section(T_z, ntot, 'photolysis', 0) ! in
      end if
      if (flag==0 .or. flag==2) then
        call get_EUV_cross_section('absorption', 0) ! in
        call get_EUV_cross_section('photoionization', 0) ! in
      end if
      
      !
      ! reset J_rate_in
      !
      J_rate_in(ix,iy,1:,1:) = 0.0_dp

      !
      ! Calculate photolysis rate
      !
      if (flag==0 .or. flag==1) then
        call photolysis_rate(alt, dalt, density, T_z, sza, & ! in
        &                    rplanet, mplanet,             & ! in
        &                    J_rate_in(ix,iy,1:,1:)        ) ! inout: note that the calculated photolysis rate is added to J_rate_in
      end if

      !
      ! Calculate photoionization rate
      !
      if (flag==0 .or. flag==2) then
        call photoionization_rate(alt, dalt, density, T_z, sza, & ! in
          &                       rplanet, mplanet,             & ! in
          &                       J_rate_in(ix,iy,1:,1:)        ) ! inout: note that the calculated photoionization rate is added to J_rate_in
      end if

    end do 
    end do 

    !
    ! Calculate reaction rate
    !
    prod = 0.0_dp
    loss = 0.0_dp
    kout = 0.0_dp

    call p__PROTEUS_source(species_var, species_fix, & ! in:  name of variable species and fixed species
      &                    nsp_v, nsp_f,             & ! in:  number of variable species and fixed species
      &                    nx, ny, nz,               & ! in:  number of grids in x, y, z direction
      &                    T_n, T_i, T_e,            & ! in:  temperature of neutrals, ions and electrons
      &                    v_var, v_fix,             & ! in:  three dimensional velocity vector of variable and fixed species
      &                    n_var, n_fix,             & ! in:  number density of variable and fixed species
      &                    J_rate_in,                & ! in:  photoionization / dissociation reaction list
      &                    prod, loss, kout          ) ! out: production and loss rates for variable species

  end subroutine wrapper_PROTEUS__exe


end module wrapper_PROTEUS