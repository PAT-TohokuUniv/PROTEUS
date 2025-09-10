module p__photolysis_rate

  implicit none
  integer(4), parameter :: sp = 4, dp = 8


  private

  ! Module-level variables ---------------------------------------------------------------------------
  integer,                         private :: nch_photolysis, nwl, nsp, nz, nch
  character(len=256), allocatable, private :: species(:)
  real(dp),           allocatable, private :: sigma_dat(:,:,:,:,:), sigma_a(:,:,:), sigma_d(:,:,:)
  real(dp),           allocatable, private :: lambda(:), dlambda(:), solar_flux(:), mass(:)
  integer,            allocatable, private :: reaction_type(:), product_list(:,:), reactant_list(:,:)
  real(dp),           allocatable, private :: Tdep_list(:)

  ! For exceptional cases of cross section calculation
  real(dp),           allocatable, private :: larr(:,:)
  real(dp),           allocatable, private :: qy_H(:), qy_CO_300K(:), qy_CO_T(:), a_300K(:), a_T(:)
  real(dp),           allocatable, private :: qy_O1D(:)

  ! For optical depth calculation 
  real(dp),           allocatable, private :: cln(:,:), tau(:,:), I_z(:,:)

  ! Public interfaces --------------------------------------------------------------------------------
  public :: p__photolysis_rate__ini, p__photolysis_rate__fin, &
    &       load_cross_section_dat, get_cross_section, photolysis_rate

contains


  !=======================================================================================================================
  !
  !  Allocate memory for the cross-section data array.
  !
  !=======================================================================================================================
  subroutine p__photolysis_rate__ini(nwl_in, nz_in, nsp_in, nch_in,                                            & ! in
    &                                species_in, mass_in, reaction_type_in, reactant_list_in, product_list_in, & ! in
    &                                lambda_in, dlambda_in, solar_flux_in                                      ) ! in
    implicit none
    integer,  intent(in) :: nwl_in, nsp_in, nz_in, nch_in
    integer,  intent(in) :: reaction_type_in(1:), reactant_list_in(1:,0:), product_list_in(1:,0:)
    real(dp), intent(in) :: lambda_in(1:), dlambda_in(1:), solar_flux_in(1:), mass_in(1:)
    character(len=*), intent(in) :: species_in(1:)
    integer ich, nch_photolysis_priv

    nz  = nz_in
    nwl = nwl_in
    nsp = nsp_in
    nch = nch_in

    ! Allocate arrays for reactions and species
    allocate(reaction_type(nch), product_list(nch, 0:20), reactant_list(nch, 0:20))
    allocate(species(1:nsp), mass(1:nsp))

    species(1:nsp) = species_in(1:nsp)
    mass(1:nsp) = mass_in(1:nsp)
    reaction_type(1:nch) = reaction_type_in(1:nch)
    reactant_list(1:nch, 0:20) = reactant_list_in(1:nch, 0:20)
    product_list(1:nch, 0:20) = product_list_in(1:nch, 0:20)

    nch_photolysis_priv = 0
    do ich = 1, nch
      if (reaction_type(ich) == 2) then
        nch_photolysis_priv = nch_photolysis_priv + 1
      end if
    end do

    ! Allocate the module-level array with the required dimensions
    allocate(sigma_dat(-2:nwl, 1:nsp+nch_photolysis_priv, 0:20, 0:10, 0:2))
    allocate(sigma_a(-2:nwl, 1:nz, 1:nsp))
    allocate(sigma_d(-2:nwl, 1:nz, 1:nch_photolysis_priv))
    allocate(lambda(1:nwl), dlambda(1:nwl), solar_flux(1:nwl))
    allocate(cln(nz, nsp), tau(nwl, nz), I_z(nwl, nz))
    allocate(Tdep_list(0:100))

    ! Allocate arrays for exceptional cases
    allocate(larr(nwl, 10))
    allocate(qy_H(nwl), qy_CO_300K(nwl), qy_CO_T(nwl), a_300K(nwl), a_T(nwl))
    allocate(qy_O1D(nwl))

    ! Initialization
    sigma_dat   = 0.0_dp
    sigma_a     = 0.0_dp
    sigma_d     = 0.0_dp

    lambda(1:nwl) = lambda_in(1:nwl)
    dlambda(1:nwl) = dlambda_in(1:nwl)
    solar_flux(1:nwl) = solar_flux_in(1:nwl)

  end subroutine p__photolysis_rate__ini

  
  !=======================================================================================================================
  !
  !  Deallocate private arrays to free memory
  !
  !=======================================================================================================================
  subroutine p__photolysis_rate__fin()
    if (allocated(sigma_dat))     deallocate(sigma_dat)
    if (allocated(sigma_a))       deallocate(sigma_a)
    if (allocated(sigma_d))       deallocate(sigma_d)
    if (allocated(lambda))        deallocate(lambda)
    if (allocated(dlambda))       deallocate(dlambda)
    if (allocated(solar_flux))    deallocate(solar_flux)
    if (allocated(reaction_type)) deallocate(reaction_type)
    if (allocated(product_list))  deallocate(product_list)
    if (allocated(reactant_list)) deallocate(reactant_list)
    if (allocated(species))       deallocate(species)
    if (allocated(mass))          deallocate(mass)
    if (allocated(larr))          deallocate(larr)
    if (allocated(Tdep_list))     deallocate(Tdep_list)
    if (allocated(qy_H))          deallocate(qy_H)
    if (allocated(qy_CO_300K))    deallocate(qy_CO_300K)
    if (allocated(qy_CO_T))       deallocate(qy_CO_T)
    if (allocated(a_300K))        deallocate(a_300K)
    if (allocated(a_T))           deallocate(a_T)
    if (allocated(qy_O1D))        deallocate(qy_O1D)
  end subroutine p__photolysis_rate__fin


  !=======================================================================================================================
  !
  !  Load photo-absorption and dissociation cross section data from files and bin them to the simulation wavelength grid
  !
  !=======================================================================================================================
  subroutine load_cross_section_dat(dirname) ! in
    implicit none
    character(len=*), intent(in) :: dirname
    integer isp, ich
    real(dp) dummy(1), dummy2(20)
    character(len=256) fname

    ! initialization
    nch_photolysis = 0
    Tdep_list = 1.0_dp 
    sigma_dat = 0.0_dp
    sigma_dat(-2,:,:,1,1) = dble(nwl)
    sigma_dat(-1,:,:,1,1) = dble(1.0_dp)
    dummy(1) = 1.0_dp ! dummy variable for temperature array
    dummy2(1:20) = 1.0_dp ! dummy variable for temperature array

    ! Template --------------------------------------------------------
    !  isp = sp_index(spl,'CO2') 
    !------------------------------------------------------------------
    !
    !  fname = trim(ADJUSTL(dirname))//'/sigma_*.dat'
    !  call read_binning_data(fname, Tdep_list,         & ! fixed, do not modify these three
    !    &                    bin1, bin2, unit1, unit2, & ! You should modify these three inputs
    !    &                    sigma_dat, type, index    )
    !   n_T_bin: number of cross section data column, e.g.) 12 (if there are cross sections for 12 temperatures)
    !   unit1: unit of column 1 of the datafile, 'nm', 'cm', or 'cm-1'.
    !   unit2: unit of column 2 of the datafile, 'cm2', 'm2', or ''.
    !          if data is a quantum yeild, it should be ''.

    ! CO2 data --------------------------------------------------------
    isp = get_species_index('CO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_a_12C16O16O.dat'
    Tdep_list(0) = 12.0_dp
    Tdep_list(1:12) = (/120.0_dp, 145.0_dp, 170.0_dp, 195.0_dp, 220.0_dp, 245.0_dp, &
      &                 270.0_dp, 295.0_dp, 320.0_dp, 345.0_dp, 370.0_dp, 395.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 12, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_d_CO2_CO+O1D.dat'
    ich = get_reaction_index(['CO2'], ['CO   ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)
    
    fname = trim(ADJUSTL(dirname))//'/COx/qy_CO2_CO+O.dat'
    ich = get_reaction_index(['CO2'], ['CO','O '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! ^13CO2 data -----------------------------------------------------
    isp = get_species_index('^13CO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_a_13C16O16O.dat'
    Tdep_list(0) = 12.0_dp
    Tdep_list(1:12) = (/120.0_dp, 145.0_dp, 170.0_dp, 195.0_dp, 220.0_dp, 245.0_dp, &
      &                 270.0_dp, 295.0_dp, 320.0_dp, 345.0_dp, 370.0_dp, 395.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 12, 'nm', 'cm^2', &
      &                    'a', isp)

    ! H2O data --------------------------------------------------------
    isp = get_species_index('H2O')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_a_H2O.dat'
    call read_binning_data(fname, Tdep_list, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/HOx/qy_H2O_H+OH.dat'
    ich = get_reaction_index(['H2O'], ['H ','OH'])
    call read_binning_data(fname, Tdep_list, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/HOx/qy_H2O_H2+O1D.dat'
    ich = get_reaction_index(['H2O'], ['H2   ','O(1D)'])
    call read_binning_data(fname, Tdep_list, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! O2 data ---------------------------------------------------------
    isp = get_species_index('O2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O2.dat'
    Tdep_list(0) = 2.0_dp
    Tdep_list(1:2) = (/90.0_dp, 295.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 2, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O2_Schumann-Runge_bands_130-190K.dat'
    call read_binning_data(fname, dummy2, &
      &                    3, 7, 'cm^-1', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O2_Schumann-Runge_bands_190-280K.dat'
    call read_binning_data(fname, dummy2, &
      &                    8, 12, 'cm^-1', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O2_Schumann-Runge_bands_280-500K.dat'
    call read_binning_data(fname, dummy2, &
      &                    13, 17, 'cm^-1', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/qy_O2_O+O.dat'
    ich = get_reaction_index(['O2'], ['O','O'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/Ox/qy_O2_O+O1D.dat'
    ich = get_reaction_index(['O2'], ['O    ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! H2O2 data -------------------------------------------------------
    isp = get_species_index('H2O2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_a_H2O2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/HOx/qy_H2O2_OH+OH.dat'
    ich = get_reaction_index(['H2O2'], ['OH','OH'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/HOx/qy_H2O2_HO2+H.dat'
    ich = get_reaction_index(['H2O2'], ['HO2','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! H2 data ---------------------------------------------------------
    isp = get_species_index('H2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_a_H2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_d_H2_H+H.dat'
    ich = get_reaction_index(['H2'], ['H','H'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! OH data ---------------------------------------------------------
    isp = get_species_index('OH')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_a_OH.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_d_OH_O+H.dat'
    ich = get_reaction_index(['OH'], ['O','H'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_d_OH_O1D+H.dat'
    ich = get_reaction_index(['OH'], ['O(1D)','H    '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! HO2 data --------------------------------------------------------
    isp = get_species_index('HO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/HOx/sigma_a_HO2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/HOx/qy_HO2_OH+O.dat'
    ich = get_reaction_index(['HO2'], ['OH','O '])

    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! O3 data ---------------------------------------------------------
    isp = get_species_index('O3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O3.dat'
    Tdep_list(0) = 2.0_dp
    Tdep_list(1:11) = (/193.0_dp, 203.0_dp, 213.0_dp, 223.0_dp, 233.0_dp, 243.0_dp, &
      &                 253.0_dp, 263.0_dp, 273.0_dp, 283.0_dp, 293.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 11, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/qy_O3_O1D_dummy.dat'
    ich = get_reaction_index(['O3'], ['O2   ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy-excep', ich)

    fname = trim(ADJUSTL(dirname))//'/Ox/qy_O3_O_dummy.dat'
    ich = get_reaction_index(['O3'], ['O2','O '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy-excep', ich)

    ! H2CO data -------------------------------------------------------
    isp = get_species_index('H2CO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_H2CO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_gamma_H2CO.dat'
    call read_binning_data(fname, dummy, &
      &                    2, 2, 'nm', 'cm^2', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_H2CO_H2+COat300K1atm.dat'
    ich = get_reaction_index(['H2CO'], ['H2','CO'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_H2CO_HCO+H_dummy.dat'
    ich = get_reaction_index(['H2CO'], ['HCO','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy-excep', ich)

    ! HCO data --------------------------------------------------------
    isp = get_species_index('HCO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_HCO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_HCO_H+CO.dat'
    ich = get_reaction_index(['HCO'], ['H ','CO'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! N2 data ---------------------------------------------------------
    isp = get_species_index('N2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_N2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_d_N2_N+N.dat'
    ich = get_reaction_index(['N2'], ['N    ','N(2D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! NO data ---------------------------------------------------------
    isp = get_species_index('NO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_NO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_d_NO_N+O.dat'
    ich = get_reaction_index(['NO'], ['N','O'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! NO2 data --------------------------------------------------------
    isp = get_species_index('NO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_NO2.dat'
    Tdep_list(0) = 2.0_dp
    Tdep_list(1:2) = (/220.0_dp, 294.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 2, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_NO2_NO+O.dat'
    Tdep_list(0) = 2.0_dp
    Tdep_list(1:2) = (/248.0_dp, 298.0_dp/)
    ich = get_reaction_index(['NO2'], ['NO','O '])
    call read_binning_data(fname, Tdep_list, &
      &                    1, 2, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_d_NO2_NO+O(1D).dat'
    ich = get_reaction_index(['NO2'], ['NO   ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! NO3 data --------------------------------------------------------
    isp = get_species_index('NO3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_NO3.dat'
    Tdep_list(0) = 2.0_dp
    Tdep_list(1:2) = (/220.0_dp, 298.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 2, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_NO3_NO+O2.dat'
    Tdep_list(0) = 3.0_dp
    Tdep_list(1:3) = (/190.0_dp, 230.0_dp, 298.0_dp/)
    ich = get_reaction_index(['NO3'], ['NO','O2'])
    call read_binning_data(fname, Tdep_list, &
      &                    1, 3, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_NO3_NO2+O.dat'
    Tdep_list(0) = 3.0_dp
    Tdep_list(1:3) = (/190.0_dp, 230.0_dp, 298.0_dp/)
    ich = get_reaction_index(['NO3'], ['NO2','O  '])
    call read_binning_data(fname, Tdep_list, &
      &                    1, 3, 'nm', '', &
      &                    'qy', ich)

    ! N2O data --------------------------------------------------------
    isp = get_species_index('N2O')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_N2O.dat'
    Tdep_list(0) = 5.0_dp
    Tdep_list(1:5) = (/194.0_dp, 225.0_dp, 243.0_dp, 263.0_dp, 302.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 5, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_N2O_N2+O(1D).dat'
    ich = get_reaction_index(['N2O'], ['N2   ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! N2O5 data -------------------------------------------------------
    isp = get_species_index('N2O5')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_N2O5.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_N2O5_Tdep.dat'
    call read_binning_data(fname, dummy2, &
      &                    2, 3, 'nm', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_N2O5_NO3+NO2.dat'
    ich = get_reaction_index(['N2O5'], ['NO3','NO2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_N2O5_NO3+NO+O.dat'
    ich = get_reaction_index(['N2O5'], ['NO3','NO ','O  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HONO data -------------------------------------------------------
    isp = get_species_index('HONO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HONO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HONO_OH+NO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HNO3 data -------------------------------------------------------
    isp = get_species_index('HNO3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HNO3.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HNO3_Tdep.dat'
    call read_binning_data(fname, dummy, &
      &                    2, 2, 'nm', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HNO3_OH+NO2.dat'
    ich = get_reaction_index(['HNO3'], ['OH ','NO2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HNO3_HONO+O.dat'
    ich = get_reaction_index(['HNO3'], ['HONO','O   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HNO3_HONO+O(1D).dat'
    ich = get_reaction_index(['HNO3'], ['HONO ','O(1D)'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HO2NO2 data -----------------------------------------------------
    isp = get_species_index('HO2NO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HO2NO2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HO2NO2_Tdep.dat'
    call read_binning_data(fname, dummy2, &
      &                    2, 3, 'nm', '', &
      &                    'a-excep', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HO2NO2_OH+NO3.dat'
    ich = get_reaction_index(['HO2NO2'], ['OH ','NO3'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HO2NO2_HO2+NO2.dat'
    ich = get_reaction_index(['HO2NO2'], ['HO2','NO2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH4 data --------------------------------------------------------
    isp = get_species_index('CH4')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_CH4.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_d_CH4_1CH2+H2.dat'
    ich = get_reaction_index(['CH4'], ['^1CH2','H2   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_d_CH4_3CH2+H.dat'
    ich = get_reaction_index(['CH4'], ['^3CH2','H    '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_d_CH4_CH3+H.dat'
    ich = get_reaction_index(['CH4'], ['CH3 ','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! C2H6 data -------------------------------------------------------
    isp = get_species_index('C2H6')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H6.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_3CH2+H2.dat'
    ich = get_reaction_index(['C2H6'], ['^3CH2','H2   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_CH4+1CH2.dat'
    ich = get_reaction_index(['C2H6'], ['CH4  ','^1CH2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_C2H2+H2.dat'
    ich = get_reaction_index(['C2H6'], ['C2H2 ','H2   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_C2H4+H.dat'
    ich = get_reaction_index(['C2H6'], ['C2H4','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_C2H4+H2.dat'
    ich = get_reaction_index(['C2H6'], ['C2H4','H2  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_C2H6_CH3+CH3.dat'
    ich = get_reaction_index(['C2H6'], ['CH3','CH3'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HNO data --------------------------------------------------------
    isp = get_species_index('HNO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HNO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HNO_NO+H.dat'
    ich = get_reaction_index(['HNO'], ['NO ','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HNO2 data -------------------------------------------------------
    isp = get_species_index('HNO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NOy/sigma_a_HNO2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NOy/qy_HNO2_NO+OH.dat'
    ich = get_reaction_index(['HNO2'], ['NO','OH'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH3 data --------------------------------------------------------
    isp = get_species_index('CH3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_CH3.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_CH3_1CH2+H.dat'
    ich = get_reaction_index(['CH3'], ['^1CH2','H    '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! NH3 data --------------------------------------------------------
    isp = get_species_index('NH3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NHx/sigma_a_NH3.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NHx/qy_NH3_NH2+H.dat'
    ich = get_reaction_index(['NH3'], ['NH2','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)
      
    ! N2H4 data -------------------------------------------------------
    isp = get_species_index('N2H4')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NHx/sigma_a_N2H4.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)


    fname = trim(ADJUSTL(dirname))//'/NHx/qy_N2H4_N2H3+H.dat'
    ich = get_reaction_index(['N2H4'], ['N2H3','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! NH data ---------------------------------------------------------
    isp = get_species_index('NH')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NHx/sigma_a_NH.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NHx/qy_NH_N+H.dat'
    ich = get_reaction_index(['NH'], ['N  ','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! NH2 data --------------------------------------------------------
    isp = get_species_index('NH2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/NHx/sigma_a_NH2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/NHx/qy_NH2_NH+H.dat'
    ich = get_reaction_index(['NH2'], ['NH ','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! C2H2 data -------------------------------------------------------
    isp = get_species_index('C2H2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_d_C2H2_C2H+H.dat'
    ich = get_reaction_index(['C2H2'], ['C2H','H  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_d_C2H2_C2+H2.dat'
    ich = get_reaction_index(['C2H2'], ['C2','H2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! C2H4 data -------------------------------------------------------
    isp = get_species_index('C2H4')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H4.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_d_C2H4_C2H2+H2.dat'
    ich = get_reaction_index(['C2H4'], ['C2H2','H2  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_d_C2H4_C2H2+H+H.dat'
    ich = get_reaction_index(['C2H4'], ['C2H2','H   ','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'd', ich)

    ! C3H8 data -------------------------------------------------------
    isp = get_species_index('C3H8')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_C3H8.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H8_C3H6+H2.dat'
    ich = get_reaction_index(['C3H8'], ['C3H6','H2  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H8_C2H6+1CH2.dat'
    ich = get_reaction_index(['C3H8'], ['C2H6 ','^1CH2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H8_C2H4+CH4.dat'
    ich = get_reaction_index(['C3H8'], ['C2H4','CH4 '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H8_C2H5+CH3.dat'
    ich = get_reaction_index(['C3H8'], ['C2H5','CH3 '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! C3H6 data -------------------------------------------------------
    isp = get_species_index('C3H6')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_C3H6.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H6_C2H2+CH3.dat'
    ich = get_reaction_index(['C3H6'], ['C2H2','CH3 '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H6_CH2CCH2+H2.dat'
    ich = get_reaction_index(['C3H6'], ['CH2CCH2','H2     '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H6_C2H4+3CH2.dat'
    ich = get_reaction_index(['C3H6'], ['C2H4 ','^3CH2'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H6_C2H+CH4.dat'
    ich = get_reaction_index(['C3H6'], ['C2H','CH4'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH data ---------------------------------------------------------
    isp = get_species_index('CH')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_CH.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_CH_C+H.dat'
    ich = get_reaction_index(['CH'], ['C','H'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH2CO data ------------------------------------------------------
    isp = get_species_index('CH2CO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_CH2CO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_CH2CO_CH2+CO.dat'
    ich = get_reaction_index(['CH2CO'], ['CH2','CO '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)
    
    ! CH3CHO data -----------------------------------------------------
    isp = get_species_index('CH3CHO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_CH3CHO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_CH3CHO_CH3+HCO.dat'
    ich = get_reaction_index(['CH3CHO'], ['CH3','HCO'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/qy_CH3CHO_CH4+CO.dat'
    ich = get_reaction_index(['CH3CHO'], ['CH4','CO '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! C2H5CHO data ----------------------------------------------------
    isp = get_species_index('C2H5CHO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_C2H5CHO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C2H5CHO_C2H5+HCO.dat'
    ich = get_reaction_index(['C2H5CHO'], ['C2H5','HCO '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! C3H3 data -------------------------------------------------------
    isp = get_species_index('C3H3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_C3H3.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_C3H3_C3H2+H.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH3C2H data -----------------------------------------------------
    isp = get_species_index('CH3C2H')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_CH3C2H.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH3C2H_C3H3+H.dat'
    ich = get_reaction_index(['CH3C2H'], ['C3H3','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH3C2H_C3H2+H2.dat'
    ich = get_reaction_index(['CH3C2H'], ['C3H2','H2  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH3C2H_CH3+C2H.dat'
    ich = get_reaction_index(['CH3C2H'], ['CH3','C2H'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! CH2CCH2 data ----------------------------------------------------
    isp = get_species_index('CH2CCH2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C3/sigma_a_CH2CCH2.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH2CCH2_C3H3+H.dat'
    ich = get_reaction_index(['CH2CCH2'], ['C3H3','H   '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH2CCH2_C3H2+H2.dat'
    ich = get_reaction_index(['CH2CCH2'], ['C3H2','H2  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    fname = trim(ADJUSTL(dirname))//'/C3/qy_CH2CCH2_C2H2+CH2.dat'
    ich = get_reaction_index(['CH2CCH2'], ['C2H2','CH2 '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HCN data --------------------------------------------------------
    isp = get_species_index('HCN')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_HCN.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_HCN_H+CN.dat'
    ich = get_reaction_index(['HCN'], ['H ','CN'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! HNCO data -------------------------------------------------------
    isp = get_species_index('HNCO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_HNCO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_HNCO_H+NCO.dat'
    ich = get_reaction_index(['HNCO'], ['H  ','NCO'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! NCO data --------------------------------------------------------
    isp = get_species_index('NCO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_NCO.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/qy_NCO_N+CO.dat'
    ich = get_reaction_index(['NCO'], ['N ','CO'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! SO data ---------------------------------------------------------
    isp = get_species_index('SO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/SOx/sigma_a_SO_Marcq.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/SOx/qy_SO_S+O.dat'
    ich = get_reaction_index(['SO'], ['S','O'])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! SO2 data --------------------------------------------------------
    isp = get_species_index('SO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/SOx/sigma_a_SO2.dat'
    Tdep_list(0) = 3.0_dp
    Tdep_list(1:3) = (/200.0_dp, 295.0_dp, 400.0_dp/)
    call read_binning_data(fname, Tdep_list, &
      &                    1, 3, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/SOx/qy_SO2_SO+O.dat'
    ich = get_reaction_index(['SO2'], ['SO','O '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! SO3 data --------------------------------------------------------
    isp = get_species_index('SO3')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/SOx/sigma_a_SO3.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/SOx/qy_SO3_SO2+O.dat'
    ich = get_reaction_index(['SO3'], ['SO2','O  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)

    ! OCS data --------------------------------------------------------
    isp = get_species_index('OCS')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/OCS/sigma_a_OCS.dat'
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/OCS/qy_OCS_CO+S.dat'
    ich = get_reaction_index(['OCS'], ['CO ','S  '])
    call read_binning_data(fname, dummy, &
      &                    1, 1, 'nm', '', &
      &                    'qy', ich)


    ! finish reading cross sections. -----------------------------------
    print *, 'Finished reading cross section data.'


  end subroutine load_cross_section_dat


  !=======================================================================================================================
  !
  !     Get cross section for absorption (sigma_a) or photolysis (sigma_d) 
  !
  !=======================================================================================================================
  subroutine get_cross_section(T, flag, outflag)
    implicit none       
    integer,    intent(in) :: outflag
    real(dp),   intent(in) :: T(1:)
    character(len=*), intent(in) :: flag
    integer, parameter :: photolysis = 2 ! reaction type index for photolysis
    integer iwl, iz, isp, ksp, ich, jch, kch, swl, ewl
    character(len=256) fname, dirname, command
    real(dp), parameter :: pi = dacos(-1.0_dp)

    !----------------------------------------------------------------------------------------------------
    ! UV photoabsorption cross section

    if (flag == 'absorption') then

       sigma_a(0,:,:) = 0 ! reset cross section existence flag
      
       do isp = 1, nsp

        if (nint(sigma_dat(0,isp,0,0,0)) == 0) cycle ! skip if no data

        ! For normal species ------------------------------------------

        call calc_sigma_a(isp, T) ! in

        ! For exceptional species -------------------------------------

        ! H2O2 for 260 ~ 340 nm range
        if (species(isp) == 'H2O2') then 
          !call calc_sigma_a_H2O2_analytic(isp, T) ! in
        end if

        ! O2 Schumann-Runge bands
        if (species(isp) == 'O2') then 
          call calc_sigma_a_O2_Schumann_Runge_bands(isp, T) ! in
          call calc_sigma_a_O2_Herzberg_continuum(isp) ! in
        end if

        ! H2CO temperature dependence
        if (species(isp) == 'H2CO') then
          call calc_sigma_a_H2CO_Tdependent(isp, T) ! in
        end if

        ! N2O5 temperature dependence
        if (species(isp) == 'N2O5') then
          call calc_sigma_a_N2O5_Tdependent(isp, T) ! in
        end if

        ! HNO3 temperature dependence
        if (species(isp) == 'HNO3') then
          call calc_sigma_a_HNO3_Tdependent(isp, T) ! in
        end if

        ! HO2NO2 temperature dependence
        if (species(isp) == 'HO2NO2') then
          call calc_sigma_a_HO2NO2_Tdependent(isp, T) ! in
        end if

       end do !isp


      ! if cross sections are less than 0, they are set to 0 ------
      where (sigma_a(:,:,:) < 0.0_dp)
        sigma_a(:,:,:) = 0.0_dp
      end where

      ! output cross sections ------
      if (outflag == 1) then 
        dirname = './UV/xsect/absorption'
        write(command,*) 'mkdir -p ', trim(dirname)
        call system(command)
        iz = 1 ! select altitude grid
        do isp = 1, nsp
          if (nint(sigma_a(0,1,isp)) == 1) then 
            fname = './UV/xsect/absorption/'//trim(ADJUSTL(species(isp)))//'.dat'
            open(11, file = fname, status = 'replace' )
            swl = nint(sigma_a(-2,1,isp))
            ewl = nint(sigma_a(-1,1,isp))
            !print *, isp, trim(species(isp)), swl, ewl
            if (swl >= 1) then 
              do iwl = swl, ewl
                write(11, *) lambda(iwl), sigma_a(iwl,1,isp)
              end do
            end if
            close(11)
          end if
        end do 
      end if

    end if ! end if flag == absorption

    !----------------------------------------------------------------------------------------------------
    ! UV photodissociation cross section

    if (flag == 'photolysis') then

      do ich = 1, nch_photolysis

        ! For normal photolysis reactions -----------------------------

        jch = nint(sigma_dat(0,nsp+ich,0,0,0)) ! reaction index
        isp = reactant_list(jch,1) ! reactant species index
        call calc_sigma_d(isp, ich, T) ! inout

        ! For exceptional photolysis reactions ------------------------

        ! O3 quantum yield
        ksp = get_species_index('O3')
        if (ksp==isp) then 
          kch = get_reaction_index(['O3'], ['O2','O '])
          if (kch==jch) call calc_sigma_d_O3(isp, ich, T, 'O') ! in
          kch = get_reaction_index(['O3'], ['O2   ','O(1D)'])
          if (kch==jch) call calc_sigma_d_O3(isp, ich, T, 'O(1D)') ! in
        end if

        ! H2CO quantum yield
        ksp = get_species_index('H2CO')
        if (ksp==isp) then 
          kch = get_reaction_index(['H2CO'], ['H2','CO'])
          if (kch==jch) call calc_sigma_d_H2CO(isp, ich, T, 'CO') ! in
          kch = get_reaction_index(['H2CO'], ['HCO','H  '])
          if (kch==jch) call calc_sigma_d_H2CO(isp, ich, T, 'H') ! in
        end if

      end do !ich


      ! if cross sections are less than 0, they are set to 0 ------
      where (sigma_d(:,:,:)<0.0_dp)
        sigma_d(:,:,:) = 0.0_dp
      end where

      ! output cross sections ------
      if (outflag == 1) then 
        dirname = './UV/xsect/dissociation'
        write(command,*) 'mkdir -p ', trim(dirname)
        call system(command)
        iz = 1 ! select altitude grid
        do jch = 1, nch_photolysis
          ich = nint(sigma_dat(0,nsp+jch,0,0,0)) ! reaction index
          if (reaction_type(ich) == 2) then
            if (nint(sigma_d(0,1,jch)) == 1) then 
              if (product_list(ich,0)==1) then 
                fname = './UV/xsect/dissociation/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'.dat'
              else if (product_list(ich,0)==2) then 
                fname = './UV/xsect/dissociation/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'.dat'
              else if (product_list(ich,0)==3) then 
                fname = './UV/xsect/dissociation/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,3))))//'.dat'
              else if (product_list(ich,0)==4) then 
                fname = './UV/xsect/dissociation/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,3))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,4))))//'.dat'
              end if
              open(11, file = fname, status = 'replace' )
              swl = nint(sigma_d(-2,1,jch))
              ewl = nint(sigma_d(-1,1,jch))
              !print *, ich, swl, ewl, trim(fname)
              if (swl >= 1) then 
                do iwl = swl, ewl
                  write(11, *) lambda(iwl), sigma_d(iwl,1,jch)
                end do
              end if
              close(11)
            end if
          end if
        end do 
      end if

    end if ! end if flag == photolysis


  end subroutine get_cross_section


  !=======================================================================================================================
  !
  !    Calculate optical depth and photolysis rate 
  !
  !=======================================================================================================================
  subroutine photolysis_rate(alt, dalt, density, T, sza, & ! in
    &                        rplanet, mplanet,           & ! in
    &                        J_rate                      ) ! out
    implicit none
    real(dp), intent(in)    :: density(1:,1:), alt(1:), dalt(1:), T(1:)
    real(dp), intent(in)    :: sza, rplanet, mplanet
    real(dp), intent(inout) :: J_rate(1:,1:) ! 1:nz, 1:nch
    real(dp) hmin
    integer isp, ich, jch, iz, swl, ewl

    ! for solar zenith angle near and greater than 90deg [Smith et al., 1972]
    real(dp) Hz, yz, Xz, chiz, Chfunc, cln_Ch, g
    ! ap : parameter in approximating exp(x^2)*erfc(x) recommended by Ren and MacKenzie [2007]
    real(dp), parameter :: ap = 2.7889_dp 
    ! Physical constant
    real(dp), parameter :: pi = dacos(-1.0_dp)
    real(dp), parameter :: k_B = 1.38064852e-23_dp ! Boltzmann constant in J/K
    real(dp), parameter :: Grav = 6.67430e-11_dp ! Gravitational constant in m^3 kg^-1 s^-2

    ! column density ----------------------------
    cln = 0.0_dp
    do isp = 1, nsp
      if (species(isp) == 'M') cycle
      cln(nz,isp) = density(nz,isp)*dalt(nz)
      do iz = nz, 2, -1
        cln(iz-1,isp) = cln(iz,isp) &
          &           + 0.5_dp * (density(iz-1,isp) + density(iz,isp)) * dalt(iz-1)
      end do
    end do

    ! radiative transfer ------------------------
    tau = 0.0_dp ! optical depth
    I_z = 0.0_dp ! solar flux at each altitude

    do isp = 1, nsp

      if (nint(sigma_a(0,1,isp)) == 0) cycle ! skip if no data
      swl = nint(sigma_a(-2,1,isp))
      ewl = nint(sigma_a(-1,1,isp))

      do iz = 1, nz
        g = Grav * mplanet / (rplanet + alt(iz))**2
        chiz = sza
        Hz   = k_B * T(iz) / mass(isp) / g
        Xz   = (rplanet + alt(iz)) / Hz
        yz   = dsqrt(0.5_dp * Xz) * dabs(dcos(chiz))

        ! exp(x^2)*erfc(x) is approximated by using the formula of Ren & MacKenzie [2007]
        if (chiz <= pi / 2.0_dp) then 
          hmin = - 1.0e10_dp
          Chfunc = dsqrt(pi/2.0_dp*Xz) &
            &        * ap / ((ap-1.0_dp)*dsqrt(pi*yz*yz) + dsqrt(pi*yz*yz + ap*ap))
        else if (chiz > pi / 2.0_dp) then 
          hmin = rplanet * (1.0_dp / dcos(chiz-pi/2.0_dp) - 1.0_dp)
          Chfunc = dsqrt(2.0_dp*pi*Xz) &
            &        * ( dsqrt(dsin(chiz)) * dexp( Xz*(1.0_dp - dsin(chiz)) ) &
            &          - 0.5_dp * ap / ((ap-1.0_dp)*dsqrt(pi*yz*yz) + dsqrt(pi*yz*yz + ap*ap)) )
        end if

        ! upper limit of Chfunc is 10^100 in order not to cause infinity tau
        if (Chfunc > 1.0e100_dp) then
          Chfunc = 1.0e100_dp
        end if

        cln_Ch = cln(iz,isp) * Chfunc
        if (alt(iz) > hmin) then 
          tau(swl:ewl,iz) = tau(swl:ewl,iz) &
            &             + cln_Ch * sigma_a(swl:ewl,iz,isp)
        else if (alt(iz) <= hmin) then
          tau(1:nwl,iz) = 1.0e10_dp
        end if
      end do ! iz

    end do ! isp

    do iz  = 1, nz
      I_z(1:nwl,iz) = solar_flux(1:nwl) * dlambda(1:nwl) * dexp(-tau(1:nwl,iz)) 
    end do

    ! photolysis rate ---------------------------
    do ich = 1, nch_photolysis

      if (nint(sigma_d(0,1,ich)) == 0) cycle ! skip if no data
      swl = nint(sigma_d(-2,1,ich))
      ewl = nint(sigma_d(-1,1,ich))
      jch = nint(sigma_dat(0,nsp+ich,0,0,0)) ! reaction index

      do iz = 1, nz
        J_rate(iz,jch) = dot_product(I_z(swl:ewl,iz), sigma_d(swl:ewl,iz,ich))
      end do

    end do ! ich


  end subroutine photolysis_rate


  !----------------------------------------------------------------------------------------------
  ! Flux, cross sections and quantum yield data are automatically adapted to the model bins.
  !   if the data bin is larger than the model bin, the data is interpolated.
  !   if the data bin is smaller than the model bin, the data is binned.
  !----------------------------------------------------------------------------------------------
  subroutine binning_cross_section(nl, idata,      & ! in
    &                              odata, swl, ewl ) ! out
    implicit none
    integer,    intent(in)  :: nl
    real(dp),   intent(in)  :: idata(nl,2)
    real(dp),   intent(out) :: odata(nwl)
    integer,    intent(out) :: swl, ewl
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
    swl = 9999999
    ewl = 0
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
          if (swl > iwl) swl = iwl
          if (ewl < iwl) ewl = iwl
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

    if (swl == 9999999) then 
      swl = nwl
      ewl = nwl
    end if
    if (ewl > nwl) then 
      ewl = nwl
    end if
    if (ewl == 0) then 
      swl = 1
      ewl = 1
    end if
    if (swl > ewl) then 
      print *, 'error in binning.'
      stop
    end if

    !print *, nl, swl, ewl
  
  end subroutine binning_cross_section


  !------------------------------------------------------------
  ! For the use of reading cross section data easily.
  !------------------------------------------------------------
  subroutine read_binning_data(fname, Tdep_list,           & ! inout, in
    &                          sdata, edata, unit1, unit2, & ! in
    &                          dtype, id_in                ) ! inout, in, in
    implicit none
    real(dp),         intent(in)    :: Tdep_list(0:)
    integer,          intent(in)    :: id_in
    integer,          intent(in)    :: sdata, edata
    character(len=*), intent(in)    :: unit1, unit2, fname, dtype
    integer i, j, nh, il, nl, swl, ewl, idtype, ndata, id, id_flag
    real(dp), allocatable :: idata(:,:), odata(:), rdata(:,:)

    if (id_in >= 1) then 

      ndata = edata - sdata + 1

      if (dtype == 'a') then 
        idtype = 2**0 ! absorption
        id = id_in
      else if (dtype == 'a-excep') then 
        idtype = 2**3 ! exception for absorption
        id = id_in
      else if (dtype == 'd') then 
        idtype = 2**6 ! dissociation
      else if (dtype == 'd-excep') then 
        idtype = 2**9 ! exception for dissociation
      else if (dtype == 'qy') then 
        idtype = 2**12 ! quantum yield
      else if (dtype == 'qy-excep') then 
        idtype = 2**15 ! exception
      else 
        print *, ''
        print *, 'error!'
        print *, '  The dtype "'//trim(adjustl(dtype))//'" is not recognized.'
        stop
      end if

      id_flag = 0
      if (idtype >= 10) then
        do i = nsp+1, nsp+nch_photolysis
          if (nint(sigma_dat(0,i,0,0,0)) == id_in) then 
            id = i
            id_flag = 1
          end if
        end do 
      end if
      if (idtype >= 10 .and. id_flag == 0) then 
        nch_photolysis = nch_photolysis + 1
        id = nsp + nch_photolysis
      end if

      write(*,'(a)',advance='no')  '  Reading datafile: '//trim(ADJUSTL(fname))//'...'

      allocate(odata(nwl))

      if (ndata == 1) then 

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
            if (unit2/='m^2' .and. unit2/='cm^2' &
            &  .and. unit2/='/cm^2/s/nm' .and. unit2/='/m^2/s/nm' .and. unit2/='W/m^2/nm' .and. unit2/='') then 
              print *, ''
              print *, 'error!'
              print *, '  The unit of the datafile "'//trim(adjustl(fname))//'" column 2 is not recognized.'
              print *, '  The unit of column 2 should be ["m^2", "cm^2"] for cross sections, '
              print *, '  "" for unitless data, and ["/m^2/s/nm", "/cm^2/s/nm", "W/m^2/nm"] for solar flux.'
              stop
            end if
            if (unit2 == 'm^2') then 
              if (idata(il,2)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "m" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
            else if (unit2 == 'cm^2') then 
              if (idata(il,2)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "cm" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
              idata(il,2) = idata(il,2) * 1.0e-4_dp ! [cm^2 -> m^2]
            end if

          end do
        close(11)

        call binning_cross_section(nl, idata,      &
          &                        odata, swl, ewl )
        if (sigma_dat(-2,id,1,1,1) > dble(swl)) sigma_dat(-2,id,1,1,1) = dble(swl) ! start wavelength
        if (sigma_dat(-1,id,1,1,1) < dble(ewl)) sigma_dat(-1,id,1,1,1) = dble(ewl) ! end wavelength
        sigma_dat(0,id,1,1,1) = 1.0_dp ! existence flag
        if (dtype=='a' .or. dtype=='d' .or. dtype=='qy')sigma_dat(1,id,0,1,1) = sigma_dat(1,id,0,1,1) + dble(ndata)
        sigma_dat(0,id,1,1,0) = sigma_dat(0,id,1,1,0) + dble(idtype)
        sigma_dat(swl:ewl,id,sdata,1,2) = odata(swl:ewl)
        sigma_dat(0,id,0,0,0) = dble(id_in) ! reaction index

        deallocate(idata)

      else if (ndata >= 2) then 

        call get_header_line_number(nh, fname)
        call get_data_line_number(nl, fname)
        allocate(idata(nl,2), rdata(nl,ndata))

        open(11, file = fname, status = 'old')
          do i = 1, nh; read(11,*); end do
          do il = 1, nl
            read(11,*) idata(il,1), (rdata(il,j), j = 1, ndata) 

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
            if (unit2/='m^2' .and. unit2/='cm^2' .and. unit2/='/cm^2/s/nm' .and. unit2/='W/m^2/nm' .and. unit2/='') then 
              print *, ''
              print *, 'error!'
              print *, '  The unit of the datafile "'//trim(adjustl(fname))//'" column 2 is not recognized.'
              print *, '  The unit of column 2 should be ["m^2", "cm^2"] for cross sections, '
              print *, '  "" for unitless data, and ["/m^2/s/nm", "/cm^2/s/nm", "W/m^2/nm"] for solar flux.'
              stop
            end if
            if (unit2 == 'm^2') then 
              if (rdata(il,1)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "m" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
            else if (unit2 == 'cm^2') then 
              if (rdata(il,1)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "cm" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
              rdata(il,:) = rdata(il,:) * 1.0e-4_dp ! [cm^2 -> m^2]
            end if

          end do
        close(11)

        sigma_dat(0,id,1,1,1) = 1.0_dp ! existence flag
        sigma_dat(0,id,1,1,0) = sigma_dat(0,id,1,1,0) + dble(idtype)
        if (dtype=='a' .or. dtype=='d' .or. dtype=='qy') sigma_dat(1,id,0,1,1) = sigma_dat(1,id,0,1,1) + dble(ndata) * 1000.0_dp
        do i = sdata, edata
          idata(:,2) = rdata(:,i-sdata+1) ! copy data to idata
          call binning_cross_section(nl, idata,      &
            &                        odata, swl, ewl )
          if (sigma_dat(-2,id,i,1,1) > dble(swl)) sigma_dat(-2,id,i,1,1) = dble(swl) ! start wavelength
          if (sigma_dat(-1,id,i,1,1) < dble(ewl)) sigma_dat(-1,id,i,1,1) = dble(ewl) ! end wavelength
          sigma_dat(1,id,i,1,1) = Tdep_list(i-sdata+1) ! temperature
          sigma_dat(swl:ewl,id,i,1,2) = odata(swl:ewl)
        end do 
        sigma_dat(0,id,0,0,0) = dble(id_in) ! reaction index

        deallocate(idata, rdata)

      end if

      deallocate(odata)

      write(*,*) 'done.'

    end if

  end subroutine read_binning_data


  !------------------------------------------------------------
  ! Absorption cross section data are automatically adapted.
  !------------------------------------------------------------
  subroutine calc_sigma_a(isp, T) ! inout
    implicit none
    real(dp),              intent(in)    :: T(1:)
    integer,               intent(in)    :: isp
    integer i, iz, Tlabel, Tlabel1, ndata, swl, ewl, dtype
    real(dp) T0, T1, Tfrac, Tclamp, one_Tfrac

    sigma_a(0,:,isp) = sigma_dat(0,isp,1,1,1)

    swl = nint(sigma_dat(-2,isp,1,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,1,1,1)) ! end wavelength
    sigma_a(-2,:,isp) = dble(swl)
    sigma_a(-1,:,isp) = dble(ewl)
    ndata = nint(sigma_dat(1,isp,0,1,1)) ! number of data columns
    dtype = nint(sigma_dat(0,isp,1,1,0)) ! data type

    do iz = 1, nz

      if (mod(dtype,2**3) /= 0) then 
        if (mod(ndata,1000) == 1) then 
          sigma_a(swl:ewl,iz,isp) = sigma_dat(swl:ewl,isp,1,1,2)
        else if (ndata >= 1000) then 
          Tclamp = T(iz)
          if (T(iz)<=sigma_dat(1,isp,1,1,1)) then 
            Tlabel = 1
            Tlabel1 = 1
            Tclamp = sigma_dat(1,isp,1,1,1)
            Tfrac  = 0.0_dp
          end if
          do i = 1, ndata/1000-1
            T0 = sigma_dat(1,isp,i,1,1)
            T1 = sigma_dat(1,isp,i+1,1,1)
            if (T(iz)>T0 .and. T(iz)<=T1) then
              Tlabel = i
              Tlabel1 = i+1
              Tfrac  = (Tclamp-T0)/(T1-T0)
            end if
          end do
          if (T(iz)>sigma_dat(1,isp,ndata/1000,1,1)) then 
            Tlabel = ndata/1000
            Tlabel1 = ndata/1000
            Tclamp = sigma_dat(1,isp,ndata/1000,1,1)
            Tfrac  = 1.0_dp
          end if
          one_Tfrac = 1.0_dp - Tfrac

          sigma_a(swl:ewl,iz,isp) = one_Tfrac*sigma_dat(swl:ewl,isp,Tlabel, 1,2) &
            &                     +     Tfrac*sigma_dat(swl:ewl,isp,Tlabel1,1,2)
        end if
      end if

    end do 

  end subroutine calc_sigma_a


  !------------------------------------------------------------
  ! Absorption cross section data are automatically adapted.
  !------------------------------------------------------------
  subroutine calc_sigma_d(isp, ich, T) ! inout
    implicit none
    real(dp),              intent(in)    :: T(1:)
    integer,               intent(in)    :: ich, isp
    integer i, iz, Tlabel, Tlabel1, ndata, swl, ewl, dtype
    real(dp) T0, T1, Tfrac, Tclamp, one_Tfrac

    sigma_d(0,:,ich) = sigma_dat(0,nsp+ich,1,1,1)

    swl = nint(sigma_dat(-2,nsp+ich,1,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,nsp+ich,1,1,1)) ! end wavelength
    sigma_d(-2,:,ich) = dble(swl)
    sigma_d(-1,:,ich) = dble(ewl)

    ndata = nint(sigma_dat(1,nsp+ich,0,1,1)) ! number of data columns
    dtype = nint(sigma_dat(0,nsp+ich,1,1,0)) ! data type

    do iz = 1, nz

      ! Use dissociation cross section data
      if (mod(dtype,2**9) /= 0) then 

        if (mod(ndata,1000) == 1) then 
          sigma_d(swl:ewl,iz,ich) = sigma_dat(swl:ewl,nsp+ich,1,1,2)
        else if (ndata >= 1000) then 
          Tclamp = T(iz)
          if (T(iz)<=sigma_dat(1,nsp+ich,1,1,1)) then 
            Tlabel = 1
            Tlabel1 = 1
            Tclamp = sigma_dat(1,nsp+ich,1,1,1)
            Tfrac  = 0.0_dp
          end if
          do i = 1, ndata/1000-1
            T0 = sigma_dat(1,nsp+ich,i,1,1)
            T1 = sigma_dat(1,nsp+ich,i+1,1,1)
            if (T(iz)>T0 .and. T(iz)<=T1) then
              Tlabel = i
              Tlabel1 = i+1
              Tfrac  = (Tclamp-T0)/(T1-T0)
            end if
          end do
          if (T(iz)>sigma_dat(1,nsp+ich,ndata/1000,1,1)) then 
            Tlabel = ndata/1000
            Tlabel1 = ndata/1000
            Tclamp = sigma_dat(1,nsp+ich,ndata/1000,1,1)
            Tfrac  = 1.0_dp
          end if
          one_Tfrac = 1.0_dp - Tfrac

          sigma_d(swl:ewl,iz,ich) = one_Tfrac*sigma_dat(swl:ewl,nsp+ich,Tlabel, 1,2) &
            &                     +     Tfrac*sigma_dat(swl:ewl,nsp+ich,Tlabel1,1,2)
        end if

      ! use quantum yield and absorption cross section data
      else if (mod(dtype,2**15) /= 0) then

        if (mod(ndata,1000) == 1) then 
          sigma_d(swl:ewl,iz,ich) = sigma_a(swl:ewl,iz,isp) * sigma_dat(swl:ewl,nsp+ich,1,1,2)
        else if (ndata >= 1000) then 
          Tclamp = T(iz)
          if (T(iz)<=sigma_dat(1,nsp+ich,1,1,1)) then 
            Tlabel = 1
            Tlabel1 = 1
            Tclamp = sigma_dat(1,nsp+ich,1,1,1)
            Tfrac  = 0.0_dp
          end if
          do i = 1, ndata/1000-1
            T0 = sigma_dat(1,nsp+ich,i,1,1)
            T1 = sigma_dat(1,nsp+ich,i+1,1,1)
            if (T(iz)>T0 .and. T(iz)<=T1) then
              Tlabel = i
              Tlabel1 = i+1
              Tfrac  = (Tclamp-T0)/(T1-T0)
            end if
          end do
          if (T(iz)>sigma_dat(1,nsp+ich,ndata/1000,1,1)) then 
            Tlabel = ndata/1000
            Tlabel1 = ndata/1000
            Tclamp = sigma_dat(1,nsp+ich,ndata/1000,1,1)
            Tfrac  = 1.0_dp
          end if
          one_Tfrac = 1.0_dp - Tfrac

          sigma_d(swl:ewl,iz,ich) = one_Tfrac*sigma_a(swl:ewl,iz,isp)*sigma_dat(swl:ewl,nsp+ich,Tlabel, 1,2) &
            &                     +     Tfrac*sigma_a(swl:ewl,iz,isp)*sigma_dat(swl:ewl,nsp+ich,Tlabel1,1,2)
        end if

      end if

    end do 

  end subroutine calc_sigma_d


  !-----------------------------------------------------------------
  ! Absorption cross sections for O2 Schumann-Runge bands @ 130-500K
  !-----------------------------------------------------------------
  subroutine calc_sigma_a_O2_Schumann_Runge_bands(isp, T)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer iz, swl, ewl
    real(dp) Tz, delta, delta2

    swl = nint(sigma_dat(-2,isp,3,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,3,1,1)) ! end wavelength
    if (sigma_a(-2,1,1) > dble(swl)) sigma_a(-2,:,isp) = dble(swl)
    if (sigma_a(-1,1,1) < dble(ewl)) sigma_a(-1,:,isp) = dble(ewl)

    ! Schumann-Runge bands according to Minschwaner et al. (1992)
    do iz = 1, nz
      Tz = T(iz)
      if (T(iz)<=130.0_dp) then 
        Tz = 130.0_dp
      else if (T(iz)>500.0_dp) then 
        Tz = 500.0_dp
      end if

      delta = ((Tz-100.0_dp)/10.0_dp)**2.0_dp
      delta2 = delta*delta

      if (Tz >= 130.0_dp .and. Tz < 190.0_dp) then 
        sigma_a(swl:ewl,iz,isp) = 1.0e-24_dp * ( sigma_dat(swl:ewl,isp,3,1,2) * delta2 &
          &                            +         sigma_dat(swl:ewl,isp,4,1,2) * delta &
          &                            +         sigma_dat(swl:ewl,isp,5,1,2) )
      else if (Tz >= 190.0_dp .and. Tz < 280.0_dp) then 
        sigma_a(swl:ewl,iz,isp) = 1.0e-24_dp * ( sigma_dat(swl:ewl,isp,8,1,2) * delta2 &
          &                            +         sigma_dat(swl:ewl,isp,9,1,2) * delta &
          &                            +         sigma_dat(swl:ewl,isp,10,1,2) )
      else if (Tz >= 280.0_dp .and. Tz <= 500.0_dp) then 
        sigma_a(swl:ewl,iz,isp) = 1.0e-24_dp * ( sigma_dat(swl:ewl,isp,13,1,2) * delta2 &
          &                            +         sigma_dat(swl:ewl,isp,14,1,2) * delta &
          &                            +         sigma_dat(swl:ewl,isp,15,1,2) )
      end if

    end do


  end subroutine calc_sigma_a_O2_Schumann_Runge_bands


  !------------------------------------------------------
  ! Absorption cross sectinos for O2 Herzberg continuum
  !------------------------------------------------------
  subroutine calc_sigma_a_O2_Herzberg_continuum(isp) ! inout
    implicit none
    integer,               intent(in)    :: isp
    integer iz, iwl, swl, ewl
    real(dp) l

    sigma_a(0,:,isp) = 1.0_dp ! existence flag
    larr = 0.0_dp

    do iwl = 1, nwl-1
      if (     lambda(iwl) <= 193.0_dp .and. 193.0_dp <  lambda(iwl+1)) then 
        swl = iwl
      else if (lambda(iwl) <  245.0_dp .and. 245.0_dp <= lambda(iwl+1)) then 
        ewl = iwl
      end if

      l = lambda(iwl)
      larr(iwl,1) = l
      larr(iwl,2) = l*l
      larr(iwl,3) = larr(iwl,2)*l
      larr(iwl,4) = larr(iwl,3)*l
    end do 

    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,:,isp) = dble(swl)
    if (sigma_a(-2,1,isp) < dble(ewl)) sigma_a(-1,:,isp) = dble(ewl)

    ! Yoshino et al. (1992)
    sigma_a(swl:ewl,1,isp) = sigma_a(swl:ewl,1,isp) &
        & + 1.0e-28_dp*(-2.3837947e4_dp &
        &               +4.1973085e2_dp * larr(swl:ewl,1) &
        &               -2.7640139e0_dp * larr(swl:ewl,2) &
        &               +8.0723193e-3_dp* larr(swl:ewl,3) &
        &               -8.8255447e-6_dp* larr(swl:ewl,4) )
    do iz = 2, nz
      sigma_a(swl:ewl,iz,isp) = sigma_a(swl:ewl,1,isp) 
    end do


  end subroutine calc_sigma_a_O2_Herzberg_continuum


  !------------------------------------------------------
  ! Absorption cross sectinos for H2O2 > 260 nm by JPL 2015
  !------------------------------------------------------
  subroutine calc_sigma_a_H2O2_analytic(isp, T)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer i, iz, iwl, swl, ewl
    real(dp) l, A(0:7), B(0:4), x, one_x

    do iwl = 1, nwl-1
      if (     lambda(iwl) <  260.0_dp .and. 260.0_dp <= lambda(iwl+1)) then 
        swl = iwl+1
      else if (lambda(iwl) <  340.0_dp .and. 340.0_dp <= lambda(iwl+1)) then 
        ewl = iwl
      end if
    end do 

    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,1,isp) = dble(swl)
    if (sigma_a(-2,1,isp) < dble(ewl)) sigma_a(-1,1,isp) = dble(ewl)

    sigma_a(0,1,isp)  = 1.0_dp ! existence flag

    A(0) = 6.4761e4_dp
    A(1) = -9.2170972e2_dp
    A(2) = 4.535649_dp
    A(3) = -4.4589016e-3_dp
    A(4) = -4.035101e-5_dp
    A(5) = 1.6878206e-7_dp
    A(6) = -2.652014e-10_dp
    A(7) = 1.5534675e-13_dp

    B(0) = 6.8123e3_dp
    B(1) = -5.1351e1_dp
    B(2) = 1.1522e-1_dp
    B(3) = -3.0493e-5_dp
    B(4) = -1.0924e-7_dp

    larr = 0.0_dp

    do iwl = swl, ewl
      l = lambda(iwl)
      do i = 0, 7
        larr(iwl,1) = larr(iwl,1) + A(i) * l**dble(i) 
      end do
      larr(iwl,1) = larr(iwl,1) * 1.0e-25_dp
      do i = 0, 4
        larr(iwl,2) = larr(iwl,2) + B(i) * l**dble(i)
      end do
      larr(iwl,2) = larr(iwl,2) * 1.0e-25_dp
    end do

    do iz = 1, nz
      x = 1.0_dp / (1.0_dp + dexp(-1265.0_dp/T(iz)))
      one_x = 1.0_dp - x
      sigma_a(swl:ewl,iz,isp) = (x*larr(swl:ewl,1) + one_x*larr(swl:ewl,2))
    end do

  end subroutine calc_sigma_a_H2O2_analytic


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for H2CO
  !-----------------------------------------------------------------
  subroutine calc_sigma_a_H2CO_Tdependent(isp, T) ! inout
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer iz, swl, ewl
    real(dp) Tz, Tz_298, sigma0(nwl), sigma1(nwl)

    sigma_a(0,:,isp) = 1.0_dp ! existence flag

    ! Temperature gradient by Meller and Moortgat(2000)
    swl = nint(sigma_dat(-2,isp,1,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,1,1,1)) ! end wavelength
    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,1,isp) = dble(swl)
    if (sigma_a(-1,1,isp) < dble(ewl)) sigma_a(-1,1,isp) = dble(ewl)

    do iz = 1, nz
      Tz = T(iz)
      if (T(iz)<=223.0_dp) then 
        Tz = 223.0_dp
      else if (T(iz)>323.0_dp) then 
        Tz = 323.0_dp
      end if
      Tz_298 = Tz - 298.0_dp

      sigma0(swl:ewl) = sigma_dat(swl:ewl,isp,1,1,2)
      sigma1(swl:ewl) = sigma_dat(swl:ewl,isp,2,1,2)
      sigma_a(swl:ewl,iz,isp) = sigma0(swl:ewl) + sigma1(swl:ewl) * Tz_298
    end do

  end subroutine calc_sigma_a_H2CO_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for N2O5
  !-----------------------------------------------------------------
  subroutine calc_sigma_a_N2O5_Tdependent(isp, T)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer  iz, swl, ewl
    real(dp) Tz, invTz1000

    sigma_a(0,isp,:) = 1.0_dp ! existence flag

    ! Harwood et al., 1993
    swl = nint(sigma_dat(-2,isp,2,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,2,1,1)) ! end wavelength
    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,1,isp) = dble(swl)
    if (sigma_a(-1,1,isp) < dble(ewl)) sigma_a(-1,1,isp) = dble(ewl)

    do iz = 1, nz
      Tz = T(iz)
      if (T(iz)<=223.0_dp) then 
        Tz = 223.0_dp
      else if (T(iz)>295.0_dp) then 
        Tz = 295.0_dp
      end if
      invTz1000 = 1000.0_dp / Tz

      sigma_a(swl:ewl,iz,isp) = 10.0_dp**(sigma_dat(swl:ewl,isp,2,1,2) &
          &                             + sigma_dat(swl:ewl,isp,3,1,2) * invTz1000) &
          &                             * 1.0e-4_dp
    end do

  end subroutine calc_sigma_a_N2O5_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for HNO3
  !-----------------------------------------------------------------
  subroutine calc_sigma_a_HNO3_Tdependent(isp, T)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer iz, swl, ewl
    real(dp) Tz

    sigma_a(0,isp,:) = 1.0_dp ! existence flag

    ! Burkholder et al., 1993
    swl = nint(sigma_dat(-2,isp,2,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,2,1,1)) ! end wavelength
    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,1,isp) = dble(swl)
    if (sigma_a(-1,1,isp) < dble(ewl)) sigma_a(-1,1,isp) = dble(ewl)

    do iz = 1, nz
      Tz = T(iz)
      sigma_a(swl:ewl,iz,isp) = sigma_a(swl:ewl,iz,isp) &
          &                   * dexp(sigma_dat(swl:ewl,isp,2,1,2) * 1.0e-3_dp * (Tz-298.0_dp))
    end do

  end subroutine calc_sigma_a_HNO3_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for HO2NO2
  !-----------------------------------------------------------------
  subroutine calc_sigma_a_HO2NO2_Tdependent(isp, T)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(1:)
    integer iz, swl, ewl
    real(dp) Tz, invQ

    sigma_a(0,isp,:) = 1.0_dp ! existence flag

    ! Knight et al. 2002
    swl = nint(sigma_dat(-2,isp,2,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,2,1,1)) ! end wavelength
    if (sigma_a(-2,1,isp) > dble(swl)) sigma_a(-2,1,isp) = dble(swl)
    if (sigma_a(-1,1,isp) < dble(ewl)) sigma_a(-1,1,isp) = dble(ewl)

    do iz = 1, nz
      Tz = T(iz)
      if (Tz < 273.0_dp) Tz = 273.0_dp
      if (Tz > 343.0_dp) Tz = 343.0_dp
      invQ = 1.0_dp + dexp(-988.0_dp/(0.69_dp*Tz))

      sigma_a(swl:ewl,iz,isp) = ( sigma_dat(swl:ewl,isp,2,1,2) * invQ &
          &                     + sigma_dat(swl:ewl,isp,3,1,2) * (1.0_dp - invQ)) &
          &                   * 1.0e-24_dp
    end do

  end subroutine calc_sigma_a_HO2NO2_Tdependent


  !------------------------------------------------------------------------------------
  ! Analytic expression of quantum yield for O3 photodissociation by JPL 2015
  !------------------------------------------------------------------------------------
  subroutine calc_sigma_d_O3(isp, ich, T, OorO1D) ! inout
    implicit none
    integer,               intent(in)    :: isp, ich
    real(dp),              intent(in)    :: T(1:)
    character(len=*),      intent(in)    :: OorO1D 
    integer iz, iwl, l193, l225, l306, l329, l340, swl, ewl
    real(dp) Tclamp
    real(dp) X1, X2, X3, w1, w2, w3, A1, A2, A3, v1, v2, c, R, q1, q2, q3

    X1 = 304.225_dp; X2 = 314.957_dp; X3 = 310.737
    w1 = 5.576_dp  ; w2 = 6.601_dp  ; w3 = 2.187_dp
    A1 = 0.8036_dp ; A2 = 8.9061_dp ; A3 = 0.1192_dp
    v1 = 0.0_dp    ; v2 = 825.518_dp
    c  = 0.0765_dp
    R  = 0.0695_dp

    do iwl = 1, nwl-1
      if (     lambda(iwl) <  193.0_dp .and. 193.0_dp <= lambda(iwl+1)) then 
        l193 = iwl+1
      else if (lambda(iwl) <= 225.0_dp .and. 225.0_dp < lambda(iwl+1)) then 
        l225 = iwl
      else if (lambda(iwl) <= 306.0_dp .and. 306.0_dp < lambda(iwl+1)) then 
        l306 = iwl
      else if (lambda(iwl) <= 329.0_dp .and. 329.0_dp < lambda(iwl+1)) then 
        l329 = iwl
      else if (lambda(iwl) <= 340.0_dp .and. 340.0_dp < lambda(iwl+1)) then 
        l340 = iwl
      end if
    end do 

    swl = nint(sigma_a(-2,1,isp)) ! start wavelength
    ewl = nint(sigma_a(-1,1,isp)) ! end wavelength
    if (sigma_d(-2,1,ich) > dble(swl)) sigma_d(-2,:,ich) = dble(swl)
    if (sigma_d(-1,1,ich) < dble(ewl)) sigma_d(-1,:,ich) = dble(ewl)

    !   Evaluation: J. B. Burkholder, S. P. Sander, J. Abbatt, J. R. Barker, R. E. Huie, C. E. Kolb, 
    !               M. J. Kurylo, V. L. Orkin, D. M. Wilmouth, and P. H. Wine 
    !               "Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 18," 
    !               JPL Publication 15-10, Jet Propulsion Laboratory, Pasadena, 2015 http://jpldataeval.jpl.nasa.gov.
    
    do iz = 1, nz

      Tclamp = T(iz)
      if (Tclamp < 200.0_dp) Tclamp = 200.0_dp
      if (Tclamp > 320.0_dp) Tclamp = 320.0_dp

      qy_O1D(l193:l225) = 1.37e-2_dp * lambda(l193:l225) - 2.16_dp
      if (OorO1D == 'O(1D)') sigma_d(l193:l225,iz,ich) = sigma_a(l193:l225,iz,isp) * qy_O1D(l193:l225)
      if (OorO1D == 'O')     sigma_d(l193:l225,iz,ich) = sigma_a(l193:l225,iz,isp) * (1.0_dp-qy_O1D(l193:l225))

      qy_O1D(l225+1:l306-1) = 0.9_dp
      if (OorO1D == 'O(1D)') sigma_d(l225+1:l306-1,iz,ich) = sigma_a(l225+1:l306-1,iz,isp) * qy_O1D(l225+1:l306-1)
      if (OorO1D == 'O')     sigma_d(l225+1:l306-1,iz,ich) = sigma_a(l225+1:l306-1,iz,isp) * (1.0_dp-qy_O1D(l225+1:l306-1))

      ! Analytic expression of O3 -> O1D at wavelength range 306-328 nm and at temperature range 200-320 K
      !   based on the review of Matsumi et al. (2002) doi:10.1029/2001JD000510.
      q1 = dexp(-v1/(R*Tclamp))
      q2 = dexp(-v2/(R*Tclamp))
      q3 = Tclamp/300.0_dp
      qy_O1D(l306:l329-1) = q1/(q1+q2)              * A1 * dexp(-((X1-lambda(l306:l329-1))/w1)**4.0_dp) &
        &                 + q2/(q1+q2) * q3**2.0_dp * A2 * dexp(-((X2-lambda(l306:l329-1))/w2)**2.0_dp) &
        &                 +              q3**1.5_dp * A3 * dexp(-((X3-lambda(l306:l329-1))/w3)**2.0_dp) &
        &                 + c
      if (OorO1D == 'O(1D)') sigma_d(l306:l329-1,iz,ich) = sigma_a(l306:l329-1,iz,isp) * qy_O1D(l306:l329-1)
      if (OorO1D == 'O')     sigma_d(l306:l329-1,iz,ich) = sigma_a(l306:l329-1,iz,isp) * (1.0_dp-qy_O1D(l306:l329-1))

      qy_O1D(l329:l340) = 0.08_dp
      if (OorO1D == 'O(1D)') sigma_d(l329:l340,iz,ich) = sigma_a(l329:l340,iz,isp) * qy_O1D(l329:l340)
      if (OorO1D == 'O')     sigma_d(l329:l340,iz,ich) = sigma_a(l329:l340,iz,isp) * (1.0_dp-qy_O1D(l329:l340))

      qy_O1D(l340+1:nwl) = 0.0_dp
      if (OorO1D == 'O(1D)') sigma_d(l340+1:nwl,iz,ich) = sigma_a(l340+1:nwl,iz,isp) * qy_O1D(l340+1:nwl)
      if (OorO1D == 'O')     sigma_d(l340+1:nwl,iz,ich) = sigma_a(l340+1:nwl,iz,isp) * (1.0_dp-qy_O1D(l340+1:nwl))

    end do

  end subroutine calc_sigma_d_O3

  
  !------------------------------------------------------------------------------------
  ! Analytic expression of quantum yield for O3 photodissociation by JPL 2015
  !------------------------------------------------------------------------------------
  subroutine calc_sigma_d_H2CO(isp, ich, T, HorCO)
    implicit none
    integer,               intent(in)    :: isp, ich
    real(dp),              intent(in)    :: T(1:)
    character(len=*),      intent(in)    :: HorCO
    integer i, iz, iwl, l250, l330, l338, l360, swl, ewl
    real(dp) l, Tclamp, P, kB
    real(dp) a(0:4)

    sigma_d(0,:,ich) = 1.0_dp ! existence flag

    kB = 1.38064852e-23_dp

    swl = nint(sigma_a(-2,1,isp)) ! start wavelength
    ewl = nint(sigma_a(-1,1,isp)) ! end wavelength
    if (sigma_d(-2,1,isp) > dble(swl)) sigma_d(-2,1,isp) = dble(swl)
    if (sigma_d(-1,1,isp) < dble(ewl)) sigma_d(-1,1,isp) = dble(ewl)

    a(0) = 557.95835182_dp
    a(1) = -7.31994058026_dp
    a(2) = 0.03553521598_dp
    a(3) = -7.54849718e-5_dp
    a(4) = 5.91001021e-8_dp

    do iwl = 1, nwl-1
      if (     lambda(iwl) <  250.0_dp .and. 250.0_dp <= lambda(iwl+1)) then 
        l250 = iwl+1
      else if (lambda(iwl) <= 330.0_dp .and. 330.0_dp < lambda(iwl+1)) then 
        l330 = iwl
      else if (lambda(iwl) <= 338.0_dp .and. 338.0_dp < lambda(iwl+1)) then 
        l338 = iwl
      else if (lambda(iwl) <= 360.0_dp .and. 360.0_dp < lambda(iwl+1)) then 
        l360 = iwl
      end if
    end do 

    do iwl = l250, l360
      l = lambda(iwl)
      if (iwl <= l338) then 
        qy_H(iwl) = 0.0_dp
        do i = 0, 4
          qy_H(iwl) = qy_H(iwl) + a(i) * l**dble(i)
        end do 
      else if (iwl > l338) then 
        qy_H(iwl) = 0.0_dp
      end if 
    end do 

    !   Evaluation: J. B. Burkholder, S. P. Sander, J. Abbatt, J. R. Barker, R. E. Huie, C. E. Kolb, 
    !               M. J. Kurylo, V. L. Orkin, D. M. Wilmouth, and P. H. Wine 
    !               "Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 18," 
    !               JPL Publication 15-10, Jet Propulsion Laboratory, Pasadena, 2015 http://jpldataeval.jpl.nasa.gov.
    
    do iz = 1, nz

      Tclamp = T(iz)
      if (Tclamp < 220.0_dp) Tclamp = 220.0_dp
      if (Tclamp > 300.0_dp) Tclamp = 300.0_dp

      P = 0.006_dp!var%n_tot(iz) * kB * T(iz) / 101325.0_dp ! [atm]

      qy_CO_300K(l250:l360) = sigma_dat(l250:l360,ich,1,1,2) * 1.0e-24_dp ! quantum yield at 300K

      ! Temperature and pressure dependent yield of H2 + CO
      a_300K(l250:l360) = 1.0_dp / qy_CO_300K(l250:l360) - 1.0_dp / (1.0_dp - qy_H(l250:l360))
      a_T(l250:l360)    = a_300K(l250:l360) * ( 1.0_dp + 0.05_dp * (lambda(l250:l360)-329.0_dp)*((300.0_dp-Tclamp)/80.0_dp) )
      qy_CO_T(l250:l360) = 1.0_dp / ( 1.0_dp / (1.0_dp-qy_H(l250:l360)) + a_T(l250:l360)*P )

      if (HorCO == 'H')  sigma_d(l250:l360,iz,ich) = sigma_a(l250:l360,iz,isp) * qy_H(l250:l360)
      if (HorCO == 'CO') sigma_d(l250:l360,iz,ich) = sigma_a(l250:l360,iz,isp) * qy_CO_T(l250:l360)

      !print *, l, qy_H, qy_CO_300K, qy_CO_T
      
    end do

  end subroutine calc_sigma_d_H2CO


  function get_species_index(species_in)
    implicit none
    character(len=*), intent(in) :: species_in
    integer get_species_index
    integer isp

    get_species_index = 0
    loop: do isp = 1, nsp
      if( trim(ADJUSTL(species(isp))) == trim(ADJUSTL(species_in)) ) then
        get_species_index = isp
        exit loop
      end if
    end do loop

  end function get_species_index


  function get_reaction_index(in_reactant, in_product)
    !
    ! this function should not be used for searching reactions other than photodissociation.
    !
    implicit none
    character(len=*), intent(in) :: in_reactant(:), in_product(:)
    integer get_reaction_index
    integer i, ich
    character(len=256) reactants(20),  products(20)

    get_reaction_index = 0

    do ich = 1, nch
      ! check if the reaction type matches the given type

      if (reaction_type(ich) == 2) then
        reactants = ''
        do i = 1, reactant_list(ich,0)
          reactants(i) = trim(species(reactant_list(ich,i)))
        end do
        products = ''
        do i = 1, product_list(ich,0)
          products(i) = trim(species(product_list(ich,i)))
        end do

        ! if the number of products is 1: 1 case
        if (product_list(ich,0) == 1 .and. size(in_product) == 1) then 
          if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1)))) then 
              get_reaction_index = ich; exit
          end if
        ! if the number of products is 2: 2 cases
        else if (product_list(ich,0) == 2 .and. size(in_product) == 2) then 
          if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1)))) then 
              get_reaction_index = ich; exit
          end if
        ! if the number of products is 3: 6 cases
        else if (product_list(ich,0) == 3 .and. size(in_product) == 3) then 
          if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2)))) then 
              get_reaction_index = ich; exit
          end if

        ! if the number of products is 4: 24 cases
        else if (product_list(ich,0) == 4 .and. size(in_product) == 4) then 
          if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then 
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(4)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(3)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(2)))) then
              get_reaction_index = ich; exit
          else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
            & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(4))) &
            & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
            & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2))) &
            & .and. trim(ADJUSTL(products(4)))==trim(ADJUSTL(in_product(1)))) then
              get_reaction_index = ich; exit
          end if

        end if

      end if

    end do 

  end function get_reaction_index


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


end module p__photolysis_rate
