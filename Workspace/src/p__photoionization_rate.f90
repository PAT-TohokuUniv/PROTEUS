 module p__photoionization_rate

  implicit none
  integer(4), parameter :: sp = 4, dp = 8


  private

  ! Module-level variables ---------------------------------------------------------------------------
  integer,                         private :: nch_photoionization, nwl, nsp, nz, nch
  character(len=256), allocatable, private :: species(:)
  real(dp),           allocatable, private :: sigma_dat(:,:,:,:,:), sigma_a(:,:), sigma_i(:,:)
  real(dp),           allocatable, private :: lambda(:), dlambda(:), solar_flux(:), mass(:)
  integer,            allocatable, private :: reaction_type(:), product_list(:,:), reactant_list(:,:)

  ! For optical depth calculation 
  real(dp),           allocatable, private :: cln(:,:), tau(:,:), I_z(:,:)

  ! Public interfaces --------------------------------------------------------------------------------
  public :: p__photoionization_rate__ini, p__photoionization_rate__fin, &
    &       load_cross_section_dat, get_cross_section, photoionization_rate

contains


  !=======================================================================================================================
  !
  !  Allocate memory for the cross-section data array.
  !
  !=======================================================================================================================
  subroutine p__photoionization_rate__ini(nwl_in, nz_in, nsp_in, nch_in,                                            & ! in
    &                                     species_in, mass_in, reaction_type_in, reactant_list_in, product_list_in, & ! in
    &                                     lambda_in, dlambda_in, solar_flux_in                                      ) ! in
    implicit none
    integer,  intent(in) :: nwl_in, nsp_in, nz_in, nch_in
    integer,  intent(in) :: reaction_type_in(1:), reactant_list_in(1:,0:), product_list_in(1:,0:)
    real(dp), intent(in) :: lambda_in(1:), dlambda_in(1:), solar_flux_in(1:), mass_in(1:)
    character(len=*), intent(in) :: species_in(1:)
    integer ich, nch_photoionization_priv

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

    nch_photoionization_priv = 0
    do ich = 1, nch
      if (reaction_type(ich) == 1) then
        nch_photoionization_priv = nch_photoionization_priv + 1
      end if
    end do

    ! Allocate the module-level array with the required dimensions
    allocate(sigma_dat(-2:nwl, 1:nsp+nch_photoionization_priv, 0:1, 0:10, 0:2))
    allocate(sigma_a(-2:nwl, 1:nsp))
    allocate(sigma_i(-2:nwl, 1:nch_photoionization_priv))
    allocate(lambda(1:nwl), dlambda(1:nwl), solar_flux(1:nwl))
    allocate(cln(nz, nsp), tau(nwl, nz), I_z(nwl, nz))

    ! Initialization
    sigma_dat   = 0.0_dp
    sigma_a     = 0.0_dp
    sigma_i     = 0.0_dp

    lambda(1:nwl) = lambda_in(1:nwl)
    dlambda(1:nwl) = dlambda_in(1:nwl)
    solar_flux(1:nwl) = solar_flux_in(1:nwl)

  end subroutine p__photoionization_rate__ini

  
  !=======================================================================================================================
  !
  !  Deallocate the module-level array to free memory.
  !
  !=======================================================================================================================
  subroutine p__photoionization_rate__fin()
    if (allocated(sigma_dat))     deallocate(sigma_dat)
    if (allocated(sigma_a))       deallocate(sigma_a)
    if (allocated(sigma_i))       deallocate(sigma_i)
    if (allocated(lambda))        deallocate(lambda)
    if (allocated(dlambda))       deallocate(dlambda)
    if (allocated(solar_flux))    deallocate(solar_flux)
    if (allocated(reaction_type)) deallocate(reaction_type)
    if (allocated(product_list))  deallocate(product_list)
    if (allocated(reactant_list)) deallocate(reactant_list)
    if (allocated(species))       deallocate(species)
    if (allocated(mass))          deallocate(mass)
  end subroutine p__photoionization_rate__fin


  !=======================================================================================================================
  !
  !  Load photo-absorption and ionization cross section data from files and bin them to the simulation wavelength grid
  !
  !=======================================================================================================================
  subroutine load_cross_section_dat(dirname) ! in
    implicit none
    character(len=*), intent(in) :: dirname
    integer isp, ich
    character(len=256) fname

    ! initialization
    nch_photoionization = 0
    sigma_dat = 0.0_dp
    sigma_dat(-2,:,:,1,1) = dble(nwl)
    sigma_dat(-1,:,:,1,1) = dble(1.0_dp)

    ! Template --------------------------------------------------------
    !  isp = get_species_index('CO2') 
    !------------------------------------------------------------------
    !
    !  fname = trim(ADJUSTL(dirname))//'/sigma_*.dat'
    !  call read_binning_data(fname,        & ! fixed, do not modify these three
    !    &                    unit1, unit2, & ! Units of column 1 and 2 of the datafile
    !    &                    flag, index   ) ! flag: 'a' for absorption, 'i' for ionization, index: species or reaction index, respectively
    !   unit1: unit of column 1 of the datafile, 'nm', 'cm', 'eV', 'A', or 'cm-1'.
    !   unit2: unit of column 2 of the datafile, 'cm2', 'm2', or 'Mb'.
    
    ! H2 data ---------------------------------------------------------
    isp = get_species_index('H2')
    !------------------------------------------------------------------
    
    fname = trim(ADJUSTL(dirname))//'/Hx/sigma_a_H2_eV.dat'
    call read_binning_data(fname, &
      &                    'eV', 'Mb', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Hx/sigma_i_H2_H2+_Kossmann1989.dat'
    ich = get_reaction_index(['H2'], ['H2+', 'e- '])
    call read_binning_data(fname, &
      &                    'eV', 'Mb', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/Hx/sigma_i_H2_H+.dat'
    ich = get_reaction_index(['H2'], ['H+', 'e-', 'H '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! H data ----------------------------------------------------------
    isp = get_species_index('H')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/Hx/sigma_a_H.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Hx/sigma_i_H_H+.dat'
    ich = get_reaction_index(['H'], ['H+', 'e-'])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! He data ---------------------------------------------------------
    isp = get_species_index('He')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/He/sigma_a_He.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/He/sigma_i_He_He+.dat'
    ich = get_reaction_index(['He'], ['He+', 'e- '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! CH4 data --------------------------------------------------------
    isp = get_species_index('CH4')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_a_CH4.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_i_CH4_CH4+.dat'
    ich = get_reaction_index(['CH4'], ['CH4+', 'e-  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_i_CH4_CH3+.dat'
    ich = get_reaction_index(['CH4'], ['CH3+', 'e-  ', 'H   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_i_CH4_CH2+.dat'
    ich = get_reaction_index(['CH4'], ['CH2+', 'e-  ', 'H2  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_i_CH4_CH+.dat'
    ich = get_reaction_index(['CH4'], ['CH+', 'e- ', 'H2 ', 'H  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C1/sigma_i_CH4_H+.dat'
    ich = get_reaction_index(['CH4'], ['H+ ', 'e- ', 'CH3'])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! C2H2 data -------------------------------------------------------
    isp = get_species_index('C2H2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H2.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)
      
    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H2_C2H2+.dat'
    ich = get_reaction_index(['C2H2'], ['C2H2+', 'e-   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H2_C2H+.dat'
    ich = get_reaction_index(['C2H2'], ['C2H+', 'e-  ', 'H   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! C2H4 data -------------------------------------------------------
    isp = get_species_index('C2H4')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H4.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H4_C2H4+.dat'
    ich = get_reaction_index(['C2H4'], ['C2H4+', 'e-   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H4_C2H3+.dat'
    ich = get_reaction_index(['C2H4'], ['C2H3+', 'e-   ', 'H    '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H4_C2H2+.dat'
    ich = get_reaction_index(['C2H4'], ['C2H2+', 'e-   ', 'H2   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! C2H6 data -------------------------------------------------------
    isp = get_species_index('C2H6')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_a_C2H6.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/C2/sigma_i_C2H6_C2H6+.dat'
    ich = get_reaction_index(['C2H6'], ['C2H6+', 'e-   '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich) 

    ! CO2 data --------------------------------------------------------
    isp = get_species_index('CO2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_a_CO2.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO2_CO2+.dat'
    ich = get_reaction_index(['CO2'], ['CO2+', 'e-  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO2_CO+.dat'
    ich = get_reaction_index(['CO2'], ['CO+', 'e- ', 'O  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO2_O+.dat'
    ich = get_reaction_index(['CO2'], ['O+(4S)', 'e-    ', 'CO    '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO2_C+.dat'
    ich = get_reaction_index(['CO2'], ['C+', 'e-', 'O2'])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! CO data ---------------------------------------------------------
    isp = get_species_index('CO')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_a_CO.dat'
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'a', isp)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO_CO+.dat'
    ich = get_reaction_index(['CO'], ['CO+', 'e- ', 'O  '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO_C+.dat'
    ich = get_reaction_index(['CO'], ['C+', 'e-', 'O '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/COx/sigma_i_CO_O+.dat'
    ich = get_reaction_index(['CO'], ['O+(4S)', 'e-    ', 'C     '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! N2 data ---------------------------------------------------------
    isp = get_species_index('N2')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/Nx/sigma_a_N2.dat'
    call read_binning_data(fname, &
      &                   'A', 'cm^2', &
      &                   'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Nx/sigma_i_N2_N2+.dat'
    ich = get_reaction_index(['N2'], ['N2+', 'e- '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    fname = trim(ADJUSTL(dirname))//'/Nx/sigma_i_N2_N+.dat'
    ich = get_reaction_index(['N2'], ['N+   ', 'e-   ', 'N(2D)'])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

    ! O data ----------------------------------------------------------
    isp = get_species_index('O')
    !------------------------------------------------------------------

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_a_O.dat'
    call read_binning_data(fname, &
      &                   'A', 'cm^2', &
      &                   'a', isp)

    fname = trim(ADJUSTL(dirname))//'/Ox/sigma_i_O_O+.dat'
    ich = get_reaction_index(['O'], ['O+(4S)', 'e-    ', 'O     '])
    call read_binning_data(fname, &
      &                    'A', 'cm^2', &
      &                    'i', ich)

  end subroutine load_cross_section_dat


  !=======================================================================================================================
  !
  !     Get cross section for absorption (sigma_a) or photoionization (sigma_i) 
  !
  !=======================================================================================================================
  subroutine get_cross_section(flag, outflag)
    implicit none       
    integer,    intent(in) :: outflag
    character(len=*), intent(in) :: flag
    integer, parameter :: photoionization = 2 ! reaction type index for photoionization
    integer iwl, iz, isp, ich, jch, swl, ewl
    character(len=256) fname, dirname, command
    real(dp), parameter :: pi = dacos(-1.0_dp)

    !----------------------------------------------------------------------------------------------------
    ! UV photoabsorption cross section

    if (flag == 'absorption') then

       sigma_a(0,:) = 0 ! reset cross section existence flag
      
       do isp = 1, nsp

        if (nint(sigma_dat(0,isp,0,0,0)) == 0) cycle ! skip if no data

        ! For normal species ------------------------------------------

        call calc_sigma_a(isp) ! in

        ! For exceptional species -------------------------------------

       end do !isp


      ! if cross sections are less than 0, they are set to 0 ------
      where (sigma_a(:,:) < 0.0_dp)
        sigma_a(:,:) = 0.0_dp
      end where

      ! output cross sections ------
      if (outflag == 1) then 
        iz = 1 ! select altitude grid
        do isp = 1, nsp
          if (nint(sigma_a(0,isp)) == 1) then 
            fname = './EUV/xsect/absorption/'//trim(ADJUSTL(species(isp)))//'.dat'
            open(11, file = fname, status = 'replace' )
            swl = nint(sigma_a(-2,isp))
            ewl = nint(sigma_a(-1,isp))
            !print *, isp, trim(species(isp)), swl, ewl
            if (swl >= 1) then 
              do iwl = swl, ewl
                write(11, *) lambda(iwl), sigma_a(iwl,isp)
              end do
            end if
            close(11)
          end if
        end do 
      end if

    end if ! end if flag == absorption

    !----------------------------------------------------------------------------------------------------
    ! UV photoionization cross section

    if (flag == 'photoionization') then

      do ich = 1, nch_photoionization

        ! For normal photoionization reactions -----------------------------

        call calc_sigma_i(ich) ! in

        ! For exceptional photoionization reactions ------------------------

      end do !ich


      ! if cross sections are less than 0, they are set to 0 ------
      where (sigma_i(:,:)<0.0_dp)
        sigma_i(:,:) = 0.0_dp
      end where

      ! output cross sections ------
      if (outflag == 1) then 
        iz = 1 ! select altitude grid
        do jch = 1, nch_photoionization
          ich = nint(sigma_dat(0,nsp+jch,0,0,0)) ! reaction index
          if (reaction_type(ich) == 1) then
            if (nint(sigma_i(0,jch)) == 1) then 
              if (product_list(ich,0)==1) then 
                fname = './EUV/xsect/ionization/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'.dat'
              else if (product_list(ich,0)==2) then 
                fname = './EUV/xsect/ionization/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'.dat'
              else if (product_list(ich,0)==3) then 
                fname = './EUV/xsect/ionization/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,3))))//'.dat'
              else if (product_list(ich,0)==4) then 
                fname = './EUV/xsect/ionization/'&
                  &     //trim(ADJUSTL(species(reactant_list(ich,1))))//'_plus_hv_to_'&
                  &     //trim(ADJUSTL(species(product_list(ich,1))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,2))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,3))))//'_plus_'&
                  &     //trim(ADJUSTL(species(product_list(ich,4))))//'.dat'
              end if
              open(11, file = fname, status = 'replace' )
              swl = nint(sigma_i(-2,jch))
              ewl = nint(sigma_i(-1,jch))
              !print *, ich, swl, ewl, trim(fname)
              if (swl >= 1) then 
                do iwl = swl, ewl
                  write(11, *) lambda(iwl), sigma_i(iwl,jch)
                end do
              end if
              close(11)
            end if
          end if
        end do 
      end if

    end if ! end if flag == photoionization


  end subroutine get_cross_section


  !=======================================================================================================================
  !
  !    Calculate optical depth and photoionization rate 
  !
  !=======================================================================================================================
  subroutine photoionization_rate(alt, dalt, density, T, sza, & ! in
    &                             rplanet, mplanet,           & ! in
    &                             J_rate                      ) ! out
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

      if (nint(sigma_a(0,isp)) == 0) cycle ! skip if no data
      swl = nint(sigma_a(-2,isp))
      ewl = nint(sigma_a(-1,isp))

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
            &             + cln_Ch * sigma_a(swl:ewl,isp)
        else if (alt(iz) <= hmin) then
          tau(1:nwl,iz) = 1.0e10_dp
        end if
      end do ! iz

    end do ! isp

    do iz  = 1, nz
      I_z(1:nwl,iz) = solar_flux(1:nwl) * dlambda(1:nwl) * dexp(-tau(1:nwl,iz)) 
    end do

    ! photoionization rate ---------------------------
    do ich = 1, nch_photoionization

      if (nint(sigma_i(0,ich)) == 0) cycle ! skip if no data
      swl = nint(sigma_i(-2,ich))
      ewl = nint(sigma_i(-1,ich))
      jch = nint(sigma_dat(0,nsp+ich,0,0,0)) ! reaction index

      do iz = 1, nz
        J_rate(iz,jch) = dot_product(I_z(swl:ewl,iz), sigma_i(swl:ewl,ich))
      end do

    end do ! ich


  end subroutine photoionization_rate


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
  subroutine read_binning_data(fname,        & ! inout, in
    &                          unit1, unit2, & ! in
    &                          dtype, id_in  ) ! inout, in, in
    implicit none
    integer,          intent(in)    :: id_in
    character(len=*), intent(in)    :: unit1, unit2, fname, dtype
    integer i, nh, il, nl, swl, ewl, idtype, ndata, id, id_flag
    real(dp), allocatable :: idata(:,:), odata(:)

    if (id_in >= 1) then 

      ndata = 1

      if (dtype == 'a') then 
        idtype = 2**0 ! absorption
        id = id_in
      else if (dtype == 'a-excep') then 
        idtype = 2**3 ! exception for absorption
        id = id_in
      else if (dtype == 'i') then 
        idtype = 2**6 ! ionization
      else if (dtype == 'i-excep') then 
        idtype = 2**9 ! exception for ionization
      else 
        print *, ''
        print *, 'error!'
        print *, '  The dtype "'//trim(adjustl(dtype))//'" is not recognized.'
        stop
      end if

      id_flag = 0
      if (idtype >= 10) then
        do i = nsp+1, nsp+nch_photoionization
          if (nint(sigma_dat(0,i,0,0,0)) == id_in) then 
            id = i
            id_flag = 1
          end if
        end do 
      end if
      if (idtype >= 10 .and. id_flag == 0) then 
        nch_photoionization = nch_photoionization + 1
        id = nsp + nch_photoionization
      end if
      
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
          if (unit1/='nm' .and. unit1/='A' .and. unit1/='/cm' .and. unit1/='cm^-1' .and. unit1/='eV') then 
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
          else if (unit1 == 'eV') then 
            idata(il,1) = 1239.8_dp/idata(il,1) ! [eV -> nm]
          end if

          ! Unit conversion: column 2
          if (unit2/='m^2' .and. unit2/='cm^2' &
          &  .and. unit2/='/cm^2/s/nm' .and. unit2/='/m^2/s/nm' .and. unit2/='W/m^2/nm' .and. unit2/='Mb') then 
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
          else if (unit2 == 'Mb') then 
            idata(il,2) = idata(il,2) * 1.0e-22_dp ! [Mb -> m^2]
          end if

        end do
      close(11)

      call binning_cross_section(nl, idata,      &
        &                        odata, swl, ewl )
      if (sigma_dat(-2,id,1,1,1) > dble(swl)) sigma_dat(-2,id,1,1,1) = dble(swl) ! start wavelength
      if (sigma_dat(-1,id,1,1,1) < dble(ewl)) sigma_dat(-1,id,1,1,1) = dble(ewl) ! end wavelength
      sigma_dat(0,id,1,1,1) = 1.0_dp ! existence flag
      if (dtype=='a' .or. dtype=='i' .or. dtype=='qy')sigma_dat(1,id,0,1,1) = sigma_dat(1,id,0,1,1) + dble(ndata)
      sigma_dat(0,id,1,1,0) = sigma_dat(0,id,1,1,0) + dble(idtype)
      sigma_dat(swl:ewl,id,1,1,2) = odata(swl:ewl)
      sigma_dat(0,id,0,0,0) = dble(id_in) ! reaction index

      deallocate(idata)
      deallocate(odata)

      write(*,*) 'done.'

    end if

  end subroutine read_binning_data


  !------------------------------------------------------------
  ! Absorption cross section data are automatically adapted.
  !------------------------------------------------------------
  subroutine calc_sigma_a(isp) ! inout
    implicit none
    integer, intent(in) :: isp
    integer ndata, swl, ewl, dtype

    sigma_a(0,isp) = sigma_dat(0,isp,1,1,1)

    swl = nint(sigma_dat(-2,isp,1,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,isp,1,1,1)) ! end wavelength
    sigma_a(-2,isp) = dble(swl)
    sigma_a(-1,isp) = dble(ewl)
    ndata = nint(sigma_dat(1,isp,0,1,1)) ! number of data columns
    dtype = nint(sigma_dat(0,isp,1,1,0)) ! data type

    if (mod(dtype,2**3) /= 0) then 
      if (mod(ndata,1000) == 1) then 
        sigma_a(swl:ewl,isp) = sigma_dat(swl:ewl,isp,1,1,2)
      end if
    end if

  end subroutine calc_sigma_a


  !------------------------------------------------------------
  ! Absorption cross section data are automatically adapted.
  !------------------------------------------------------------
  subroutine calc_sigma_i(ich) ! inout
    implicit none
    integer, intent(in) :: ich
    integer ndata, swl, ewl, dtype

    sigma_i(0,ich) = sigma_dat(0,nsp+ich,1,1,1)

    swl = nint(sigma_dat(-2,nsp+ich,1,1,1)) ! start wavelength
    ewl = nint(sigma_dat(-1,nsp+ich,1,1,1)) ! end wavelength
    sigma_i(-2,ich) = dble(swl)
    sigma_i(-1,ich) = dble(ewl)

    ndata = nint(sigma_dat(1,nsp+ich,0,1,1)) ! number of data columns
    dtype = nint(sigma_dat(0,nsp+ich,1,1,0)) ! data type

    if (mod(dtype,2**9) /= 0) then 
      if (mod(ndata,1000) == 1) then 
        sigma_i(swl:ewl,ich) = sigma_dat(swl:ewl,nsp+ich,1,1,2)
      end if
    end if

  end subroutine calc_sigma_i


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
    ! this function should not be used for searching reactions other than photoionization.
    implicit none
    character(len=*), intent(in) :: in_reactant(:), in_product(:)
    integer get_reaction_index
    integer i, ich
    character(len=256) reactants(20),  products(20)

    get_reaction_index = 0

    do ich = 1, nch
      ! check if the reaction type matches the given type

      if (reaction_type(ich) == 1) then
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


end module p__photoionization_rate
