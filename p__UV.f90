module p__UV
  use v__tdec,    only : var_, grd_, cst_, xct_, spl_, flx_
  use c__prm,     only : c__prm__ini
  use p__search,  only : p__search_reactant, p__search_product, ch_identify, sp_index


  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__UV_flux, p__UV_cross_section_dat, p__UV_cross_section_exe, &
    &       p__UV_EUV_xct_convergence_exe

contains


  !-------------------------------------------------
  ! solar UV model: Woods et al., 2016
  !-------------------------------------------------
  subroutine p__UV_flux(spl, grd, cst, & ! in
    &                   var, xct, flx  ) ! inout
    implicit none
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(cst_),           intent(in)    :: cst
    type(var_),           intent(inout) :: var
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer iwl, i, nwl, nwli, n_wl_cases, il, nh, nl, label_solarflux
    real(dp), allocatable :: dlambda(:,:), rdata(:,:)
    character(len=256) fname

    label_solarflux = 1

    !---------------------------------------------------
    ! Defining UV wavelength bin
    !---------------------------------------------------
    ! bin definition is given by the following file.
    fname = './UV/bin/1nm.dat'
    !fname = './UV/bin/CO2_isotope.dat'

    ! The bin datafile must have a following format.
    !
    !   swl1  ewl1  bin1
    !   swl2  ewl2  bin2
    !   swl3  ewl3  bin3
    !    :     :     :
    !
    !   swl*: start of wavelangth range [nm]
    !   ewl*: end of wavelength range [nm]
    !   bin*: wavelength bin corresponding to the above wavelength range [nm]
    ! 
    ! ** Each wavelength range must not overlap. **


    !---------------------------------------------------
    ! wavelength bins are automatically adapted
    !---------------------------------------------------
    call p__count_header(nh, fname)
    call p__count_lines(nl, fname)
    n_wl_cases = nl
    allocate(dlambda(n_wl_cases,3))
    open(11, file = fname, status = 'unknown')
      do il = 1, nh; read(11,*); end do
      do il = 1, nl
        read(11,*) dlambda(il,1), dlambda(il,2), dlambda(il,3)
      end do 
    close(11)

    flx%nwl_UV = 0
    do i = 1, n_wl_cases
      nwli = nint( (dlambda(i,2)-dlambda(i,1))/(dlambda(i,3)) ) + 1
      flx%nwl_UV = flx%nwl_UV + nwli 
    end do

    allocate(flx%solar_UV(flx%nwl_UV), flx%lambda_UV(flx%nwl_UV), flx%dlambda_UV(flx%nwl_UV))
    allocate(var%tau_UV(flx%nwl_UV,grd%nz),var%tau_UV_subsolar(flx%nwl_UV,grd%nz),var%tau_RS(flx%nwl_UV,grd%nz))
    allocate(var%I_UV(flx%nwl_UV,grd%nz))
    allocate(xct%sigma_a_UV(flx%nwl_UV,grd%nz,spl%nsp), xct%sigma_d_UV(flx%nwl_UV,grd%nz,spl%nch))
    allocate(xct%wlrange_a_UV(2,spl%nsp), xct%wlrange_d_UV(2,spl%nch))
    allocate(xct%sigma_a_UV_EUV(flx%nwl_UV,spl%nsp))
    allocate(xct%label_sigma_a_UV(spl%nsp), xct%label_sigma_a_UV_EUV_conv(spl%nsp))
    allocate(xct%sigma_a_RS(flx%nwl_UV,spl%nsp), xct%label_sigma_a_RS(spl%nsp))

    nwl = 0
    do i = 1, n_wl_cases
      nwli = nint( (dlambda(i,2)-dlambda(i,1))/(dlambda(i,3)) ) + 1
      do iwl = 1, nwli
        flx%dlambda_UV(nwl+iwl) = dlambda(i,3)
        flx%lambda_UV(nwl+iwl)  = dlambda(i,1) + dlambda(i,3)*dble(iwl-1) ![nm]
      end do
      nwl = nwl + nwli
    end do

    deallocate(dlambda)

    !---------------------------------------------------
    ! solar flux 
    !---------------------------------------------------
    !units: photons m^-2 s^-1 nm^-1

    !fname = './UV/flux_photon_1-1000nm.txt'
    !call p__UV_read_data(xct, flx, fname, label_solarflux, &
    !  &                  rdata, 1, '', 1)
    !flx%solar_UV(:) = rdata(:,2) / (flx%orbit * flx%orbit) 

    !fname = './UV/Claire_2012_Sun3.8Ga.txt'
    !write(*,fmt='(a)') '  Reading solar photon flux..'
    !call p__UV_read_data(xct, flx, fname, label_solarflux, &
    !  &                  rdata, 1, '')
    !do iwl = 1, flx%nwl_UV
    !  flx%solar_UV(iwl) = rdata(iwl,2) / (flx%orbit * flx%orbit) * 1.0e4_dp ![/cm^2/s] -> [/m^2/s]
    !end do

    fname = './UV/ref_solar_irradiance_whi-2008_ver2.dat'
    write(*,fmt='(a)') '  Reading solar photon flux..'
    call p__UV_read_data(xct, flx, fname, label_solarflux, &
      &                  rdata, 4, '')
    do iwl = 1, flx%nwl_UV
      flx%solar_UV(iwl) = rdata(iwl,4) * flx%lambda_UV(iwl) * 1.0e-9_dp / cst%planck / cst%c &
        &               / (flx%orbit * flx%orbit) 
    end do 

  end subroutine p__UV_flux


  !-------------------------------------------------
  ! cross section data for UV 
  !-------------------------------------------------
  subroutine p__UV_cross_section_dat(spl,    & ! in
    &                                xct, flx) ! inout
    implicit none
    type(spl_),           intent(in)    :: spl
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer i, iwl, isp
    character(len=256) fname

    ! Template --------------------------------------------------------
    !  isp = sp_index(spl,'CO2') 
    !------------------------------------------------------------------
    !
    !  fname = './UV/xct_data/sigma_*.dat' 
    !  print *, 'Reading absorption cross section for *..'
    !  allocate(xct%T_*(N)) ! it should be needed if there are N (more than 2) columns (temperature dependence) 
    !                         in the cross section data. 
    !                         Otherwise, skip this line.
    !  call p__UV_read_data(xct, flx, fname, isp, & ! fixed, do not modify these three
    !    &                  a, b, c ) ! You should modify only these three inputs.
    !                         a: xct%variable_name, e.g.) xct%sigma_a_CO2 
    !                              : cross section data of CO2 will be stored to this variable after binning.
    !                         b: number of cross section data column, e.g.) 12 (if there are cross sections for 12 temperatures)
    !                         c: unit in 'cm2' or in 'm2', e.g.) 'cm2' (if the cross sections are in units of 'cm^2')
    !                            if data is a solar flux or a quantum yeild, it should be ''.
    !                            if wavelength is in 'cm^-1' and cross sections are in 'cm^2', it should be 'cm-1, cm2'
    !                            if wavelength is in 'cm^-1' and cross sections are in 'm^2', it should be 'cm-1, m2'
    !                            if wavelength is in 'cm^-1' and data is no unit, it should be 'cm-1'

    ! CO2 data --------------------------------------------------------
    isp = sp_index(spl,'CO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_12C16O16O.dat'
    allocate(xct%T_x_CO2(12))
    xct%T_x_CO2(:) = (/120.0_dp, 145.0_dp, 170.0_dp, 195.0_dp, 220.0_dp, 245.0_dp, &
      &                270.0_dp, 295.0_dp, 320.0_dp, 345.0_dp, 370.0_dp, 395.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CO2, 12, 'cm2')

    fname = './UV/xct_data/sigma_d_CO2_CO+O1D.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_CO2_O1D, 1, 'cm2')
    
    fname = './UV/xct_data/qy_CO2_CO+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_CO2_O, 1, '')

    ! ^13CO2 data -----------------------------------------------------
    isp = sp_index(spl,'^13CO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_13C16O16O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_u13CO2, 12, 'cm2')

    ! H2O2 data -------------------------------------------------------
    isp = sp_index(spl,'H2O2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_H2O2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_H2O2, 1, 'cm2')

    fname = './UV/xct_data/qy_H2O2_OH+OH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_H2O2_OH, 1, '')

    fname = './UV/xct_data/qy_H2O2_HO2+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_H2O2_HO2, 1, '')

    ! H2O data --------------------------------------------------------
    isp = sp_index(spl,'H2O')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_H2O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_H2O, 1, 'cm2')

    fname = './UV/xct_data/qy_H2O_H+OH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_H2O_OH, 1, '')

    fname = './UV/xct_data/qy_H2O_H2+O1D.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_H2O_O1D, 1, '')

    ! H2 data ---------------------------------------------------------
    isp = sp_index(spl,'H2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_H2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_H2_H+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_H2_H, 1, 'cm2')

    ! OH data ---------------------------------------------------------
    isp = sp_index(spl,'OH')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_OH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_OH, 1, 'cm2')

    fname = './UV/xct_data/qy_OH_O+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_OH_O, 1, '')

    fname = './UV/xct_data/sigma_d_OH_O+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_OH_O, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_OH_O1D+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_OH_O1D, 1, 'cm2')

    ! HO2 data --------------------------------------------------------
    isp = sp_index(spl,'HO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HO2, 1, 'cm2')

    fname = './UV/xct_data/qy_HO2_OH+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HO2_OH, 1, '')

    ! O3 data ---------------------------------------------------------
    isp = sp_index(spl,'O3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_O3.dat'
    allocate(xct%T_x_O3(11))
    xct%T_x_O3(:) = (/193.0_dp, 203.0_dp, 213.0_dp, 223.0_dp, 233.0_dp, 243.0_dp, &
      &               253.0_dp, 263.0_dp, 273.0_dp, 283.0_dp, 293.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_O3, 11, 'cm2')

    ! O2 data ---------------------------------------------------------
    isp = sp_index(spl,'O2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_O2.dat'
    allocate(xct%T_x_O2(2))
    xct%T_x_O2(:) = (/90.0_dp, 295.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_O2, 2, 'cm2')

    fname = './UV/xct_data/sigma_a_O2_Schumann-Runge_bands_130-190K.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_O2_SR_bands_130K, 5, 'cm-1')

    fname = './UV/xct_data/sigma_a_O2_Schumann-Runge_bands_190-280K.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_O2_SR_bands_190K, 5, 'cm-1')

    fname = './UV/xct_data/sigma_a_O2_Schumann-Runge_bands_280-500K.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_O2_SR_bands_280K, 5, 'cm-1')

    fname = './UV/xct_data/qy_O2_O+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_O2_O, 1, '')

    fname = './UV/xct_data/qy_O2_O+O1D.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_O2_O1D, 1, '')

    ! H2CO data --------------------------------------------------------
    isp = sp_index(spl,'H2CO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_H2CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_H2CO, 1, 'cm2')

    fname = './UV/xct_data/sigma_gamma_H2CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_gamma_H2CO, 1, 'cm2')

    fname = './UV/xct_data/qy_H2CO_H2+COat300K1atm.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_H2CO_CO_300K, 1, '')

    ! N2 data ---------------------------------------------------------
    isp = sp_index(spl,'N2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_N2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_N2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_N2_N+N.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_N2_N, 1, 'cm2')

    ! NO data ---------------------------------------------------------
    isp = sp_index(spl,'NO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NO, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_NO_N+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_NO_N, 1, 'cm2')

    ! NO2 data ---------------------------------------------------------
    isp = sp_index(spl,'NO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NO2.dat'
    allocate(xct%T_x_NO2(2))
    xct%T_x_NO2(:) = (/220.0_dp, 294.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NO2, 2, 'cm2')

    fname = './UV/xct_data/qy_NO2_NO+O.dat'
    allocate(xct%T_qy_NO2(2))
    xct%T_qy_NO2(:) = (/248.0_dp, 298.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_NO2_O, 2, '')

    fname = './UV/xct_data/sigma_d_NO2_NO+O(1D).dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_d_NO2_O1D, 1, 'cm2')

    ! NO3 data ---------------------------------------------------------
    isp = sp_index(spl,'NO3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NO3.dat'
    allocate(xct%T_x_NO3(2))
    xct%T_x_NO3(:) = (/220.0_dp, 298.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NO3, 2, 'cm2')

    fname = './UV/xct_data/qy_NO3_NO+O2.dat'
    allocate(xct%T_qy_NO3(3))
    xct%T_qy_NO3(:) = (/190.0_dp, 230.0_dp, 298.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_NO3_NO, 3, '')

    fname = './UV/xct_data/qy_NO3_NO2+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_NO3_NO2, 3, '')

    ! N2O data ---------------------------------------------------------
    isp = sp_index(spl,'N2O')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_N2O.dat'
    allocate(xct%T_x_N2O(5))
    xct%T_x_N2O(:) = (/194.0_dp, 225.0_dp, 243.0_dp, 263.0_dp, 302.0_dp/)
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_N2O, 5, 'cm2')

    fname = './UV/xct_data/qy_N2O_N2+O(1D).dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_N2O_O1D, 1, '')

    ! N2O5 data ---------------------------------------------------------
    isp = sp_index(spl,'N2O5')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_N2O5.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_N2O5, 1, 'cm2')

    fname = './UV/xct_data/sigma_a_N2O5_Tdep.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_N2O5_above260nm, 2, '')

    fname = './UV/xct_data/qy_N2O5_NO3+NO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_N2O5_NO2, 1, '')

    fname = './UV/xct_data/qy_N2O5_NO3+NO+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_N2O5_O, 1, '')

    ! HONO data ---------------------------------------------------------
    isp = sp_index(spl,'HONO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HONO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HONO, 1, 'cm2')

    fname = './UV/xct_data/qy_HONO_OH+NO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HONO_OH, 1, '')

    ! HNO3 data ---------------------------------------------------------
    isp = sp_index(spl,'HNO3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HNO3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HNO3, 1, 'cm2')

    fname = './UV/xct_data/sigma_a_HNO3_Tdep.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HNO3_Tdep, 1, '')

    fname = './UV/xct_data/qy_HNO3_OH+NO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HNO3_OH, 1, '')

    fname = './UV/xct_data/qy_HNO3_HONO+O.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HNO3_O, 1, '')

    fname = './UV/xct_data/qy_HNO3_HONO+O(1D).dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HNO3_O1D, 1, '')

    ! HO2NO2 data ---------------------------------------------------------
    isp = sp_index(spl,'HO2NO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HO2NO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HO2NO2, 1, 'cm2')

    fname = './UV/xct_data/sigma_a_HO2NO2_Tdep.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HO2NO2_Tdep, 2, '')

    fname = './UV/xct_data/qy_HO2NO2_OH+NO3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HO2NO2_OH, 1, '')

    fname = './UV/xct_data/qy_HO2NO2_HO2+NO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HO2NO2_HO2, 1, '')

    ! HCO data --------------------------------------------------------
    isp = sp_index(spl,'HCO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HCO, 1, 'cm2')

    fname = './UV/xct_data/qy_HCO_H+CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%qy_HCO_H, 1, '')

    ! CH4 data --------------------------------------------------------
    isp = sp_index(spl,'CH4')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH4.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH4, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_CH4_1CH2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_CH4_1CH2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_CH4_3CH2+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_CH4_3CH2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_CH4_CH3+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_CH4_CH3, 1, 'cm2')

    ! C2H6 data --------------------------------------------------------
    isp = sp_index(spl,'C2H6')
    !------------------------------------------------------------------
    fname = './UV/xct_data/sigma_a_C2H6.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C2H6, 1, 'cm2')

    fname = './UV/xct_data/qy_C2H6_3CH2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_3CH2, 1, '')

    fname = './UV/xct_data/qy_C2H6_CH4+1CH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_1CH2, 1, '')

    fname = './UV/xct_data/qy_C2H6_C2H2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_C2H2, 1, '')

    fname = './UV/xct_data/qy_C2H6_C2H4+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_C2H4_H, 1, '')

    fname = './UV/xct_data/qy_C2H6_C2H4+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_C2H4_H2, 1, '')

    fname = './UV/xct_data/qy_C2H6_CH3+CH3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H6_CH3, 1, '')

    ! HNO data --------------------------------------------------------
    isp = sp_index(spl,'HNO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HNO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HNO, 1, 'cm2')

    fname = './UV/xct_data/qy_HNO_NO+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_HNO_NO, 1, '')

    ! HNO2 data --------------------------------------------------------
    isp = sp_index(spl,'HNO2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HNO2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HNO2, 1, 'cm2')

    fname = './UV/xct_data/qy_HNO2_NO+OH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_HNO2_NO, 1, '')

    ! CH3 data --------------------------------------------------------
    isp = sp_index(spl,'CH3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH3, 1, 'cm2')

    fname = './UV/xct_data/qy_CH3_1CH2+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3_1CH2, 1, '')

    ! NH3 data --------------------------------------------------------
    isp = sp_index(spl,'NH3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NH3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NH3, 1, 'cm2')

    fname = './UV/xct_data/qy_NH3_NH2+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_NH3_NH2, 1, '')

    ! N2H4 data --------------------------------------------------------
    isp = sp_index(spl,'N2H4')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_N2H4.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_N2H4, 1, 'cm2')

    fname = './UV/xct_data/qy_N2H4_N2H3+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_N2H4_N2H3, 1, '')

    ! NH data --------------------------------------------------------
    isp = sp_index(spl,'NH')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NH, 1, 'cm2')

    fname = './UV/xct_data/qy_NH_N+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_NH_N, 1, '')

    ! NH2 data --------------------------------------------------------
    isp = sp_index(spl,'NH2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NH2, 1, 'cm2')

    fname = './UV/xct_data/qy_NH2_NH+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_NH2_NH, 1, '')

    ! C2H2 data --------------------------------------------------------
    isp = sp_index(spl,'C2H2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C2H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C2H2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_C2H2_C2H+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_C2H2_C2H, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_C2H2_C2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_C2H2_C2, 1, 'cm2')

    ! C2H4 data --------------------------------------------------------
    isp = sp_index(spl,'C2H4')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C2H4.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C2H4, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_C2H4_C2H2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_C2H4_C2H2_H2, 1, 'cm2')

    fname = './UV/xct_data/sigma_d_C2H4_C2H2+H+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%sigma_d_C2H4_C2H2_H, 1, 'cm2')

    ! C3H8 data --------------------------------------------------------
    isp = sp_index(spl,'C3H8')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C3H8.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C3H8, 1, 'cm2')

    fname = './UV/xct_data/qy_C3H8_C3H6+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H8_C3H6, 1, '')

    fname = './UV/xct_data/qy_C3H8_C2H6+1CH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H8_C2H6, 1, '')

    fname = './UV/xct_data/qy_C3H8_C2H4+CH4.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H8_C2H4, 1, '')

    fname = './UV/xct_data/qy_C3H8_C2H5+CH3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H8_C2H5, 1, '')

    ! C3H6 data --------------------------------------------------------
    isp = sp_index(spl,'C3H6')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C3H6.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C3H6, 1, 'cm2')

    fname = './UV/xct_data/qy_C3H6_C2H2+CH3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H6_C2H2, 1, '')

    fname = './UV/xct_data/qy_C3H6_CH2CCH2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H6_CH2CCH2, 1, '')

    fname = './UV/xct_data/qy_C3H6_C2H4+3CH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H6_C2H4, 1, '')

    fname = './UV/xct_data/qy_C3H6_C2H+CH4.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H6_C2H, 1, '')

    ! CH data --------------------------------------------------------
    isp = sp_index(spl,'CH')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH, 1, 'cm2')

    fname = './UV/xct_data/qy_CH_C+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH_C, 1, '')

    ! CH2CO data --------------------------------------------------------
    isp = sp_index(spl,'CH2CO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH2CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH2CO, 1, 'cm2')

    fname = './UV/xct_data/qy_CH2CO_CH2+CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH2CO_CH2, 1, '')
    
    ! CH3CHO data --------------------------------------------------------
    isp = sp_index(spl,'CH3CHO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH3CHO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH3CHO, 1, 'cm2')

    fname = './UV/xct_data/qy_CH3CHO_CH3+HCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3CHO_CH3, 1, '')

    fname = './UV/xct_data/qy_CH3CHO_CH4+CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3CHO_CH4, 1, '')

    ! C2H5CHO data --------------------------------------------------------
    isp = sp_index(spl,'C2H5CHO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C2H5CHO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C2H5CHO, 1, 'cm2')

    fname = './UV/xct_data/qy_C2H5CHO_C2H5+HCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C2H5CHO_C2H5, 1, '')

    ! C3H3 data --------------------------------------------------------
    isp = sp_index(spl,'C3H3')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_C3H3.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_C3H3, 1, 'cm2')

    fname = './UV/xct_data/qy_C3H3_C3H2+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_C3H3_C3H2, 1, '')

    ! CH3C2H data --------------------------------------------------------
    isp = sp_index(spl,'CH3C2H')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH3C2H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH3C2H, 1, 'cm2')

    fname = './UV/xct_data/qy_CH3C2H_C3H3+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3C2H_C3H3, 1, '')

    fname = './UV/xct_data/qy_CH3C2H_C3H2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3C2H_C3H2, 1, '')

    fname = './UV/xct_data/qy_CH3C2H_CH3+C2H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH3C2H_CH3, 1, '')

    ! CH2CCH2 data --------------------------------------------------------
    isp = sp_index(spl,'CH2CCH2')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_CH2CCH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_CH2CCH2, 1, 'cm2')

    fname = './UV/xct_data/qy_CH2CCH2_C3H3+H.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH2CCH2_C3H3, 1, '')

    fname = './UV/xct_data/qy_CH2CCH2_C3H2+H2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH2CCH2_C3H2, 1, '')

    fname = './UV/xct_data/qy_CH2CCH2_C2H2+CH2.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_CH2CCH2_C2H2, 1, '')

    ! HCN data --------------------------------------------------------
    isp = sp_index(spl,'HCN')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HCN.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HCN, 1, 'cm2')

    fname = './UV/xct_data/qy_HCN_H+CN.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_HCN_H, 1, '')

    ! HNCO data --------------------------------------------------------
    isp = sp_index(spl,'HNCO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_HNCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_HNCO, 1, 'cm2')

    fname = './UV/xct_data/qy_HNCO_H+NCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_HNCO_H, 1, '')

    ! NCO data --------------------------------------------------------
    isp = sp_index(spl,'NCO')
    !------------------------------------------------------------------

    fname = './UV/xct_data/sigma_a_NCO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
      &                  xct%sigma_a_NCO, 1, 'cm2')

    fname = './UV/xct_data/qy_NCO_N+CO.dat'
    call p__UV_read_data(xct, flx, fname, isp, &
         &                  xct%qy_NCO_N, 1, '')              

    ! finish reading cross sections. -----------------------------------
    print *, 'Finished reading cross section data.'

  end subroutine p__UV_cross_section_dat


  !-------------------------------------------------
  ! cross section data for UV
  !-------------------------------------------------
  subroutine p__UV_cross_section_exe(T, nz, spl, flx, var, & ! in
    &                                xct                   ) ! inout
    implicit none
    integer,    intent(in)    :: nz
    real(dp),   intent(in)    :: T(nz)
    type(spl_), intent(in)    :: spl
    type(flx_), intent(in)    :: flx
    type(var_), intent(in)    :: var
    type(xct_), intent(inout) :: xct
    integer iwl, i, iz, isp, ich, swl, ewl, R_label
    character(len=256) reactants(10),  products(10), fname
    real(dp)  alpha_R         !! polarilizability for Rayleigh scattering
    real(dp)  pi2             !! circle ratio
    real(dp) dummy(1)

    pi2 = 3.1415926535_dp

    !----------------------------------------------------------------------------------------------------
    ! UV photoabsorption cross section

    if (xct%type == 'absorption') then

       xct%label_sigma_a_UV = 0
       xct%label_sigma_a_RS = 0
       xct%wlrange_a_UV(1,:) = flx%nwl_UV
       xct%wlrange_a_UV(2,:) = 1
       
       do isp = 1, spl%nsp
        ! Template --------------------------------------------------------
        !  if ( spl%species(isp) == 'CO2') then
        !    call p__UV_sigma_a_adapt(xct, T, isp, & ! fixed, do not modify these three.
        !      &                  a, b) ! You should modify only a and b.
        !           a: xct%variable_name e.g.) xct%sigma_a_CO2
        !           b: temperature array corresponding to the column of read data 
        !              e.g.) xct%T_x_CO2
        !              if there is only 1 column of cross section data, it should be 'dummy', which has only 1 element.
        !  end if

        ! CO2 ------------------------------------------------------------
          if (     spl%species(isp) == 'CO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CO2, xct%T_x_CO2)
             
          else if (spl%species(isp) == '^13CO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_u13CO2, xct%T_x_CO2)

        ! H2O2 ------------------------------------------------------------
          else if (spl%species(isp) == 'H2O2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_H2O2, dummy)
             call p__UV_sigma_a_H2O2_analytic(xct, flx, T, isp)
             
        ! H2O ------------------------------------------------------------
          else if (spl%species(isp) == 'H2O') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_H2O, dummy)
             
        ! HO2 ------------------------------------------------------------
          else if (spl%species(isp) == 'HO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HO2, dummy)

        ! H2 ------------------------------------------------------------
          else if (spl%species(isp) == 'H2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_H2, dummy)

        ! OH ------------------------------------------------------------
          else if (spl%species(isp) == 'OH') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_OH, dummy)
             
        ! O3 ------------------------------------------------------------
          else if (spl%species(isp) == 'O3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_O3, xct%T_x_O3)

        ! O2 ------------------------------------------------------------
          else if (spl%species(isp) == 'O2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_O2, xct%T_x_O2)
             call p__UV_sigma_a_O2_Schumann_Runge_bands(xct, T, isp)
             call p__UV_sigma_a_O2_Herzberg_continuum(xct, flx, T, isp)

        ! H2CO ------------------------------------------------------------
          else if (spl%species(isp) == 'H2CO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_H2CO, dummy)
             call p__UV_sigma_a_H2CO_Tdependent(xct, T, isp)
             
        ! N2 ------------------------------------------------------------
          else if (spl%species(isp) == 'N2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_N2, dummy)
             
        ! NO ------------------------------------------------------------
          else if (spl%species(isp) == 'NO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NO, dummy)

        ! NO2 ------------------------------------------------------------
          else if (spl%species(isp) == 'NO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NO2, xct%T_x_NO2)

        ! NO3 ------------------------------------------------------------
          else if (spl%species(isp) == 'NO3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NO3, xct%T_x_NO3)

        ! N2O ------------------------------------------------------------
          else if (spl%species(isp) == 'N2O') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_N2O, xct%T_x_N2O)

        ! N2O5 ------------------------------------------------------------
          else if (spl%species(isp) == 'N2O5') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_N2O5, dummy)
             call p__UV_sigma_a_N2O5_Tdependent(xct, T, isp)

        ! HONO ------------------------------------------------------------
          else if (spl%species(isp) == 'HONO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HONO, dummy)

        ! HNO3 ------------------------------------------------------------
          else if (spl%species(isp) == 'HNO3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HNO3, dummy)
             call p__UV_sigma_a_HNO3_Tdependent(xct, T, isp)

        ! HO2NO2 ------------------------------------------------------------
          else if (spl%species(isp) == 'HO2NO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HO2NO2, dummy)
             call p__UV_sigma_a_HO2NO2_Tdependent(xct, T, isp)
          
        ! HCO ------------------------------------------------------------ 
          else if (spl%species(isp) == 'HCO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HCO, dummy)

        ! CH4 ------------------------------------------------------------ 
          else if (spl%species(isp) == 'CH4') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH4, dummy)

        ! C2H6 ------------------------------------------------------------
          else if (spl%species(isp) == 'C2H6') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C2H6, dummy)

        ! HNO2 ------------------------------------------------------------
          else if (spl%species(isp) == 'HNO2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HNO2, dummy)

        ! CH3 ------------------------------------------------------------
          else if (spl%species(isp) == 'CH3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH3, dummy)

        ! HNO ------------------------------------------------------------
          else if (spl%species(isp) == 'HNO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HNO, dummy)

        ! NH3 ------------------------------------------------------------
          else if (spl%species(isp) == 'NH3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NH3, dummy)

        ! N2H4 ------------------------------------------------------------
          else if (spl%species(isp) == 'N2H4') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_N2H4, dummy)

        ! NH ------------------------------------------------------------
          else if (spl%species(isp) == 'NH') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NH, dummy)

        ! NH2 ------------------------------------------------------------
          else if (spl%species(isp) == 'NH2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NH2, dummy)

        ! C2H2 ------------------------------------------------------------
          else if (spl%species(isp) == 'C2H2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C2H2, dummy)

        ! C2H4 ------------------------------------------------------------ 
          else if (spl%species(isp) == 'C2H4') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C2H4, dummy)

        ! C3H8 ------------------------------------------------------------ 
          else if (spl%species(isp) == 'C3H8') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C3H8, dummy)

        ! C3H6 ------------------------------------------------------------
          else if (spl%species(isp) == 'C3H6') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C3H6, dummy)

        ! CH ------------------------------------------------------------ 
          else if (spl%species(isp) == 'CH') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH, dummy)

        ! CH2CO ------------------------------------------------------------
          else if (spl%species(isp) == 'CH2CO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH2CO, dummy)          

        ! CH3CHO ------------------------------------------------------------ 
          else if (spl%species(isp) == 'CH3CHO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH3CHO, dummy)

        ! C2H5CHO ------------------------------------------------------------ 
          else if (spl%species(isp) == 'C2H5CHO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C2H5CHO, dummy)

        ! C3H3 ------------------------------------------------------------
          else if (spl%species(isp) == 'C3H3') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_C3H3, dummy)

        ! CH3C2H ------------------------------------------------------------ 
          else if (spl%species(isp) == 'CH3C2H') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH3C2H, dummy)

        ! CH2CCH2 ------------------------------------------------------------ 
          else if (spl%species(isp) == 'CH2CCH2') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_CH2CCH2, dummy)

        ! HCN ------------------------------------------------------------
          else if (spl%species(isp) == 'HCN') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HCN, dummy)

        ! HNCO ------------------------------------------------------------
          else if (spl%species(isp) == 'HNCO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_HNCO, dummy)

        ! NCO ------------------------------------------------------------
          else if (spl%species(isp) == 'NCO') then
             call p__UV_sigma_a_adapt(xct, T, isp, &
                  &                      xct%sigma_a_NCO, dummy)          
          
          end if ! end if species

       end do !isp

      ! Rayleigh scattering for major species (Liou 2002) -----------------
        ! Polarizabilities (alpha) are taken from 
        ! Computational Chemistry Comparison and Benchmark DataBase
        ! alpha  -> nm^3
        ! lambda -> nm
      do isp = 1, spl%nsp

!        swl = xct%wlrange_a_UV(1,isp)
!        ewl = xct%wlrange_a_UV(2,isp)

        do iwl = 1, flx%nwl_UV

          if (spl%species(isp) == 'CO2') then
            xct%label_sigma_a_RS(isp) = 1
            alpha_R = 2.507e-3_dp
            xct%sigma_a_RS(iwl,isp) = 128.0_dp * pi2**5   &
            &                       * alpha_R * alpha_R      &
            &                       / 3.0_dp / flx%lambda_UV(iwl)**4  &
            &                       * 1.0e-18_dp

          else if (spl%species(isp) == 'H2') then
            xct%label_sigma_a_RS(isp) = 1
            alpha_R = 0.787e-3_dp
            xct%sigma_a_RS(iwl,isp) = 128.0_dp * pi2**5   &
            &                       * alpha_R * alpha_R      &
            &                       / 3.0_dp / flx%lambda_UV(iwl)**4  &
            &                       * 1.0e-18_dp

          else if (spl%species(isp) == 'CH4') then
            xct%label_sigma_a_RS(isp) = 1
            alpha_R = 2.448e-3_dp
            xct%sigma_a_RS(iwl,isp) = 128.0_dp * pi2**5   &
            &                       * alpha_R * alpha_R      &
            &                       / 3.0_dp / flx%lambda_UV(iwl)**4  &
            &                       * 1.0e-18_dp

          else if (spl%species(isp) == 'H2O') then
            xct%label_sigma_a_RS(isp) = 1
            alpha_R = 1.501e-3_dp
            xct%sigma_a_RS(iwl,isp) = 128.0_dp * pi2**5   &
            &                       * alpha_R * alpha_R      &
            &                       / 3.0_dp / flx%lambda_UV(iwl)**4  &
            &                       * 1.0e-18_dp

          else if (spl%species(isp) == 'N2') then
            xct%label_sigma_a_RS(isp) = 1
            alpha_R = 1.710e-3_dp
            xct%sigma_a_RS(iwl,isp) = 128.0_dp * pi2**5   &
            &                       * alpha_R * alpha_R      &
            &                       / 3.0_dp / flx%lambda_UV(iwl)**4  &
            &                       * 1.0e-18_dp

          else
            xct%sigma_a_RS(iwl,isp) = 0.0_dp

          end if
        
        end do

      end do  


      !! if cross sections are less than 0, they are set to 0 ------
      !where (xct%sigma_a_UV(:,:,:)<0.0_dp)
      !  xct%sigma_a_UV(:,:,:) = 0.0_dp
      !end where

      ! output cross sections
      if (var%istep == 0) then 
        iz = 1 ! select altitude grid
        do isp = 1, spl%nsp
          if (xct%label_sigma_a_UV(isp) == 1) then 
            fname = './UV/xct/absorption/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
            open(11, file = fname, status = 'unknown' )
            swl = xct%wlrange_a_UV(1,isp)
            ewl = xct%wlrange_a_UV(2,isp)
            do iwl = swl, ewl
              write(11, *) flx%lambda_UV(iwl), xct%sigma_a_UV(iwl,iz,isp)
            end do
            close(11)
          end if
        end do 
      end if

      if (var%istep == 0) then 
        do isp = 1, spl%nsp
          if (xct%label_sigma_a_UV(isp) == 1) then 
            fname = './UV/xct/RS/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
            open(71, file = fname, status = 'unknown' )
            do iwl = 1,flx%nwl_UV
              write(71, *) flx%lambda_UV(iwl), xct%sigma_a_RS(iwl,isp)
            end do
            close(71)
          end if
        end do 
      end if

    end if ! end if label == absorption

    !----------------------------------------------------------------------------------------------------
    ! UV photodissociation cross section

    if (xct%type == 'dissociation') then

      xct%wlrange_d_UV(1,:) = flx%nwl_UV
      xct%wlrange_d_UV(2,:) = 1

      do ich = 1, spl%nch

        reactants = ''
        products = ''

        if (spl%reaction_type_char(ich) == 'photodissociation') then
          do i = 2, spl%reactant_list(ich,1)+1
            reactants(i-1) = trim(spl%species(spl%reactant_list(ich,i)))
          end do
          do i = 2, spl%product_list(ich,1)+1
            products(i-1) = trim(spl%species(spl%product_list(ich,i)))
          end do
          isp = spl%reactant_list(ich,2)

          ! Template --------------------------------------------------------
          !  
          ! R_label = ch_identify(spl, ich, 'photodissociation', ['reactant'], ['product1','product2',,,])
          !        'reactant' : e.g.) 'CO2'
          !        'product*' : e.g.) 'CO', 'O ' 
          !          ** It shoule be noted that the number of characters for each product must be the same!!!
          !             Please add some spaces in order to make the number of characters for each product the same.
          !          e.g.) 'H2O2', 'O   ' will work fine. Each product has the same number of characters (4 characters).
          !             If the number of characters of products differs each other, compile error occurs if you use gfortran. 
          !             (It will work fine without any errors if you use Intel fortran (ifort) even if the number of characters 
          !              is different for each product, but it is not a normal usage of fortran language.)
          !
          ! ! if you apply a quantum yield for photodissociation reaction, 
          ! call p__UV_qy_adapt(xct, T, isp, ich, R_label, & ! fixed, do not modify these five.
          !                     a, b) ! you should modify a and b.
          !         a: xct%variable_name e.g.) xct%qy_CO2_O
          !         b: temperature array corresponding to the column of read data 
          !              e.g.) xct%T_qy_CO2
          !              if there is only 1 column of cross section data, it should be 'dummy', which has only 1 element.
          !
          ! ! if you apply a photodissociation cross section, 
          ! call p__UV_sigma_d_adapt(xct, T, ich, R_label, & ! fixed, do not modify these five.
          !                          a, b) ! you should modify a and b.
          !         a: xct%variable_name e.g.) xct%sigma_d_CO2_O1D
          !         b: temperature array corresponding to the column of read data 
          !              e.g.) xct%T_x_CO2_O1D
          !              if there is only 1 column of cross section data, it should be 'dummy', which has only 1 element.
          !    

          ! CO2 -------------------------------------------
            ! CO2 + hv -> CO + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['CO2'], ['CO','O '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_CO2_O, dummy)

            ! CO2 + hv -> CO + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['CO2'], ['CO   ','O(1D)'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_CO2_O1D, dummy)

            ! ^13CO2 + hv -> ^13CO + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['^13CO2'], ['^13CO','O    '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_CO2_O, dummy)

            ! ^13CO2 + hv -> ^13CO + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['^13CO2'], ['^13CO','O(1D)'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_CO2_O1D, dummy)
            
          ! H2O -------------------------------------------
            ! H2O + hv -> H + OH
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2O'], ['H ','OH'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_H2O_OH, dummy)

            ! H2O + hv -> H2 + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2O'], ['H2   ','O(1D)'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_H2O_O1D, dummy)

            ! H2O + hv -> H + H + O
              ! no reactions?

          ! H2O2 -------------------------------------------
            ! H2O2 + hv -> OH + OH
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2O2'], ['OH','OH'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_H2O2_OH, dummy)

            ! H2O2 + hv -> HO2 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2O2'], ['HO2','H  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_H2O2_HO2, dummy)

            ! H2O2 + hv -> H2O + O(1D)
              ! no reactions?

          ! O3 -------------------------------------------
            ! O3 + hv -> O2 + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['O3'], ['O2','O '])
            call p__UV_qy_O3(xct, flx, T, isp, ich, R_label, 'O')

            ! O3 + hv -> O2 + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['O3'], ['O2   ','O(1D)'])
            call p__UV_qy_O3(xct, flx, T, isp, ich, R_label, 'O(1D)')

          ! O2 -------------------------------------------
            ! O2 + hv -> O + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['O2'], ['O','O'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_O2_O, dummy)

            ! O2 + hv -> O + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['O2'], ['O    ','O(1D)'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_O2_O1D, dummy)

          ! H2 -------------------------------------------
            ! H2 + hv -> H + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2'], ['H','H'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_H2_H, dummy)

          ! OH -------------------------------------------
            ! OH + hv -> O + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['OH'], ['O','H'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_OH_O, dummy)

            ! OH + hv -> O(1D) + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['OH'], ['O(1D)','H    '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_OH_O1D, dummy)

          ! HO2 -------------------------------------------
            ! HO2 + hv -> OH + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['HO2'], ['OH','O '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HO2_OH, dummy)
          
          ! H2CO -------------------------------------------
            ! H2CO + hv -> H + HCO
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2CO'], ['H  ','HCO'])
            call p__UV_qy_H2CO(var, xct, flx, T, isp, ich, R_label, 'H')

            ! H2CO + hv -> H2 + CO
            R_label = ch_identify(spl, ich, 'photodissociation', ['H2CO'], ['H2','CO'])
            call p__UV_qy_H2CO(var, xct, flx, T, isp, ich, R_label, 'CO')

          ! N2 -------------------------------------------
            ! N2 + hv -> N + N(2D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['N2'], ['N    ','N(2D)'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_N2_N, dummy)

          ! NO -------------------------------------------
            ! NO + hv -> N + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['NO'], ['N','O'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_NO_N, dummy)

          ! NO2 -------------------------------------------
            ! NO2 + hv -> NO + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['NO2'], ['NO','O '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_NO2_O, xct%T_qy_NO2)

            ! NO2 + hv -> NO + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['NO2'], ['NO   ','O(1D)'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_NO2_O1D, dummy)

          ! NO3 -------------------------------------------
            ! NO3 + hv -> NO + O2
            R_label = ch_identify(spl, ich, 'photodissociation', ['NO3'], ['NO','O2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_NO3_NO, xct%T_qy_NO3)

            ! NO3 + hv -> NO2 + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['NO3'], ['NO2','O  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_NO3_NO2, xct%T_qy_NO3)

          ! N2O -------------------------------------------
            ! N2O + hv -> N2 + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['N2O'], ['N2   ','O(1D)'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_N2O_O1D, dummy)

          ! N2O5 -------------------------------------------
            ! N2O5 + hv -> NO3 + NO2
            R_label = ch_identify(spl, ich, 'photodissociation', ['N2O5'], ['NO3','NO2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_N2O5_NO2, dummy)

            ! N2O5 + hv -> NO3 + NO + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['N2O5'], ['NO3','NO ','O  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_N2O5_O, dummy)

          ! HONO -------------------------------------------
            ! HONO + hv -> OH + NO
            R_label = ch_identify(spl, ich, 'photodissociation', ['HONO'], ['OH','NO'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HONO_OH, dummy)

          ! HNO3 -------------------------------------------
            ! HNO3 + hv -> OH + NO2
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNO3'], ['OH ','NO2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HNO3_OH, dummy)

            ! HNO3 + hv -> HONO + O
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNO3'], ['HONO','O   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HNO3_O, dummy)

            ! HNO3 + hv -> HONO + O(1D)
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNO3'], ['HONO ','O(1D)'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HNO3_O1D, dummy)

          ! HO2NO2 -------------------------------------------
            ! HO2NO2 + hv -> OH + NO3
            R_label = ch_identify(spl, ich, 'photodissociation', ['HO2NO2'], ['OH ','NO3'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HO2NO2_OH, dummy)

            ! HO2NO2 + hv -> HO2 + NO2
            R_label = ch_identify(spl, ich, 'photodissociation', ['HO2NO2'], ['HO2','NO2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
              &                 xct%qy_HO2NO2_HO2, dummy)
            
          ! HCO -------------------------------------------
            ! HCO + hv -> H + CO
            R_label = ch_identify(spl, ich, 'photodissociation', ['HCO'], ['H ','CO'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_HCO_H, dummy)

          ! CH4 -------------------------------------------
            ! CH4 + hv -> ^1CH2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH4'], ['^1CH2','H2   '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
                 &                      xct%sigma_d_CH4_1CH2, dummy)

            ! CH4 + hv -> ^3CH2 + H + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH4'], ['^3CH2','H    ','H    '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_CH4_3CH2, dummy)

            ! CH4 + hv -> CH3 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH4'], ['CH3','H  '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
              &                      xct%sigma_d_CH4_CH3, dummy)
            
          ! C2H6 -------------------------------------------
            ! C2H6 + hv -> ^3CH2 + ^3CH2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['^3CH2','^3CH2', 'H2   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_3CH2, dummy)

            ! C2H6 + hv -> CH4 + ^1CH2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['CH4  ','^1CH2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_1CH2, dummy)
            
            ! C2H6 + hv -> C2H2 + H2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['C2H2','H2  ', 'H2  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_C2H2, dummy)

            ! C2H6 + hv -> C2H4 + H + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['C2H4','H   ', 'H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_C2H4_H, dummy)

            ! C2H6 + hv -> C2H4 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['C2H4','H2  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_C2H4_H2, dummy)

            ! C2H6 + hv -> CH3 + CH3
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H6'], ['CH3','CH3'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H6_CH3, dummy)

          ! HNO2 -------------------------------------------
            ! HNO2 + hv -> NO + OH
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNO2'], ['NO','OH'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_HNO2_NO, dummy)

          ! CH3 -------------------------------------------
            ! CH3 + hv -> ^1CH2 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3'], ['^1CH2','H    '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3_1CH2, dummy)

          ! HNO -------------------------------------------
            ! HNO + hv -> NO + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNO'], ['NO','H '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_HNO_NO, dummy)

          ! NH3 -------------------------------------------
            ! NH3 + hv -> NH2 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['NH3'], ['NH2','H  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_NH3_NH2, dummy)

          ! N2H4 -------------------------------------------
            ! N2H4 + hv -> N2H3 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['N2H4'], ['N2H3','H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_N2H4_N2H3, dummy)

          ! NH -------------------------------------------
            ! NH + hv -> N + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['NH'], ['N','H'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_NH_N, dummy)

          ! NH2 -------------------------------------------
            ! NH2 + hv -> NH + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['NH2'], ['NH','H '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_NH2_NH, dummy)

          ! C2H2 -------------------------------------------
            ! C2H2 + hv -> C2H + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H2'], ['C2H','H  '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
                 &                      xct%sigma_d_C2H2_C2H, dummy)

            ! C2H2 + hv -> C2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H2'], ['C2','H2'])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
                 &                      xct%sigma_d_C2H2_C2, dummy)

          ! C2H4 -------------------------------------------
            ! C2H4 + hv -> C2H2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H4'], ['C2H2','H2  '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
                 &                      xct%sigma_d_C2H4_C2H2_H2, dummy)

            ! C2H4 + hv -> C2H2 + H + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H4'], ['C2H2','H   ','H   '])
            call p__UV_sigma_d_adapt(xct, T, ich, R_label, &
                 &                      xct%sigma_d_C2H4_C2H2_H, dummy)

          ! C3H8 -------------------------------------------
            ! C3H8 + hv -> C3H6 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H8'], ['C3H6','H2  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H8_C3H6, dummy)

            ! C3H8 + hv -> C2H6 + ^1CH2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H8'], ['C2H6 ','^1CH2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H8_C2H6, dummy)

            ! C3H8 + hv -> C2H4 + CH4
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H8'], ['C2H4','CH4 '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H8_C2H4, dummy)

            ! C3H8 + hv -> C2H5 + CH3
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H8'], ['C2H5','CH3 '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H8_C2H5, dummy)
            
          ! C3H6 -------------------------------------------
            ! C3H6 + hv -> C2H2 + CH3 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H6'], ['C2H2','CH3 ','H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H6_C2H2, dummy)

            ! C3H6 + hv -> CH2CCH2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H6'], ['C2H2   ','CH2CCH2','H2     '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H6_CH2CCH2, dummy)

            ! C3H6 + hv -> C2H4 + ^3CH2
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H6'], ['C2H4 ','^3CH2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H6_C2H4, dummy)

            ! C3H6 + hv -> C2H + CH4 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H6'], ['C2H','CH4','H  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H6_C2H, dummy)

          ! CH -------------------------------------------
            ! CH + hv -> C + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH'], ['C','H'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH_C, dummy)
            
          ! CH2CO -------------------------------------------
            ! CH2CO + hv -> ^3CH2 + CO
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH2CO'], ['^3CH2','CO   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH2CO_CH2, dummy)
            
          ! CH3CHO -------------------------------------------
            ! CH3CHO + hv -> CH3 + HCO
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3CHO'], ['CH3','HCO'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3CHO_CH3, dummy)

            ! CH3CHO + hv -> CH4 + CO
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3CHO'], ['CH4','CO '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3CHO_CH4, dummy)

          ! C2H5CHO -------------------------------------------
            ! C2H5CHO + hv -> C2H5 + HCO
            R_label = ch_identify(spl, ich, 'photodissociation', ['C2H5CHO'], ['C2H5','HCO '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C2H5CHO_C2H5, dummy)

          ! C3H3 -------------------------------------------
            ! C3H3 + hv -> C3H2 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['C3H3'], ['C3H2','H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_C3H3_C3H2, dummy)

          ! CH3C2H -------------------------------------------
            ! CH3C2H + hv -> C3H3 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3C2H'], ['C3H3','H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3C2H_C3H3, dummy)

            ! CH3C2H + hv -> C3H2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3C2H'], ['C3H2','H2  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3C2H_C3H2, dummy)

            ! CH3C2H + hv -> CH3 + C2H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH3C2H'], ['CH3','C2H'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH3C2H_CH3, dummy)

          ! CH2CCH2 -------------------------------------------
            ! CH2CCH2 + hv -> C3H3 + H
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH2CCH2'], ['C3H3','H   '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH2CCH2_C3H3, dummy)

            ! CH2CCH2 + hv -> C3H2 + H2
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH2CCH2'], ['C3H2','H2  '])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH2CCH2_C3H2, dummy)

            ! CH2CCH2 + hv -> C2H2 + ^3CH2
            R_label = ch_identify(spl, ich, 'photodissociation', ['CH2CCH2'], ['C2H2 ','^3CH2'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_CH2CCH2_C2H2, dummy)

          ! HCN -------------------------------------------
            ! HCN + hv -> H + CN
            R_label = ch_identify(spl, ich, 'photodissociation', ['HCN'], ['H ','CN'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_HCN_H, dummy)

          ! HNCO -------------------------------------------
            ! HNCO + hv -> H + NCO
            R_label = ch_identify(spl, ich, 'photodissociation', ['HNCO'], ['H  ','NCO'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_HNCO_H, dummy)

          ! NCO -------------------------------------------
            ! NCO + hv -> N + CO
            R_label = ch_identify(spl, ich, 'photodissociation', ['NCO'], ['N ','CO'])
            call p__UV_qy_adapt(xct, T, isp, ich, R_label, &
                 &                 xct%qy_NCO_N, dummy)            
              
          ! end list ---------------------------------------------------------------------------------


          !! if cross sections are less than 0, they are set to 0 ------
          !where (xct%sigma_d_UV(:,:,ich)<0.0_dp)
          !  xct%sigma_d_UV(:,:,ich) = 0.0_dp
          !end where

          ! output cross sections
          if (var%istep == 0) then 
            iz = 1 ! select altitude grid
            if (spl%product_list(ich,1)==1) then 
              fname = './UV/xct/dissociation/'//trim(ADJUSTL(reactants(1)))//' + hv -> '&
                &     //trim(ADJUSTL(products(1)))//'.dat'
            else if (spl%product_list(ich,1)==2) then 
              fname = './UV/xct/dissociation/'//trim(ADJUSTL(reactants(1)))//' + hv -> '&
                &     //trim(ADJUSTL(products(1)))//' + '//trim(ADJUSTL(products(2)))//'.dat'
            else if (spl%product_list(ich,1)==3) then 
              fname = './UV/xct/dissociation/'//trim(ADJUSTL(reactants(1)))//' + hv -> '&
                &     //trim(ADJUSTL(products(1)))//' + '//trim(ADJUSTL(products(2)))//' + '//trim(ADJUSTL(products(3)))//'.dat'
            end if
            open(11, file = fname, status = 'unknown' )
            swl = xct%wlrange_d_UV(1,ich)
            ewl = xct%wlrange_d_UV(2,ich)
            do iwl = swl, ewl
              write(11, *) flx%lambda_UV(iwl), xct%sigma_d_UV(iwl,iz,ich)
            end do
            close(11)
          end if

        end if ! end if reaction == photodissociation

      end do !ich

    end if ! end if label == dissociation


  end subroutine p__UV_cross_section_exe



  !=======================================================================================================================
  !
  !                Cross sections are automatically adapted by using following subroutines.
  !
  !=======================================================================================================================

  !----------------------------------------------------------------------------------------------
  ! Flux, cross sections and quantum yield data are automatically adapted to the model bins.
  !   if the data bin is larger than the model bin, the data is interpolated.
  !   if the data bin is smaller than the model bin, the data is binned.
  !----------------------------------------------------------------------------------------------
  subroutine p__UV_binning(flx, idata, nl, odata, nwl, swl, ewl)
    implicit none
    type(flx_), intent(in)  :: flx
    integer,    intent(in)  :: nl, nwl
    real(8),    intent(in)  :: idata(nl,2)
    real(8),    intent(out) :: odata(nwl)
    integer,    intent(out) :: swl, ewl
    integer iwl, il, label, il0, ndata
    real(8) idata_tmp(nl,2)
    real(dp) i0, ip, im, dip, dim, dipm, o0, op, om, dop, dom, dopm, sum_data, sum_wl, tmp

    if (idata(1,1) > idata(nl,1)) then 
      do il = 1, nl
        idata_tmp(il,1) = idata(nl-il+1,1)
        idata_tmp(il,2) = idata(nl-il+1,2)
      end do 
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
          i0 = idata(il,1)
          ip = idata(il+1,1)
          dip = ip - i0
          dim = dip
        else if (il > 1 .and. il < nl) then 
          im = idata(il-1,1)
          i0 = idata(il,1)
          ip = idata(il+1,1)
          dim = i0 - im
          dip = ip - i0
        else if (il == nl) then 
          im = idata(il-1,1)
          i0 = idata(il,1)
          dim = i0 - im
          dip = dim
        end if

        ! model bin
        if (iwl == 1) then 
          o0 = flx%lambda_UV(1)
          op = flx%lambda_UV(2)
          dop = op - o0
          dom = o0
        else if (iwl > 1 .and. iwl < nwl) then 
          om = flx%lambda_UV(iwl-1)
          o0 = flx%lambda_UV(iwl)
          op = flx%lambda_UV(iwl+1)
          dom = o0 - om
          dop = op - o0
        else if (iwl == nwl) then 
          om = flx%lambda_UV(iwl-1)
          o0 = flx%lambda_UV(iwl)
          dom = o0 - om
          dop = 1.0e10_dp ! to take large enough value
        end if 
        
        ! for binning
        if (o0-dom*0.5d0 <= i0 .and. i0 <= o0+dop*0.5d0) then 
          sum_data = sum_data + idata(il,2)*(dim+dip)
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
          tmp = (idata(il,2)*(ip-o0) + idata(il+1,2)*(o0-i0))/(ip-i0)
          if (label /= 1) label = 0
        end if

        ! to escape from this loop
        if (label == 1 .and. i0 > o0+dop*0.5d0) then 
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
  
  end subroutine p__UV_binning


  !------------------------------------------------------------
  ! For the use of reading cross section data easily.
  !------------------------------------------------------------
  subroutine p__UV_read_data(xct, flx, fname, isp, &
    &                        xctvar, ndata, unit   )
    implicit none
    type(xct_),            intent(inout) :: xct
    type(flx_),            intent(inout) :: flx
    real(dp), allocatable, intent(inout) :: xctvar(:,:)
    integer,               intent(in)    :: ndata, isp
    character(len=*),      intent(in)    :: unit, fname
    integer iwl, jnu, i, j, nh, il, nl, swl, ewl, swl1, ewl1
    real(dp), allocatable :: idata(:,:), odata(:), rdata(:,:)

    if (isp >= 1) then 

      write(*,'(a)',advance='no')  '  Reading datafile: '//trim(ADJUSTL(fname))//'...'

      allocate(odata(flx%nwl_UV))

      if (ndata == 1) then 

        call p__count_header(nh, fname)
        call p__count_lines(nl, fname)
        allocate(idata(nl,2))

        open(11, file = fname, status = 'unknown')
          do i = 1, nh; read(11,*); end do
          do il = 1, nl
            read(11,*) idata(il,1), idata(il,2) 
            if (unit == 'm2') then 
              if (idata(il,2)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "m" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
            end if
            if (unit == 'cm2') then 
              if (idata(il,2)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "cm" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
              idata(il,2) = idata(il,2) * 1.0e-4_dp ! [cm2 -> m2]
            else if (unit == 'cm-1, cm2') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
              idata(il,2) = idata(il,2) * 1.0e-4_dp ! [cm2 -> m2]
            else if (unit == 'cm-1, m2') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
            else if (unit == 'cm-1') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
            end if
          end do
        close(11)
        call p__UV_binning(flx, idata, nl, odata, flx%nwl_UV, swl, ewl)
        allocate(xctvar(ewl-swl+1,ndata+1))
        xctvar(:,1) = (/ (i,i=swl,ewl) /)
        xctvar(:,2) = odata(swl:ewl)
        deallocate(idata)

      else if (ndata >= 2) then 

        call p__count_header(nh, fname)
        call p__count_lines(nl, fname)
        allocate(idata(nl,2), rdata(nl,ndata))

        open(11, file = fname, status = 'unknown')
          do i = 1, nh; read(11,*); end do
          do il = 1, nl
            read(11,*) idata(il,1), (rdata(il,j), j = 1, ndata) 
            if (unit == 'm2') then 
              if (rdata(il,1)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "m" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
            else if (unit == 'cm2') then 
              if (rdata(il,1)>1.0e-3_dp) then 
                print *, ''
                print *, 'error!'
                print *, '  Is the datafile "'//trim(adjustl(fname))//'" quantum yield ?'
                print *, '  You indicated the unit of the data is in "cm" but it is probably quantum yield data..'
                print *, "  If it is quantum yield data, the unit should be ''."
                stop
              end if 
              rdata(il,:) = rdata(il,:) * 1.0e-4_dp ! [cm2 -> m2]
            else if (unit == 'cm-1, cm2') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
              rdata(il,:) = rdata(il,:) * 1.0e-4_dp ! [cm2 -> m2]
            else if (unit == 'cm-1, m2') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
            else if (unit == 'cm-1') then 
              idata(il,1) = 1.0e7_dp / idata(il,1) ! [cm-1 -> nm]
            end if
          end do
        close(11)

        idata(:,2) = rdata(:,1)
        call p__UV_binning(flx, idata, nl, odata, flx%nwl_UV, swl, ewl)
        allocate(xctvar(ewl-swl+1,ndata+1))
        xctvar(:,1) = (/ (i,i=swl,ewl) /)
        xctvar(:,2) = odata(swl:ewl)

        do i = 2, ndata
          idata(:,2) = rdata(:,i)
          call p__UV_binning(flx, idata, nl, odata, flx%nwl_UV, swl, ewl)
          xctvar(:,i+1) = odata(swl:ewl)
        end do 

        deallocate(idata, rdata)

      end if

      deallocate(odata)

      write(*,*) 'done.'

    end if

  end subroutine p__UV_read_data


  !------------------------------------------------------------
  ! Absorption cross section data are automatically adapted.
  !------------------------------------------------------------
  subroutine p__UV_sigma_a_adapt(xct, T, isp, &
    &                            xctvar, Tdata)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    real(dp),              intent(in)    :: xctvar(:,:), Tdata(:)
    integer i, iz, iwl, jwl, Tlabel, Tlabel1, nz, ndata, nwl, swl, swl_1, ewl
    real(dp) Tfrac, Tclamp, one_Tfrac

    xct%label_sigma_a_UV(isp) = 1

    ndata = size(Tdata)
    nwl   = size(xctvar)/(ndata+1) 
    nz    = size(T)

    if (ndata == 1) then 

      swl = nint(xctvar(1,1)); swl_1 = swl-1
      ewl = nint(xctvar(size(xctvar)/2,1))
      if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
      if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

      do iz = 1, nz
        !do iwl = swl, ewl
        !  jwl = iwl-swl_1
        !  xct%sigma_a_UV(iwl,iz,isp) = xctvar(jwl,2)
        !end do 
        xct%sigma_a_UV(swl:ewl,iz,isp) = xctvar(1:ewl-swl+1,2)
      end do

    else if (ndata >= 2) then 

      swl = nint(xctvar(1,1)); swl_1 = swl-1
      ewl = nint(xctvar(size(xctvar)/(ndata+1),1))
      if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
      if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

      do iz = 1, nz
        Tclamp = T(iz)
        if (T(iz)<=Tdata(1)) then 
          Tlabel = 2
          Tclamp = Tdata(1)
          Tfrac  = 0.0_dp
        end if
        do i = 1, size(Tdata)-1
          if (T(iz)>Tdata(i) .and. T(iz)<=Tdata(i+1)) then
            Tlabel = i+1
            Tfrac  = (Tclamp-Tdata(i))/(Tdata(i+1)-Tdata(i))
          end if
        end do
        if (T(iz)>Tdata(size(Tdata))) then 
          Tlabel = size(Tdata)
          Tclamp = Tdata(size(Tdata))
          Tfrac  = 1.0_dp
        end if
        one_Tfrac = 1.0_dp - Tfrac
        Tlabel1 = Tlabel+1
        
        !do iwl = swl, ewl
        !  jwl = iwl-swl_1
        !  xct%sigma_a_UV(iwl,iz,isp) = one_Tfrac*xctvar(jwl,Tlabel) + Tfrac*xctvar(jwl,Tlabel1)
        !end do 
        xct%sigma_a_UV(swl:ewl,iz,isp) = one_Tfrac*xctvar(1:ewl-swl+1,Tlabel) + Tfrac*xctvar(1:ewl-swl+1,Tlabel1)

      end do

    end if

    

  end subroutine p__UV_sigma_a_adapt


  !-----------------------------------------------------------------
  ! Absorption cross sections for O2 Schumann-Runge bands @ 130-500K
  !-----------------------------------------------------------------
  subroutine p__UV_sigma_a_O2_Schumann_Runge_bands(xct, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, jwl, nz, ndata, swl, ewl
    real(dp) Tz, delta, delta2

    xct%label_sigma_a_UV(isp) = 1

    nz    = size(T)

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

      swl = nint(xct%sigma_a_O2_SR_bands_130K(1,1))
      ewl = nint(xct%sigma_a_O2_SR_bands_130K(size(xct%sigma_a_O2_SR_bands_130K)/6,1))
      if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
      if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

      if (Tz >= 130.0_dp .and. Tz < 190.0_dp) then 
        do iwl = swl, ewl
          jwl = iwl - swl + 1
          xct%sigma_a_UV(iwl,iz,isp) = 1.0e-24_dp * ( xct%sigma_a_O2_SR_bands_130K(jwl,2) * delta2 &
            &                        +                xct%sigma_a_O2_SR_bands_130K(jwl,3) * delta &
            &                        +                xct%sigma_a_O2_SR_bands_130K(jwl,4) )
        end do
      else if (Tz >= 190.0_dp .and. Tz < 280.0_dp) then 
        do iwl = swl, ewl
          jwl = iwl - swl + 1
          xct%sigma_a_UV(iwl,iz,isp) = 1.0e-24_dp * ( xct%sigma_a_O2_SR_bands_190K(jwl,2) * delta2 &
            &                        +                xct%sigma_a_O2_SR_bands_190K(jwl,3) * delta &
            &                        +                xct%sigma_a_O2_SR_bands_190K(jwl,4) )
        end do
      else if (Tz >= 280.0_dp .and. Tz <= 500.0_dp) then 
        do iwl = swl, ewl
          jwl = iwl - swl + 1
          xct%sigma_a_UV(iwl,iz,isp) = 1.0e-24_dp * ( xct%sigma_a_O2_SR_bands_280K(jwl,2) * delta2 &
            &                        +                xct%sigma_a_O2_SR_bands_280K(jwl,3) * delta &
            &                        +                xct%sigma_a_O2_SR_bands_280K(jwl,4) )
        end do
      end if

    end do


  end subroutine p__UV_sigma_a_O2_Schumann_Runge_bands


  !------------------------------------------------------
  ! Absorption cross sectinos for O2 Herzberg continuum
  !------------------------------------------------------
  subroutine p__UV_sigma_a_O2_Herzberg_continuum(xct, flx, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(flx_),            intent(in)    :: flx
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, swl, ewl, nz
    real(dp) l, larr(flx%nwl_UV,4)

    xct%label_sigma_a_UV(isp) = 1

    nz = size(T)

    do iwl = 1, flx%nwl_UV-1
      if (     flx%lambda_UV(iwl) <= 193.0_dp .and. 193.0_dp <  flx%lambda_UV(iwl+1)) then 
        swl = iwl
      else if (flx%lambda_UV(iwl) <  245.0_dp .and. 245.0_dp <= flx%lambda_UV(iwl+1)) then 
        ewl = iwl
      end if

      l = flx%lambda_UV(iwl)
      larr(iwl,1) = l
      larr(iwl,2) = l*l
      larr(iwl,3) = larr(iwl,2)*l
      larr(iwl,4) = larr(iwl,3)*l
    end do 
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    ! Yoshino et al. (1992)
    iz = 1
      !do iwl = swl, ewl
      !  l = flx%lambda_UV(iwl)
      !  xct%sigma_a_UV(iwl,iz,isp) = xct%sigma_a_UV(iwl,iz,isp) &
      !    & + 1.0e-28_dp*(-2.3837947e4_dp &
      !    &               +4.1973085e2_dp * larr(iwl,1) &
      !    &               -2.7640139e0_dp * larr(iwl,2) &
      !    &               +8.0723193e-3_dp* larr(iwl,3) &
      !    &               -8.8255447e-6_dp* larr(iwl,4) )
      !end do

    xct%sigma_a_UV(swl:ewl,iz,isp) = xct%sigma_a_UV(swl:ewl,iz,isp) &
        & + 1.0e-28_dp*(-2.3837947e4_dp &
        &               +4.1973085e2_dp * larr(swl:ewl,1) &
        &               -2.7640139e0_dp * larr(swl:ewl,2) &
        &               +8.0723193e-3_dp* larr(swl:ewl,3) &
        &               -8.8255447e-6_dp* larr(swl:ewl,4) )

    do iz = 2, nz
      xct%sigma_a_UV(swl:ewl,iz,isp) = xct%sigma_a_UV(swl:ewl,1,isp) 
    end do


  end subroutine p__UV_sigma_a_O2_Herzberg_continuum


  !------------------------------------------------------
  ! Absorption cross sectinos for H2O2 > 260 nm by JPL 2015
  !------------------------------------------------------
  subroutine p__UV_sigma_a_H2O2_analytic(xct, flx, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(flx_),            intent(in)    :: flx
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, swl, ewl, nz
    real(dp) l, tmpA(flx%nwl_UV), tmpB(flx%nwl_UV), A(0:7), B(0:4), x, one_x

    xct%label_sigma_a_UV(isp) = 1

    nz = size(T)

    do iwl = 1, flx%nwl_UV-1
      if (     flx%lambda_UV(iwl) <  260.0_dp .and. 260.0_dp <= flx%lambda_UV(iwl+1)) then 
        swl = iwl+1
      else if (flx%lambda_UV(iwl) <  340.0_dp .and. 340.0_dp <= flx%lambda_UV(iwl+1)) then 
        ewl = iwl
      end if
    end do 
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    A(0) = 6.4761e-4_dp
    A(1) = -9.2170972e2_dp
    A(2) = 4.535649
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

    tmpA = 0.0_dp
    tmpB = 0.0_dp
    do iwl = swl, ewl
      l = flx%lambda_UV(iwl)
      do i = 0, 7
        tmpA(iwl) = tmpA(iwl) + A(i) * l**dble(i) 
      end do
      tmpA(iwl) = tmpA(iwl) * 1.0e-25_dp
      do i = 0, 4
        tmpB(iwl) = tmpB(iwl) + B(i) * l**dble(i)
      end do
      tmpB(iwl) = tmpB(iwl) * 1.0e-25_dp
    end do

    do iz = 1, nz
      x = 1.0_dp / (1.0_dp + dexp(-1265.0_dp/T(iz)))
      one_x = 1.0_dp - x
      !do iwl = swl, ewl
      !  xct%sigma_a_UV(iwl,iz,isp) = (x*tmpA(iwl) + one_x*tmpB(iwl))
      !end do
      xct%sigma_a_UV(swl:ewl,iz,isp) = (x*tmpA(swl:ewl) + one_x*tmpB(swl:ewl))
    end do

  end subroutine p__UV_sigma_a_H2O2_analytic


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for H2CO
  !-----------------------------------------------------------------
  subroutine p__UV_sigma_a_H2CO_Tdependent(xct, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, jwl, nz, ndata, swl, swl_1, ewl
    real(dp) Tz, Tz_298

    xct%label_sigma_a_UV(isp) = 1

    nz    = size(T)

    ! Temperature gradient by Meller and Moortgat(2000)
    swl = nint(xct%sigma_gamma_H2CO(1,1)); swl_1 = swl-1
    ewl = nint(xct%sigma_gamma_H2CO(size(xct%sigma_gamma_H2CO)/2,1))
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    do iz = 1, nz
      Tz = T(iz)
      if (T(iz)<=223.0_dp) then 
        Tz = 223.0_dp
      else if (T(iz)>323.0_dp) then 
        Tz = 323.0_dp
      end if
      Tz_298 = Tz - 298.0_dp

      !do iwl = swl, ewl
      !  jwl = iwl-swl_1
      !  xct%sigma_a_UV(iwl,iz,isp) = xct%sigma_a_UV(iwl,iz,isp) + xct%sigma_gamma_H2CO(jwl,2) * Tz_298
      !end do 
      xct%sigma_a_UV(swl:ewl,iz,isp) = xct%sigma_a_UV(swl:ewl,iz,isp) + xct%sigma_gamma_H2CO(1:ewl-swl+1,2) * Tz_298
    end do

  end subroutine p__UV_sigma_a_H2CO_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for N2O5
  !-----------------------------------------------------------------
  subroutine p__UV_sigma_a_N2O5_Tdependent(xct, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, jwl, nz, ndata, swl, swl_1, ewl
    real(dp) Tz, invTz1000

    xct%label_sigma_a_UV(isp) = 1

    nz    = size(T)

    ! Harwood et al., 1993
    swl = nint(xct%sigma_a_N2O5_above260nm(1,1)); swl_1 = swl-1
    ewl = nint(xct%sigma_a_N2O5_above260nm(size(xct%sigma_a_N2O5_above260nm)/3,1))
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    do iz = 1, nz
      Tz = T(iz)
      if (T(iz)<=223.0_dp) then 
        Tz = 223.0_dp
      else if (T(iz)>295.0_dp) then 
        Tz = 295.0_dp
      end if
      invTz1000 = 1000.0_dp / Tz

      !do iwl = swl, ewl
      !  jwl = iwl-swl_1
      !  xct%sigma_a_UV(iwl,iz,isp) = 10.0_dp**(xct%sigma_a_N2O5_above260nm(jwl,2) &
      !    &                                  + xct%sigma_a_N2O5_above260nm(jwl,3) * invTz1000) &
      !    &                        * 1.0e-4_dp
      !end do 
      xct%sigma_a_UV(swl:ewl,iz,isp) = 10.0_dp**(xct%sigma_a_N2O5_above260nm(1:ewl-swl+1,2) &
          &                                    + xct%sigma_a_N2O5_above260nm(1:ewl-swl+1,3) * invTz1000) &
          &                          * 1.0e-4_dp
    end do

  end subroutine p__UV_sigma_a_N2O5_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for HNO3
  !-----------------------------------------------------------------
  subroutine p__UV_sigma_a_HNO3_Tdependent(xct, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, jwl, nz, ndata, swl, swl_1, ewl
    real(dp) Tz

    xct%label_sigma_a_UV(isp) = 1

    nz    = size(T)

    ! Burkholder et al., 1993
    swl = nint(xct%sigma_a_HNO3_Tdep(1,1)); swl_1 = swl-1
    ewl = nint(xct%sigma_a_HNO3_Tdep(size(xct%sigma_a_HNO3_Tdep)/2,1))
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    do iz = 1, nz
      Tz = T(iz)
      !do iwl = swl, ewl
      !  jwl = iwl-swl_1
      !  xct%sigma_a_UV(iwl,iz,isp) = xct%sigma_a_UV(iwl,iz,isp) &
      !    &                        * dexp(xct%sigma_a_HNO3_Tdep(jwl,2) * 1.0e-3_dp * (Tz-298.0_dp))
      !end do 
      xct%sigma_a_UV(swl:ewl,iz,isp) = xct%sigma_a_UV(swl:ewl,iz,isp) &
          &                          * dexp(xct%sigma_a_HNO3_Tdep(1:ewl-swl+1,2) * 1.0e-3_dp * (Tz-298.0_dp))
    end do


  end subroutine p__UV_sigma_a_HNO3_Tdependent


  !-----------------------------------------------------------------
  ! Temperature dependent of absorption cross sections for HO2NO2
  !-----------------------------------------------------------------
  subroutine p__UV_sigma_a_HO2NO2_Tdependent(xct, T, isp)
    implicit none
    integer,               intent(in)    :: isp
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    integer i, iz, iwl, jwl, nz, ndata, swl, swl_1, ewl
    real(dp) Tz, invQ

    xct%label_sigma_a_UV(isp) = 1

    nz    = size(T)

    ! Knight et al. 2002
    swl = nint(xct%sigma_a_HO2NO2_Tdep(1,1)); swl_1 = swl-1
    ewl = nint(xct%sigma_a_HO2NO2_Tdep(size(xct%sigma_a_HO2NO2_Tdep)/3,1))
    if (swl<xct%wlrange_a_UV(1,isp)) xct%wlrange_a_UV(1,isp) = swl
    if (ewl>xct%wlrange_a_UV(2,isp)) xct%wlrange_a_UV(2,isp) = ewl

    do iz = 1, nz
      Tz = T(iz)
      if (Tz < 273.0_dp) Tz = 273.0_dp
      if (Tz > 343.0_dp) Tz = 343.0_dp
      invQ = 1.0_dp + dexp(-988.0_dp/(0.69_dp*Tz))
      !do iwl = swl, ewl
      !  jwl = iwl-swl_1
      !  xct%sigma_a_UV(iwl,iz,isp) = ( xct%sigma_a_HO2NO2_Tdep(jwl,2) * invQ &
      !    &                          + xct%sigma_a_HO2NO2_Tdep(jwl,3) * (1.0_dp - invQ)) &
      !    &                        * 1.0e-24_dp
      !end do 
      xct%sigma_a_UV(swl:ewl,iz,isp) = ( xct%sigma_a_HO2NO2_Tdep(1:ewl-swl+1,2) * invQ &
          &                            + xct%sigma_a_HO2NO2_Tdep(1:ewl-swl+1,3) * (1.0_dp - invQ)) &
          &                          * 1.0e-24_dp
    end do

  end subroutine p__UV_sigma_a_HO2NO2_Tdependent


  !-------------------------------------------------------------
  ! Photodissociation cross sections are automatically adapted.
  !-------------------------------------------------------------
  subroutine p__UV_sigma_d_adapt(xct, T, ich, R_label, &
    &                            xctvar, Tdata)
    implicit none
    integer,               intent(in)    :: ich, R_label
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    real(dp),              intent(in)    :: xctvar(:,:), Tdata(:)
    integer i, iz, iwl, jwl, Tlabel, Tlabel1, nz, ndata, nwl, swl, swl_1, ewl
    real(dp) Tfrac, one_Tfrac, Tclamp

    if (R_label == 1) then 

      ndata = size(Tdata)
      nwl   = size(xctvar)/(ndata+1) 
      nz    = size(T)

      if (ndata == 1) then 

        swl = nint(xctvar(1,1)); swl_1 = swl-1
        ewl = nint(xctvar(size(xctvar)/2,1))
        if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
        if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

        do iz = 1, nz
          !do iwl = swl, ewl
          !  jwl = iwl-swl_1
          !  xct%sigma_d_UV(iwl,iz,ich) = xctvar(jwl,2)
          !end do 
          xct%sigma_d_UV(swl:ewl,iz,ich) = xctvar(1:ewl-swl+1,2)
        end do

      else if (ndata >= 2) then 

        swl = nint(xctvar(1,1)); swl_1 = swl-1
        ewl = nint(xctvar(size(xctvar)/(ndata+1),1))
        if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
        if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

        do iz = 1, nz
          Tclamp = T(iz)
          if (T(iz)<=Tdata(1)) then 
            Tlabel = 2
            Tclamp = Tdata(1)
            Tfrac  = 0.0_dp
          end if
          do i = 1, size(Tdata)-1
            if (T(iz)>Tdata(i) .and. T(iz)<=Tdata(i+1)) then
              Tlabel = i+1
              Tfrac  = (Tclamp-Tdata(i))/(Tdata(i+1)-Tdata(i))
            end if
          end do
          if (T(iz)>Tdata(size(Tdata))) then 
            Tlabel = size(Tdata)
            Tclamp = Tdata(size(Tdata))
            Tfrac  = 1.0_dp
          end if
          one_Tfrac = 1.0_dp - Tfrac
          Tlabel1 = Tlabel+1

          !do iwl = swl, ewl
          !  jwl = iwl-swl_1
          !  xct%sigma_d_UV(iwl,iz,ich) = one_Tfrac*xctvar(jwl,Tlabel) + Tfrac*xctvar(jwl,Tlabel1)
          !end do 
          xct%sigma_d_UV(swl:ewl,iz,ich) = one_Tfrac*xctvar(1:ewl-swl+1,Tlabel) + Tfrac*xctvar(1:ewl-swl+1,Tlabel1)

        end do

      end if

    end if
    

  end subroutine p__UV_sigma_d_adapt


  !-------------------------------------------------
  ! Quantum yield
  !-------------------------------------------------
  subroutine p__UV_qy_adapt(xct, T, isp, ich, R_label, &
    &                       qyvar, Tdata)
    implicit none
    integer,               intent(in)    :: isp, ich, R_label
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    real(dp),              intent(in)    :: qyvar(:,:), Tdata(:)
    integer i, iz, iwl, jwl, Tlabel, Tlabel1, nz, ndata, nwl, swl, swl_1, ewl
    real(dp) Tfrac, one_Tfrac, Tclamp, qytmp

    if (R_label == 1) then 

      ndata = size(Tdata)
      nwl   = size(qyvar)/(ndata+1) 
      nz    = size(T)

      if (ndata == 1) then 

        swl = nint(qyvar(1,1)); swl_1 = swl-1
        ewl = nint(qyvar(size(qyvar)/2,1))
        if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
        if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

        do iz = 1, nz
          !do iwl = swl, ewl
          !  jwl = iwl-swl_1
          !  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qyvar(jwl,2)
          !end do 
          xct%sigma_d_UV(swl:ewl,iz,ich) = xct%sigma_a_UV(swl:ewl,iz,isp) * qyvar(1:ewl-swl+1,2)
        end do

      else if (ndata >= 2) then 

        swl = nint(qyvar(1,1)); swl_1 = swl-1
        ewl = nint(qyvar(size(qyvar)/(ndata+1),1))
        if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
        if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

        do iz = 1, nz
          Tclamp = T(iz)
          if (T(iz)<=Tdata(1)) then 
            Tlabel = 2
            Tclamp = Tdata(1)
            Tfrac  = 0.0_dp
          end if
          do i = 1, size(Tdata)-1
            if (T(iz)>Tdata(i) .and. T(iz)<=Tdata(i+1)) then
              Tlabel = i+1
              Tfrac  = (Tclamp-Tdata(i))/(Tdata(i+1)-Tdata(i))
            end if
          end do
          if (T(iz)>Tdata(size(Tdata))) then 
            Tlabel = size(Tdata)
            Tclamp = Tdata(size(Tdata))
            Tfrac  = 1.0_dp
          end if
          one_Tfrac = 1.0_dp - Tfrac
          Tlabel1 = Tlabel+1

          !do iwl = swl, ewl
          !  jwl = iwl-swl_1
          !  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) &
          !    &                        * (one_Tfrac*qyvar(jwl,Tlabel) + Tfrac*qyvar(jwl,Tlabel1))
          !end do 
          xct%sigma_d_UV(swl:ewl,iz,ich) = xct%sigma_a_UV(swl:ewl,iz,isp) &
            &                            * (one_Tfrac*qyvar(1:ewl-swl+1,Tlabel) + Tfrac*qyvar(1:ewl-swl+1,Tlabel1))

        end do

      end if

    end if

  end subroutine p__UV_qy_adapt


  !------------------------------------------------------------------------------------
  ! Analytic expression of quantum yield for O3 photodissociation by JPL 2015
  !------------------------------------------------------------------------------------
  subroutine p__UV_qy_O3(xct, flx, T, isp, ich, R_label, OorO1D)
    implicit none
    integer,               intent(in)    :: isp, ich, R_label
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    type(flx_),            intent(in)    :: flx
    character(len=*),      intent(in)    :: OorO1D 
    integer i, iz, iwl, nz, l193, l225, l306, l329, l340, swl, ewl
    real(dp) l, Tclamp, qy_O1D
    real(dp) X1, X2, X3, w1, w2, w3, A1, A2, A3, v1, v2, c, R, q1, q2, q3

    if (R_label == 1) then 

      X1 = 304.225_dp; X2 = 314.957_dp; X3 = 310.737
      w1 = 5.576_dp  ; w2 = 6.601_dp  ; w3 = 2.187_dp
      A1 = 0.8036_dp ; A2 = 8.9061_dp ; A3 = 0.1192_dp
      v1 = 0.0_dp    ; v2 = 825.518_dp
      c  = 0.0765_dp
      R  = 0.0695_dp

      nz = size(T)

      do iwl = 1, flx%nwl_UV-1
        if (     flx%lambda_UV(iwl) <  193.0_dp .and. 193.0_dp <= flx%lambda_UV(iwl+1)) then 
          l193 = iwl+1
        else if (flx%lambda_UV(iwl) <= 225.0_dp .and. 225.0_dp < flx%lambda_UV(iwl+1)) then 
          l225 = iwl
        else if (flx%lambda_UV(iwl) <= 306.0_dp .and. 306.0_dp < flx%lambda_UV(iwl+1)) then 
          l306 = iwl
        else if (flx%lambda_UV(iwl) <= 329.0_dp .and. 329.0_dp < flx%lambda_UV(iwl+1)) then 
          l329 = iwl
        else if (flx%lambda_UV(iwl) <= 340.0_dp .and. 340.0_dp < flx%lambda_UV(iwl+1)) then 
          l340 = iwl
        end if
      end do 

      swl = xct%wlrange_a_UV(1,isp)
      ewl = xct%wlrange_a_UV(2,isp)
      if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
      if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

      !   Evaluation: J. B. Burkholder, S. P. Sander, J. Abbatt, J. R. Barker, R. E. Huie, C. E. Kolb, 
      !               M. J. Kurylo, V. L. Orkin, D. M. Wilmouth, and P. H. Wine 
      !               "Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies, Evaluation No. 18," 
      !               JPL Publication 15-10, Jet Propulsion Laboratory, Pasadena, 2015 http://jpldataeval.jpl.nasa.gov.
      
      do iz = 1, nz

        Tclamp = T(iz)
        if (Tclamp < 200.0_dp) Tclamp = 200.0_dp
        if (Tclamp > 320.0_dp) Tclamp = 320.0_dp

        do iwl = l193, l225
          l = flx%lambda_UV(iwl)
          qy_O1D = 1.37e-2_dp * l - 2.16_dp
          if (OorO1D == 'O(1D)') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_O1D
          if (OorO1D == 'O')     xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp-qy_O1D)
        end do

        do iwl = l225+1, l306-1
          qy_O1D = 0.9_dp
          if (OorO1D == 'O(1D)') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_O1D
          if (OorO1D == 'O')     xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp-qy_O1D)
        end do

        ! Analytic expression of O3 -> O1D at wavelength range 306-328 nm and at temperature range 200-320 K
        !   based on the review of Matsumi et al. (2002) doi:10.1029/2001JD000510.
        q1 = dexp(-v1/(R*Tclamp))
        q2 = dexp(-v2/(R*Tclamp))
        q3 = Tclamp/300.0_dp
        do iwl = l306, l329-1
          l = flx%lambda_UV(iwl)
          qy_O1D = q1/(q1+q2)              * A1 * dexp(-((X1-l)/w1)**4.0_dp) &
            &    + q2/(q1+q2) * q3**2.0_dp * A2 * dexp(-((X2-l)/w2)**2.0_dp) &
            &    +              q3**1.5_dp * A3 * dexp(-((X3-l)/w3)**2.0_dp) &
            &    + c
          if (OorO1D == 'O(1D)') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_O1D
          if (OorO1D == 'O')     xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp-qy_O1D)
        end do

        do iwl = l329, l340
          qy_O1D = 0.08_dp
          if (OorO1D == 'O(1D)') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_O1D
          if (OorO1D == 'O')     xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp-qy_O1D)
        end do

        do iwl = l340+1, flx%nwl_UV
          qy_O1D = 0.0_dp
          if (OorO1D == 'O(1D)') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_O1D
          if (OorO1D == 'O')     xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * (1.0_dp-qy_O1D)
        end do

      end do

    end if

  end subroutine p__UV_qy_O3


  !------------------------------------------------------------------------------------
  ! Analytic expression of quantum yield for O3 photodissociation by JPL 2015
  !------------------------------------------------------------------------------------
  subroutine p__UV_qy_H2CO(var, xct, flx, T, isp, ich, R_label, HorCO)
    implicit none
    integer,               intent(in)    :: isp, ich, R_label
    real(dp),              intent(in)    :: T(:)
    type(xct_),            intent(inout) :: xct
    type(flx_),            intent(in)    :: flx
    type(var_),            intent(in)    :: var
    character(len=*),      intent(in)    :: HorCO
    integer i, iz, iwl, jwl, nz, l250, l330, l338, l360, swl, ewl
    real(dp) l, Tclamp, qy_H(flx%nwl_UV), qy_CO_300K, qy_CO_T, a_300K, a_T, P, kB
    real(dp) a(0:4)

    kB = 1.38064852e-23_dp

    if (R_label == 1) then 

      swl = xct%wlrange_a_UV(1,isp)
      ewl = xct%wlrange_a_UV(2,isp)
      if (swl<xct%wlrange_d_UV(1,ich)) xct%wlrange_d_UV(1,ich) = swl
      if (ewl>xct%wlrange_d_UV(2,ich)) xct%wlrange_d_UV(2,ich) = ewl

      a(0) = 557.95835182_dp
      a(1) = -7.31994058026_dp
      a(2) = 0.03553521598_dp
      a(3) = -7.54849718e-5_dp
      a(4) = 5.91001021e-8_dp

      nz = size(T)

      do iwl = 1, flx%nwl_UV-1
        if (     flx%lambda_UV(iwl) <  250.0_dp .and. 250.0_dp <= flx%lambda_UV(iwl+1)) then 
          l250 = iwl+1
        else if (flx%lambda_UV(iwl) <= 330.0_dp .and. 330.0_dp < flx%lambda_UV(iwl+1)) then 
          l330 = iwl
        else if (flx%lambda_UV(iwl) <= 338.0_dp .and. 338.0_dp < flx%lambda_UV(iwl+1)) then 
          l338 = iwl
        else if (flx%lambda_UV(iwl) <= 360.0_dp .and. 360.0_dp < flx%lambda_UV(iwl+1)) then 
          l360 = iwl
        end if
      end do 

      do iwl = l250, l360
        l = flx%lambda_UV(iwl)
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

        P = var%n_tot(iz) * kB * T(iz) / 101325.0_dp ! [atm]

        do iwl = l250, l360

          jwl = iwl - nint(xct%qy_H2CO_CO_300K(1,1)) + 1
          qy_CO_300K = xct%qy_H2CO_CO_300K(jwl,2)

          ! Temperature and pressure dependent yield of H2 + CO
          a_300K = 1.0_dp / qy_CO_300K - 1.0_dp / (1.0_dp - qy_H(iwl))
          a_T    = a_300K * ( 1.0_dp + 0.05_dp * (l-329.0_dp)*((300.0_dp-Tclamp)/80.0_dp) )
          qy_CO_T = 1.0_dp / ( 1.0_dp / (1.0_dp-qy_H(iwl)) + a_T*P )

          if (HorCO == 'H')  xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_H(iwl)
          if (HorCO == 'CO') xct%sigma_d_UV(iwl,iz,ich) = xct%sigma_a_UV(iwl,iz,isp) * qy_CO_T

          !print *, l, qy_H, qy_CO_300K, qy_CO_T
          
        end do

      end do
      !stop

    end if

  end subroutine p__UV_qy_H2CO


  !-------------------------------------------------
  ! cross section data for UV
  !-------------------------------------------------
  subroutine p__UV_EUV_xct_convergence_exe(nz, spl, var, flx, & ! in
    &                                      xct                ) ! inout
    implicit none
    integer,    intent(in)    :: nz
    type(spl_), intent(in)    :: spl
    type(var_), intent(in)    :: var
    type(flx_), intent(in)    :: flx
    type(xct_), intent(inout) :: xct
    integer iwl, jwl, i, iz, isp, nl, swl, ewl
    character(len=256) fname, num
    real(dp), allocatable :: idata(:,:), odata(:)

    nl = 37

    xct%label_sigma_a_UV_EUV_conv = 0

    allocate(odata(flx%nwl_UV))

    !--------------------------------
    ! UV absorption cross section
    !--------------------------------
    if (xct%type == 'absorption') then
      do isp = 1, spl%nsp

        allocate(idata(nl,2))
        idata(:,1) = flx%lambda_EUV(:)
        idata(:,2) = xct%sigma_a_EUV(:,isp)
        call p__UV_binning(flx, idata, nl, odata, flx%nwl_UV, swl, ewl)
        do iwl = swl, ewl
          if (xct%sigma_a_UV_EUV(iwl,isp) == 0.0_dp) then 
            xct%label_sigma_a_UV_EUV_conv(isp) = 1
            xct%sigma_a_UV_EUV(iwl,isp) = odata(iwl)
          end if
        end do
        deallocate(idata)

      end do
    end if

    deallocate(odata)

  end subroutine p__UV_EUV_xct_convergence_exe


end module p__UV
