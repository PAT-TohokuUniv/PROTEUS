module p__EUV
  
  use v__tdec,                 only : var_, grd_, cst_, xct_, spl_, flx_, set_
  use c__prm,                  only : c__prm__ini
  use p__photoionization_rate, only : p__photoionization_rate__ini, load_cross_section_dat, get_cross_section
  use p__search,               only : p__search_reactant, p__search_product, ch_identify, sp_index


  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__EUV_flux, p__EUV_cross_section, p__EUVAC_photoionization_rate

contains


  !-------------------------------------------------
  ! solar EUV flux model: EUVAC
  !-------------------------------------------------
  subroutine p__EUV_flux(spl, grd, set, & ! in
    &                    var, xct, flx  ) ! inout
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer i, iwl, IDIM, KMAX, nwl_tmp, nwli, nwl
    real(dp) F107, F107A
    real(dp), allocatable :: indata(:,:), dl(:), l(:)
    real, allocatable :: BINLAM(:), XS_OUT(:), BIN_ANGS(:)

    flx%multiplying_factor_EUV = set%euv_factor

    if (set%euv_input == "HEUVAC") then

      IDIM = 10550
      KMAX = IDIM-1
      allocate(BINLAM(IDIM),XS_OUT(IDIM), BIN_ANGS(IDIM))
      allocate(indata(IDIM,2))
      do iwl = 1, IDIM
        BINLAM(iwl) = dble(iwl) / 10.0d0
      end do
      F107 = set%F107 ! F10.7 at Earth
      F107A = set%F107

      call HEUVAC(IDIM, &  !.. in: Array dimensions
        &         F107, &  !.. in: daily 10.7 cm flux index. 
        &        F107A, &  !.. in: 81 day average of daily F10.7 
        &         KMAX, &  !.. in: number of bins for the flux
        &       BINLAM, &  !.. in: the wavelength (Angstrom) bin boundaries 
        &     BIN_ANGS, &  !.. out: fluxes (photons/cm2/s) in the bins
        &       XS_OUT)    !.. out: Weighted cross section in the bins
      
      do iwl = 1, IDIM-1
        indata(iwl,1) = dble(BINLAM(iwl)) / 10.0d0
        indata(iwl,2) = dble(BIN_ANGS(iwl)) / dble(BINLAM(iwl+1)-BINLAM(iwl)) * 10.0d0
      end do
      indata(IDIM,1) = dble(BINLAM(IDIM)) / 10.0d0
      indata(IDIM,1) = dble(BIN_ANGS(IDIM)) / dble(BINLAM(IDIM)-BINLAM(IDIM-1)) * 10.0d0

      nwl_tmp = 0
      do i = 1, set%n_wl_bin
        nwli = nint( (set%wl_bin(i,2) - set%wl_bin(i,1)) / set%wl_bin(i,3) )
        nwl_tmp = nwl_tmp + nwli 
      end do 

      allocate(dl(nwl_tmp), l(nwl_tmp))

      nwl = 0
      loop: do i = 1, set%n_wl_bin
        nwli = nint( (set%wl_bin(i,2) - set%wl_bin(i,1)) / set%wl_bin(i,3) )
        do iwl = 1, nwli
          dl(nwl+iwl) = set%wl_bin(i,3)
          l(nwl+iwl)  = set%wl_bin(i,1) + set%wl_bin(i,3) / 2.0_dp + set%wl_bin(i,3) * dble(iwl-1)
          if (l(nwl+iwl) >= 105.5_dp) then 
            flx%nwl_EUV = nwl+iwl
            exit loop
          end if
        end do 
        nwl = nwl + nwli
      end do loop

      allocate(flx%lambda_EUV(flx%nwl_EUV), flx%dlambda_EUV(flx%nwl_EUV))
      allocate(flx%solar_EUV(flx%nwl_EUV))
      allocate(var%tau_EUV(flx%nwl_EUV,grd%nz))
      allocate(var%I_EUV(flx%nwl_EUV,grd%nz))
      allocate(xct%sigma_a_EUV(flx%nwl_EUV,0:spl%nsp), xct%sigma_i_EUV(flx%nwl_EUV,0:spl%nch))
      allocate(xct%label_sigma_a_EUV(0:spl%nsp))

      do iwl = 1, flx%nwl_EUV
        flx%lambda_EUV(iwl) = l(iwl)
        flx%dlambda_EUV(iwl) = dl(iwl)
      end do

      call p__HEUVAC_binning(flx, indata, IDIM, flx%solar_EUV, flx%nwl_EUV)

      ! photon flux -----------------------------------------------------
      do iwl = 1, flx%nwl_EUV
        flx%solar_EUV(iwl) = flx%multiplying_factor_EUV * flx%solar_EUV(iwl) * 1.0e4_dp / (flx%orbit * flx%orbit)
      end do

      open(10, file=trim(ADJUSTL(set%dir_name))//'/input/EUV_flux_out.dat', status='unknown')
        do iwl = 1, flx%nwl_EUV
          write(10,*) flx%lambda_EUV(iwl), flx%solar_EUV(iwl)
        end do
      close(10)

      deallocate(indata, BINLAM, XS_OUT, BIN_ANGS)

      do i = 1, flx%nwl_UV
        if (flx%lambda_UV(i)<105.5_dp) then 
          flx%solar_UV(i) = flx%solar_EUV(i)
        end if
      end do 

    else if (set%euv_input == "EUV-DATA") then

      nwl_tmp = 0
      do i = 1, set%n_wl_bin
        nwli = nint( (set%wl_bin(i,2) - set%wl_bin(i,1)) / set%wl_bin(i,3) )
        nwl_tmp = nwl_tmp + nwli 
      end do 

      allocate(dl(nwl_tmp), l(nwl_tmp))

      nwl = 0
      loop2: do i = 1, set%n_wl_bin
        nwli = nint( (set%wl_bin(i,2) - set%wl_bin(i,1)) / set%wl_bin(i,3) )
        do iwl = 1, nwli
          dl(nwl+iwl) = set%wl_bin(i,3)
          l(nwl+iwl)  = set%wl_bin(i,1) + set%wl_bin(i,3) / 2.0_dp + set%wl_bin(i,3) * dble(iwl-1)
          if (l(nwl+iwl) >= 105.5_dp) then 
            flx%nwl_EUV = nwl+iwl
            exit loop2
          end if
        end do 
        nwl = nwl + nwli
      end do loop2

      allocate(flx%lambda_EUV(flx%nwl_EUV), flx%dlambda_EUV(flx%nwl_EUV))
      allocate(flx%solar_EUV(flx%nwl_EUV))
      allocate(var%tau_EUV(flx%nwl_EUV,grd%nz))
      allocate(var%I_EUV(flx%nwl_EUV,grd%nz))
      allocate(xct%sigma_a_EUV(flx%nwl_EUV,0:spl%nsp), xct%sigma_i_EUV(flx%nwl_EUV,0:spl%nch))
      allocate(xct%label_sigma_a_EUV(0:spl%nsp))

      do iwl = 1, flx%nwl_EUV
        flx%lambda_EUV(iwl) = l(iwl)
        flx%dlambda_EUV(iwl) = dl(iwl)
        flx%solar_EUV(iwl) = flx%solar_UV(iwl)
      end do

      deallocate(dl, l)

    end if

  end subroutine p__EUV_flux

  !-------------------------------------------------
  ! cross section data for EUV
  !-------------------------------------------------
  subroutine p__EUV_cross_section(spl, grd, flx, var) ! inout
    type(spl_),   intent(in)     :: spl
    type(grd_),   intent(inout)  :: grd
    type(flx_),   intent(inout)  :: flx
    type(var_),   intent(inout)  :: var

    call p__photoionization_rate__ini(flx%nwl_EUV, grd%nz, spl%nsp, spl%nch,                     & ! in
      &                               spl%species(1:), var%m(1:), spl%reaction_type_list,        & ! in
      &                               spl%reactant_list, spl%product_list,                       & ! in
      &                               flx%lambda_EUV(1:), flx%dlambda_EUV(1:), flx%solar_EUV(1:) ) ! in
    call load_cross_section_dat('./EUV/xsect_data') ! in
    call get_cross_section('absorption', 1) ! in
    call get_cross_section('photoionization', 1) ! in

  end subroutine p__EUV_cross_section

  !-------------------------------------------------
  ! Photoionization rate calculation
  !-------------------------------------------------
  subroutine p__EUVAC_photoionization_rate(cst, grd, xct, flx, spl, & ! in
    &                                      var, J_rate              ) ! inout
    implicit none
    type(cst_), intent(in)     :: cst
    type(grd_), intent(in)     :: grd
    type(xct_), intent(in)     :: xct
    type(flx_), intent(in)     :: flx
    type(spl_), intent(in)     :: spl
    type(var_), intent(inout)  :: var
    real(dp),   intent(inout)  :: J_rate(grd%nz,spl%nch)

    integer isp, ich, jch, iz, iwl, swl, ewl

    ! for solar zenith angle near and greater than 90deg [Smith et al., 1972]
    real(dp) Hz, yz, Xz, chiz, Chfunc, cln_Ch, g
    ! ap : parameter in approximating exp(x^2)*erfc(x) recommended by Ren and MacKenzie [2007]
    real(dp), parameter :: ap = 2.7889_dp 
    ! Physical constant
    real(dp), parameter :: pi = dacos(-1.0_dp)
    real(dp), parameter :: k_B = 1.38064852e-23_dp ! Boltzmann constant in J/K
    real(dp), parameter :: Grav = 6.67430e-11_dp ! Gravitational constant in m^3 kg^-1 s^-2

    ! radiative transfer ------------------------
    var%tau_EUV = 0.0_dp 
    var%I_EUV = 0.0_dp 

    do isp = 1, spl%nsp

      if (xct%label_sigma_a_EUV(isp) /= 1) cycle ! skip if no data

      do iz = 1, grd%nz
        g = Grav * cst%Mplanet / (cst%R + grd%alt(iz))**2
        chiz = grd%sza(grd%ix,grd%iy)
        Hz   = k_B * var%Tn(iz) / var%m(isp) / g
        Xz   = (cst%R + grd%alt(iz)) / Hz
        yz   = dsqrt(0.5_dp * Xz) * dabs(dcos(chiz))

        ! exp(x^2)*erfc(x) is approximated by using the formula of Ren & MacKenzie [2007]
        if (chiz <= pi / 2.0_dp) then 
          Chfunc = dsqrt(pi/2.0_dp*Xz) &
            &        * ap / ((ap-1.0_dp)*dsqrt(pi*yz*yz) + dsqrt(pi*yz*yz + ap*ap))
        else if (chiz > pi / 2.0_dp) then 
          Chfunc = dsqrt(2.0_dp*pi*Xz) &
            &        * ( dsqrt(dsin(chiz)) * dexp( Xz*(1.0_dp - dsin(chiz)) ) &
            &          - 0.5_dp * ap / ((ap-1.0_dp)*dsqrt(pi*yz*yz) + dsqrt(pi*yz*yz + ap*ap)) )
        end if

        ! upper limit of Chfunc is 10^10 in order not to cause infinity tau
        if (Chfunc > 1.0e10_dp) then
          Chfunc = 1.0e10_dp
        end if

        cln_Ch = var%clm_ni(iz,isp) * Chfunc
        var%tau_EUV(1:flx%nwl_EUV,iz) = var%tau_EUV(1:flx%nwl_EUV,iz) &
          &             + cln_Ch * xct%sigma_a_EUV(1:flx%nwl_EUV,isp)
      end do ! iz

    end do ! isp

    do iz  = 1, grd%nz
      var%I_EUV(1:flx%nwl_EUV,iz) = flx%solar_EUV(1:flx%nwl_EUV) * dexp(-var%tau_EUV(1:flx%nwl_EUV,iz)) 
    end do

    ! photoionization rate ---------------------------
    do ich = 1, spl%nch

      if (spl%reaction_type_char(ich) /= 'photoionization') cycle ! skip if not photoionization

      do iz = 1, grd%nz
        J_rate(iz,ich) = dot_product(var%I_EUV(1:flx%nwl_EUV,iz), xct%sigma_i_EUV(1:flx%nwl_EUV,ich))
      end do

    end do ! ich

  end subroutine p__EUVAC_photoionization_rate


  !----------------------------------------------------------------------------------------------
  ! Flux, cross sections and quantum yield data are automatically adapted to the model bins.
  !   if the data bin is larger than the model bin, the data is interpolated.
  !   if the data bin is smaller than the model bin, the data is binned.
  !----------------------------------------------------------------------------------------------
  subroutine p__HEUVAC_binning(flx, idata, nl, odata, nwl)
    implicit none
    type(flx_), intent(in)  :: flx
    integer,    intent(in)  :: nl, nwl
    real(8),    intent(in)  :: idata(nl,2)
    real(8),    intent(out) :: odata(nwl)
    integer iwl, il, label, il0, ndata
    real(8) idata_tmp(nl,2)
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
          o0 = flx%lambda_EUV(1)
          dop = flx%dlambda_EUV(1) * 0.5d0
          dom = flx%dlambda_EUV(1) * 0.5d0
        else if (iwl > 1 .and. iwl < nwl) then 
          o0 = flx%lambda_EUV(iwl)
          dom = flx%dlambda_EUV(iwl) * 0.5d0
          dop = flx%dlambda_EUV(iwl) * 0.5d0
        else if (iwl == nwl) then 
          o0 = flx%lambda_EUV(iwl)
          dom = flx%dlambda_EUV(iwl) * 0.5d0
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
  
  end subroutine p__HEUVAC_binning


end module p__EUV