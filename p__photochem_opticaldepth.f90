! p__photochem_opticaldepth.f90
!
!
!
!

module p__photochem_opticaldepth
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  use p__EUVAC,       only : p__EUVAC_cross_section
  use p__UV,          only : p__UV_cross_section_exe, p__UV_EUV_xct_convergence_exe

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_opticaldepth__exe

contains


  subroutine p__photochem_opticaldepth__exe(spl, cst, flx, grd, set, & ! in
    &                                       xct, var                 ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(flx_),           intent(in)    :: flx
    type(set_),           intent(in)    :: set
    type(xct_),           intent(inout) :: xct
    type(var_),           intent(inout) :: var

    integer isp, jsp, iz, iwl, swl, ewl
    real(dp) tmp1, tmp2

    ! for solar zenith angle near and greater than 90deg [Smith et al., 1972]
    real(dp) Hz, yz, Xz, chiz, Chfunc(grd%nz), cln_Ch
    ! ap : parameter in approximating exp(x^2)*erfc(x) recommended by Ren and MacKenzie [2007]
    real(dp), parameter :: ap = 2.7889_dp 

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    ! column density -------------------------------------------------------------------------------
      var%clm_ni = 0.0_dp
      do isp = 1, spl%nsp
        var%clm_ni(isp,grd%nz) = var%ni(isp,grd%nz)*grd%dalt(grd%nz)
        do iz = grd%nz, 2, -1
          var%clm_ni(isp,iz-1) = var%clm_ni(isp,iz) &
            &                   + (var%ni(isp,iz-1) + var%ni(isp,iz)) &
            &                   / 2.0_dp * grd%dalt(iz-1)
        end do
      end do

    ! 'M' density : sum of all ---------------------------------------------------------
      isp = sp_index(spl, 'M')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        var%ni(isp,:) = 0.0_dp
        do jsp = 1, spl%nsp
          if (isp /= jsp) then
            do iz = 1, grd%nz
              var%ni(isp,iz) &
                &  = var%ni(isp,iz) + var%ni(jsp,iz)
            end do 
          end if
        end do
      end if

    ! total density, mean mass ---------------------------------------------------------
      do iz = 1, grd%nz
        tmp1 = 0.0_dp
        tmp2 = 0.0_dp
        do isp = 1, spl%nsp
          if (spl%species(isp) /= 'M' .and. spl%species(isp) /= 'products' .and. spl%species(isp) /= 'hv') then
            tmp1 = tmp1 + var%ni(isp,iz) * var%m(isp)
            tmp2 = tmp2 + var%ni(isp,iz)
          end if
        end do
        var%m_mean(iz) = tmp1 / tmp2
        var%n_tot(iz) = tmp2
      end do

    ! Ch*(X,chi) function (Chfunc) in Smith et al., 1972: treatment of solar zenith angle near terminator ---

      do iz = 1, grd%nz
        chiz = grd%sza(grd%ix, grd%iy)
        Hz   = cst%k_B * var%Tn(iz) / var%m_mean(iz) / cst%g
        Xz   = (cst%R + grd%alt(iz)) / Hz
        yz   = dsqrt(0.5_dp * Xz) * dabs(dcos(chiz))
 
        ! exp(x^2)*erfc(x) is approximated by using the formula of Ren & MacKenzie [2007]
        if (chiz <= cst%pi / 2.0_dp) then 
          Chfunc(iz) = dsqrt(cst%pi/2.0_dp*Xz) &
            &        * ap / ((ap-1.0_dp)*dsqrt(cst%pi*yz*yz) + dsqrt(cst%pi*yz*yz + ap*ap))
        else if (chiz > cst%pi / 2.0_dp) then 
          Chfunc(iz) = dsqrt(2.0_dp*cst%pi*Xz) &
            &        * ( dsqrt(dsin(chiz)) * dexp( Xz*(1.0_dp - dsin(chiz)) ) &
            &          - 0.5_dp * ap / ((ap-1.0_dp)*dsqrt(cst%pi*yz*yz) + dsqrt(cst%pi*yz*yz + ap*ap)) )
        end if

        ! upper limit of Chfunc is 10^10 in order not to cause infinity tau
        if (Chfunc(iz) > 1.0e10_dp) then
          Chfunc(iz) = 1.0e10_dp
        end if
        
      end do

    ! solar flux at each altitude ------------------------------------------------------------------
      var%tau_EUV = 0.0_dp
      var%tau_EUV_subsolar = 0.0_dp
      var%tau_UV  = 0.0_dp
      var%tau_UV_subsolar = 0.0_dp
      var%tau_RS = 0.0_dp
      xct%type       = 'absorption'
      !xct%sigma_a_EUV = 0.0_dp
      !xct%sigma_a_UV  = 0.0_dp

      !call p__EUVAC_cross_section(spl, & ! in
      !  &                         xct  ) ! inout
      !call p__UV_cross_section_exe(var%Tn, grd%nz, spl, flx, var, & ! in
      !  &                          xct                            ) ! inout

      do isp = 1, spl%nsp

        if (xct%label_sigma_a_EUV(isp) == 1) then 
          do iz = 1, grd%nz
            cln_Ch = var%clm_ni(isp,iz) * Chfunc(iz)
            do iwl = 1, flx%nwl_EUV
              var%tau_EUV(iwl,iz) = var%tau_EUV(iwl,iz) &
                &                 + cln_Ch * xct%sigma_a_EUV(iwl,isp) 
              var%tau_EUV_subsolar(iwl,iz) = var%tau_EUV_subsolar(iwl,iz) &
                &                          + var%clm_ni(isp,iz) * xct%sigma_a_EUV(iwl,isp) 
              if ( var%tau_EUV(iwl,iz) > 100.0_dp ) then
                var%tau_EUV(iwl,iz) = 100.0_dp
              end if
            end do
          end do
        end if

        if (xct%label_sigma_a_UV(isp) == 1) then 
          swl = xct%wlrange_a_UV(1,isp)
          ewl = xct%wlrange_a_UV(2,isp)
          do iz = 1, grd%nz
            cln_Ch = var%clm_ni(isp,iz) * Chfunc(iz)
            do iwl = swl, ewl
              if (xct%sigma_a_UV(iwl,iz,isp) > 0.0_dp) then 
                var%tau_UV(iwl,iz) = var%tau_UV(iwl,iz) &
                  &                + cln_Ch * xct%sigma_a_UV(iwl,iz,isp)
              end if
            end do
          end do
        end if

        !if (xct%label_sigma_a_RS(isp) == 1) then 
        !  do iz = 1, grd%nz
        !    cln_Ch = var%clm_ni(isp,iz) * Chfunc(iz)
        !    do iwl = 1, flx%nwl_UV
        !      var%tau_RS(iwl,iz) = var%tau_RS(iwl,iz) &
        !        &                 + cln_Ch * xct%sigma_a_RS(iwl,isp) 
        !    end do
        !  end do
        !end if

      end do ! isp

      ! cross section convergence
      !   EUVAC cross sections are added to UV cross sections
      !call p__UV_EUV_xct_convergence_exe(grd%nz, spl, var, flx, & ! in
      !  &                                xct                    ) ! inout
      do isp = 1, spl%nsp
        if (xct%label_sigma_a_UV_EUV_conv(isp) == 1) then 
          do iz = 1, grd%nz
            cln_Ch = var%clm_ni(isp,iz) * Chfunc(iz)
            do iwl = 1, flx%nwl_UV
              var%tau_UV(iwl,iz) = var%tau_UV(iwl,iz) &
                &                + cln_Ch * xct%sigma_a_UV_EUV(iwl,isp) 
              if ( var%tau_UV(iwl,iz) > 100.0_dp ) then
                var%tau_UV(iwl,iz) = 100.0_dp
              end if
            end do
          end do
        end if
      end do

      do iz  = 1, grd%nz
        do iwl = 1, flx%nwl_EUV
          var%I_EUV(iwl,iz) = flx%solar_EUV(iwl) &
            &               * dexp( -var%tau_EUV(iwl,iz) ) * flx%mode_factor !& 
            ! ionization by soft electron?, photoelectron? : 1% of subsolar value
            !&               + flx%solar_EUV(iwl) * 0.01_dp &
            !&               * dexp( -var%tau_EUV_subsolar(iwl,iz) ) * flx%mode_factor
        end do
      end do

      do iz  = 1, grd%nz
        do iwl = 1, flx%nwl_UV
          var%I_UV(iwl,iz) = flx%solar_UV(iwl) * flx%dlambda_UV(iwl) &
          &                * dexp( -var%tau_UV(iwl,iz) -var%tau_RS(iwl,iz)) * flx%mode_factor 
        end do
      end do
      !stop

  end subroutine p__photochem_opticaldepth__exe


end module p__photochem_opticaldepth
