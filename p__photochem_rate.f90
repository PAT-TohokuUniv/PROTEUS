! p__photochem_rate.f90
!
!
!
!

module p__photochem_rate
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__EUVAC,       only : p__EUVAC_cross_section
  use p__UV,          only : p__UV_cross_section_exe
  use p__io,          only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                        p__io_progress

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_rate__exe

contains


  subroutine p__photochem_rate__exe(spl, cst, flx, grd, set, & ! in
    &                               xct, var                 ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(flx_),           intent(in)    :: flx
    type(set_),           intent(in)    :: set
    type(xct_),           intent(inout) :: xct
    type(var_),           intent(inout) :: var

    integer isp, jsp, ksp, ich, jch, iz, iwl, swl, ewl, icase
    integer i, j, k, ijch, ilist, ii, jj
    integer bij
    real(dp) tau_factor
    real(dp) n_M
    integer pop, label
    real(dp) stack(spl%nrpn), token
    real(dp) a0, b0, ni, nj, nk, kij, aij
    real(dp) a, b, c, d, dn(spl%nsp_i, grd%nz), T
    real(dp) kall(3), k0, kinf, k2, k3
    real(dp) tmp, tmp1, tmp2, tmpzarr(grd%nz)

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    ! reaction rate -------------------------------------------------------------------------------
      var%ki(1:spl%nch,1:grd%nz) = 0.0_dp

      ! photoionization
      xct%type        = 'ionization'
      !xct%sigma_i_EUV = 0.0_dp
      !call p__EUVAC_cross_section(spl, & ! in
      !      &                     xct  ) ! inout
      do ich = 1, spl%nch
        if (spl%reaction_type_char(ich) == 'photoionization') then
          do iz = 1, grd%nz
            do iwl = 1, flx%nwl_EUV
              var%ki(ich,iz) = var%ki(ich,iz) &
                &            + var%I_EUV(iwl,iz) * xct%sigma_i_EUV(iwl,ich)
            end do
          end do
        end if
      end do

      ! photodissociation
      xct%type       = 'dissociation'
      !xct%sigma_d_UV = 0.0_dp
      !call p__UV_cross_section_exe(var%Tn, grd%nz, spl, flx, var, & ! in
      !      &                      xct                            ) ! inout
      do ich = 1, spl%nch
        if (spl%reaction_type_char(ich) == 'photodissociation') then
          swl = xct%wlrange_d_UV(1,ich)
          ewl = xct%wlrange_d_UV(2,ich)
          do iz = 1, grd%nz
            do iwl = swl, ewl
              if (xct%sigma_d_UV(iwl,iz,ich) > 0.0_dp) then 
                var%ki(ich,iz) = var%ki(ich,iz) &
                  &            + xct%sigma_d_UV(iwl,iz,ich) * var%I_UV(iwl,iz) 
              end if
            end do
          end do
        end if
      end do

      ! chemical reaction rate
      do iz = 1, grd%nz
        do ich = 1, spl%nch

          ! Temperature range (not for presssure dependent 3 body reactions)
          if (spl%reaction_type_list(ich) == 0) then
              do icase = 1, spl%rate_cases(ich)

                if (      spl%T_range(ich,icase,1) > 0.0_dp - 0.1_dp &
                  & .and. spl%T_range(ich,icase,1) < 0.0_dp + 0.1_dp) then
                  T = var%Tn(iz)
                end if
                if (      spl%T_range(ich,icase,1) > 1.0_dp - 0.1_dp &
                  & .and. spl%T_range(ich,icase,1) < 1.0_dp + 0.1_dp) then
                  T = var%Tn(iz)
                end if
                if (      spl%T_range(ich,icase,1) > 2.0_dp - 0.1_dp &
                  & .and. spl%T_range(ich,icase,1) < 2.0_dp + 0.1_dp) then
                  T = var%Te(iz)
                end if
                if (      spl%T_range(ich,icase,1) > 3.0_dp - 0.1_dp &
                  & .and. spl%T_range(ich,icase,1) < 3.0_dp + 0.1_dp) then
                  T = var%Ti(iz)
                end if

                ! reaction rate calculation using Reversed Polish Notation
                if (T >= spl%T_range(ich,icase,2) .and. T < spl%T_range(ich,icase,3)) then
                  pop = 0
                  do i = 2, nint(spl%rate_rpn_token(ich,icase,1))+1
                    token = spl%rate_rpn_token(ich,icase,i)
                    label = spl%rate_rpn_label(ich,icase,i)
                    ! values
                    if (label == 0) then
                      pop = pop + 1
                      stack(pop) = dble(token)
                    ! operators
                    else if (label == 1) then
                      if (nint(token) == 1) then
                        stack(pop-1) = stack(pop-1) + stack(pop)
                        stack(pop) = 0.0_dp
                        pop = pop - 1
                      else if (nint(token) == 2) then
                        stack(pop-1) = stack(pop-1) - stack(pop)
                        stack(pop) = 0.0_dp
                        pop = pop - 1
                      else if (nint(token) == 3) then
                        stack(pop-1) = stack(pop-1) * stack(pop)
                        stack(pop) = 0.0_dp
                        pop = pop - 1
                      else if (nint(token) == 4) then
                        stack(pop-1) = stack(pop-1) / stack(pop)
                        stack(pop) = 0.0_dp
                        pop = pop - 1
                      else if (nint(token) == 6) then
                        stack(pop-1) = stack(pop-1) ** stack(pop)
                        stack(pop) = 0.0_dp
                        pop = pop - 1
                      else if (nint(token) == 7) then
                        stack(pop) = dexp(stack(pop))
                      else if (nint(token) == 8) then
                        stack(pop) = dsqrt(stack(pop))
                      end if
                    ! temperatures
                    else if (label == 2) then
                      pop = pop + 1
                      if (nint(token) == 1) then
                        stack(pop) = var%Tn(iz)
                      else if (nint(token) == 2) then
                        stack(pop) = var%Ti(iz)
                      else if (nint(token) == 3) then
                        stack(pop) = var%Te(iz)
                      end if
                    ! densities
                    else if (label == 3) then
                      pop = pop + 1
                      if (nint(token) == 1) then
                        stack(pop) = var%n_tot(iz) 
                      end if
                    ! altitude dependent
                    else if (label == 4) then
                      pop = pop + 1
                      if (nint(token) == 1) then
                        stack(pop) = grd%alt(iz)/1000.0_dp 
                      end if
                    end if
                  end do
                end if

                var%ki(ich,iz) = stack(1)
                if (spl%reactant_list(ich,1) == 2) then
                  var%ki(ich,iz) = var%ki(ich,iz) * 1.0e-6_dp
                else if (spl%reactant_list(ich,1) == 3) then
                  var%ki(ich,iz) = var%ki(ich,iz) * 1.0e-12_dp
                end if

              end do

          end if

          ! pressure dependent three-body reaction
          if (     spl%reaction_type_char(ich) == 'pressure_dependent_3body' &
            & .or. spl%reaction_type_char(ich) == 'pressure_dependent_3bodyM') then
            do icase = 1, 2

              if (      spl%T_range(ich,icase,1) > 0.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 0.0_dp + 0.1_dp) then
                T = var%Tn(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 1.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 1.0_dp + 0.1_dp) then
                T = var%Tn(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 2.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 2.0_dp + 0.1_dp) then
                T = var%Te(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 3.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 3.0_dp + 0.1_dp) then
                T = var%Ti(iz)
              end if

              ! reaction rate calculation using Reversed Polish Notation
              pop = 0
              do i = 2, nint(spl%rate_rpn_token(ich,icase,1))+1
                token = spl%rate_rpn_token(ich,icase,i)
                label = spl%rate_rpn_label(ich,icase,i)
                if (label == 0) then
                  pop = pop + 1
                  stack(pop) = dble(token)
                else if (label == 1) then
                  if (nint(token) == 1) then
                    stack(pop-1) = stack(pop-1) + stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 2) then
                    stack(pop-1) = stack(pop-1) - stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 3) then
                    stack(pop-1) = stack(pop-1) * stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 4) then
                    stack(pop-1) = stack(pop-1) / stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 6) then
                    stack(pop-1) = stack(pop-1) ** stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 7) then
                    stack(pop) = dexp(stack(pop))
                  else if (nint(token) == 8) then
                    stack(pop) = dsqrt(stack(pop))
                  end if
                else if (label == 2) then
                  pop = pop + 1
                  if (nint(token) == 1) then
                    stack(pop) = var%Tn(iz)
                  else if (nint(token) == 2) then
                    stack(pop) = var%Ti(iz)
                  else if (nint(token) == 3) then
                    stack(pop) = var%Te(iz)
                  end if
                end if
              end do

              kall(icase) = stack(1)

            end do

            !var%ki(ich,iz) = k0 / ( 1_dp + k0/kinf * n_M )

            ! Chaffin et al., 2017 supplement document
            if (     spl%reaction_type_char(ich) == 'pressure_dependent_3body' ) then
              k0   = kall(1) * 1.0e-12_dp
              kinf = kall(2) * 1.0e-6_dp
              n_M = var%n_tot(iz)
              var%ki(ich,iz) = k0 / ( 1.0_dp + k0*n_M/kinf ) &
                &             * 0.6_dp ** ( 1.0_dp / (1.0_dp + (dlog10(k0*n_M/kinf) )**2.0_dp) )
            else if (spl%reaction_type_char(ich) == 'pressure_dependent_3bodyM' ) then
              k0   = kall(1) * 1.0e-6_dp
              kinf = kall(2)
              n_M = var%n_tot(iz)
              var%ki(ich,iz) = k0 / n_M / ( 1.0_dp + k0*n_M/kinf ) &
                &             * 0.6_dp ** ( 1.0_dp / (1.0_dp + (dlog10(k0*n_M/kinf) )**2.0_dp) )
            end if

          end if

          ! pressure dependent three-body reaction: Lindemann-Hinshelwood expression
          if (     spl%reaction_type_char(ich) == 'Lindemann-Hinshelwood' ) then
            do icase = 1, 3

              if (      spl%T_range(ich,icase,1) > 0.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 0.0_dp + 0.1_dp) then
                T = var%Tn(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 1.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 1.0_dp + 0.1_dp) then
                T = var%Tn(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 2.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 2.0_dp + 0.1_dp) then
                T = var%Te(iz)
              end if
              if (      spl%T_range(ich,icase,1) > 3.0_dp - 0.1_dp &
                & .and. spl%T_range(ich,icase,1) < 3.0_dp + 0.1_dp) then
                T = var%Ti(iz)
              end if

              ! reaction rate calculation using Reversed Polish Notation
              pop = 0
              do i = 2, nint(spl%rate_rpn_token(ich,icase,1))+1
                token = spl%rate_rpn_token(ich,icase,i)
                label = spl%rate_rpn_label(ich,icase,i)
                if (label == 0) then
                  pop = pop + 1
                  stack(pop) = dble(token)
                else if (label == 1) then
                  if (nint(token) == 1) then
                    stack(pop-1) = stack(pop-1) + stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 2) then
                    stack(pop-1) = stack(pop-1) - stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 3) then
                    stack(pop-1) = stack(pop-1) * stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 4) then
                    stack(pop-1) = stack(pop-1) / stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 6) then
                    stack(pop-1) = stack(pop-1) ** stack(pop)
                    stack(pop) = 0.0_dp
                    pop = pop - 1
                  else if (nint(token) == 7) then
                    stack(pop) = dexp(stack(pop))
                  else if (nint(token) == 8) then
                    stack(pop) = dsqrt(stack(pop))
                  end if
                else if (label == 2) then
                  pop = pop + 1
                  if (nint(token) == 1) then
                    stack(pop) = var%Tn(iz)
                  else if (nint(token) == 2) then
                    stack(pop) = var%Ti(iz)
                  else if (nint(token) == 3) then
                    stack(pop) = var%Te(iz)
                  end if
                end if
              end do

              kall(icase) = stack(1)

            end do

            ! JPL document
            k0 = kall(1) * 1.0e-6_dp
            k2 = kall(2) * 1.0e-6_dp
            k3 = kall(3) * 1.0e-12_dp
            n_M = var%n_tot(iz)
            var%ki(ich,iz) = k0 + k3*n_M / (1.0_dp + k3*n_M/k2)
          end if

        end do

      end do ! z

      ! Special reaction rate case!
      ! if you input production or loss rate explicitly, please input as var%ki_special in p__planet__exe in p__planet.f90
      do iz = 1, grd%nz
        var%ni(0,iz) = 1.0_dp ! if there are no reactant like meteoroid ablation, P = k * ni(0,iz) = k
        do jch = 1, var%nspecial
          ich = var%ich_special(jch)
          if (    spl%reaction_type_char(ich) == 'electron impact' &
          &  .or. spl%reaction_type_char(ich) == 'proton impact' &
          &  .or. spl%reaction_type_char(ich) == 'H impact'  &
          &  .or. spl%reaction_type_char(ich) == 'Meteoroid ablation' &
          &  .or. spl%reaction_type_char(ich) == 'Rainout' ) then
            var%ki(ich,iz) = var%ki_special(jch,grd%ix,grd%iy,iz)
          end if
        end do
      end do

    ! Production rate
      do iz = 1, grd%nz
        do isp = 1, spl%nsp_i
          i = isp
          var%Pi(i,iz) = 0.0_dp
          do ich = 1, spl%Prod_list(i,1)
            a = dble(spl%Prod_list(i,2*ich+1))
            jch = spl%Prod_list(i,2*ich)
            tmp = a * var%ki(jch,iz)
            do jsp = 1, spl%reactant_list(jch,1)
              ksp = spl%reactant_list(jch,jsp+1)
              tmp = tmp * var%ni(ksp,iz)
            end do
            var%Pi(i,iz) = var%Pi(i,iz) + tmp
            var%Pij(i,iz,jch) = tmp
          end do
        end do
      end do ! z

    ! Loss rate
      do iz = 1, grd%nz
        do isp = 1, spl%nsp_i
          i = isp
          var%Li(i,iz) = 0.0_dp
          do ich = 1, spl%Loss_list(i,1)
            a = dble(spl%Loss_list(i,2*ich+1))
            jch = spl%Loss_list(i,2*ich)
            tmp = a * var%ki(jch,iz)
            do jsp = 1, spl%reactant_list(jch,1)
              ksp = spl%reactant_list(jch,jsp+1)
              tmp = tmp * var%ni(ksp,iz)
            end do
            var%Li(i,iz) = var%Li(i,iz) + tmp
          end do
        end do
      end do ! z

      do iz = 1, grd%nz
        do ich = 1, spl%nch
          tmp = var%ki(ich,iz)
          do jsp = 1, spl%reactant_list(ich,1)
            ksp = spl%reactant_list(ich,jsp+1)
            tmp = tmp * var%ni(ksp,iz)
          end do
          var%rate(ich,iz) = tmp
        end do
      end do ! z


  end subroutine p__photochem_rate__exe



end module p__photochem_rate
