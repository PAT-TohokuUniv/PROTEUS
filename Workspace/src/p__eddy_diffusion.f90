module p__eddy_diffusion

  use v__tdec,        only : spl_, var_, grd_, cst_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  
  implicit none
  integer(4), parameter :: sp = 4, dp = 8
  
  private
  public :: p__eddy_diffusion__exe
  
contains
  

  subroutine p__eddy_diffusion__exe(spl, cst, grd, & ! in
    &                               var            ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(var_),           intent(inout) :: var
    integer iz
    real(dp) K0, K1, T, P_loc, P_tro
    real(dp) tmp
    logical :: flag

    !--------------------------------------------------------------------------------------
    !
    !                                      Earth
    !
    !--------------------------------------------------------------------------------------
    ! var%n_tot(iz) is total neutral density in the unit of [/m^3]
    if (spl%planet == 'Earth') then
      
      ! preset eddy diffusion for Earth atmosphere
      ! Hu et al., 2012; Yoshida et al., 2024
      P_tro = 0.1e5_dp
      do iz = 1, grd%nz
        P_loc = var%n_tot(iz) * cst%k_B * var%Tn(iz)
        if (P_loc >= P_tro) then
          var%K_eddy(iz) = 10.0_dp
        else
          var%K_eddy(iz) = min(10.0_dp, 0.3_dp*dsqrt(P_tro/P_loc))
        end if
      end do

    end if

    !--------------------------------------------------------------------------------------
    !
    !                                       Mars
    !
    !--------------------------------------------------------------------------------------
    ! var%n_tot(iz) is total neutral density in the unit of [/m^3]
    if (spl%planet == 'Mars') then
      
      ! Chaffin et al., 2017 case
      do iz = 1, grd%nz
          var%K_eddy(iz) = 1.0d2 ![m^2/s]
          tmp = 2.0e9_dp / dsqrt(var%n_tot(iz)*1.0e-6_dp) ![m^2/s]
          if (tmp > 1.0e2_dp) var%K_eddy(iz) = tmp
      end do

      ! For early Mars case
      !do iz = 1, grd%nz
      !  var%K_eddy(iz) = 1.0d1 ![m^2/s]
      !  tmp = 2.0e9_dp / dsqrt(var%n_tot(iz)*1.0e-6_dp) ![m^2/s]
      !  if (tmp > 1.0e1_dp) var%K_eddy(iz) = tmp
      !  if (tmp > 5.0e1_dp) var%K_eddy(iz) = 5.0e1_dp
      !end do

    end if

    !--------------------------------------------------------------------------------------
    !
    !                                      Jupiter
    !
    !--------------------------------------------------------------------------------------
    if (spl%planet == 'Jupiter') then

      var%K_eddy = 2.0d2 ![m^2/s]
      
    end if


    !--------------------------------------------------------------------------------------
    !
    !                                      Venus
    !
    !--------------------------------------------------------------------------------------
    if (spl%planet == 'Venus') then

      ! no preset is available now.
      var%K_eddy = 2.2d-1 

    end if
  end subroutine p__eddy_diffusion__exe


end module p__eddy_diffusion