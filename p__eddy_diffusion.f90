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
    real(dp) K0, K1, T
    real(dp) tmp

    !--------------------------------------------------------------------------------------
    !
    !                                      Earth
    !
    !--------------------------------------------------------------------------------------
    ! var%n_tot(iz) is total neutral density in the unit of [/m^3]
    if (spl%planet == 'Earth') then
      
      do iz = 1, grd%nz
        if (grd%alt(iz) <= 1.0e4_dp) var%K_eddy(iz) = 1.0d0 ![m^2/s]
        if (grd%alt(iz) > 1.0e4_dp .and. grd%alt(iz) <= 4.0e4) var%K_eddy(iz) = 1.0d-1 ![m^2/s]
        if (grd%alt(iz) > 4.0e4_dp) var%K_eddy(iz) = 1.0d0 ![m^2/s]
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

      !! for Masunaga et al.
      !k0 = 100.0_dp ![m^2/s]
      !k1 = 400.0_dp ![m^2/s]
      !T  = 6.0 * 86400.0_dp ! sec
      !var%K_eddy = k0 
      !var%K_eddy = k0 + k1 * (1.0d0 - dcos(2.0d0*cst%pi*var%sum_time/T))/2.0d0

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

! eddy diffusion from Krasnopolsky 2007
      do iz = 1, 30
        var%K_eddy(iz) = 2.2d-1 
      enddo
      do iz = 31, 48
        var%K_eddy(iz) = 10 ** ( 0.0387_dp * real(iz) - 1.819_dp )
      enddo

      var%K_eddy(:) = 0.8_dp * var%K_eddy(:)
      
    end if
  end subroutine p__eddy_diffusion__exe


end module p__eddy_diffusion