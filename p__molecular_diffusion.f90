module p__molecular_diffusion
  use v__tdec,        only : spl_, var_, grd_, cst_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  
  implicit none
  integer(4), parameter :: sp = 4, dp = 8
  
  private
  public :: p__molecular_diffusion__exe
  
contains
  

  subroutine p__molecular_diffusion__exe(spl, cst, grd, & ! in
    &                                    var            ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(var_),           intent(inout) :: var
    integer iz, isp
    real(dp) Tn(grd%nz), Ti(grd%nz), Te(grd%nz)
    real(dp) nu(spl%nsp, grd%nz)
    real(dp) D_max, sigma_atom
    real(dp) tmp, tmp1, tmp2

    !--------------------------------------------------------------------------------------
    !
    !                                       Earth
    !
    !--------------------------------------------------------------------------------------
    ! var%n_tot(iz) is total neutral density in the unit of [/m^3]
    if (spl%planet == 'Earth') then
      do isp = 1, spl%nsp
        if ( spl%label_fix(isp) == 0 .and. spl%species(isp) /= 'M' ) then

          do iz = 1, grd%nz

            Tn(iz) = var%Tn(iz)! neutral temperature
            Ti(iz) = var%Ti(iz)! ion temperature
            Te(iz) = var%Te(iz)! electron temperature

            ! all species other than H, H2
            var%D_mol(isp,iz) = 1.0e17_dp*Tn(iz)**(0.75_dp)/(var%n_tot(iz)/1.0e6_dp)! [cm^2/s]
            var%D_mol(isp,iz) = var%D_mol(isp,iz)*1.0e-4_dp ! [cm^2/s] -> [m^2/s]

          end do

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
      do isp = 1, spl%nsp
        if ( spl%label_fix(isp) == 0 .and. spl%species(isp) /= 'M' ) then

          do iz = 1, grd%nz

            Tn(iz) = var%Tn(iz)! neutral temperature
            Ti(iz) = var%Ti(iz)! ion temperature
            Te(iz) = var%Te(iz)! electron temperature

            ! all species other than H, H2
            var%D_mol(isp,iz) = 1.0e17_dp*Tn(iz)**(0.75_dp)/(var%n_tot(iz)/1.0e6_dp)! [cm^2/s]
            var%D_mol(isp,iz) = var%D_mol(isp,iz)*1.0e-4_dp ! [cm^2/s] -> [m^2/s]

            ! H
            if (spl%species(isp)=='H') then
              var%D_mol(isp,iz) = 8.4_dp*1.0e17_dp*Tn(iz)**(0.597_dp)/(var%n_tot(iz)/1.0e6_dp)! [cm^2/s]
              var%D_mol(isp,iz) = var%D_mol(isp,iz)*1.0e-4_dp ! [cm^2/s] -> [m^2/s]

            ! H2
            else if (spl%species(isp)=='H2') then 
              var%D_mol(isp,iz) = 2.23_dp*1.0e17_dp*Tn(iz)**(0.75_dp)/(var%n_tot(iz)/1.0e6_dp)! [cm^2/s]
              var%D_mol(isp,iz) = var%D_mol(isp,iz)*1.0e-4_dp ! [cm^2/s] -> [m^2/s]

            end if 
            
          end do

        end if
      end do
    end if

    !--------------------------------------------------------------------------------------
    !
    !                                      Jupiter
    !
    !--------------------------------------------------------------------------------------
    if (spl%planet == 'Jupiter') then
      sigma_atom = 2.7d-10 ! H2 molecule radius
      do isp = 1, spl%nsp
        if (  spl%species(isp) /= 'M' ) then
          do iz = 1, grd%nz

            Tn(iz) = var%Tn(iz)
            Ti(iz) = var%Ti(iz)

            ! Collision Frequency
            tmp1 = 2.0_dp * var%n_tot(iz) * sigma_atom**2.0_dp
            tmp2 = 2.0_dp * cst%pi * cst%k_B * Tn(iz) / var%m_mean(iz)
            nu(isp,iz) = tmp1 * ( tmp2 * ( 1.0_dp  + var%m_mean(iz) / var%m(isp) ) )**0.5_dp

            ! Diffusion Coefficient limit
            D_max = 1.0d18
            tmp = cst%k_B * Ti(iz)
            var%D_mol(isp,iz) = tmp / var%m(isp) / nu(isp,iz)
            if ( var%D_mol(isp,iz) > D_max ) then
              var%D_mol(isp,iz) = D_max
            end if

          end do

        end if
      end do
    end if



  end subroutine p__molecular_diffusion__exe


end module p__molecular_diffusion