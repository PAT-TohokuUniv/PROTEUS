module p__EUVAC
  
  use v__tdec,    only : var_, grd_, cst_, xct_, spl_, flx_, set_
  use c__prm,     only : c__prm__ini
  use p__search,  only : p__search_reactant, p__search_product, ch_identify, sp_index


  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__EUVAC_flux, p__EUVAC_cross_section

contains


  !-------------------------------------------------
  ! solar EUV flux model: EUVAC
  !-------------------------------------------------
  subroutine p__EUVAC_flux(spl, grd, set, & ! in
    &                      var, xct, flx  ) ! inout
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer iwl
    real(dp) F107A

    flx%nwl_EUV = 37
    allocate(flx%lambda_EUV(flx%nwl_EUV),flx%dlambda_EUV(flx%nwl_EUV))
    allocate(flx%F74113(flx%nwl_EUV),flx%Ai(flx%nwl_EUV),flx%solar_EUV(flx%nwl_EUV))
    allocate(var%tau_EUV(flx%nwl_EUV,grd%nz))
    allocate(var%I_EUV(flx%nwl_EUV,grd%nz))
    allocate(xct%sigma_a_EUV(flx%nwl_EUV,spl%nsp), xct%sigma_i_EUV(flx%nwl_EUV,spl%nch))
    allocate(xct%label_sigma_a_EUV(spl%nsp))

    flx%multiplying_factor_EUV = set%euv_factor
    flx%F107 = set%F107 ! F10.7 at Earth
    F107A = flx%F107

    ! [nm]
    flx%lambda_EUV(1)  =    75.0d-1
    flx%lambda_EUV(2)  =   125.0d-1
    flx%lambda_EUV(3)  =   175.0d-1
    flx%lambda_EUV(4)  =   225.0d-1
    flx%lambda_EUV(5)  =   256.32d-1
    flx%lambda_EUV(6)  =   284.15d-1
    flx%lambda_EUV(7)  =   275.0d-1
    flx%lambda_EUV(8)  =   303.31d-1
    flx%lambda_EUV(9)  =   303.78d-1
    flx%lambda_EUV(10) =   325.0d-1
    flx%lambda_EUV(11) =   368.07d-1
    flx%lambda_EUV(12) =   375.0d-1
    flx%lambda_EUV(13) =   425.0d-1
    flx%lambda_EUV(14) =   465.22d-1
    flx%lambda_EUV(15) =   475.0d-1
    flx%lambda_EUV(16) =   525.0d-1
    flx%lambda_EUV(17) =   554.37d-1
    flx%lambda_EUV(18) =   584.33d-1
    flx%lambda_EUV(19) =   575.0d-1
    flx%lambda_EUV(20) =   609.76d-1
    flx%lambda_EUV(21) =   629.73d-1
    flx%lambda_EUV(22) =   625.0d-1
    flx%lambda_EUV(23) =   650.0d-1
    flx%lambda_EUV(24) =   703.36d-1
    flx%lambda_EUV(25) =   725.0d-1
    flx%lambda_EUV(26) =   765.15d-1
    flx%lambda_EUV(27) =   770.41d-1
    flx%lambda_EUV(28) =   789.36d-1
    flx%lambda_EUV(29) =   775.0d-1
    flx%lambda_EUV(30) =   825.0d-1
    flx%lambda_EUV(31) =   875.0d-1
    flx%lambda_EUV(32) =   925.0d-1
    flx%lambda_EUV(33) =   977.02d-1
    flx%lambda_EUV(34) =   975.0d-1
    flx%lambda_EUV(35) =  1025.72d-1
    flx%lambda_EUV(36) =  1031.91d-1
    flx%lambda_EUV(37) =  1025.0d-1
    flx%dlambda_EUV(:) =  1.0_dp

    !--------------------------------------------------------------- F74113 -----
     flx%F74113(1)  = 1.2d0  !/ 3.0_dp
     flx%F74113(2)  = 4.5d-1 !/ 3.0_dp
     flx%F74113(3)  = 4.8d0  !/ 2.0_dp
     flx%F74113(4)  = 3.1d0  !/ 2.0_dp
     flx%F74113(5)  = 4.6d-1
     flx%F74113(6)  = 2.1d-1
     flx%F74113(7)  = 1.679d0
     flx%F74113(8)  = 8.0d-1
     flx%F74113(9)  = 6.9d0
     flx%F74113(10) = 9.65d-1
     flx%F74113(11) = 6.5d-1
     flx%F74113(12) = 3.14d-1
     flx%F74113(13) = 3.83d-1
     flx%F74113(14) = 2.9d-1
     flx%F74113(15) = 2.85d-1
     flx%F74113(16) = 4.52d-1
     flx%F74113(17) = 7.20d-1
     flx%F74113(18) = 1.27d0
     flx%F74113(19) = 3.57d-1
     flx%F74113(20) = 5.3d-1
     flx%F74113(21) = 1.59d0
     flx%F74113(22) = 3.42d-1
     flx%F74113(23) = 2.3d-1
     flx%F74113(24) = 3.6d-1
     flx%F74113(25) = 1.41d-1
     flx%F74113(26) = 1.7d-1
     flx%F74113(27) = 2.6d-1
     flx%F74113(28) = 7.02d-1
     flx%F74113(29) = 7.58d-1
     flx%F74113(30) = 1.625d0
     flx%F74113(31) = 3.537d0
     flx%F74113(32) = 3.0d0
     flx%F74113(33) = 4.4d0
     flx%F74113(34) = 1.475d0
     flx%F74113(35) = 3.5d0
     flx%F74113(36) = 2.1d0
     flx%F74113(37) = 2.467d0
    !----------------------------------------------------------------------------

    !------------------------------------------------------------------- Ai -----
     flx%Ai(1)  = 1.0017d-2
     flx%Ai(2)  = 7.125d-3
     flx%Ai(3)  = 1.3375d-2
     flx%Ai(4)  = 1.945d-2
     flx%Ai(5)  = 2.775d-3
     flx%Ai(6)  = 1.3768d-1
     flx%Ai(7)  = 2.6467d-2
     flx%Ai(8)  = 2.5d-2
     flx%Ai(9)  = 3.3333d-3
     flx%Ai(10) = 2.245d-2
     flx%Ai(11) = 6.5917d-3
     flx%Ai(12) = 3.6542d-2
     flx%Ai(13) = 7.4083d-3
     flx%Ai(14) = 7.4917d-3
     flx%Ai(15) = 2.0225d-2
     flx%Ai(16) = 8.7583d-3
     flx%Ai(17) = 3.2667d-3
     flx%Ai(18) = 5.1583d-3
     flx%Ai(19) = 3.6583d-3
     flx%Ai(20) = 1.6175d-2
     flx%Ai(21) = 3.325d-3
     flx%Ai(22) = 1.18d-2
     flx%Ai(23) = 4.2667d-3
     flx%Ai(24) = 3.0417d-3
     flx%Ai(25) = 4.75d-3
     flx%Ai(26) = 3.85d-3
     flx%Ai(27) = 1.2808d-2
     flx%Ai(28) = 3.2750d-3
     flx%Ai(29) = 4.7667d-3
     flx%Ai(30) = 4.8167d-3
     flx%Ai(31) = 5.675d-3
     flx%Ai(32) = 4.9833d-3
     flx%Ai(33) = 3.9417d-3
     flx%Ai(34) = 4.4167d-3
     flx%Ai(35) = 5.1833d-3
     flx%Ai(36) = 5.2833d-3
     flx%Ai(37) = 4.3750d-3
    !----------------------------------------------------------------------------

    !------------------------------------------------------- flux of photon -----
    do iwl = 1, flx%nwl_EUV
      flx%solar_EUV(iwl) = flx%F74113(iwl) * (1.0_dp + flx%Ai(iwl) * (0.5_dp*(flx%F107 + F107A) - 80.0_dp)) * 1.0e+13_dp &![m^-2 s^-1]
        &                /  (flx%orbit * flx%orbit) ! ........................ for orbit of Planet
      flx%solar_EUV(iwl) = flx%multiplying_factor_EUV * flx%solar_EUV(iwl)
    end do
    !----------------------------------------------------------------------------

     return

   end subroutine p__EUVAC_flux



  !-------------------------------------------------
  ! cross sections: EUVAC
  !-------------------------------------------------
  subroutine p__EUVAC_cross_section(spl, & ! in
    &                               xct  ) ! out
    type(spl_),           intent(in)    :: spl
    type(xct_),           intent(inout) :: xct
    character(len=256) reactants(10),  products(10)
    integer iwl, isp, ich, i

    !----------------------------------------------------------------------------------------------------
    ! Photoabsorption cross sections [m^2]

    if (xct%type == 'absorption') then

      xct%label_sigma_a_EUV = 0

      do isp = 1, spl%nsp

        if (spl%species(isp) == 'O2') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  O2  ---------
          xct%sigma_a_EUV(1,isp)  =  1.316d-22
          xct%sigma_a_EUV(2,isp)  =  3.806d-22
          xct%sigma_a_EUV(3,isp)  =  7.509d-22
          xct%sigma_a_EUV(4,isp)  = 10.900d-22
          xct%sigma_a_EUV(5,isp)  = 13.370d-22
          xct%sigma_a_EUV(6,isp)  = 15.790d-22
          xct%sigma_a_EUV(7,isp)  = 14.387d-22
          xct%sigma_a_EUV(8,isp)  = 16.800d-22
          xct%sigma_a_EUV(9,isp)  = 16.810d-22
          xct%sigma_a_EUV(10,isp) = 17.438d-22
          xct%sigma_a_EUV(11,isp) = 18.320d-22
          xct%sigma_a_EUV(12,isp) = 18.118d-22
          xct%sigma_a_EUV(13,isp) = 20.310d-22
          xct%sigma_a_EUV(14,isp) = 21.910d-22
          xct%sigma_a_EUV(15,isp) = 23.101d-22
          xct%sigma_a_EUV(16,isp) = 24.606d-22
          xct%sigma_a_EUV(17,isp) = 26.040d-22
          xct%sigma_a_EUV(18,isp) = 22.720d-22
          xct%sigma_a_EUV(19,isp) = 26.610d-22
          xct%sigma_a_EUV(20,isp) = 28.070d-22
          xct%sigma_a_EUV(21,isp) = 32.060d-22
          xct%sigma_a_EUV(22,isp) = 26.017d-22
          xct%sigma_a_EUV(23,isp) = 21.919d-22
          xct%sigma_a_EUV(24,isp) = 27.440d-22
          xct%sigma_a_EUV(25,isp) = 28.535d-22
          xct%sigma_a_EUV(26,isp) = 20.800d-22
          xct%sigma_a_EUV(27,isp) = 18.910d-22
          xct%sigma_a_EUV(28,isp) = 26.668d-22
          xct%sigma_a_EUV(29,isp) = 22.145d-22
          xct%sigma_a_EUV(30,isp) = 16.631d-22
          xct%sigma_a_EUV(31,isp) =  8.562d-22
          xct%sigma_a_EUV(32,isp) = 12.817d-22
          xct%sigma_a_EUV(33,isp) = 18.730d-22
          xct%sigma_a_EUV(34,isp) = 21.108d-22
          xct%sigma_a_EUV(35,isp) =  1.630d-22
          xct%sigma_a_EUV(36,isp) =  1.050d-22
          xct%sigma_a_EUV(37,isp) =  1.346d-22
        end if

        if (spl%species(isp) == 'O') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  O  ---------
          xct%sigma_a_EUV(1,isp)  =  0.730d-22
          xct%sigma_a_EUV(2,isp)  =  1.839d-22
          xct%sigma_a_EUV(3,isp)  =  3.732d-22
          xct%sigma_a_EUV(4,isp)  =  5.202d-22
          xct%sigma_a_EUV(5,isp)  =  6.050d-22
          xct%sigma_a_EUV(6,isp)  =  7.080d-22
          xct%sigma_a_EUV(7,isp)  =  6.461d-22
          xct%sigma_a_EUV(8,isp)  =  7.680d-22
          xct%sigma_a_EUV(9,isp)  =  7.700d-22
          xct%sigma_a_EUV(10,isp) =  8.693d-22
          xct%sigma_a_EUV(11,isp) =  9.840d-22
          xct%sigma_a_EUV(12,isp) =  9.687d-22
          xct%sigma_a_EUV(13,isp) = 11.496d-22
          xct%sigma_a_EUV(14,isp) = 11.930d-22
          xct%sigma_a_EUV(15,isp) = 12.127d-22
          xct%sigma_a_EUV(16,isp) = 12.059d-22
          xct%sigma_a_EUV(17,isp) = 12.590d-22
          xct%sigma_a_EUV(18,isp) = 13.090d-22
          xct%sigma_a_EUV(19,isp) = 13.024d-22
          xct%sigma_a_EUV(20,isp) = 13.400d-22
          xct%sigma_a_EUV(21,isp) = 13.400d-22
          xct%sigma_a_EUV(22,isp) = 13.365d-22
          xct%sigma_a_EUV(23,isp) = 17.245d-22
          xct%sigma_a_EUV(24,isp) = 11.460d-22
          xct%sigma_a_EUV(25,isp) = 10.736d-22
          xct%sigma_a_EUV(26,isp) =  4.000d-22
          xct%sigma_a_EUV(27,isp) =  3.890d-22
          xct%sigma_a_EUV(28,isp) =  3.749d-22
          xct%sigma_a_EUV(29,isp) =  5.091d-22
          xct%sigma_a_EUV(30,isp) =  3.498d-22
          xct%sigma_a_EUV(31,isp) =  4.554d-22
          xct%sigma_a_EUV(32,isp) =  1.315d-22
          xct%sigma_a_EUV(33,isp) =  0.000d-22
          xct%sigma_a_EUV(34,isp) =  0.000d-22
          xct%sigma_a_EUV(35,isp) =  0.000d-22
          xct%sigma_a_EUV(36,isp) =  0.000d-22
          xct%sigma_a_EUV(37,isp) =  0.000d-22
        end if

        if (spl%species(isp) == 'CO2') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  CO2  ---------
          xct%sigma_a_EUV(1,isp)  =  1.550d-22
          xct%sigma_a_EUV(2,isp)  =  4.616d-22
          xct%sigma_a_EUV(3,isp)  =  9.089d-22
          xct%sigma_a_EUV(4,isp)  = 14.361d-22
          xct%sigma_a_EUV(5,isp)  = 16.505d-22
          xct%sigma_a_EUV(6,isp)  = 19.016d-22
          xct%sigma_a_EUV(7,isp)  = 17.518d-22
          xct%sigma_a_EUV(8,isp)  = 21.492d-22
          xct%sigma_a_EUV(9,isp)  = 21.594d-22
          xct%sigma_a_EUV(10,isp) = 23.574d-22
          xct%sigma_a_EUV(11,isp) = 25.269d-22
          xct%sigma_a_EUV(12,isp) = 24.871d-22
          xct%sigma_a_EUV(13,isp) = 28.271d-22
          xct%sigma_a_EUV(14,isp) = 29.526d-22
          xct%sigma_a_EUV(15,isp) = 30.254d-22
          xct%sigma_a_EUV(16,isp) = 31.491d-22
          xct%sigma_a_EUV(17,isp) = 33.202d-22
          xct%sigma_a_EUV(18,isp) = 34.200d-22
          xct%sigma_a_EUV(19,isp) = 34.913d-22
          xct%sigma_a_EUV(20,isp) = 35.303d-22
          xct%sigma_a_EUV(21,isp) = 34.300d-22
          xct%sigma_a_EUV(22,isp) = 34.447d-22
          xct%sigma_a_EUV(23,isp) = 33.699d-22
          xct%sigma_a_EUV(24,isp) = 23.518d-22
          xct%sigma_a_EUV(25,isp) = 32.832d-22
          xct%sigma_a_EUV(26,isp) = 93.839d-22
          xct%sigma_a_EUV(27,isp) = 61.939d-22
          xct%sigma_a_EUV(28,isp) = 26.493d-22
          xct%sigma_a_EUV(29,isp) = 39.831d-22
          xct%sigma_a_EUV(30,isp) = 13.980d-22
          xct%sigma_a_EUV(31,isp) = 44.673d-22
          xct%sigma_a_EUV(32,isp) = 52.081d-22
          xct%sigma_a_EUV(33,isp) = 42.869d-22
          xct%sigma_a_EUV(34,isp) = 50.311d-22
          xct%sigma_a_EUV(35,isp) = 15.100d-22
          xct%sigma_a_EUV(36,isp) = 14.200d-22
          xct%sigma_a_EUV(37,isp) = 18.241d-22
        end if

        if (spl%species(isp) == 'CO') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  CO  ---------
          xct%sigma_a_EUV(1,isp)  =  0.866d-22
          xct%sigma_a_EUV(2,isp)  =  2.391d-22
          xct%sigma_a_EUV(3,isp)  =  4.671d-22
          xct%sigma_a_EUV(4,isp)  =  7.011d-22
          xct%sigma_a_EUV(5,isp)  =  8.614d-22
          xct%sigma_a_EUV(6,isp)  = 10.541d-22
          xct%sigma_a_EUV(7,isp)  =  9.424d-22
          xct%sigma_a_EUV(8,isp)  = 11.867d-22
          xct%sigma_a_EUV(9,isp)  = 11.900d-22
          xct%sigma_a_EUV(10,isp) = 13.441d-22
          xct%sigma_a_EUV(11,isp) = 15.259d-22
          xct%sigma_a_EUV(12,isp) = 14.956d-22
          xct%sigma_a_EUV(13,isp) = 17.956d-22
          xct%sigma_a_EUV(14,isp) = 20.173d-22
          xct%sigma_a_EUV(15,isp) = 20.574d-22
          xct%sigma_a_EUV(16,isp) = 21.085d-22
          xct%sigma_a_EUV(17,isp) = 21.624d-22
          xct%sigma_a_EUV(18,isp) = 22.000d-22
          xct%sigma_a_EUV(19,isp) = 21.910d-22
          xct%sigma_a_EUV(20,isp) = 22.100d-22
          xct%sigma_a_EUV(21,isp) = 22.025d-22
          xct%sigma_a_EUV(22,isp) = 21.915d-22
          xct%sigma_a_EUV(23,isp) = 21.036d-22
          xct%sigma_a_EUV(24,isp) = 23.853d-22
          xct%sigma_a_EUV(25,isp) = 25.501d-22
          xct%sigma_a_EUV(26,isp) = 26.276d-22
          xct%sigma_a_EUV(27,isp) = 15.262d-22
          xct%sigma_a_EUV(28,isp) = 33.132d-22
          xct%sigma_a_EUV(29,isp) = 20.535d-22
          xct%sigma_a_EUV(30,isp) = 22.608d-22
          xct%sigma_a_EUV(31,isp) = 36.976d-22
          xct%sigma_a_EUV(32,isp) = 50.318d-22
          xct%sigma_a_EUV(33,isp) = 28.500d-22
          xct%sigma_a_EUV(34,isp) = 52.827d-22
          xct%sigma_a_EUV(35,isp) =  1.388d-22
          xct%sigma_a_EUV(36,isp) =  1.388d-22
          xct%sigma_a_EUV(37,isp) =  8.568d-22
        end if

        if (spl%species(isp) == 'N2') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  N2  ---------
          xct%sigma_a_EUV(1,isp)  =   0.720d-22
          xct%sigma_a_EUV(2,isp)  =   2.261d-22
          xct%sigma_a_EUV(3,isp)  =   4.958d-22
          xct%sigma_a_EUV(4,isp)  =   8.392d-22
          xct%sigma_a_EUV(5,isp)  =  10.210d-22
          xct%sigma_a_EUV(6,isp)  =  10.900d-22
          xct%sigma_a_EUV(7,isp)  =  10.493d-22
          xct%sigma_a_EUV(8,isp)  =  11.670d-22
          xct%sigma_a_EUV(9,isp)  =  11.700d-22
          xct%sigma_a_EUV(10,isp) =  13.857d-22
          xct%sigma_a_EUV(11,isp) =  16.910d-22
          xct%sigma_a_EUV(12,isp) =  16.395d-22
          xct%sigma_a_EUV(13,isp) =  21.675d-22
          xct%sigma_a_EUV(14,isp) =  23.160d-22
          xct%sigma_a_EUV(15,isp) =  23.471d-22
          xct%sigma_a_EUV(16,isp) =  24.501d-22
          xct%sigma_a_EUV(17,isp) =  24.130d-22
          xct%sigma_a_EUV(18,isp) =  22.400d-22
          xct%sigma_a_EUV(19,isp) =  22.787d-22
          xct%sigma_a_EUV(20,isp) =  22.790d-22
          xct%sigma_a_EUV(21,isp) =  23.370d-22
          xct%sigma_a_EUV(22,isp) =  23.339d-22
          xct%sigma_a_EUV(23,isp) =  31.755d-22
          xct%sigma_a_EUV(24,isp) =  26.540d-22
          xct%sigma_a_EUV(25,isp) =  24.662d-22
          xct%sigma_a_EUV(26,isp) = 120.490d-22
          xct%sigma_a_EUV(27,isp) =  14.180d-22
          xct%sigma_a_EUV(28,isp) =  16.487d-22
          xct%sigma_a_EUV(29,isp) =  33.578d-22
          xct%sigma_a_EUV(30,isp) =  16.992d-22
          xct%sigma_a_EUV(31,isp) =  20.249d-22
          xct%sigma_a_EUV(32,isp) =   9.680d-22
          xct%sigma_a_EUV(33,isp) =   2.240d-22
          xct%sigma_a_EUV(34,isp) =  50.988d-22
          xct%sigma_a_EUV(35,isp) =   0.000d-22
          xct%sigma_a_EUV(36,isp) =   0.000d-22
          xct%sigma_a_EUV(37,isp) =   0.000d-22
        end if

        if (spl%species(isp) == 'H') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  H  ---------
          xct%sigma_a_EUV(1,isp)  = 0.0024d-22
          xct%sigma_a_EUV(2,isp)  = 0.0169d-22
          xct%sigma_a_EUV(3,isp)  = 0.0483d-22
          xct%sigma_a_EUV(4,isp)  = 0.1007d-22
          xct%sigma_a_EUV(5,isp)  = 0.1405d-22
          xct%sigma_a_EUV(6,isp)  = 0.1913d-22
          xct%sigma_a_EUV(7,isp)  = 0.1676d-22
          xct%sigma_a_EUV(8,isp)  = 0.2324d-22
          xct%sigma_a_EUV(9,isp)  = 0.2334d-22
          xct%sigma_a_EUV(10,isp) = 0.3077d-22
          xct%sigma_a_EUV(11,isp) = 0.4152d-22
          xct%sigma_a_EUV(12,isp) = 0.3984d-22
          xct%sigma_a_EUV(13,isp) = 0.6163d-22
          xct%sigma_a_EUV(14,isp) = 0.8387d-22
          xct%sigma_a_EUV(15,isp) = 0.9739d-22
          xct%sigma_a_EUV(16,isp) = 1.1990d-22
          xct%sigma_a_EUV(17,isp) = 1.4190d-22
          xct%sigma_a_EUV(18,isp) = 1.6620d-22
          xct%sigma_a_EUV(19,isp) = 1.6200d-22
          xct%sigma_a_EUV(20,isp) = 1.8880d-22
          xct%sigma_a_EUV(21,isp) = 2.0790d-22
          xct%sigma_a_EUV(22,isp) = 2.0760d-22
          xct%sigma_a_EUV(23,isp) = 2.6410d-22
          xct%sigma_a_EUV(24,isp) = 2.8970d-22
          xct%sigma_a_EUV(25,isp) = 3.1730d-22
          xct%sigma_a_EUV(26,isp) = 3.7300d-22
          xct%sigma_a_EUV(27,isp) = 3.8070d-22
          xct%sigma_a_EUV(28,isp) = 4.0930d-22
          xct%sigma_a_EUV(29,isp) = 3.8680d-22
          xct%sigma_a_EUV(30,isp) = 4.7840d-22
          xct%sigma_a_EUV(31,isp) = 5.6700d-22
          xct%sigma_a_EUV(32,isp) = 3.4690d-22
          xct%sigma_a_EUV(33,isp) = 0.0000d-22
          xct%sigma_a_EUV(34,isp) = 0.0000d-22
          xct%sigma_a_EUV(35,isp) = 0.0000d-22
          xct%sigma_a_EUV(36,isp) = 0.0000d-22
          xct%sigma_a_EUV(37,isp) = 0.0000d-22
        end if

        if (spl%species(isp) == 'H2') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  H2  --------
          xct%sigma_a_EUV(1,isp)  =  0.0108d-22
          xct%sigma_a_EUV(2,isp)  =  0.0798d-22
          xct%sigma_a_EUV(3,isp)  =  0.2085d-22
          xct%sigma_a_EUV(4,isp)  =  0.4333d-22
          xct%sigma_a_EUV(5,isp)  =  0.6037d-22
          xct%sigma_a_EUV(6,isp)  =  0.8388d-22
          xct%sigma_a_EUV(7,isp)  =  0.7296d-22
          xct%sigma_a_EUV(8,isp)  =  1.0180d-22
          xct%sigma_a_EUV(9,isp)  =  1.0220d-22
          xct%sigma_a_EUV(10,isp) =  1.4170d-22
          xct%sigma_a_EUV(11,isp) =  1.9420d-22
          xct%sigma_a_EUV(12,isp) =  1.9010d-22
          xct%sigma_a_EUV(13,isp) =  3.0250d-22
          xct%sigma_a_EUV(14,isp) =  3.8700d-22
          xct%sigma_a_EUV(15,isp) =  4.5020d-22
          xct%sigma_a_EUV(16,isp) =  5.3560d-22
          xct%sigma_a_EUV(17,isp) =  6.1680d-22
          xct%sigma_a_EUV(18,isp) =  7.0210d-22
          xct%sigma_a_EUV(19,isp) =  6.8640d-22
          xct%sigma_a_EUV(20,isp) =  7.8110d-22
          xct%sigma_a_EUV(21,isp) =  8.4640d-22
          xct%sigma_a_EUV(22,isp) =  8.4450d-22
          xct%sigma_a_EUV(23,isp) =  9.9000d-22
          xct%sigma_a_EUV(24,isp) = 10.7310d-22
          xct%sigma_a_EUV(25,isp) = 11.3720d-22
          xct%sigma_a_EUV(26,isp) = 10.7550d-22
          xct%sigma_a_EUV(27,isp) =  8.6400d-22
          xct%sigma_a_EUV(28,isp) =  7.3390d-22
          xct%sigma_a_EUV(29,isp) =  8.7480d-22
          xct%sigma_a_EUV(30,isp) =  8.2530d-22
          xct%sigma_a_EUV(31,isp) =  0.4763d-22
          xct%sigma_a_EUV(32,isp) =  0.1853d-22
          xct%sigma_a_EUV(33,isp) =  0.0000d-22
          xct%sigma_a_EUV(34,isp) =  0.0456d-22
          xct%sigma_a_EUV(35,isp) =  0.0000d-22
          xct%sigma_a_EUV(36,isp) =  0.0000d-22
          xct%sigma_a_EUV(37,isp) =  0.0000d-22
        end if

        if (spl%species(isp) == 'He') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  He  --------
          xct%sigma_a_EUV(1,isp)  =  0.1441d-22
          xct%sigma_a_EUV(2,isp)  =  0.4785d-22
          xct%sigma_a_EUV(3,isp)  =  1.1571d-22
          xct%sigma_a_EUV(4,isp)  =  1.6008d-22
          xct%sigma_a_EUV(5,isp)  =  2.1212d-22
          xct%sigma_a_EUV(6,isp)  =  2.5947d-22
          xct%sigma_a_EUV(7,isp)  =  2.3205d-22
          xct%sigma_a_EUV(8,isp)  =  2.9529d-22
          xct%sigma_a_EUV(9,isp)  =  2.9618d-22
          xct%sigma_a_EUV(10,isp) =  3.5437d-22
          xct%sigma_a_EUV(11,isp) =  4.2675d-22
          xct%sigma_a_EUV(12,isp) =  4.1424d-22
          xct%sigma_a_EUV(13,isp) =  5.4466d-22
          xct%sigma_a_EUV(14,isp) =  6.5631d-22
          xct%sigma_a_EUV(15,isp) =  7.2084d-22
          xct%sigma_a_EUV(16,isp) =  0.9581d-22
          xct%sigma_a_EUV(17,isp) =  0.0000d-22
          xct%sigma_a_EUV(18,isp) =  0.0000d-22
          xct%sigma_a_EUV(19,isp) =  0.0000d-22
          xct%sigma_a_EUV(20,isp) =  0.0000d-22
          xct%sigma_a_EUV(21,isp) =  0.0000d-22
          xct%sigma_a_EUV(22,isp) =  0.0000d-22
          xct%sigma_a_EUV(23,isp) =  0.0000d-22
          xct%sigma_a_EUV(24,isp) =  0.0000d-22
          xct%sigma_a_EUV(25,isp) =  0.0000d-22
          xct%sigma_a_EUV(26,isp) =  0.0000d-22
          xct%sigma_a_EUV(27,isp) =  0.0000d-22
          xct%sigma_a_EUV(28,isp) =  0.0000d-22
          xct%sigma_a_EUV(29,isp) =  0.0000d-22
          xct%sigma_a_EUV(30,isp) =  0.0000d-22
          xct%sigma_a_EUV(31,isp) =  0.0000d-22
          xct%sigma_a_EUV(32,isp) =  0.0000d-22
          xct%sigma_a_EUV(33,isp) =  0.0000d-22
          xct%sigma_a_EUV(34,isp) =  0.0000d-22
          xct%sigma_a_EUV(35,isp) =  0.0000d-22
          xct%sigma_a_EUV(36,isp) =  0.0000d-22
          xct%sigma_a_EUV(37,isp) =  0.0000d-22
        end if

        if (spl%species(isp) == 'CH4') then
          xct%label_sigma_a_EUV(isp) = 1
          !-----------------------------------  CH4  --------
          xct%sigma_a_EUV(1,isp)  =  0.204d-22
          xct%sigma_a_EUV(2,isp)  =  0.593d-22
          xct%sigma_a_EUV(3,isp)  =  1.496d-22
          xct%sigma_a_EUV(4,isp)  =  2.794d-22
          xct%sigma_a_EUV(5,isp)  =  3.857d-22
          xct%sigma_a_EUV(6,isp)  =  5.053d-22
          xct%sigma_a_EUV(7,isp)  =  4.360d-22
          xct%sigma_a_EUV(8,isp)  =  6.033d-22
          xct%sigma_a_EUV(9,isp)  =  6.059d-22
          xct%sigma_a_EUV(10,isp) =  7.829d-22
          xct%sigma_a_EUV(11,isp) = 10.165d-22
          xct%sigma_a_EUV(12,isp) =  9.776d-22
          xct%sigma_a_EUV(13,isp) = 14.701d-22
          xct%sigma_a_EUV(14,isp) = 18.770d-22
          xct%sigma_a_EUV(15,isp) = 21.449d-22
          xct%sigma_a_EUV(16,isp) = 24.644d-22
          xct%sigma_a_EUV(17,isp) = 27.924d-22
          xct%sigma_a_EUV(18,isp) = 31.052d-22
          xct%sigma_a_EUV(19,isp) = 30.697d-22
          xct%sigma_a_EUV(20,isp) = 33.178d-22
          xct%sigma_a_EUV(21,isp) = 35.276d-22
          xct%sigma_a_EUV(22,isp) = 34.990d-22
          xct%sigma_a_EUV(23,isp) = 39.280d-22
          xct%sigma_a_EUV(24,isp) = 41.069d-22
          xct%sigma_a_EUV(25,isp) = 42.927d-22
          xct%sigma_a_EUV(26,isp) = 45.458d-22
          xct%sigma_a_EUV(27,isp) = 45.716d-22
          xct%sigma_a_EUV(28,isp) = 46.472d-22
          xct%sigma_a_EUV(29,isp) = 45.921d-22
          xct%sigma_a_EUV(30,isp) = 48.327d-22
          xct%sigma_a_EUV(31,isp) = 48.968d-22
          xct%sigma_a_EUV(32,isp) = 48.001d-22
          xct%sigma_a_EUV(33,isp) = 41.154d-22
          xct%sigma_a_EUV(34,isp) = 38.192d-22
          xct%sigma_a_EUV(35,isp) = 32.700d-22
          xct%sigma_a_EUV(36,isp) = 30.121d-22
          xct%sigma_a_EUV(37,isp) = 29.108d-22
        end if

        if (spl%species(isp) == 'C2H2') then
          xct%label_sigma_a_EUV(isp) = 1
          !----------------------------------  C2H2  --------
          xct%sigma_a_EUV(1,isp)  =  0.0000d-22
          xct%sigma_a_EUV(2,isp)  =  0.1700d-22
          xct%sigma_a_EUV(3,isp)  =  1.3100d-22
          xct%sigma_a_EUV(4,isp)  =  3.5620d-22
          xct%sigma_a_EUV(5,isp)  =  5.2829d-22
          xct%sigma_a_EUV(6,isp)  =  7.2894d-22
          xct%sigma_a_EUV(7,isp)  =  6.7353d-22
          xct%sigma_a_EUV(8,isp)  =  8.8883d-22
          xct%sigma_a_EUV(9,isp)  =  8.9250d-22
          xct%sigma_a_EUV(10,isp) = 10.5676d-22
          xct%sigma_a_EUV(11,isp) = 13.8224d-22
          xct%sigma_a_EUV(12,isp) = 14.2143d-22
          xct%sigma_a_EUV(13,isp) = 18.2193d-22
          xct%sigma_a_EUV(14,isp) = 21.5703d-22
          xct%sigma_a_EUV(15,isp) = 22.1713d-22
          xct%sigma_a_EUV(16,isp) = 25.6502d-22
          xct%sigma_a_EUV(17,isp) = 27.7643d-22
          xct%sigma_a_EUV(18,isp) = 29.6222d-22
          xct%sigma_a_EUV(19,isp) = 28.9982d-22
          xct%sigma_a_EUV(20,isp) = 31.0559d-22
          xct%sigma_a_EUV(21,isp) = 32.7357d-22
          xct%sigma_a_EUV(22,isp) = 32.7562d-22
          xct%sigma_a_EUV(23,isp) = 37.4962d-22
          xct%sigma_a_EUV(24,isp) = 40.3110d-22
          xct%sigma_a_EUV(25,isp) = 42.1700d-22
          xct%sigma_a_EUV(26,isp) = 45.5431d-22
          xct%sigma_a_EUV(27,isp) = 45.9497d-22
          xct%sigma_a_EUV(28,isp) = 44.4371d-22
          xct%sigma_a_EUV(29,isp) = 43.7200d-22
          xct%sigma_a_EUV(30,isp) = 40.1824d-22
          xct%sigma_a_EUV(31,isp) = 34.8990d-22
          xct%sigma_a_EUV(32,isp) = 30.9561d-22
          xct%sigma_a_EUV(33,isp) = 27.8087d-22
          xct%sigma_a_EUV(34,isp) = 27.9595d-22
          xct%sigma_a_EUV(35,isp) = 25.1875d-22
          xct%sigma_a_EUV(36,isp) = 24.8720d-22
          xct%sigma_a_EUV(37,isp) = 26.1462d-22
        end if

        if (spl%species(isp) == 'C2H4') then
          xct%label_sigma_a_EUV(isp) = 1
          !----------------------------------  C2H4  --------
          ! C2H4: Ibuki+ 1989
          xct%sigma_a_EUV(1,isp)  =  0.00d-22
          xct%sigma_a_EUV(2,isp)  =  0.67d-22
          xct%sigma_a_EUV(3,isp)  =  1.71d-22
          xct%sigma_a_EUV(4,isp)  =  3.43d-22
          xct%sigma_a_EUV(5,isp)  =  6.43d-22
          xct%sigma_a_EUV(6,isp)  =  8.19d-22
          xct%sigma_a_EUV(7,isp)  =  6.04d-22
          xct%sigma_a_EUV(8,isp)  =  9.69d-22
          xct%sigma_a_EUV(9,isp)  =  9.72d-22
          xct%sigma_a_EUV(10,isp) =  9.48d-22
          xct%sigma_a_EUV(11,isp) = 15.17d-22
          xct%sigma_a_EUV(12,isp) = 13.63d-22
          xct%sigma_a_EUV(13,isp) = 18.98d-22
          xct%sigma_a_EUV(14,isp) = 28.53d-22
          xct%sigma_a_EUV(15,isp) = 25.75d-22
          xct%sigma_a_EUV(16,isp) = 34.17d-22
          xct%sigma_a_EUV(17,isp) = 41.89d-22
          xct%sigma_a_EUV(18,isp) = 44.87d-22
          xct%sigma_a_EUV(19,isp) = 41.39d-22
          xct%sigma_a_EUV(20,isp) = 47.57d-22
          xct%sigma_a_EUV(21,isp) = 49.66d-22
          xct%sigma_a_EUV(22,isp) = 46.64d-22
          xct%sigma_a_EUV(23,isp) = 52.23d-22
          xct%sigma_a_EUV(24,isp) = 56.79d-22
          xct%sigma_a_EUV(25,isp) = 56.56d-22
          xct%sigma_a_EUV(26,isp) = 54.62d-22
          xct%sigma_a_EUV(27,isp) = 54.24d-22
          xct%sigma_a_EUV(28,isp) = 54.61d-22
          xct%sigma_a_EUV(29,isp) = 55.72d-22
          xct%sigma_a_EUV(30,isp) = 55.11d-22
          xct%sigma_a_EUV(31,isp) = 54.86d-22
          xct%sigma_a_EUV(32,isp) = 47.18d-22
          xct%sigma_a_EUV(33,isp) = 45.80d-22
          xct%sigma_a_EUV(34,isp) = 43.64d-22
          xct%sigma_a_EUV(35,isp) = 44.13d-22
          xct%sigma_a_EUV(36,isp) = 43.57d-22
          xct%sigma_a_EUV(37,isp) = 46.48d-22
        end if

        if (spl%species(isp) == 'C2H6') then
          xct%label_sigma_a_EUV(isp) = 1
          !----------------------------------  C2H6  --------
          ! C2H6: Au+ 1993, assumed to be total ionization
          xct%sigma_a_EUV(1,isp)  =  0.00d-22
          xct%sigma_a_EUV(2,isp)  =  0.00d-22
          xct%sigma_a_EUV(3,isp)  =  0.00d-22
          xct%sigma_a_EUV(4,isp)  =  2.61d-22
          xct%sigma_a_EUV(5,isp)  =  5.77d-22
          xct%sigma_a_EUV(6,isp)  =  7.74d-22
          xct%sigma_a_EUV(7,isp)  =  5.36d-22
          xct%sigma_a_EUV(8,isp)  =  9.41d-22
          xct%sigma_a_EUV(9,isp)  =  9.47d-22
          xct%sigma_a_EUV(10,isp) =  9.07d-22
          xct%sigma_a_EUV(11,isp) = 16.72d-22
          xct%sigma_a_EUV(12,isp) = 14.44d-22
          xct%sigma_a_EUV(13,isp) = 20.87d-22
          xct%sigma_a_EUV(14,isp) = 32.18d-22
          xct%sigma_a_EUV(15,isp) = 30.08d-22
          xct%sigma_a_EUV(16,isp) = 38.26d-22
          xct%sigma_a_EUV(17,isp) = 47.77d-22
          xct%sigma_a_EUV(18,isp) = 53.03d-22
          xct%sigma_a_EUV(19,isp) = 46.87d-22
          xct%sigma_a_EUV(20,isp) = 57.42d-22
          xct%sigma_a_EUV(21,isp) = 62.04d-22
          xct%sigma_a_EUV(22,isp) = 55.91d-22
          xct%sigma_a_EUV(23,isp) = 65.42d-22
          xct%sigma_a_EUV(24,isp) = 69.73d-22
          xct%sigma_a_EUV(25,isp) = 69.81d-22
          xct%sigma_a_EUV(26,isp) = 76.21d-22
          xct%sigma_a_EUV(27,isp) = 76.38d-22
          xct%sigma_a_EUV(28,isp) = 75.88d-22
          xct%sigma_a_EUV(29,isp) = 75.47d-22
          xct%sigma_a_EUV(30,isp) = 75.39d-22
          xct%sigma_a_EUV(31,isp) = 65.67d-22
          xct%sigma_a_EUV(32,isp) = 60.14d-22
          xct%sigma_a_EUV(33,isp) = 38.89d-22
          xct%sigma_a_EUV(34,isp) = 50.63d-22
          xct%sigma_a_EUV(35,isp) = 19.96d-22
          xct%sigma_a_EUV(36,isp) = 17.77d-22
          xct%sigma_a_EUV(37,isp) = 29.07d-22
        end if

        if (spl%species(isp) == 'Na') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  Na  ---------
          xct%sigma_a_EUV(1,isp)  = 8.384d-23
          xct%sigma_a_EUV(2,isp)  = 3.717d-22
          xct%sigma_a_EUV(3,isp)  = 7.048d-22
          xct%sigma_a_EUV(4,isp)  = 8.230d-22
          xct%sigma_a_EUV(5,isp)  = 7.036d-22
          xct%sigma_a_EUV(6,isp)  = 5.689d-22
          xct%sigma_a_EUV(7,isp)  = 7.313d-22
          xct%sigma_a_EUV(8,isp)  = 4.960d-22
          xct%sigma_a_EUV(9,isp)  = 4.948d-22
          xct%sigma_a_EUV(10,isp) = 5.056d-22
          xct%sigma_a_EUV(11,isp) = 6.129d-24
          xct%sigma_a_EUV(12,isp) = 5.773d-24
          xct%sigma_a_EUV(13,isp) = 6.746d-24
          xct%sigma_a_EUV(14,isp) = 7.943d-24
          xct%sigma_a_EUV(15,isp) = 7.673d-24
          xct%sigma_a_EUV(16,isp) = 8.538d-24
          xct%sigma_a_EUV(17,isp) = 9.395d-24
          xct%sigma_a_EUV(18,isp) = 9.826d-24
          xct%sigma_a_EUV(19,isp) = 9.330d-24
          xct%sigma_a_EUV(20,isp) = 1.016d-23
          xct%sigma_a_EUV(21,isp) = 1.041d-23
          xct%sigma_a_EUV(22,isp) = 1.003d-23
          xct%sigma_a_EUV(23,isp) = 1.065d-23
          xct%sigma_a_EUV(24,isp) = 1.120d-23
          xct%sigma_a_EUV(25,isp) = 1.117d-23
          xct%sigma_a_EUV(26,isp) = 1.170d-23
          xct%sigma_a_EUV(27,isp) = 1.173d-23
          xct%sigma_a_EUV(28,isp) = 1.185d-23
          xct%sigma_a_EUV(29,isp) = 1.159d-23
          xct%sigma_a_EUV(30,isp) = 1.191d-23
          xct%sigma_a_EUV(31,isp) = 1.212d-23
          xct%sigma_a_EUV(32,isp) = 1.222d-23
          xct%sigma_a_EUV(33,isp) = 1.218d-23
          xct%sigma_a_EUV(34,isp) = 1.223d-23
          xct%sigma_a_EUV(35,isp) = 1.203d-23
          xct%sigma_a_EUV(36,isp) = 1.201d-23
          xct%sigma_a_EUV(37,isp) = 1.212d-23
        end if

        if (spl%species(isp) == 'Mg') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  Mg  ---------
          xct%sigma_a_EUV(1,isp)  = 1.209d-22
          xct%sigma_a_EUV(2,isp)  = 4.936d-22
          xct%sigma_a_EUV(3,isp)  = 6.757d-22
          xct%sigma_a_EUV(4,isp)  = 4.220d-22
          xct%sigma_a_EUV(5,isp)  = 1.579d-23
          xct%sigma_a_EUV(6,isp)  = 1.755d-23
          xct%sigma_a_EUV(7,isp)  = 1.538d-23
          xct%sigma_a_EUV(8,isp)  = 1.870d-23
          xct%sigma_a_EUV(9,isp)  = 1.872d-23
          xct%sigma_a_EUV(10,isp) = 1.850d-23
          xct%sigma_a_EUV(11,isp) = 2.208d-23
          xct%sigma_a_EUV(12,isp) = 2.122d-23
          xct%sigma_a_EUV(13,isp) = 2.344d-23
          xct%sigma_a_EUV(14,isp) = 2.553d-23
          xct%sigma_a_EUV(15,isp) = 2.513d-23
          xct%sigma_a_EUV(16,isp) = 2.626d-23
          xct%sigma_a_EUV(17,isp) = 2.683d-23
          xct%sigma_a_EUV(18,isp) = 2.686d-23
          xct%sigma_a_EUV(19,isp) = 2.681d-23
          xct%sigma_a_EUV(20,isp) = 2.673d-23
          xct%sigma_a_EUV(21,isp) = 2.654d-23
          xct%sigma_a_EUV(22,isp) = 2.680d-23
          xct%sigma_a_EUV(23,isp) = 2.625d-23
          xct%sigma_a_EUV(24,isp) = 2.511d-23
          xct%sigma_a_EUV(25,isp) = 2.520d-23
          xct%sigma_a_EUV(26,isp) = 2.314d-23
          xct%sigma_a_EUV(27,isp) = 2.295d-23
          xct%sigma_a_EUV(28,isp) = 2.221d-23
          xct%sigma_a_EUV(29,isp) = 2.368d-23
          xct%sigma_a_EUV(30,isp) = 2.177d-23
          xct%sigma_a_EUV(31,isp) = 1.952d-23
          xct%sigma_a_EUV(32,isp) = 1.702d-23
          xct%sigma_a_EUV(33,isp) = 1.288d-23
          xct%sigma_a_EUV(34,isp) = 1.436d-23
          xct%sigma_a_EUV(35,isp) = 1.022d-23
          xct%sigma_a_EUV(36,isp) = 9.891d-24
          xct%sigma_a_EUV(37,isp) = 1.162d-23
        end if

        if (spl%species(isp) == 'Fe') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  Fe  ---------
          xct%sigma_a_EUV(1,isp)  = 1.656d-22
          xct%sigma_a_EUV(2,isp)  = 4.475d-22
          xct%sigma_a_EUV(3,isp)  = 6.524d-22
          xct%sigma_a_EUV(4,isp)  = 7.745d-22
          xct%sigma_a_EUV(5,isp)  = 8.748d-22
          xct%sigma_a_EUV(6,isp)  = 8.891d-22
          xct%sigma_a_EUV(7,isp)  = 8.686d-22
          xct%sigma_a_EUV(8,isp)  = 8.881d-22
          xct%sigma_a_EUV(9,isp)  = 8.879d-22
          xct%sigma_a_EUV(10,isp) = 8.888d-22
          xct%sigma_a_EUV(11,isp) = 8.369d-22
          xct%sigma_a_EUV(12,isp) = 8.571d-22
          xct%sigma_a_EUV(13,isp) = 7.936d-22
          xct%sigma_a_EUV(14,isp) = 6.893d-22
          xct%sigma_a_EUV(15,isp) = 7.145d-22
          xct%sigma_a_EUV(16,isp) = 6.318d-22
          xct%sigma_a_EUV(17,isp) = 5.478d-22
          xct%sigma_a_EUV(18,isp) = 5.069d-22
          xct%sigma_a_EUV(19,isp) = 5.541d-22
          xct%sigma_a_EUV(20,isp) = 4.760d-22
          xct%sigma_a_EUV(21,isp) = 4.545d-22
          xct%sigma_a_EUV(22,isp) = 4.874d-22
          xct%sigma_a_EUV(23,isp) = 4.354d-22
          xct%sigma_a_EUV(24,isp) = 3.988d-22
          xct%sigma_a_EUV(25,isp) = 4.005d-22
          xct%sigma_a_EUV(26,isp) = 3.826d-22
          xct%sigma_a_EUV(27,isp) = 3.826d-22
          xct%sigma_a_EUV(28,isp) = 3.842d-22
          xct%sigma_a_EUV(29,isp) = 3.839d-22
          xct%sigma_a_EUV(30,isp) = 3.863d-22
          xct%sigma_a_EUV(31,isp) = 4.074d-22
          xct%sigma_a_EUV(32,isp) = 4.469d-22
          xct%sigma_a_EUV(33,isp) = 1.608d-25
          xct%sigma_a_EUV(34,isp) = 5.041d-22
          xct%sigma_a_EUV(35,isp) = 6.817d-25
          xct%sigma_a_EUV(36,isp) = 8.101d-25
          xct%sigma_a_EUV(37,isp) = 3.029d-25
        end if

        if (spl%species(isp) == 'Si') then
          xct%label_sigma_a_EUV(isp) = 1
          !------------------------------------  Si  ---------
          xct%sigma_a_EUV(1,isp)  = 2.209d-22
          xct%sigma_a_EUV(2,isp)  = 5.039d-22
          xct%sigma_a_EUV(3,isp)  = 4.122d-23
          xct%sigma_a_EUV(4,isp)  = 5.743d-23
          xct%sigma_a_EUV(5,isp)  = 6.910d-23
          xct%sigma_a_EUV(6,isp)  = 7.209d-23
          xct%sigma_a_EUV(7,isp)  = 6.817d-23
          xct%sigma_a_EUV(8,isp)  = 7.314d-23
          xct%sigma_a_EUV(9,isp)  = 7.316d-23
          xct%sigma_a_EUV(10,isp) = 7.302d-23
          xct%sigma_a_EUV(11,isp) = 7.177d-23
          xct%sigma_a_EUV(12,isp) = 7.280d-23
          xct%sigma_a_EUV(13,isp) = 6.912d-23
          xct%sigma_a_EUV(14,isp) = 6.249d-23
          xct%sigma_a_EUV(15,isp) = 6.400d-23
          xct%sigma_a_EUV(16,isp) = 5.970d-23
          xct%sigma_a_EUV(17,isp) = 5.876d-23
          xct%sigma_a_EUV(18,isp) = 6.100d-23
          xct%sigma_a_EUV(19,isp) = 5.862d-23
          xct%sigma_a_EUV(20,isp) = 6.493d-23
          xct%sigma_a_EUV(21,isp) = 6.956d-23
          xct%sigma_a_EUV(22,isp) = 6.318d-23
          xct%sigma_a_EUV(23,isp) = 7.580d-23
          xct%sigma_a_EUV(24,isp) = 1.008d-22
          xct%sigma_a_EUV(25,isp) = 9.888d-23
          xct%sigma_a_EUV(26,isp) = 1.484d-22
          xct%sigma_a_EUV(27,isp) = 1.535d-22
          xct%sigma_a_EUV(28,isp) = 1.734d-22
          xct%sigma_a_EUV(29,isp) = 1.347d-22
          xct%sigma_a_EUV(30,isp) = 1.856d-22
          xct%sigma_a_EUV(31,isp) = 2.538d-22
          xct%sigma_a_EUV(32,isp) = 3.413d-22
          xct%sigma_a_EUV(33,isp) = 5.125d-22
          xct%sigma_a_EUV(34,isp) = 4.470d-22
          xct%sigma_a_EUV(35,isp) = 6.456d-22
          xct%sigma_a_EUV(36,isp) = 6.640d-22
          xct%sigma_a_EUV(37,isp) = 5.728d-22
        end if

      end do!isp

    end if


    !----------------------------------------------------------------------------------------------------
    ! Photoionization cross sections [m^2]

    if (xct%type == 'ionization') then

      do ich = 1, spl%nch
        if (spl%reaction_type_char(ich) == 'photoionization') then
          do i = 1, spl%reactant_list(ich,0)
            reactants(i) = trim(spl%species(spl%reactant_list(ich,i)))
          end do
          do i = 1, spl%product_list(ich,0)
            products(i) = trim(spl%species(spl%product_list(ich,i)))
          end do
          isp = spl%reactant_list(ich,1)

          if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'CO2+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO2 + hv -> CO2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.447d-22
            xct%sigma_i_EUV(2,ich)  =  2.083d-22
            xct%sigma_i_EUV(3,ich)  =  4.960d-22
            xct%sigma_i_EUV(4,ich)  =  8.515d-22
            xct%sigma_i_EUV(5,ich)  = 11.113d-22
            xct%sigma_i_EUV(6,ich)  = 13.004d-22
            xct%sigma_i_EUV(7,ich)  = 11.906d-22
            xct%sigma_i_EUV(8,ich)  = 14.390d-22
            xct%sigma_i_EUV(9,ich)  = 14.414d-22
            xct%sigma_i_EUV(10,ich) = 15.954d-22
            xct%sigma_i_EUV(11,ich) = 18.271d-22
            xct%sigma_i_EUV(12,ich) = 17.982d-22
            xct%sigma_i_EUV(13,ich) = 21.082d-22
            xct%sigma_i_EUV(14,ich) = 24.378d-22
            xct%sigma_i_EUV(15,ich) = 27.163d-22
            xct%sigma_i_EUV(16,ich) = 30.138d-22
            xct%sigma_i_EUV(17,ich) = 31.451d-22
            xct%sigma_i_EUV(18,ich) = 32.382d-22
            xct%sigma_i_EUV(19,ich) = 33.482d-22
            xct%sigma_i_EUV(20,ich) = 34.318d-22
            xct%sigma_i_EUV(21,ich) = 33.795d-22
            xct%sigma_i_EUV(22,ich) = 34.003d-22
            xct%sigma_i_EUV(23,ich) = 32.287d-22
            xct%sigma_i_EUV(24,ich) = 20.856d-22
            xct%sigma_i_EUV(25,ich) = 27.490d-22
            xct%sigma_i_EUV(26,ich) = 86.317d-22
            xct%sigma_i_EUV(27,ich) = 51.765d-22
            xct%sigma_i_EUV(28,ich) = 21.676d-22
            xct%sigma_i_EUV(29,ich) = 34.094d-22
            xct%sigma_i_EUV(30,ich) = 10.930d-22
            xct%sigma_i_EUV(31,ich) =  7.135d-22
            xct%sigma_i_EUV(32,ich) =  0.000d-22
            xct%sigma_i_EUV(33,ich) =  0.000d-22
            xct%sigma_i_EUV(34,ich) =  0.000d-22
            xct%sigma_i_EUV(35,ich) =  0.000d-22
            xct%sigma_i_EUV(36,ich) =  0.000d-22
            xct%sigma_i_EUV(37,ich) =  0.000d-22
          end if

          if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'CO+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO2 + hv -> CO+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.163d-22
            xct%sigma_i_EUV(2,ich)  = 0.510d-22
            xct%sigma_i_EUV(3,ich)  = 1.052d-22
            xct%sigma_i_EUV(4,ich)  = 1.618d-22
            xct%sigma_i_EUV(5,ich)  = 1.467d-22
            xct%sigma_i_EUV(6,ich)  = 1.640d-22
            xct%sigma_i_EUV(7,ich)  = 1.539d-22
            xct%sigma_i_EUV(8,ich)  = 1.959d-22
            xct%sigma_i_EUV(9,ich)  = 1.968d-22
            xct%sigma_i_EUV(10,ich) = 2.442d-22
            xct%sigma_i_EUV(11,ich) = 3.040d-22
            xct%sigma_i_EUV(12,ich) = 2.995d-22
            xct%sigma_i_EUV(13,ich) = 3.369d-22
            xct%sigma_i_EUV(14,ich) = 2.247d-22
            xct%sigma_i_EUV(15,ich) = 1.504d-22
            xct%sigma_i_EUV(16,ich) = 0.820d-22
            xct%sigma_i_EUV(17,ich) = 0.409d-22
            xct%sigma_i_EUV(18,ich) = 0.305d-22
            xct%sigma_i_EUV(19,ich) = 0.306d-22
            xct%sigma_i_EUV(20,ich) = 0.135d-22
            xct%sigma_i_EUV(21,ich) = 0.037d-22
            xct%sigma_i_EUV(22,ich) = 0.043d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'O+(4S)' .and. products(2) == 'e-' ) then
            !------------------------------------  CO2 + hv -> O+(4S) + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.626d-22
            xct%sigma_i_EUV(2,ich)  = 1.320d-22
            xct%sigma_i_EUV(3,ich)  = 1.929d-22
            xct%sigma_i_EUV(4,ich)  = 2.622d-22
            xct%sigma_i_EUV(5,ich)  = 2.260d-22
            xct%sigma_i_EUV(6,ich)  = 2.572d-22
            xct%sigma_i_EUV(7,ich)  = 2.382d-22
            xct%sigma_i_EUV(8,ich)  = 3.271d-22
            xct%sigma_i_EUV(9,ich)  = 3.280d-22
            xct%sigma_i_EUV(10,ich) = 3.426d-22
            xct%sigma_i_EUV(11,ich) = 3.128d-22
            xct%sigma_i_EUV(12,ich) = 3.224d-22
            xct%sigma_i_EUV(13,ich) = 2.597d-22
            xct%sigma_i_EUV(14,ich) = 2.130d-22
            xct%sigma_i_EUV(15,ich) = 1.911d-22
            xct%sigma_i_EUV(16,ich) = 1.636d-22
            xct%sigma_i_EUV(17,ich) = 1.351d-22
            xct%sigma_i_EUV(18,ich) = 1.170d-22
            xct%sigma_i_EUV(19,ich) = 1.171d-22
            xct%sigma_i_EUV(20,ich) = 0.850d-22
            xct%sigma_i_EUV(21,ich) = 0.468d-22
            xct%sigma_i_EUV(22,ich) = 0.527d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CO2' &
          & .and. products(1) == 'C+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO2 + hv -> C+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.306d-22
            xct%sigma_i_EUV(2,ich)  = 0.658d-22
            xct%sigma_i_EUV(3,ich)  = 1.033d-22
            xct%sigma_i_EUV(4,ich)  = 1.433d-22
            xct%sigma_i_EUV(5,ich)  = 1.168d-22
            xct%sigma_i_EUV(6,ich)  = 1.287d-22
            xct%sigma_i_EUV(7,ich)  = 1.219d-22
            xct%sigma_i_EUV(8,ich)  = 1.706d-22
            xct%sigma_i_EUV(9,ich)  = 1.715d-22
            xct%sigma_i_EUV(10,ich) = 1.794d-22
            xct%sigma_i_EUV(11,ich) = 1.104d-22
            xct%sigma_i_EUV(12,ich) = 1.310d-22
            xct%sigma_i_EUV(13,ich) = 0.124d-22
            xct%sigma_i_EUV(14,ich) = 0.000d-22
            xct%sigma_i_EUV(15,ich) = 0.000d-22
            xct%sigma_i_EUV(16,ich) = 0.000d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CO' &
          & .and. products(1) == 'CO+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO + hv -> CO+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.291d-22
            xct%sigma_i_EUV(2,ich)  =  1.074d-22
            xct%sigma_i_EUV(3,ich)  =  2.459d-22
            xct%sigma_i_EUV(4,ich)  =  4.082d-22
            xct%sigma_i_EUV(5,ich)  =  5.449d-22
            xct%sigma_i_EUV(6,ich)  =  7.713d-22
            xct%sigma_i_EUV(7,ich)  =  6.361d-22
            xct%sigma_i_EUV(8,ich)  =  9.209d-22
            xct%sigma_i_EUV(9,ich)  =  9.246d-22
            xct%sigma_i_EUV(10,ich) = 11.532d-22
            xct%sigma_i_EUV(11,ich) = 13.980d-22
            xct%sigma_i_EUV(12,ich) = 13.609d-22
            xct%sigma_i_EUV(13,ich) = 16.876d-22
            xct%sigma_i_EUV(14,ich) = 19.085d-22
            xct%sigma_i_EUV(15,ich) = 19.669d-22
            xct%sigma_i_EUV(16,ich) = 20.454d-22
            xct%sigma_i_EUV(17,ich) = 21.565d-22
            xct%sigma_i_EUV(18,ich) = 22.000d-22
            xct%sigma_i_EUV(19,ich) = 21.895d-22
            xct%sigma_i_EUV(20,ich) = 21.918d-22
            xct%sigma_i_EUV(21,ich) = 22.025d-22
            xct%sigma_i_EUV(22,ich) = 21.845d-22
            xct%sigma_i_EUV(23,ich) = 20.097d-22
            xct%sigma_i_EUV(24,ich) = 22.115d-22
            xct%sigma_i_EUV(25,ich) = 21.084d-22
            xct%sigma_i_EUV(26,ich) = 13.033d-22
            xct%sigma_i_EUV(27,ich) =  9.884d-22
            xct%sigma_i_EUV(28,ich) = 17.350d-22
            xct%sigma_i_EUV(29,ich) = 11.375d-22
            xct%sigma_i_EUV(30,ich) = 17.559d-22
            xct%sigma_i_EUV(31,ich) = 11.701d-22
            xct%sigma_i_EUV(32,ich) =  0.000d-22
            xct%sigma_i_EUV(33,ich) =  0.000d-22
            xct%sigma_i_EUV(34,ich) =  0.000d-22
            xct%sigma_i_EUV(35,ich) =  0.000d-22
            xct%sigma_i_EUV(36,ich) =  0.000d-22
            xct%sigma_i_EUV(37,ich) =  0.000d-22
          end if

          if ( reactants(1) == 'CO' &
          & .and. products(1) == 'C+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO + hv -> C+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.282d-22
            xct%sigma_i_EUV(2,ich)  = 0.672d-22
            xct%sigma_i_EUV(3,ich)  = 1.156d-22
            xct%sigma_i_EUV(4,ich)  = 1.514d-22
            xct%sigma_i_EUV(5,ich)  = 1.593d-22
            xct%sigma_i_EUV(6,ich)  = 1.141d-22
            xct%sigma_i_EUV(7,ich)  = 1.502d-22
            xct%sigma_i_EUV(8,ich)  = 1.076d-22
            xct%sigma_i_EUV(9,ich)  = 1.073d-22
            xct%sigma_i_EUV(10,ich) = 0.963d-22
            xct%sigma_i_EUV(11,ich) = 0.771d-22
            xct%sigma_i_EUV(12,ich) = 0.814d-22
            xct%sigma_i_EUV(13,ich) = 0.962d-22
            xct%sigma_i_EUV(14,ich) = 1.029d-22
            xct%sigma_i_EUV(15,ich) = 0.895d-22
            xct%sigma_i_EUV(16,ich) = 0.631d-22
            xct%sigma_i_EUV(17,ich) = 0.060d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CO' &
          & .and. products(1) == 'O+' .and. products(2) == 'e-' ) then
            !------------------------------------  CO + hv -> O+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.247d-22
            xct%sigma_i_EUV(2,ich)  = 0.600d-22
            xct%sigma_i_EUV(3,ich)  = 1.029d-22
            xct%sigma_i_EUV(4,ich)  = 1.411d-22
            xct%sigma_i_EUV(5,ich)  = 1.572d-22
            xct%sigma_i_EUV(6,ich)  = 1.687d-22
            xct%sigma_i_EUV(7,ich)  = 1.561d-22
            xct%sigma_i_EUV(8,ich)  = 1.582d-22
            xct%sigma_i_EUV(9,ich)  = 1.581d-22
            xct%sigma_i_EUV(10,ich) = 0.946d-22
            xct%sigma_i_EUV(11,ich) = 0.509d-22
            xct%sigma_i_EUV(12,ich) = 0.533d-22
            xct%sigma_i_EUV(13,ich) = 0.118d-22
            xct%sigma_i_EUV(14,ich) = 0.058d-22
            xct%sigma_i_EUV(15,ich) = 0.009d-22
            xct%sigma_i_EUV(16,ich) = 0.000d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'O2' &
          & .and. products(1) == 'O2+' .and. products(2) == 'e-' ) then
            !------------------------------------  O2 + hv -> O2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  1.316d-22
            xct%sigma_i_EUV(2,ich)  =  2.346d-22
            xct%sigma_i_EUV(3,ich)  =  4.139d-22
            xct%sigma_i_EUV(4,ich)  =  6.619d-22
            xct%sigma_i_EUV(5,ich)  =  8.460d-22
            xct%sigma_i_EUV(6,ich)  =  9.890d-22
            xct%sigma_i_EUV(7,ich)  =  9.056d-22
            xct%sigma_i_EUV(8,ich)  = 10.860d-22
            xct%sigma_i_EUV(9,ich)  = 10.880d-22
            xct%sigma_i_EUV(10,ich) = 12.229d-22
            xct%sigma_i_EUV(11,ich) = 13.760d-22
            xct%sigma_i_EUV(12,ich) = 13.418d-22
            xct%sigma_i_EUV(13,ich) = 15.490d-22
            xct%sigma_i_EUV(14,ich) = 16.970d-22
            xct%sigma_i_EUV(15,ich) = 17.754d-22
            xct%sigma_i_EUV(16,ich) = 19.469d-22
            xct%sigma_i_EUV(17,ich) = 21.600d-22
            xct%sigma_i_EUV(18,ich) = 18.840d-22
            xct%sigma_i_EUV(19,ich) = 22.789d-22
            xct%sigma_i_EUV(20,ich) = 24.540d-22
            xct%sigma_i_EUV(21,ich) = 30.070d-22
            xct%sigma_i_EUV(22,ich) = 23.974d-22
            xct%sigma_i_EUV(23,ich) = 21.116d-22
            xct%sigma_i_EUV(24,ich) = 23.750d-22
            xct%sigma_i_EUV(25,ich) = 23.805d-22
            xct%sigma_i_EUV(26,ich) = 11.720d-22
            xct%sigma_i_EUV(27,ich) =  8.470d-22
            xct%sigma_i_EUV(28,ich) = 10.191d-22
            xct%sigma_i_EUV(29,ich) = 10.597d-22
            xct%sigma_i_EUV(30,ich) =  6.413d-22
            xct%sigma_i_EUV(31,ich) =  5.494d-22
            xct%sigma_i_EUV(32,ich) =  9.374d-22
            xct%sigma_i_EUV(33,ich) = 15.540d-22
            xct%sigma_i_EUV(34,ich) = 13.940d-22
            xct%sigma_i_EUV(35,ich) =  1.050d-22
            xct%sigma_i_EUV(36,ich) =  0.000d-22
            xct%sigma_i_EUV(37,ich) =  0.259d-22
          end if

          if ( reactants(1) == 'O' &
          & .and. products(1) == 'O+(4S)' .and. products(2) == 'e-' ) then ! 4S
            !------------------------------------  O + hv -> O+(4S) + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.190d-22
            xct%sigma_i_EUV(2,ich)  = 0.486d-22
            xct%sigma_i_EUV(3,ich)  = 0.952d-22
            xct%sigma_i_EUV(4,ich)  = 1.311d-22
            xct%sigma_i_EUV(5,ich)  = 1.539d-22
            xct%sigma_i_EUV(6,ich)  = 1.770d-22
            xct%sigma_i_EUV(7,ich)  = 1.628d-22
            xct%sigma_i_EUV(8,ich)  = 1.920d-22
            xct%sigma_i_EUV(9,ich)  = 1.925d-22
            xct%sigma_i_EUV(10,ich) = 2.259d-22
            xct%sigma_i_EUV(11,ich) = 2.559d-22
            xct%sigma_i_EUV(12,ich) = 2.523d-22
            xct%sigma_i_EUV(13,ich) = 3.073d-22
            xct%sigma_i_EUV(14,ich) = 3.340d-22
            xct%sigma_i_EUV(15,ich) = 3.394d-22
            xct%sigma_i_EUV(16,ich) = 3.421d-22
            xct%sigma_i_EUV(17,ich) = 3.650d-22
            xct%sigma_i_EUV(18,ich) = 3.920d-22
            xct%sigma_i_EUV(19,ich) = 3.620d-22
            xct%sigma_i_EUV(20,ich) = 3.610d-22
            xct%sigma_i_EUV(21,ich) = 3.880d-22
            xct%sigma_i_EUV(22,ich) = 4.250d-22
            xct%sigma_i_EUV(23,ich) = 5.128d-22
            xct%sigma_i_EUV(24,ich) = 4.890d-22
            xct%sigma_i_EUV(25,ich) = 6.739d-22
            xct%sigma_i_EUV(26,ich) = 4.000d-22
            xct%sigma_i_EUV(27,ich) = 3.890d-22
            xct%sigma_i_EUV(28,ich) = 3.749d-22
            xct%sigma_i_EUV(29,ich) = 5.091d-22
            xct%sigma_i_EUV(30,ich) = 3.498d-22
            xct%sigma_i_EUV(31,ich) = 4.554d-22
            xct%sigma_i_EUV(32,ich) = 1.315d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'O' &
          & .and. products(1) == 'O+(2D)' .and. products(2) == 'e-' ) then
            !------------------------------------  O + hv -> O+(2D) + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.206d-22
            xct%sigma_i_EUV(2,ich)  =  0.529d-22
            xct%sigma_i_EUV(3,ich)  =  1.171d-22
            xct%sigma_i_EUV(4,ich)  =  1.762d-22
            xct%sigma_i_EUV(5,ich)  =  2.138d-22
            xct%sigma_i_EUV(6,ich)  =  2.620d-22
            xct%sigma_i_EUV(7,ich)  =  2.325d-22
            xct%sigma_i_EUV(8,ich)  =  2.842d-22
            xct%sigma_i_EUV(9,ich)  =  2.849d-22
            xct%sigma_i_EUV(10,ich) =  3.446d-22
            xct%sigma_i_EUV(11,ich) =  3.936d-22
            xct%sigma_i_EUV(12,ich) =  3.883d-22
            xct%sigma_i_EUV(13,ich) =  4.896d-22
            xct%sigma_i_EUV(14,ich) =  5.370d-22
            xct%sigma_i_EUV(15,ich) =  5.459d-22
            xct%sigma_i_EUV(16,ich) =  5.427d-22
            xct%sigma_i_EUV(17,ich) =  5.670d-22
            xct%sigma_i_EUV(18,ich) =  6.020d-22
            xct%sigma_i_EUV(19,ich) =  5.910d-22
            xct%sigma_i_EUV(20,ich) =  6.170d-22
            xct%sigma_i_EUV(21,ich) =  6.290d-22
            xct%sigma_i_EUV(22,ich) =  6.159d-22
            xct%sigma_i_EUV(23,ich) = 11.453d-22
            xct%sigma_i_EUV(24,ich) =  6.570d-22
            xct%sigma_i_EUV(25,ich) =  3.997d-22
            xct%sigma_i_EUV(26,ich) =  0.000d-22
            xct%sigma_i_EUV(27,ich) =  0.000d-22
            xct%sigma_i_EUV(28,ich) =  0.000d-22
            xct%sigma_i_EUV(29,ich) =  0.000d-22
            xct%sigma_i_EUV(30,ich) =  0.000d-22
            xct%sigma_i_EUV(31,ich) =  0.000d-22
            xct%sigma_i_EUV(32,ich) =  0.000d-22
            xct%sigma_i_EUV(33,ich) =  0.000d-22
            xct%sigma_i_EUV(34,ich) =  0.000d-22
            xct%sigma_i_EUV(35,ich) =  0.000d-22
            xct%sigma_i_EUV(36,ich) =  0.000d-22
            xct%sigma_i_EUV(37,ich) =  0.000d-22
          end if

          if ( reactants(1) == 'O' &
          & .and. products(1) == 'O+(2P)' .and. products(2) == 'e-' ) then
            !------------------------------------  O + hv -> O+(2P) + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.134d-22
            xct%sigma_i_EUV(2,ich)  = 0.345d-22
            xct%sigma_i_EUV(3,ich)  = 0.768d-22
            xct%sigma_i_EUV(4,ich)  = 1.144d-22
            xct%sigma_i_EUV(5,ich)  = 1.363d-22
            xct%sigma_i_EUV(6,ich)  = 1.630d-22
            xct%sigma_i_EUV(7,ich)  = 1.488d-22
            xct%sigma_i_EUV(8,ich)  = 1.920d-22
            xct%sigma_i_EUV(9,ich)  = 1.925d-22
            xct%sigma_i_EUV(10,ich) = 2.173d-22
            xct%sigma_i_EUV(11,ich) = 2.558d-22
            xct%sigma_i_EUV(12,ich) = 2.422d-22
            xct%sigma_i_EUV(13,ich) = 2.986d-22
            xct%sigma_i_EUV(14,ich) = 3.220d-22
            xct%sigma_i_EUV(15,ich) = 3.274d-22
            xct%sigma_i_EUV(16,ich) = 3.211d-22
            xct%sigma_i_EUV(17,ich) = 3.270d-22
            xct%sigma_i_EUV(18,ich) = 3.150d-22
            xct%sigma_i_EUV(19,ich) = 3.494d-22
            xct%sigma_i_EUV(20,ich) = 3.620d-22
            xct%sigma_i_EUV(21,ich) = 3.230d-22
            xct%sigma_i_EUV(22,ich) = 2.956d-22
            xct%sigma_i_EUV(23,ich) = 0.664d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'O' &
          & .and. products(1) == 'O+(4P)' .and. products(2) == 'e-' ) then
            !------------------------------------  O + hv -> O+(4P) + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.062d-22
            xct%sigma_i_EUV(2,ich)  = 0.163d-22
            xct%sigma_i_EUV(3,ich)  = 0.348d-22
            xct%sigma_i_EUV(4,ich)  = 0.508d-22
            xct%sigma_i_EUV(5,ich)  = 0.598d-22
            xct%sigma_i_EUV(6,ich)  = 0.710d-22
            xct%sigma_i_EUV(7,ich)  = 0.637d-22
            xct%sigma_i_EUV(8,ich)  = 0.691d-22
            xct%sigma_i_EUV(9,ich)  = 0.693d-22
            xct%sigma_i_EUV(10,ich) = 0.815d-22
            xct%sigma_i_EUV(11,ich) = 0.787d-22
            xct%sigma_i_EUV(12,ich) = 0.859d-22
            xct%sigma_i_EUV(13,ich) = 0.541d-22
            xct%sigma_i_EUV(14,ich) = 0.000d-22
            xct%sigma_i_EUV(15,ich) = 0.000d-22
            xct%sigma_i_EUV(16,ich) = 0.000d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'N2' &
          & .and. products(1) == 'N2+' .and. products(2) == 'e-' ) then
            !------------------------------------  N2 + hv -> N2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.443d-22
            xct%sigma_i_EUV(2,ich)  =  1.479d-22
            xct%sigma_i_EUV(3,ich)  =  3.153d-22
            xct%sigma_i_EUV(4,ich)  =  5.226d-22
            xct%sigma_i_EUV(5,ich)  =  6.781d-22
            xct%sigma_i_EUV(6,ich)  =  8.100d-22
            xct%sigma_i_EUV(7,ich)  =  7.347d-22
            xct%sigma_i_EUV(8,ich)  =  9.180d-22
            xct%sigma_i_EUV(9,ich)  =  9.210d-22
            xct%sigma_i_EUV(10,ich) = 11.600d-22
            xct%sigma_i_EUV(11,ich) = 15.350d-22
            xct%sigma_i_EUV(12,ich) = 14.669d-22
            xct%sigma_i_EUV(13,ich) = 20.692d-22
            xct%sigma_i_EUV(14,ich) = 22.100d-22
            xct%sigma_i_EUV(15,ich) = 22.772d-22
            xct%sigma_i_EUV(16,ich) = 24.468d-22
            xct%sigma_i_EUV(17,ich) = 24.130d-22
            xct%sigma_i_EUV(18,ich) = 22.400d-22
            xct%sigma_i_EUV(19,ich) = 22.787d-22
            xct%sigma_i_EUV(20,ich) = 22.790d-22
            xct%sigma_i_EUV(21,ich) = 23.370d-22
            xct%sigma_i_EUV(22,ich) = 23.339d-22
            xct%sigma_i_EUV(23,ich) = 29.235d-22
            xct%sigma_i_EUV(24,ich) = 25.480d-22
            xct%sigma_i_EUV(25,ich) = 15.060d-22
            xct%sigma_i_EUV(26,ich) = 65.800d-22
            xct%sigma_i_EUV(27,ich) =  8.500d-22
            xct%sigma_i_EUV(28,ich) =  8.860d-22
            xct%sigma_i_EUV(29,ich) = 14.274d-22
            xct%sigma_i_EUV(30,ich) =  0.000d-22
            xct%sigma_i_EUV(31,ich) =  0.000d-22
            xct%sigma_i_EUV(32,ich) =  0.000d-22
            xct%sigma_i_EUV(33,ich) =  0.000d-22
            xct%sigma_i_EUV(34,ich) =  0.000d-22
            xct%sigma_i_EUV(35,ich) =  0.000d-22
            xct%sigma_i_EUV(36,ich) =  0.000d-22
            xct%sigma_i_EUV(37,ich) =  0.000d-22
          end if

          if ( reactants(1) == 'N2' &
          & .and. products(1) == 'N+' .and. products(2) == 'e-' ) then
            !------------------------------------  N2 + hv -> N+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.277d-22
            xct%sigma_i_EUV(2,ich)  = 0.782d-22
            xct%sigma_i_EUV(3,ich)  = 1.805d-22
            xct%sigma_i_EUV(4,ich)  = 3.166d-22
            xct%sigma_i_EUV(5,ich)  = 3.420d-22
            xct%sigma_i_EUV(6,ich)  = 2.800d-22
            xct%sigma_i_EUV(7,ich)  = 3.145d-22
            xct%sigma_i_EUV(8,ich)  = 2.490d-22
            xct%sigma_i_EUV(9,ich)  = 2.490d-22
            xct%sigma_i_EUV(10,ich) = 2.257d-22
            xct%sigma_i_EUV(11,ich) = 1.560d-22
            xct%sigma_i_EUV(12,ich) = 1.726d-22
            xct%sigma_i_EUV(13,ich) = 0.982d-22
            xct%sigma_i_EUV(14,ich) = 1.060d-22
            xct%sigma_i_EUV(15,ich) = 0.699d-22
            xct%sigma_i_EUV(16,ich) = 0.033d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'H' &
          & .and. products(1) == 'H+' .and. products(2) == 'e-' ) then
            !------------------------------------  H + hv -> H+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.0024d-22
            xct%sigma_i_EUV(2,ich)  = 0.0169d-22
            xct%sigma_i_EUV(3,ich)  = 0.0483d-22
            xct%sigma_i_EUV(4,ich)  = 0.1007d-22
            xct%sigma_i_EUV(5,ich)  = 0.1405d-22
            xct%sigma_i_EUV(6,ich)  = 0.1913d-22
            xct%sigma_i_EUV(7,ich)  = 0.1676d-22
            xct%sigma_i_EUV(8,ich)  = 0.2324d-22
            xct%sigma_i_EUV(9,ich)  = 0.2334d-22
            xct%sigma_i_EUV(10,ich) = 0.3077d-22
            xct%sigma_i_EUV(11,ich) = 0.4152d-22
            xct%sigma_i_EUV(12,ich) = 0.3984d-22
            xct%sigma_i_EUV(13,ich) = 0.6163d-22
            xct%sigma_i_EUV(14,ich) = 0.8387d-22
            xct%sigma_i_EUV(15,ich) = 0.9739d-22
            xct%sigma_i_EUV(16,ich) = 1.1990d-22
            xct%sigma_i_EUV(17,ich) = 1.4190d-22
            xct%sigma_i_EUV(18,ich) = 1.6620d-22
            xct%sigma_i_EUV(19,ich) = 1.6200d-22
            xct%sigma_i_EUV(20,ich) = 1.8880d-22
            xct%sigma_i_EUV(21,ich) = 2.0790d-22
            xct%sigma_i_EUV(22,ich) = 2.0760d-22
            xct%sigma_i_EUV(23,ich) = 2.6410d-22
            xct%sigma_i_EUV(24,ich) = 2.8970d-22
            xct%sigma_i_EUV(25,ich) = 3.1730d-22
            xct%sigma_i_EUV(26,ich) = 3.7300d-22
            xct%sigma_i_EUV(27,ich) = 3.8070d-22
            xct%sigma_i_EUV(28,ich) = 4.0930d-22
            xct%sigma_i_EUV(29,ich) = 3.8680d-22
            xct%sigma_i_EUV(30,ich) = 4.7840d-22
            xct%sigma_i_EUV(31,ich) = 5.6700d-22
            xct%sigma_i_EUV(32,ich) = 3.4690d-22
            xct%sigma_i_EUV(33,ich) = 0.0000d-22
            xct%sigma_i_EUV(34,ich) = 0.0000d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'H2' &
          & .and. products(1) == 'H2+' .and. products(2) == 'e-' ) then
            !------------------------------------  H2 + hv -> H2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.0097d-22
            xct%sigma_i_EUV(2,ich)  = 0.0758d-22
            xct%sigma_i_EUV(3,ich)  = 0.2009d-22
            xct%sigma_i_EUV(4,ich)  = 0.4028d-22
            xct%sigma_i_EUV(5,ich)  = 0.5509d-22
            xct%sigma_i_EUV(6,ich)  = 0.7454d-22
            xct%sigma_i_EUV(7,ich)  = 0.6538d-22
            xct%sigma_i_EUV(8,ich)  = 0.8999d-22
            xct%sigma_i_EUV(9,ich)  = 0.9041d-22
            xct%sigma_i_EUV(10,ich) = 1.2960d-22
            xct%sigma_i_EUV(11,ich) = 1.7840d-22
            xct%sigma_i_EUV(12,ich) = 1.7420d-22
            xct%sigma_i_EUV(13,ich) = 2.8900d-22
            xct%sigma_i_EUV(14,ich) = 3.7780d-22
            xct%sigma_i_EUV(15,ich) = 4.0470d-22
            xct%sigma_i_EUV(16,ich) = 5.2540d-22
            xct%sigma_i_EUV(17,ich) = 6.0500d-22
            xct%sigma_i_EUV(18,ich) = 6.9000d-22
            xct%sigma_i_EUV(19,ich) = 6.7410d-22
            xct%sigma_i_EUV(20,ich) = 7.6680d-22
            xct%sigma_i_EUV(21,ich) = 8.2990d-22
            xct%sigma_i_EUV(22,ich) = 8.2880d-22
            xct%sigma_i_EUV(23,ich) = 9.7020d-22
            xct%sigma_i_EUV(24,ich) = 10.731d-22
            xct%sigma_i_EUV(25,ich) = 9.7610d-22
            xct%sigma_i_EUV(26,ich) = 8.6240d-22
            xct%sigma_i_EUV(27,ich) = 7.0710d-22
            xct%sigma_i_EUV(28,ich) = 5.0720d-22
            xct%sigma_i_EUV(29,ich) = 6.6290d-22
            xct%sigma_i_EUV(30,ich) = 0.0889d-22
            xct%sigma_i_EUV(31,ich) = 0.0000d-22
            xct%sigma_i_EUV(32,ich) = 0.0000d-22
            xct%sigma_i_EUV(33,ich) = 0.0000d-22
            xct%sigma_i_EUV(34,ich) = 0.0000d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'H2' &
          & .and. products(1) == 'H+' .and. products(2) == 'e-' ) then
            !------------------------------------  H2 + hv -> H+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.0011d-22
            xct%sigma_i_EUV(2,ich)  = 0.0040d-22
            xct%sigma_i_EUV(3,ich)  = 0.0075d-22
            xct%sigma_i_EUV(4,ich)  = 0.0305d-22
            xct%sigma_i_EUV(5,ich)  = 0.0527d-22
            xct%sigma_i_EUV(6,ich)  = 0.0773d-22
            xct%sigma_i_EUV(7,ich)  = 0.0661d-22
            xct%sigma_i_EUV(8,ich)  = 0.1005d-22
            xct%sigma_i_EUV(9,ich)  = 0.1011d-22
            xct%sigma_i_EUV(10,ich) = 0.1200d-22
            xct%sigma_i_EUV(11,ich) = 0.1577d-22
            xct%sigma_i_EUV(12,ich) = 0.1594d-22
            xct%sigma_i_EUV(13,ich) = 0.1255d-22
            xct%sigma_i_EUV(14,ich) = 0.0925d-22
            xct%sigma_i_EUV(15,ich) = 0.0944d-22
            xct%sigma_i_EUV(16,ich) = 0.1020d-22
            xct%sigma_i_EUV(17,ich) = 0.1184d-22
            xct%sigma_i_EUV(18,ich) = 0.1208d-22
            xct%sigma_i_EUV(19,ich) = 0.1237d-22
            xct%sigma_i_EUV(20,ich) = 0.1429d-22
            xct%sigma_i_EUV(21,ich) = 0.1573d-22
            xct%sigma_i_EUV(22,ich) = 0.1524d-22
            xct%sigma_i_EUV(23,ich) = 0.0287d-22
            xct%sigma_i_EUV(24,ich) = 0.0000d-22
            xct%sigma_i_EUV(25,ich) = 0.0000d-22
            xct%sigma_i_EUV(26,ich) = 0.0000d-22
            xct%sigma_i_EUV(27,ich) = 0.0000d-22
            xct%sigma_i_EUV(28,ich) = 0.0000d-22
            xct%sigma_i_EUV(29,ich) = 0.0000d-22
            xct%sigma_i_EUV(30,ich) = 0.0000d-22
            xct%sigma_i_EUV(31,ich) = 0.0000d-22
            xct%sigma_i_EUV(32,ich) = 0.0000d-22
            xct%sigma_i_EUV(33,ich) = 0.0000d-22
            xct%sigma_i_EUV(34,ich) = 0.0000d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'He' &
          & .and. products(1) == 'He+' .and. products(2) == 'e-' ) then
            !------------------------------------  He + hv -> He+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.1441d-22
            xct%sigma_i_EUV(2,ich)  = 0.4785d-22
            xct%sigma_i_EUV(3,ich)  = 1.1571d-22
            xct%sigma_i_EUV(4,ich)  = 1.6008d-22
            xct%sigma_i_EUV(5,ich)  = 2.1212d-22
            xct%sigma_i_EUV(6,ich)  = 2.5947d-22
            xct%sigma_i_EUV(7,ich)  = 2.3205d-22
            xct%sigma_i_EUV(8,ich)  = 2.9529d-22
            xct%sigma_i_EUV(9,ich)  = 2.9618d-22
            xct%sigma_i_EUV(10,ich) = 3.5437d-22
            xct%sigma_i_EUV(11,ich) = 4.2675d-22
            xct%sigma_i_EUV(12,ich) = 4.1424d-22
            xct%sigma_i_EUV(13,ich) = 5.4466d-22
            xct%sigma_i_EUV(14,ich) = 6.5631d-22
            xct%sigma_i_EUV(15,ich) = 7.2084d-22
            xct%sigma_i_EUV(16,ich) = 0.9581d-22
            xct%sigma_i_EUV(17,ich) = 0.0000d-22
            xct%sigma_i_EUV(18,ich) = 0.0000d-22
            xct%sigma_i_EUV(19,ich) = 0.0000d-22
            xct%sigma_i_EUV(20,ich) = 0.0000d-22
            xct%sigma_i_EUV(21,ich) = 0.0000d-22
            xct%sigma_i_EUV(22,ich) = 0.0000d-22
            xct%sigma_i_EUV(23,ich) = 0.0000d-22
            xct%sigma_i_EUV(24,ich) = 0.0000d-22
            xct%sigma_i_EUV(25,ich) = 0.0000d-22
            xct%sigma_i_EUV(26,ich) = 0.0000d-22
            xct%sigma_i_EUV(27,ich) = 0.0000d-22
            xct%sigma_i_EUV(28,ich) = 0.0000d-22
            xct%sigma_i_EUV(29,ich) = 0.0000d-22
            xct%sigma_i_EUV(30,ich) = 0.0000d-22
            xct%sigma_i_EUV(31,ich) = 0.0000d-22
            xct%sigma_i_EUV(32,ich) = 0.0000d-22
            xct%sigma_i_EUV(33,ich) = 0.0000d-22
            xct%sigma_i_EUV(34,ich) = 0.0000d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'CH4' &
          & .and. products(1) == 'CH4+' .and. products(2) == 'e-' ) then
            !------------------------------------  CH4 + hv -> CH4+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.051d-22
            xct%sigma_i_EUV(2,ich)  =  0.147d-22
            xct%sigma_i_EUV(3,ich)  =  0.387d-22
            xct%sigma_i_EUV(4,ich)  =  0.839d-22
            xct%sigma_i_EUV(5,ich)  =  1.192d-22
            xct%sigma_i_EUV(6,ich)  =  1.681d-22
            xct%sigma_i_EUV(7,ich)  =  1.398d-22
            xct%sigma_i_EUV(8,ich)  =  2.095d-22
            xct%sigma_i_EUV(9,ich)  =  2.103d-22
            xct%sigma_i_EUV(10,ich) =  2.957d-22
            xct%sigma_i_EUV(11,ich) =  3.972d-22
            xct%sigma_i_EUV(12,ich) =  3.820d-22
            xct%sigma_i_EUV(13,ich) =  6.255d-22
            xct%sigma_i_EUV(14,ich) =  8.442d-22
            xct%sigma_i_EUV(15,ich) =  9.837d-22
            xct%sigma_i_EUV(16,ich) = 11.432d-22
            xct%sigma_i_EUV(17,ich) = 13.398d-22
            xct%sigma_i_EUV(18,ich) = 14.801d-22
            xct%sigma_i_EUV(19,ich) = 14.640d-22
            xct%sigma_i_EUV(20,ich) = 15.734d-22
            xct%sigma_i_EUV(21,ich) = 17.102d-22
            xct%sigma_i_EUV(22,ich) = 16.883d-22
            xct%sigma_i_EUV(23,ich) = 19.261d-22
            xct%sigma_i_EUV(24,ich) = 20.222d-22
            xct%sigma_i_EUV(25,ich) = 21.314d-22
            xct%sigma_i_EUV(26,ich) = 22.599d-22
            xct%sigma_i_EUV(27,ich) = 22.763d-22
            xct%sigma_i_EUV(28,ich) = 23.198d-22
            xct%sigma_i_EUV(29,ich) = 22.886d-22
            xct%sigma_i_EUV(30,ich) = 25.607d-22
            xct%sigma_i_EUV(31,ich) = 24.233d-22
            xct%sigma_i_EUV(32,ich) = 13.863d-22
            xct%sigma_i_EUV(33,ich) =  0.136d-22
            xct%sigma_i_EUV(34,ich) =  0.475d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'CH4' &
          & .and. products(1) == 'CH3+' .and. products(2) == 'e-' ) then
            !------------------------------------  CH4 + hv -> CH3+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.052d-22
            xct%sigma_i_EUV(2,ich)  =  0.152d-22
            xct%sigma_i_EUV(3,ich)  =  0.409d-22
            xct%sigma_i_EUV(4,ich)  =  0.884d-22
            xct%sigma_i_EUV(5,ich)  =  1.290d-22
            xct%sigma_i_EUV(6,ich)  =  1.824d-22
            xct%sigma_i_EUV(7,ich)  =  1.514d-22
            xct%sigma_i_EUV(8,ich)  =  2.287d-22
            xct%sigma_i_EUV(9,ich)  =  2.302d-22
            xct%sigma_i_EUV(10,ich) =  3.108d-22
            xct%sigma_i_EUV(11,ich) =  4.305d-22
            xct%sigma_i_EUV(12,ich) =  4.101d-22
            xct%sigma_i_EUV(13,ich) =  6.573d-22
            xct%sigma_i_EUV(14,ich) =  8.776d-22
            xct%sigma_i_EUV(15,ich) = 10.212d-22
            xct%sigma_i_EUV(16,ich) = 11.974d-22
            xct%sigma_i_EUV(17,ich) = 13.853d-22
            xct%sigma_i_EUV(18,ich) = 15.501d-22
            xct%sigma_i_EUV(19,ich) = 15.374d-22
            xct%sigma_i_EUV(20,ich) = 16.719d-22
            xct%sigma_i_EUV(21,ich) = 17.494d-22
            xct%sigma_i_EUV(22,ich) = 17.422d-22
            xct%sigma_i_EUV(23,ich) = 19.266d-22
            xct%sigma_i_EUV(24,ich) = 20.092d-22
            xct%sigma_i_EUV(25,ich) = 20.850d-22
            xct%sigma_i_EUV(26,ich) = 21.436d-22
            xct%sigma_i_EUV(27,ich) = 21.316d-22
            xct%sigma_i_EUV(28,ich) = 20.899d-22
            xct%sigma_i_EUV(29,ich) = 21.145d-22
            xct%sigma_i_EUV(30,ich) = 14.651d-22
            xct%sigma_i_EUV(31,ich) =  1.294d-22
            xct%sigma_i_EUV(32,ich) = 0.0000d-22
            xct%sigma_i_EUV(33,ich) = 0.0000d-22
            xct%sigma_i_EUV(34,ich) = 0.0000d-22
            xct%sigma_i_EUV(35,ich) = 0.0000d-22
            xct%sigma_i_EUV(36,ich) = 0.0000d-22
            xct%sigma_i_EUV(37,ich) = 0.0000d-22
          end if

          if ( reactants(1) == 'CH4' &
          & .and. products(1) == 'CH2+' .and. products(2) == 'e-' ) then
            !------------------------------------  CH4 + hv -> CH2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.033d-22
            xct%sigma_i_EUV(2,ich)  = 0.095d-22
            xct%sigma_i_EUV(3,ich)  = 0.201d-22
            xct%sigma_i_EUV(4,ich)  = 0.416d-22
            xct%sigma_i_EUV(5,ich)  = 0.576d-22
            xct%sigma_i_EUV(6,ich)  = 0.665d-22
            xct%sigma_i_EUV(7,ich)  = 0.614d-22
            xct%sigma_i_EUV(8,ich)  = 0.701d-22
            xct%sigma_i_EUV(9,ich)  = 0.701d-22
            xct%sigma_i_EUV(10,ich) = 0.781d-22
            xct%sigma_i_EUV(11,ich) = 0.867d-22
            xct%sigma_i_EUV(12,ich) = 0.852d-22
            xct%sigma_i_EUV(13,ich) = 1.074d-22
            xct%sigma_i_EUV(14,ich) = 1.097d-22
            xct%sigma_i_EUV(15,ich) = 1.014d-22
            xct%sigma_i_EUV(16,ich) = 0.926d-22
            xct%sigma_i_EUV(17,ich) = 0.652d-22
            xct%sigma_i_EUV(18,ich) = 0.750d-22
            xct%sigma_i_EUV(19,ich) = 0.683d-22
            xct%sigma_i_EUV(20,ich) = 0.726d-22
            xct%sigma_i_EUV(21,ich) = 0.680d-22
            xct%sigma_i_EUV(22,ich) = 0.685d-22
            xct%sigma_i_EUV(23,ich) = 0.754d-22
            xct%sigma_i_EUV(24,ich) = 0.755d-22
            xct%sigma_i_EUV(25,ich) = 0.764d-22
            xct%sigma_i_EUV(26,ich) = 0.765d-22
            xct%sigma_i_EUV(27,ich) = 0.717d-22
            xct%sigma_i_EUV(28,ich) = 0.510d-22
            xct%sigma_i_EUV(29,ich) = 0.662d-22
            xct%sigma_i_EUV(30,ich) = 0.025d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CH4' &
          & .and. products(1) == 'H2+' .and. products(2) == 'e-' ) then
            !------------------------------------  CH4 + hv -> H2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.005d-22
            xct%sigma_i_EUV(2,ich)  = 0.015d-22
            xct%sigma_i_EUV(3,ich)  = 0.038d-22
            xct%sigma_i_EUV(4,ich)  = 0.046d-22
            xct%sigma_i_EUV(5,ich)  = 0.058d-22
            xct%sigma_i_EUV(6,ich)  = 0.049d-22
            xct%sigma_i_EUV(7,ich)  = 0.055d-22
            xct%sigma_i_EUV(8,ich)  = 0.052d-22
            xct%sigma_i_EUV(9,ich)  = 0.052d-22
            xct%sigma_i_EUV(10,ich) = 0.055d-22
            xct%sigma_i_EUV(11,ich) = 0.053d-22
            xct%sigma_i_EUV(12,ich) = 0.054d-22
            xct%sigma_i_EUV(13,ich) = 0.018d-22
            xct%sigma_i_EUV(14,ich) = 0.000d-22
            xct%sigma_i_EUV(15,ich) = 0.000d-22
            xct%sigma_i_EUV(16,ich) = 0.000d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'CH4' &
          & .and. products(1) == 'H+' .and. products(2) == 'e-' ) then
            !------------------------------------  CH4 + hv -> H+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.047d-22
            xct%sigma_i_EUV(2,ich)  = 0.137d-22
            xct%sigma_i_EUV(3,ich)  = 0.344d-22
            xct%sigma_i_EUV(4,ich)  = 0.414d-22
            xct%sigma_i_EUV(5,ich)  = 0.492d-22
            xct%sigma_i_EUV(6,ich)  = 0.559d-22
            xct%sigma_i_EUV(7,ich)  = 0.519d-22
            xct%sigma_i_EUV(8,ich)  = 0.559d-22
            xct%sigma_i_EUV(9,ich)  = 0.561d-22
            xct%sigma_i_EUV(10,ich) = 0.552d-22
            xct%sigma_i_EUV(11,ich) = 0.521d-22
            xct%sigma_i_EUV(12,ich) = 0.527d-22
            xct%sigma_i_EUV(13,ich) = 0.362d-22
            xct%sigma_i_EUV(14,ich) = 0.238d-22
            xct%sigma_i_EUV(15,ich) = 0.225d-22
            xct%sigma_i_EUV(16,ich) = 0.181d-22
            xct%sigma_i_EUV(17,ich) = 0.000d-22
            xct%sigma_i_EUV(18,ich) = 0.000d-22
            xct%sigma_i_EUV(19,ich) = 0.000d-22
            xct%sigma_i_EUV(20,ich) = 0.000d-22
            xct%sigma_i_EUV(21,ich) = 0.000d-22
            xct%sigma_i_EUV(22,ich) = 0.000d-22
            xct%sigma_i_EUV(23,ich) = 0.000d-22
            xct%sigma_i_EUV(24,ich) = 0.000d-22
            xct%sigma_i_EUV(25,ich) = 0.000d-22
            xct%sigma_i_EUV(26,ich) = 0.000d-22
            xct%sigma_i_EUV(27,ich) = 0.000d-22
            xct%sigma_i_EUV(28,ich) = 0.000d-22
            xct%sigma_i_EUV(29,ich) = 0.000d-22
            xct%sigma_i_EUV(30,ich) = 0.000d-22
            xct%sigma_i_EUV(31,ich) = 0.000d-22
            xct%sigma_i_EUV(32,ich) = 0.000d-22
            xct%sigma_i_EUV(33,ich) = 0.000d-22
            xct%sigma_i_EUV(34,ich) = 0.000d-22
            xct%sigma_i_EUV(35,ich) = 0.000d-22
            xct%sigma_i_EUV(36,ich) = 0.000d-22
            xct%sigma_i_EUV(37,ich) = 0.000d-22
          end if

          if ( reactants(1) == 'C2H2' &
          & .and. products(1) == 'C2H2+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H2 + hv -> C2H2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.0500d-22
            xct%sigma_i_EUV(2,ich)  =  1.0500d-22
            xct%sigma_i_EUV(3,ich)  =  3.0000d-22
            xct%sigma_i_EUV(4,ich)  =  5.1250d-22
            xct%sigma_i_EUV(5,ich)  =  7.0000d-22
            xct%sigma_i_EUV(6,ich)  =  8.0000d-22
            xct%sigma_i_EUV(7,ich)  =  7.6250d-22
            xct%sigma_i_EUV(8,ich)  =  9.5000d-22
            xct%sigma_i_EUV(9,ich)  =  9.5000d-22
            xct%sigma_i_EUV(10,ich) = 10.7500d-22
            xct%sigma_i_EUV(11,ich) = 14.0000d-22
            xct%sigma_i_EUV(12,ich) = 14.0000d-22
            xct%sigma_i_EUV(13,ich) = 17.7500d-22
            xct%sigma_i_EUV(14,ich) = 22.0000d-22
            xct%sigma_i_EUV(15,ich) = 22.4000d-22
            xct%sigma_i_EUV(16,ich) = 27.4000d-22
            xct%sigma_i_EUV(17,ich) = 31.5000d-22
            xct%sigma_i_EUV(18,ich) = 33.0000d-22
            xct%sigma_i_EUV(19,ich) = 32.6000d-22
            xct%sigma_i_EUV(20,ich) = 38.0000d-22
            xct%sigma_i_EUV(21,ich) = 41.0000d-22
            xct%sigma_i_EUV(22,ich) = 38.6000d-22
            xct%sigma_i_EUV(23,ich) = 42.5000d-22
            xct%sigma_i_EUV(24,ich) = 44.0000d-22
            xct%sigma_i_EUV(25,ich) = 44.2500d-22
            xct%sigma_i_EUV(26,ich) = 47.5000d-22
            xct%sigma_i_EUV(27,ich) = 55.0000d-22
            xct%sigma_i_EUV(28,ich) = 54.0000d-22
            xct%sigma_i_EUV(29,ich) = 51.2500d-22
            xct%sigma_i_EUV(30,ich) = 49.7500d-22
            xct%sigma_i_EUV(31,ich) = 42.2500d-22
            xct%sigma_i_EUV(32,ich) = 38.5000d-22
            xct%sigma_i_EUV(33,ich) = 35.0000d-22
            xct%sigma_i_EUV(34,ich) = 36.2500d-22
            xct%sigma_i_EUV(35,ich) = 32.0000d-22
            xct%sigma_i_EUV(36,ich) = 32.5000d-22
            xct%sigma_i_EUV(37,ich) = 33.1250d-22
          end if

          if ( reactants(1) == 'C2H4' &
          & .and. products(1) == 'C2H4+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H4 + hv -> C2H4+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.00d-22
            xct%sigma_i_EUV(2,ich)  =  0.00d-22
            xct%sigma_i_EUV(3,ich)  =  0.16d-22
            xct%sigma_i_EUV(4,ich)  =  0.69d-22
            xct%sigma_i_EUV(5,ich)  =  1.53d-22
            xct%sigma_i_EUV(6,ich)  =  2.04d-22
            xct%sigma_i_EUV(7,ich)  =  1.42d-22
            xct%sigma_i_EUV(8,ich)  =  2.62d-22
            xct%sigma_i_EUV(9,ich)  =  2.63d-22
            xct%sigma_i_EUV(10,ich) =  2.52d-22
            xct%sigma_i_EUV(11,ich) =  4.61d-22
            xct%sigma_i_EUV(12,ich) =  4.04d-22
            xct%sigma_i_EUV(13,ich) =  6.06d-22
            xct%sigma_i_EUV(14,ich) =  9.35d-22
            xct%sigma_i_EUV(15,ich) =  8.48d-22
            xct%sigma_i_EUV(16,ich) = 10.96d-22
            xct%sigma_i_EUV(17,ich) = 13.53d-22
            xct%sigma_i_EUV(18,ich) = 14.91d-22
            xct%sigma_i_EUV(19,ich) = 13.36d-22
            xct%sigma_i_EUV(20,ich) = 16.01d-22
            xct%sigma_i_EUV(21,ich) = 17.24d-22
            xct%sigma_i_EUV(22,ich) = 15.56d-22
            xct%sigma_i_EUV(23,ich) = 18.18d-22
            xct%sigma_i_EUV(24,ich) = 20.43d-22
            xct%sigma_i_EUV(25,ich) = 20.26d-22
            xct%sigma_i_EUV(26,ich) = 22.08d-22
            xct%sigma_i_EUV(27,ich) = 22.07d-22
            xct%sigma_i_EUV(28,ich) = 22.92d-22
            xct%sigma_i_EUV(29,ich) = 22.06d-22
            xct%sigma_i_EUV(30,ich) = 23.53d-22
            xct%sigma_i_EUV(31,ich) = 25.43d-22
            xct%sigma_i_EUV(32,ich) = 23.81d-22
            xct%sigma_i_EUV(33,ich) = 18.06d-22
            xct%sigma_i_EUV(34,ich) = 20.49d-22
            xct%sigma_i_EUV(35,ich) = 11.81d-22
            xct%sigma_i_EUV(36,ich) = 10.92d-22
            xct%sigma_i_EUV(37,ich) = 15.49d-22
          end if

          if ( reactants(1) == 'C2H4' &
          & .and. products(1) == 'C2H3+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H4 + hv -> C2H3+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.00d-22
            xct%sigma_i_EUV(2,ich)  = 0.00d-22
            xct%sigma_i_EUV(3,ich)  = 0.15d-22
            xct%sigma_i_EUV(4,ich)  = 0.72d-22
            xct%sigma_i_EUV(5,ich)  = 1.57d-22
            xct%sigma_i_EUV(6,ich)  = 2.17d-22
            xct%sigma_i_EUV(7,ich)  = 1.46d-22
            xct%sigma_i_EUV(8,ich)  = 2.61d-22
            xct%sigma_i_EUV(9,ich)  = 2.61d-22
            xct%sigma_i_EUV(10,ich) = 2.56d-22
            xct%sigma_i_EUV(11,ich) = 4.56d-22
            xct%sigma_i_EUV(12,ich) = 3.98d-22
            xct%sigma_i_EUV(13,ich) = 6.03d-22
            xct%sigma_i_EUV(14,ich) = 0.45d-22
            xct%sigma_i_EUV(15,ich) = 9.12d-22
            xct%sigma_i_EUV(16,ich) = 2.26d-22
            xct%sigma_i_EUV(17,ich) = 5.27d-22
            xct%sigma_i_EUV(18,ich) = 6.74d-22
            xct%sigma_i_EUV(19,ich) = 5.03d-22
            xct%sigma_i_EUV(20,ich) = 7.49d-22
            xct%sigma_i_EUV(21,ich) = 8.13d-22
            xct%sigma_i_EUV(22,ich) = 7.24d-22
            xct%sigma_i_EUV(23,ich) = 8.49d-22
            xct%sigma_i_EUV(24,ich) = 9.44d-22
            xct%sigma_i_EUV(25,ich) = 9.46d-22
            xct%sigma_i_EUV(26,ich) = 5.87d-22
            xct%sigma_i_EUV(27,ich) = 5.32d-22
            xct%sigma_i_EUV(28,ich) = 2.75d-22
            xct%sigma_i_EUV(29,ich) = 7.43d-22
            xct%sigma_i_EUV(30,ich) = 1.21d-22
            xct%sigma_i_EUV(31,ich) = 6.07d-22
            xct%sigma_i_EUV(32,ich) = 1.96d-22
            xct%sigma_i_EUV(33,ich) = 0.13d-22
            xct%sigma_i_EUV(34,ich) = 0.40d-22
            xct%sigma_i_EUV(35,ich) = 0.00d-22
            xct%sigma_i_EUV(36,ich) = 0.00d-22
            xct%sigma_i_EUV(37,ich) = 0.00d-22
          end if

          if ( reactants(1) == 'C2H4' &
          & .and. products(1) == 'C2H2+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H4 + hv -> C2H2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.00d-22
            xct%sigma_i_EUV(2,ich)  =  0.00d-22
            xct%sigma_i_EUV(3,ich)  =  0.21d-22
            xct%sigma_i_EUV(4,ich)  =  0.86d-22
            xct%sigma_i_EUV(5,ich)  =  1.66d-22
            xct%sigma_i_EUV(6,ich)  =  2.11d-22
            xct%sigma_i_EUV(7,ich)  =  1.55d-22
            xct%sigma_i_EUV(8,ich)  =  2.48d-22
            xct%sigma_i_EUV(9,ich)  =  2.49d-22
            xct%sigma_i_EUV(10,ich) =  2.44d-22
            xct%sigma_i_EUV(11,ich) =  3.82d-22
            xct%sigma_i_EUV(12,ich) =  3.42d-22
            xct%sigma_i_EUV(13,ich) =  4.74d-22
            xct%sigma_i_EUV(14,ich) =  7.42d-22
            xct%sigma_i_EUV(15,ich) =  6.66d-22
            xct%sigma_i_EUV(16,ich) =  8.43d-22
            xct%sigma_i_EUV(17,ich) =  9.52d-22
            xct%sigma_i_EUV(18,ich) =  9.85d-22
            xct%sigma_i_EUV(19,ich) =  9.48d-22
            xct%sigma_i_EUV(20,ich) = 10.18d-22
            xct%sigma_i_EUV(21,ich) = 10.18d-22
            xct%sigma_i_EUV(22,ich) = 10.11d-22
            xct%sigma_i_EUV(23,ich) = 10.10d-22
            xct%sigma_i_EUV(24,ich) = 10.06d-22
            xct%sigma_i_EUV(25,ich) = 10.03d-22
            xct%sigma_i_EUV(26,ich) =  9.70d-22
            xct%sigma_i_EUV(27,ich) =  9.55d-22
            xct%sigma_i_EUV(28,ich) =  9.05d-22
            xct%sigma_i_EUV(29,ich) = 10.10d-22
            xct%sigma_i_EUV(30,ich) =  8.79d-22
            xct%sigma_i_EUV(31,ich) =  7.30d-22
            xct%sigma_i_EUV(32,ich) =  4.99d-22
            xct%sigma_i_EUV(33,ich) =  0.86d-22
            xct%sigma_i_EUV(34,ich) =  2.10d-22
            xct%sigma_i_EUV(35,ich) =  0.04d-22
            xct%sigma_i_EUV(36,ich) =  0.01d-22
            xct%sigma_i_EUV(37,ich) =  0.17d-22
          end if

          if ( reactants(1) == 'C2H4' &
          & .and. products(1) == 'C2H+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H4 + hv -> C2H+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.00d-22
            xct%sigma_i_EUV(2,ich)  = 0.00d-22
            xct%sigma_i_EUV(3,ich)  = 0.04d-22
            xct%sigma_i_EUV(4,ich)  = 0.17d-22
            xct%sigma_i_EUV(5,ich)  = 0.28d-22
            xct%sigma_i_EUV(6,ich)  = 0.33d-22
            xct%sigma_i_EUV(7,ich)  = 0.27d-22
            xct%sigma_i_EUV(8,ich)  = 0.37d-22
            xct%sigma_i_EUV(9,ich)  = 0.37d-22
            xct%sigma_i_EUV(10,ich) = 0.37d-22
            xct%sigma_i_EUV(11,ich) = 0.55d-22
            xct%sigma_i_EUV(12,ich) = 0.50d-22
            xct%sigma_i_EUV(13,ich) = 0.71d-22
            xct%sigma_i_EUV(14,ich) = 0.39d-22
            xct%sigma_i_EUV(15,ich) = 0.50d-22
            xct%sigma_i_EUV(16,ich) = 0.20d-22
            xct%sigma_i_EUV(17,ich) = 0.05d-22
            xct%sigma_i_EUV(18,ich) = 0.00d-22
            xct%sigma_i_EUV(19,ich) = 0.06d-22
            xct%sigma_i_EUV(20,ich) = 0.00d-22
            xct%sigma_i_EUV(21,ich) = 0.00d-22
            xct%sigma_i_EUV(22,ich) = 0.00d-22
            xct%sigma_i_EUV(23,ich) = 0.00d-22
            xct%sigma_i_EUV(24,ich) = 0.00d-22
            xct%sigma_i_EUV(25,ich) = 0.00d-22
            xct%sigma_i_EUV(26,ich) = 0.00d-22
            xct%sigma_i_EUV(27,ich) = 0.00d-22
            xct%sigma_i_EUV(28,ich) = 0.00d-22
            xct%sigma_i_EUV(29,ich) = 0.00d-22
            xct%sigma_i_EUV(30,ich) = 0.00d-22
            xct%sigma_i_EUV(31,ich) = 0.00d-22
            xct%sigma_i_EUV(32,ich) = 0.00d-22
            xct%sigma_i_EUV(33,ich) = 0.00d-22
            xct%sigma_i_EUV(34,ich) = 0.00d-22
            xct%sigma_i_EUV(35,ich) = 0.00d-22
            xct%sigma_i_EUV(36,ich) = 0.00d-22
            xct%sigma_i_EUV(37,ich) = 0.00d-22
          end if

          if ( reactants(1) == 'C2H6' &
          & .and. products(1) == 'C2H6+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H6 + hv -> C2H6+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.00d-22
            xct%sigma_i_EUV(2,ich)  =  0.00d-22
            xct%sigma_i_EUV(3,ich)  =  0.00d-22
            xct%sigma_i_EUV(4,ich)  =  0.21d-22
            xct%sigma_i_EUV(5,ich)  =  0.50d-22
            xct%sigma_i_EUV(6,ich)  =  0.71d-22
            xct%sigma_i_EUV(7,ich)  =  0.47d-22
            xct%sigma_i_EUV(8,ich)  =  0.89d-22
            xct%sigma_i_EUV(9,ich)  =  0.90d-22
            xct%sigma_i_EUV(10,ich) =  0.87d-22
            xct%sigma_i_EUV(11,ich) =  1.87d-22
            xct%sigma_i_EUV(12,ich) =  1.56d-22
            xct%sigma_i_EUV(13,ich) =  2.45d-22
            xct%sigma_i_EUV(14,ich) =  4.28d-22
            xct%sigma_i_EUV(15,ich) =  3.93d-22
            xct%sigma_i_EUV(16,ich) =  5.35d-22
            xct%sigma_i_EUV(17,ich) =  6.91d-22
            xct%sigma_i_EUV(18,ich) =  7.65d-22
            xct%sigma_i_EUV(19,ich) =  6.76d-22
            xct%sigma_i_EUV(20,ich) =  8.28d-22
            xct%sigma_i_EUV(21,ich) =  8.97d-22
            xct%sigma_i_EUV(22,ich) =  8.05d-22
            xct%sigma_i_EUV(23,ich) =  9.49d-22
            xct%sigma_i_EUV(24,ich) = 10.24d-22
            xct%sigma_i_EUV(25,ich) = 10.25d-22
            xct%sigma_i_EUV(26,ich) = 11.29d-22
            xct%sigma_i_EUV(27,ich) = 11.35d-22
            xct%sigma_i_EUV(28,ich) = 11.84d-22
            xct%sigma_i_EUV(29,ich) = 11.09d-22
            xct%sigma_i_EUV(30,ich) = 12.16d-22
            xct%sigma_i_EUV(31,ich) = 12.30d-22
            xct%sigma_i_EUV(32,ich) = 12.18d-22
            xct%sigma_i_EUV(33,ich) = 12.23d-22
            xct%sigma_i_EUV(34,ich) = 12.00d-22
            xct%sigma_i_EUV(35,ich) = 10.08d-22
            xct%sigma_i_EUV(36,ich) =  9.66d-22
            xct%sigma_i_EUV(37,ich) = 11.83d-22
          end if

          if ( reactants(1) == 'C2H6' &
          & .and. products(1) == 'C2H5+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H6 + hv -> C2H5+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.00d-22
            xct%sigma_i_EUV(2,ich)  =  0.00d-22
            xct%sigma_i_EUV(3,ich)  =  0.00d-22
            xct%sigma_i_EUV(4,ich)  =  0.22d-22
            xct%sigma_i_EUV(5,ich)  =  0.54d-22
            xct%sigma_i_EUV(6,ich)  =  0.77d-22
            xct%sigma_i_EUV(7,ich)  =  0.49d-22
            xct%sigma_i_EUV(8,ich)  =  0.95d-22
            xct%sigma_i_EUV(9,ich)  =  0.95d-22
            xct%sigma_i_EUV(10,ich) =  0.92d-22
            xct%sigma_i_EUV(11,ich) =  1.81d-22
            xct%sigma_i_EUV(12,ich) =  1.48d-22
            xct%sigma_i_EUV(13,ich) =  2.29d-22
            xct%sigma_i_EUV(14,ich) =  3.66d-22
            xct%sigma_i_EUV(15,ich) =  3.37d-22
            xct%sigma_i_EUV(16,ich) =  4.46d-22
            xct%sigma_i_EUV(17,ich) =  5.90d-22
            xct%sigma_i_EUV(18,ich) =  6.61d-22
            xct%sigma_i_EUV(19,ich) =  5.75d-22
            xct%sigma_i_EUV(20,ich) =  7.41d-22
            xct%sigma_i_EUV(21,ich) =  8.18d-22
            xct%sigma_i_EUV(22,ich) =  7.13d-22
            xct%sigma_i_EUV(23,ich) =  8.67d-22
            xct%sigma_i_EUV(24,ich) =  9.10d-22
            xct%sigma_i_EUV(25,ich) =  9.15d-22
            xct%sigma_i_EUV(26,ich) = 10.26d-22
            xct%sigma_i_EUV(27,ich) = 10.28d-22
            xct%sigma_i_EUV(28,ich) = 10.59d-22
            xct%sigma_i_EUV(29,ich) = 10.16d-22
            xct%sigma_i_EUV(30,ich) = 10.81d-22
            xct%sigma_i_EUV(31,ich) =  9.68d-22
            xct%sigma_i_EUV(32,ich) =  7.61d-22
            xct%sigma_i_EUV(33,ich) =  2.63d-22
            xct%sigma_i_EUV(34,ich) =  4.72d-22
            xct%sigma_i_EUV(35,ich) =  0.66d-22
            xct%sigma_i_EUV(36,ich) =  0.51d-22
            xct%sigma_i_EUV(37,ich) =  1.27d-22
          end if

          if ( reactants(1) == 'C2H6' &
          & .and. products(1) == 'C2H4+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H6 + hv -> C2H4+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  =  0.00d-22
            xct%sigma_i_EUV(2,ich)  =  0.00d-22
            xct%sigma_i_EUV(3,ich)  =  0.00d-22
            xct%sigma_i_EUV(4,ich)  =  0.98d-22
            xct%sigma_i_EUV(5,ich)  =  2.34d-22
            xct%sigma_i_EUV(6,ich)  =  3.26d-22
            xct%sigma_i_EUV(7,ich)  =  2.14d-22
            xct%sigma_i_EUV(8,ich)  =  4.12d-22
            xct%sigma_i_EUV(9,ich)  =  4.15d-22
            xct%sigma_i_EUV(10,ich) =  3.91d-22
            xct%sigma_i_EUV(11,ich) =  7.69d-22
            xct%sigma_i_EUV(12,ich) =  6.59d-22
            xct%sigma_i_EUV(13,ich) =  9.85d-22
            xct%sigma_i_EUV(14,ich) = 16.18d-22
            xct%sigma_i_EUV(15,ich) = 14.85d-22
            xct%sigma_i_EUV(16,ich) = 19.68d-22
            xct%sigma_i_EUV(17,ich) = 25.37d-22
            xct%sigma_i_EUV(18,ich) = 28.52d-22
            xct%sigma_i_EUV(19,ich) = 24.88d-22
            xct%sigma_i_EUV(20,ich) = 31.53d-22
            xct%sigma_i_EUV(21,ich) = 34.32d-22
            xct%sigma_i_EUV(22,ich) = 30.48d-22
            xct%sigma_i_EUV(23,ich) = 36.34d-22
            xct%sigma_i_EUV(24,ich) = 39.16d-22
            xct%sigma_i_EUV(25,ich) = 39.19d-22
            xct%sigma_i_EUV(26,ich) = 44.13d-22
            xct%sigma_i_EUV(27,ich) = 44.36d-22
            xct%sigma_i_EUV(28,ich) = 44.94d-22
            xct%sigma_i_EUV(29,ich) = 43.28d-22
            xct%sigma_i_EUV(30,ich) = 45.21d-22
            xct%sigma_i_EUV(31,ich) = 42.12d-22
            xct%sigma_i_EUV(32,ich) = 40.06d-22
            xct%sigma_i_EUV(33,ich) = 24.03d-22
            xct%sigma_i_EUV(34,ich) = 33.90d-22
            xct%sigma_i_EUV(35,ich) =  9.22d-22
            xct%sigma_i_EUV(36,ich) =  7.59d-22
            xct%sigma_i_EUV(37,ich) = 15.97d-22
          end if

          if ( reactants(1) == 'C2H6' &
          & .and. products(1) == 'C2H3+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H6 + hv -> C2H3+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.00d-22
            xct%sigma_i_EUV(2,ich)  = 0.00d-22
            xct%sigma_i_EUV(3,ich)  = 0.00d-22
            xct%sigma_i_EUV(4,ich)  = 0.79d-22
            xct%sigma_i_EUV(5,ich)  = 1.59d-22
            xct%sigma_i_EUV(6,ich)  = 1.99d-22
            xct%sigma_i_EUV(7,ich)  = 1.50d-22
            xct%sigma_i_EUV(8,ich)  = 2.30d-22
            xct%sigma_i_EUV(9,ich)  = 2.31d-22
            xct%sigma_i_EUV(10,ich) = 2.25d-22
            xct%sigma_i_EUV(11,ich) = 3.46d-22
            xct%sigma_i_EUV(12,ich) = 3.11d-22
            xct%sigma_i_EUV(13,ich) = 4.00d-22
            xct%sigma_i_EUV(14,ich) = 5.57d-22
            xct%sigma_i_EUV(15,ich) = 5.33d-22
            xct%sigma_i_EUV(16,ich) = 6.21d-22
            xct%sigma_i_EUV(17,ich) = 6.85d-22
            xct%sigma_i_EUV(18,ich) = 7.28d-22
            xct%sigma_i_EUV(19,ich) = 6.75d-22
            xct%sigma_i_EUV(20,ich) = 7.25d-22
            xct%sigma_i_EUV(21,ich) = 7.47d-22
            xct%sigma_i_EUV(22,ich) = 7.27d-22
            xct%sigma_i_EUV(23,ich) = 7.59d-22
            xct%sigma_i_EUV(24,ich) = 7.65d-22
            xct%sigma_i_EUV(25,ich) = 7.67d-22
            xct%sigma_i_EUV(26,ich) = 7.15d-22
            xct%sigma_i_EUV(27,ich) = 7.04d-22
            xct%sigma_i_EUV(28,ich) = 5.79d-22
            xct%sigma_i_EUV(29,ich) = 7.45d-22
            xct%sigma_i_EUV(30,ich) = 4.93d-22
            xct%sigma_i_EUV(31,ich) = 1.08d-22
            xct%sigma_i_EUV(32,ich) = 0.19d-22
            xct%sigma_i_EUV(33,ich) = 0.00d-22
            xct%sigma_i_EUV(34,ich) = 0.02d-22
            xct%sigma_i_EUV(35,ich) = 0.00d-22
            xct%sigma_i_EUV(36,ich) = 0.00d-22
            xct%sigma_i_EUV(37,ich) = 0.00d-22
          end if

          if ( reactants(1) == 'C2H6' &
          & .and. products(1) == 'C2H2+' .and. products(2) == 'e-' ) then
            !------------------------------------  C2H6 + hv -> C2H2+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 0.00d-22
            xct%sigma_i_EUV(2,ich)  = 0.00d-22
            xct%sigma_i_EUV(3,ich)  = 0.00d-22
            xct%sigma_i_EUV(4,ich)  = 0.42d-22
            xct%sigma_i_EUV(5,ich)  = 0.81d-22
            xct%sigma_i_EUV(6,ich)  = 1.02d-22
            xct%sigma_i_EUV(7,ich)  = 0.76d-22
            xct%sigma_i_EUV(8,ich)  = 1.16d-22
            xct%sigma_i_EUV(9,ich)  = 1.16d-22
            xct%sigma_i_EUV(10,ich) = 1.13d-22
            xct%sigma_i_EUV(11,ich) = 1.89d-22
            xct%sigma_i_EUV(12,ich) = 1.70d-22
            xct%sigma_i_EUV(13,ich) = 2.28d-22
            xct%sigma_i_EUV(14,ich) = 2.49d-22
            xct%sigma_i_EUV(15,ich) = 2.61d-22
            xct%sigma_i_EUV(16,ich) = 2.56d-22
            xct%sigma_i_EUV(17,ich) = 2.75d-22
            xct%sigma_i_EUV(18,ich) = 2.97d-22
            xct%sigma_i_EUV(19,ich) = 2.74d-22
            xct%sigma_i_EUV(20,ich) = 2.96d-22
            xct%sigma_i_EUV(21,ich) = 3.10d-22
            xct%sigma_i_EUV(22,ich) = 2.98d-22
            xct%sigma_i_EUV(23,ich) = 3.33d-22
            xct%sigma_i_EUV(24,ich) = 3.58d-22
            xct%sigma_i_EUV(25,ich) = 3.56d-22
            xct%sigma_i_EUV(26,ich) = 3.39d-22
            xct%sigma_i_EUV(27,ich) = 3.36d-22
            xct%sigma_i_EUV(28,ich) = 2.72d-22
            xct%sigma_i_EUV(29,ich) = 3.50d-22
            xct%sigma_i_EUV(30,ich) = 2.28d-22
            xct%sigma_i_EUV(31,ich) = 0.49d-22
            xct%sigma_i_EUV(32,ich) = 0.09d-22
            xct%sigma_i_EUV(33,ich) = 0.00d-22
            xct%sigma_i_EUV(34,ich) = 0.00d-22
            xct%sigma_i_EUV(35,ich) = 0.00d-22
            xct%sigma_i_EUV(36,ich) = 0.00d-22
            xct%sigma_i_EUV(37,ich) = 0.00d-22
          end if

          if ( reactants(1) == 'Na' &
          & .and. products(1) == 'Na+' .and. products(2) == 'e-' ) then
            !------------------------------------  Na + hv -> Na+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 8.384d-23
            xct%sigma_i_EUV(2,ich)  = 3.717d-22
            xct%sigma_i_EUV(3,ich)  = 7.048d-22
            xct%sigma_i_EUV(4,ich)  = 8.230d-22
            xct%sigma_i_EUV(5,ich)  = 7.036d-22
            xct%sigma_i_EUV(6,ich)  = 5.689d-22
            xct%sigma_i_EUV(7,ich)  = 7.313d-22
            xct%sigma_i_EUV(8,ich)  = 4.960d-22
            xct%sigma_i_EUV(9,ich)  = 4.948d-22
            xct%sigma_i_EUV(10,ich) = 5.056d-22
            xct%sigma_i_EUV(11,ich) = 6.129d-24
            xct%sigma_i_EUV(12,ich) = 5.773d-24
            xct%sigma_i_EUV(13,ich) = 6.746d-24
            xct%sigma_i_EUV(14,ich) = 7.943d-24
            xct%sigma_i_EUV(15,ich) = 7.673d-24
            xct%sigma_i_EUV(16,ich) = 8.538d-24
            xct%sigma_i_EUV(17,ich) = 9.395d-24
            xct%sigma_i_EUV(18,ich) = 9.826d-24
            xct%sigma_i_EUV(19,ich) = 9.330d-24
            xct%sigma_i_EUV(20,ich) = 1.016d-23
            xct%sigma_i_EUV(21,ich) = 1.041d-23
            xct%sigma_i_EUV(22,ich) = 1.003d-23
            xct%sigma_i_EUV(23,ich) = 1.065d-23
            xct%sigma_i_EUV(24,ich) = 1.120d-23
            xct%sigma_i_EUV(25,ich) = 1.117d-23
            xct%sigma_i_EUV(26,ich) = 1.170d-23
            xct%sigma_i_EUV(27,ich) = 1.173d-23
            xct%sigma_i_EUV(28,ich) = 1.185d-23
            xct%sigma_i_EUV(29,ich) = 1.159d-23
            xct%sigma_i_EUV(30,ich) = 1.191d-23
            xct%sigma_i_EUV(31,ich) = 1.212d-23
            xct%sigma_i_EUV(32,ich) = 1.222d-23
            xct%sigma_i_EUV(33,ich) = 1.218d-23
            xct%sigma_i_EUV(34,ich) = 1.223d-23
            xct%sigma_i_EUV(35,ich) = 1.203d-23
            xct%sigma_i_EUV(36,ich) = 1.201d-23
            xct%sigma_i_EUV(37,ich) = 1.212d-23
          end if

          if ( reactants(1) == 'Mg' &
          & .and. products(1) == 'Mg+' .and. products(2) == 'e-' ) then
            !------------------------------------  Mg + hv -> Mg+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 1.209d-22
            xct%sigma_i_EUV(2,ich)  = 4.936d-22
            xct%sigma_i_EUV(3,ich)  = 6.757d-22
            xct%sigma_i_EUV(4,ich)  = 4.220d-22
            xct%sigma_i_EUV(5,ich)  = 1.579d-23
            xct%sigma_i_EUV(6,ich)  = 1.755d-23
            xct%sigma_i_EUV(7,ich)  = 1.538d-23
            xct%sigma_i_EUV(8,ich)  = 1.870d-23
            xct%sigma_i_EUV(9,ich)  = 1.872d-23
            xct%sigma_i_EUV(10,ich) = 1.850d-23
            xct%sigma_i_EUV(11,ich) = 2.208d-23
            xct%sigma_i_EUV(12,ich) = 2.122d-23
            xct%sigma_i_EUV(13,ich) = 2.344d-23
            xct%sigma_i_EUV(14,ich) = 2.553d-23
            xct%sigma_i_EUV(15,ich) = 2.513d-23
            xct%sigma_i_EUV(16,ich) = 2.626d-23
            xct%sigma_i_EUV(17,ich) = 2.683d-23
            xct%sigma_i_EUV(18,ich) = 2.686d-23
            xct%sigma_i_EUV(19,ich) = 2.681d-23
            xct%sigma_i_EUV(20,ich) = 2.673d-23
            xct%sigma_i_EUV(21,ich) = 2.654d-23
            xct%sigma_i_EUV(22,ich) = 2.680d-23
            xct%sigma_i_EUV(23,ich) = 2.625d-23
            xct%sigma_i_EUV(24,ich) = 2.511d-23
            xct%sigma_i_EUV(25,ich) = 2.520d-23
            xct%sigma_i_EUV(26,ich) = 2.314d-23
            xct%sigma_i_EUV(27,ich) = 2.295d-23
            xct%sigma_i_EUV(28,ich) = 2.221d-23
            xct%sigma_i_EUV(29,ich) = 2.368d-23
            xct%sigma_i_EUV(30,ich) = 2.177d-23
            xct%sigma_i_EUV(31,ich) = 1.952d-23
            xct%sigma_i_EUV(32,ich) = 1.702d-23
            xct%sigma_i_EUV(33,ich) = 1.288d-23
            xct%sigma_i_EUV(34,ich) = 1.436d-23
            xct%sigma_i_EUV(35,ich) = 1.022d-23
            xct%sigma_i_EUV(36,ich) = 9.891d-24
            xct%sigma_i_EUV(37,ich) = 1.162d-23
          end if

          if ( reactants(1) == 'Fe' &
          & .and. products(1) == 'Fe+' .and. products(2) == 'e-' ) then
            !------------------------------------  Fe + hv -> Fe+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 1.656d-22
            xct%sigma_i_EUV(2,ich)  = 4.475d-22
            xct%sigma_i_EUV(3,ich)  = 6.524d-22
            xct%sigma_i_EUV(4,ich)  = 7.745d-22
            xct%sigma_i_EUV(5,ich)  = 8.748d-22
            xct%sigma_i_EUV(6,ich)  = 8.891d-22
            xct%sigma_i_EUV(7,ich)  = 8.686d-22
            xct%sigma_i_EUV(8,ich)  = 8.881d-22
            xct%sigma_i_EUV(9,ich)  = 8.879d-22
            xct%sigma_i_EUV(10,ich) = 8.888d-22
            xct%sigma_i_EUV(11,ich) = 8.369d-22
            xct%sigma_i_EUV(12,ich) = 8.571d-22
            xct%sigma_i_EUV(13,ich) = 7.936d-22
            xct%sigma_i_EUV(14,ich) = 6.893d-22
            xct%sigma_i_EUV(15,ich) = 7.145d-22
            xct%sigma_i_EUV(16,ich) = 6.318d-22
            xct%sigma_i_EUV(17,ich) = 5.478d-22
            xct%sigma_i_EUV(18,ich) = 5.069d-22
            xct%sigma_i_EUV(19,ich) = 5.541d-22
            xct%sigma_i_EUV(20,ich) = 4.760d-22
            xct%sigma_i_EUV(21,ich) = 4.545d-22
            xct%sigma_i_EUV(22,ich) = 4.874d-22
            xct%sigma_i_EUV(23,ich) = 4.354d-22
            xct%sigma_i_EUV(24,ich) = 3.988d-22
            xct%sigma_i_EUV(25,ich) = 4.005d-22
            xct%sigma_i_EUV(26,ich) = 3.826d-22
            xct%sigma_i_EUV(27,ich) = 3.826d-22
            xct%sigma_i_EUV(28,ich) = 3.842d-22
            xct%sigma_i_EUV(29,ich) = 3.839d-22
            xct%sigma_i_EUV(30,ich) = 3.863d-22
            xct%sigma_i_EUV(31,ich) = 4.074d-22
            xct%sigma_i_EUV(32,ich) = 4.469d-22
            xct%sigma_i_EUV(33,ich) = 1.608d-25
            xct%sigma_i_EUV(34,ich) = 5.041d-22
            xct%sigma_i_EUV(35,ich) = 6.817d-25
            xct%sigma_i_EUV(36,ich) = 8.101d-25
            xct%sigma_i_EUV(37,ich) = 3.029d-25
          end if

          if ( reactants(1) == 'Si' &
          & .and. products(1) == 'Si+' .and. products(2) == 'e-' ) then
            !------------------------------------  Si + hv -> Si+ + e-  ---------
            xct%sigma_i_EUV(1,ich)  = 2.209d-22
            xct%sigma_i_EUV(2,ich)  = 5.039d-22
            xct%sigma_i_EUV(3,ich)  = 4.122d-23
            xct%sigma_i_EUV(4,ich)  = 5.743d-23
            xct%sigma_i_EUV(5,ich)  = 6.910d-23
            xct%sigma_i_EUV(6,ich)  = 7.209d-23
            xct%sigma_i_EUV(7,ich)  = 6.817d-23
            xct%sigma_i_EUV(8,ich)  = 7.314d-23
            xct%sigma_i_EUV(9,ich)  = 7.316d-23
            xct%sigma_i_EUV(10,ich) = 7.302d-23
            xct%sigma_i_EUV(11,ich) = 7.177d-23
            xct%sigma_i_EUV(12,ich) = 7.280d-23
            xct%sigma_i_EUV(13,ich) = 6.912d-23
            xct%sigma_i_EUV(14,ich) = 6.249d-23
            xct%sigma_i_EUV(15,ich) = 6.400d-23
            xct%sigma_i_EUV(16,ich) = 5.970d-23
            xct%sigma_i_EUV(17,ich) = 5.876d-23
            xct%sigma_i_EUV(18,ich) = 6.100d-23
            xct%sigma_i_EUV(19,ich) = 5.862d-23
            xct%sigma_i_EUV(20,ich) = 6.493d-23
            xct%sigma_i_EUV(21,ich) = 6.956d-23
            xct%sigma_i_EUV(22,ich) = 6.318d-23
            xct%sigma_i_EUV(23,ich) = 7.580d-23
            xct%sigma_i_EUV(24,ich) = 1.008d-22
            xct%sigma_i_EUV(25,ich) = 9.888d-23
            xct%sigma_i_EUV(26,ich) = 1.484d-22
            xct%sigma_i_EUV(27,ich) = 1.535d-22
            xct%sigma_i_EUV(28,ich) = 1.734d-22
            xct%sigma_i_EUV(29,ich) = 1.347d-22
            xct%sigma_i_EUV(30,ich) = 1.856d-22
            xct%sigma_i_EUV(31,ich) = 2.538d-22
            xct%sigma_i_EUV(32,ich) = 3.413d-22
            xct%sigma_i_EUV(33,ich) = 5.125d-22
            xct%sigma_i_EUV(34,ich) = 4.470d-22
            xct%sigma_i_EUV(35,ich) = 6.456d-22
            xct%sigma_i_EUV(36,ich) = 6.640d-22
            xct%sigma_i_EUV(37,ich) = 5.728d-22
          end if

        end if

      end do! ich

    end if

  end subroutine p__EUVAC_cross_section


end module p__EUVAC
