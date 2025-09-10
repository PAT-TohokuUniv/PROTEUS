module c__prm
  
  use v__tdec, only : cst_, spl_, set_, var_, grd_, flx_

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: c__prm__ini, c__prm__planet

contains


  subroutine c__prm__planet(spl, set,    & ! in
    &                       cst, grd, flx) ! inout
    type(spl_),   intent(in)     :: spl
    type(set_),   intent(in)     :: set
    type(cst_),   intent(inout)  :: cst
    type(grd_),   intent(inout)  :: grd
    type(flx_),   intent(inout)  :: flx
    integer i, j, ix, iy, iz
    real(dp) e, phi, theta, delta, h, omega, omega0, r0, alpha
    real(dp) tmp, tmp1, tmp2, tmpzarr(grd%nz)

    real(dp) Ls  ! solar longitude of Mars
    real(dp) DOY ! day of year for Earth

    allocate(grd%lt(grd%nx), grd%lat(grd%ny), grd%sza(grd%nx,grd%ny), grd%sza_xact(grd%nx,grd%ny))

    ! adapt for the selected planet

    !------------------------------------------------------------------------------
    !
    !                                    Mars
    !
    !------------------------------------------------------------------------------
    if (spl%planet == 'Mars') then

      Ls         = set%Ls

      r0         = 1.524_dp ! [AU]
      cst%R      = 3389.5e3_dp ! planetary raduis [m]
      grd%tilt   = 25.2_dp  ! tilt angle of rotational axis [deg] !! NOT radian!!
      omega0     = 100.0_dp   ! angle between equinox and periherion
      e          = 0.093_dp ! eccentricity

      cst%g      = 3.711_dp ! gravity acceleration @ Mars
      cst%Mplanet= 0.1075_dp*5.972e24_dp ! Mass of Mars
      cst%daysec = 24.0_dp * 3600.0_dp + 39.0_dp * 60.0_dp + 35.244_dp ! 1 day [sec] of Mars

      omega      = Ls + omega0 ! angle between equinox and chosen day

    !------------------------------------------------------------------------------
    !
    !                                   Jupiter
    !
    !------------------------------------------------------------------------------
    else if (spl%planet == 'Jupiter') then

      r0         = 5.2_dp  ! [AU]
      cst%R      = 71492.0e3_dp ! planetary raduis [m]
      grd%tilt   = 0.0_dp    ! tilt angle of rotational axis [deg] !! NOT radian!!
      omega0     = 0.0_dp    ! angle between equinox and periherion
      e          = 0.0_dp    ! eccentricity

      cst%g      = 24.8_dp ! gravity acceleration @ Jupiter
      cst%Mplanet= 1.898e27_dp ! Mass of Jupiter
      cst%daysec = 35729.685_dp ! rotational period of Jupiter [sec]

      omega      = 0.0_dp    ! angle from solstice

    !------------------------------------------------------------------------------
    !
    !                                    Earth
    !
    !------------------------------------------------------------------------------
    else if (spl%planet == 'Earth') then

      DOY        = 1.0_dp!set%DOY ! day of year ( Jan 1 = 1, Dec 31 = 365 )

      r0         = 1.0_dp  ! [AU]
      cst%R      = 6378.0e3_dp ! planetary raduis [m]
      grd%tilt   = 23.0_dp   ! tilt angle of rotational axis [deg] !! NOT radian!!
      omega0     = 0.0_dp    ! angle between equinox and periherion
      e          = 0.0_dp    ! eccentricity

      cst%g      = 9.80_dp ! gravity acceleration @ Earth
      cst%Mplanet=  5.9724e24_dp ! Mass of Earth
      cst%daysec = 24.0_dp * 3600.0_dp

      omega      = (DOY - 82.0_dp)*2.0_dp*360.0_dp/365.0_dp ! angle between equinox and chosen day

    !------------------------------------------------------------------------------
    !
    !                                    Titan
    !
    !------------------------------------------------------------------------------
    else if (spl%planet == 'Titan') then

      DOY        = set%DOY ! day of year ( Jan 1 = 1, Dec 31 = 365 )

      r0         = 10.0_dp  ! [AU]
      cst%R      = 2574.0e3_dp ! planetary raduis [m]
      grd%tilt   = 25.0_dp   ! tilt angle of rotational axis [deg] !! NOT radian!! ! shoule be modified !!!!!!!!!!!!!!
      omega0     = 0.0_dp    ! angle between equinox and periherion
      e          = 0.0_dp    ! eccentricity

      cst%g      = 1.35_dp ! gravity acceleration @ Earth
      cst%daysec = 24.0_dp * 3600.0_dp ! shoule be modified !!!!!!!!!!!!!!

      omega      = 0.0_dp ! angle between equinox and chosen day ! shoule be modified !!!!!!!!!!!!!!

    !------------------------------------------------------------------------------
    !
    !                                    Venus
    !
    !------------------------------------------------------------------------------
    else if (spl%planet == 'Venus') then

      Ls         = set%Ls

      r0         = 0.723_dp ! [AU]
      cst%R      = 6051.8e3_dp ! planetary raduis [m]
      grd%tilt   = 0.0_dp  ! tilt angle of rotational axis [deg] !! NOT radian!!
      omega0     = 100.0_dp   ! angle between equinox and periherion
      e          = 0.0_dp ! eccentricity

      cst%g      = 8.87_dp ! gravity acceleration @ Venus
      cst%Mplanet= 0.815_dp*5.972e24_dp ! Mass of Venus
      cst%daysec = 24.0_dp * 3600.0_dp * 117_dp ! 1 day [sec] of Venus

      omega      = 0.0 ! angle between equinox and chosen day

    end if
    !---------------------------------------------------------------------------------

    grd%tilt = grd%tilt * cst%pi / 180.0_dp
    omega = omega * cst%pi / 180.0_dp
    omega0 = omega0 * cst%pi / 180.0_dp

    ! seasonal dependence of solar flux
    flx%orbit = r0 * (1.0_dp - e*e) / (1.0_dp + e*dcos(omega))

    ! latitude, local time grid settings
    if (grd%nx == 1) then
      grd%lt = cst%pi ! LT 1200 noon
    else if (grd%nx >= 2) then
      do ix = 1, grd%nx
        ! ix = 1, grd%nx : midnight
        ! lt = 0 : midnight, lt = pi : noon
        grd%lt(ix) = ( dble(ix) - 1.0d0 ) * cst%pi * 2.0_dp / dble(grd%nx - 1)
      end do
    end if

    if (grd%ny == 1) then
      do iy = 1, grd%ny
        grd%lat(iy) = grd%latitude * cst%pi / 2.0_dp / 90.0_dp
      end do
    else if (grd%ny >= 2) then
      do iy = 1, grd%ny
        ! iy = 1      : south pole
        ! iy = grd%ny : north pole
        grd%lat(iy) = ( dble(iy) - dble(grd%ny+1)/2.0_dp ) * cst%pi / dble(grd%ny - 1)
      end do
    end if

    ! solar zenith angle
    alpha = grd%tilt
    delta = dasin(dsin(alpha)*dsin(omega-omega0))
    do ix = 1, grd%nx
      h = grd%lt(ix) - cst%pi
      do iy = 1, grd%ny
        phi = grd%lat(iy)
        tmp = dacos( dsin(phi) * dsin(delta) + dcos(phi) * dcos(delta) * dcos(h) )
        grd%sza_xact(ix,iy) = tmp
        !if( tmp >= 1.36_dp ) then
        !  tmp = 1.36_dp
        !end if
        grd%sza(ix,iy) = tmp
      end do
    end do

  end subroutine c__prm__planet

  subroutine c__prm__ini(cst) !inout
    type(cst_),   intent(inout)  :: cst

    ! Physical constant
    cst%pi  = dacos(-1.0_dp)
    cst%NA  = 6.022e23_dp
    cst%eV  = 1.6022e-19_dp ![eV]
    cst%k_B = 1.38064852e-23_dp !Boltzmann constant [m^2 kg s^-2 K^-1]
    cst%q_e = 1.6022e-19_dp ![C]
    cst%m_u = 1.660538921e-27_dp ![kg]
    cst%c   = 2.99792458e8_dp ![m/s]
    cst%h   = 6.626e-34_dp ![J s]
    cst%BigG = 6.67e-11_dp
    cst%R_gas = 8.2057e-2_dp ![L atm/K/mol] 

  end subroutine c__prm__ini

end module c__prm
