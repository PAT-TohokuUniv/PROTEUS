!==============================================================================================================
!
!                　Photochemical and RadiatiOn Transport model for Extensive USe (PROTEUS)
!
!==============================================================================================================
!
! 23 July 2025, Yuki Nakamura
!
!         This photochemical model is designed for application to many planets like Earth, Mars, Venus, Jupiter,
!       Titan, exoplanets, etc. Graphical User Interface developed with Github Electron framework is connected 
!       to Fortran modules for the flexibility to select a planet, to add and remove chemical reactions.
!

program e__main

  use v__tdec,                   only : set_, grd_, var_, cst_, xct_, spl_, flx_
  use c__prm,                    only : c__prm__ini, c__prm__planet
  use p__search,                 only : p__search_reactant, p__search_product, ch_identify, sp_index
  use v__Earth,                  only : v__Earth__ini
  use v__Venus,                  only : v__Venus__ini
  use v__Mars,                   only : v__Mars__ini
  use v__Jupiter,                only : v__Jupiter__ini
  use v__in,                     only : v__in__ini, v__in__exe
  use p__EUVAC,                  only : p__EUVAC_flux, p__EUVAC_cross_section
  use p__EUV,                    only : p__EUV_flux, p__EUV_cross_section
  use p__UV,                     only : p__UV_flux, p__UV_cross_section
  use p__PROTEUS,                only : p__PROTEUS_source, p__PROTEUS_Jacobian, &
    &                                   get_reactant_product_list, get_reaction_type_list
  use p__vertical_diffusion,     only : p__vertical_diffusion__exe
  use p__photolysis_rate,        only : p__photolysis_rate__ini, p__photolysis_rate__fin, &
    &                                   load_cross_section_dat, get_cross_section, photolysis_rate
  use p__photochem_rate,         only : p__photochem_rate__ini, p__photochem_rate__exe
  use p__photochem_transport,    only : p__photochem_transport__ini, p__photochem_transport__exe
  use p__chemical_scheme,        only : p__chemical_scheme__exe, p__LU_decomposition, p__LU_solver
  use p__photochem_scheme,       only : p__photochem_scheme__ini, p__photochem_scheme__exe
  use p__io,                     only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                                   p__io_timeseries__ini, p__io_timeseries__exe, p__io_timeseries__fin, &
    &                                   p__io_progress

  implicit none
  integer(4), parameter  :: sp = 4, dp = 8
  type(set_) :: set ! type of calculation settings
  type(grd_) :: grd ! type of grid
  type(var_) :: var ! type of variables
  type(cst_) :: cst ! type of physical constant
  type(xct_) :: xct ! type of cross section
  type(spl_) :: spl ! type of species list and planet info
  type(flx_) :: flx ! type of solar flux
  integer i, j, k, ix, xs, iy, iz, ich, isp, jsp, is, iday, iout
  real(dp) sumt
  real(dp) tmp, tmp1, tmp2, tmp3
  character(len=256) fname, command, ci

  real(dp) t1, t2, t3, t4

  ! Call cpu time
  call cpu_time(t1)

  !----------------------------------------------------------------------------------------------------------
  !
  !                                     Initialization : automatically adapted
  !
  !----------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------
  ! read physical constant
  !-----------------------------------------------------
  call c__prm__ini(cst) ! out

  !-----------------------------------------------------
  ! Planet selection: automatically
  !-----------------------------------------------------
  call v__in__ini(spl, set) ! out
  write(*,*) '----------------------------------------------------------------------------------------------------------'
  write(*,*) '|                                                                                                        |'
  write(*,*) '|        * Welcome to PROTEUS (Photochemical and RadiatiOn Transport model for Extensive USe) *          |'
  write(*,*) '|                                                                                                        |'
  write(*,*) '----------------------------------------------------------------------------------------------------------'
  write(*,'("  Selected planet       :   ",a80)') spl%planet
  write(*,'("  Selected project dir  :   ",a80)') set%dir_name

  !-----------------------------------------------------
  ! initialize variables
  !-----------------------------------------------------
  call v__in__exe(cst,      & ! in
    &             set, spl, & ! inout
    &             var, grd  ) ! out

  set%start_rot = 0
  grd%latitude = set%latitude

  !-----------------------------------------------------
  ! apply physical constant and geometry to the chosen planet
  !-----------------------------------------------------
  call c__prm__planet(spl, set,    & ! in
    &                 cst, grd, flx) ! inout

  write(*,'("  Distance from the sun : ",f7.2," [AU]")') flx%orbit
  write(*,*) '----------------------------------------------------------------------------------------------------------'

  !-----------------------------------------------------
  ! Special treatment for selected planet
  !-----------------------------------------------------
  if (      spl%planet == 'Venus' ) then
    call v__Venus__ini(spl, cst, grd, flx, & ! in
      &                var                 ) ! inout
  else if ( spl%planet == 'Earth' ) then
    call v__Earth__ini(spl, cst, grd, flx, & ! in
      &                var                 ) ! inout
  else if ( spl%planet == 'Mars' ) then
    call v__Mars__ini(spl, cst, grd, flx, & ! in
      &               var                 ) ! inout
  else if ( spl%planet == 'Jupiter' ) then
    call v__Jupiter__ini(spl, cst, grd, flx, set, & ! in
      &                  var                      ) ! inout
  end if

  !-----------------------------------------------------
  ! Solar irradiance data
  !-----------------------------------------------------
  call p__UV_flux(spl, grd, set, & ! in
    &             var, xct, flx  ) ! inout

  !-----------------------------------------------------
  ! UV cross section data 
  !-----------------------------------------------------
  call p__UV_cross_section(spl, grd, flx, var)

  !-----------------------------------------------------
  ! Solar EUV irradiance data
  !-----------------------------------------------------
  if (set%euv_input /= 'EUVAC') then 
    call p__EUV_flux(spl, grd, set, & ! in
      &              var, xct, flx  ) ! inout
  else if (set%euv_input == 'EUVAC') then
    call p__EUVAC_flux(spl, grd, set, & ! in
      &                var, xct, flx  ) ! inout
  end if

  !-----------------------------------------------------
  ! EUV cross section data 
  !-----------------------------------------------------
  if (set%euv_input /= 'EUVAC') then 
    call p__EUV_cross_section(spl, grd, flx, var)
  else if (set%euv_input == 'EUVAC') then
    xct%type = 'absorption'
    call p__EUVAC_cross_section(spl, & ! in
      &                         xct  ) ! inout
    xct%type = 'ionization'
    call p__EUVAC_cross_section(spl, & ! in
      &                         xct  ) ! inout
  end if

  !-----------------------------------------------------
  ! Initialization of module-level variables
  !-----------------------------------------------------
  call p__photochem_rate__ini(spl, grd) ! in
  call p__photochem_transport__ini(spl, grd) ! in
  call p__photochem_scheme__ini(spl, grd) ! in

  !-----------------------------------------------------
  ! Initialization of timeseries output
  !-----------------------------------------------------
  call p__io_timeseries__ini(set, grd, spl, var)


  ! Call cpu time
  call cpu_time(t2)

  write(*,*) '----------------------------------------------------------------------------------------------------------'
  write(*,*) 'Start calculation!'


  !----------------------------------------------------------------------------------------------------------
  !
  !                          Start photochemical calculation: 1D and 2D stable solution
  !
  !----------------------------------------------------------------------------------------------------------
  if (set%mode == '1D' .or. set%mode == '2DLAT') set%calc_stable = 1
  if (set%calc_stable == 1) then

    call cpu_time(var%t1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( set%fix_sza == 1 ) then
      grd%sza = set%sza * cst%pi / 180.0_dp
      grd%sza_xact = set%sza * cst%pi / 180.0_dp
    end if
    flx%irradiance_factor = 1.0_dp 
    if (set%diurnal_ave == 1) then
      flx%irradiance_factor = 0.5_dp 
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    grd%ix = (grd%nx-1)/2+1  ! 12:00 LT

    do isp = 1, spl%nsp
      do iy = 1, grd%ny
        var%ni_stable(iy,1:grd%nz,isp) = var%ni(1:grd%nz,isp)
      end do
    end do

    do iy = 1, grd%ny
      grd%iy = iy

      do isp = 1, spl%nsp
        var%ni(1:grd%nz,isp) = var%ni_stable(iy,1:grd%nz,isp)
      end do

      var%dtime = 1.0e-8_dp
      var%sum_time = 0.0_dp
      iout = 0
      call cpu_time(var%t2)
      call p__io_progress(spl, var, grd, set) ! in

      loopt: do is = 1, set%nstep
        var%istep = is

        ! dt upper limit
        if ( var%dtime  >= set%dtime_limit ) then
          var%dtime = set%dtime_limit ! dt DO NOT excess set%dtime_limit 
        end if

        !-----------------------------------------------------
        !            Photochemical calculation
        !-----------------------------------------------------
        call p__photochem_rate__exe(spl, cst, grd, flx, set, & ! in
          &                         xct, var                 ) ! inout
        call p__photochem_transport__exe(spl, cst, grd, set, & ! in
          &                              var                 ) ! inout
        call p__photochem_scheme__exe(spl, grd, set, & ! in
          &                           var            ) ! inout
          
        if ( var%sum_time > set%fin_sec ) then
          exit loopt
        end if

        if (set%mode == '1D') then 
          write(*,'("  time step = ",i6,"  dt = ",e9.3," s  t = ",e9.3," s  : max(dN/N) = ",e9.3,"  @ zgrid = ",i3,"  by ",a10)') &
            &             var%istep, &
            &             var%dtime, var%sum_time, &
            &             var%max_dn_n(3), nint(var%max_dn_n(2)), spl%species(nint(var%max_dn_n(1)))
        end if

        call p__io_timeseries__exe(set, grd, spl, var, iout)

      end do loopt ! end of time step

      if (set%mode == '2DLAT') then
        write(*,*) grd%iy, '/', grd%ny
      end if

      do isp = 1, spl%nsp
        var%ni_stable(iy,1:grd%nz,isp) = var%ni(1:grd%nz,isp)
      end do

      call cpu_time(var%t2)

      call p__io_progress(spl, var, grd, set) ! in

    end do ! end y

    !----------------------------------------------------------------------------------------------------------
    !
    !                                           Output stable results
    !
    !----------------------------------------------------------------------------------------------------------

    call cpu_time(t3)

    call p__io_stable__fin(set, spl, var, grd) ! in
    call p__io_timeseries__fin(set, grd, spl, iout) ! in

  end if



  !----------------------------------------------------------------------------------------------------------
  !
  !            Start photochemical calculation: 2D and 3D rotation (not exactly the global model)
  !
  !                  !!! YOU MUST CALCULATE OR INPUT STABLE DENSITY BEFORE THIS MODE !!!
  !
  !----------------------------------------------------------------------------------------------------------
  set%start_rot   = 1

  if (set%mode == '2DROT' .or. set%mode == '3DROT') then

    if (set%mode == '2DROT') then
      if (set%use_1d == 1) then
        xs = (grd%nx-1)/2+1
        do isp = 1, spl%nsp
          fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/num/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'unknown' )
            do iz = 1, grd%nz
              read(11, *) tmp, var%ni_stable(1,iz,isp) 
            end do
          close(11)
        end do
      end if
    else if (set%mode == '3DROT') then
      if (set%use_2d == 1) then
        xs = (grd%nx-1)/2+1
        do isp = 1, spl%nsp
          fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/2DLAT/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'unknown' )
            do iz = 1, grd%nz
              read(11, *) tmp, (var%ni_stable(iy,iz,isp), iy = 1, grd%ny) 
            end do
          close(11)
        end do
      end if
    end if

    var%ni_3d = 0.0_dp
    ix = (grd%nx-1)/2+1
    do isp = 1, spl%nsp
      do iz = 1, grd%nz
        var%ni_3d(ix,1:grd%ny,iz,isp) = var%ni_stable(1:grd%ny,iz,isp)
      end do
    end do

    call cpu_time(var%t3)

    flx%irradiance_factor = 1.0_dp

    var%dtime = cst%daysec / grd%nx
    set%dtime_limit = cst%daysec / grd%nx

    write(*,*) 'latitude index, rotation index, local time index'

    do iy = 1, grd%ny
      grd%iy = iy

      call cpu_time(var%t4)
      call p__io_progress(spl, var, grd, set) ! in

      do iday = 1, set%nday
        grd%iday = iday

        if ( iday == 1 ) then
          xs = (grd%nx-1)/2+1
        else if (iday >= 2) then
          xs = 1
        end if

        do ix = xs, grd%nx
          grd%ix = ix
          
          if ( ix == xs ) then
            do isp = 1, spl%nsp
              var%ni(1:grd%nz,isp) = var%ni_3d(ix,iy,1:grd%nz,isp)
            end do
          else if ( ix >= xs+1 ) then
            do isp = 1, spl%nsp
              var%ni(1:grd%nz,isp) = var%ni_3d(ix-1,iy,1:grd%nz,isp)
            end do
          end if

          ! dt upper limit
          if ( var%dtime  >= set%dtime_limit ) then
            var%dtime = set%dtime_limit ! dt DO NOT excess set%dtime_limit 
          end if

          !-----------------------------------------------------
          !            Photochemical calculation
          !-----------------------------------------------------
          call p__photochem_rate__exe(spl, cst, grd, flx, set, & ! in
            &                         xct, var                 ) ! inout
          call p__photochem_transport__exe(spl, cst, grd, set, & ! in
            &                              var                 ) ! inout
          call p__photochem_scheme__exe(spl, grd, set, & ! in
            &                           var            ) ! inout

          do isp = 1, spl%nsp
            var%ni_3d(ix,iy,1:grd%nz,isp) = var%ni(1:grd%nz,isp)
          end do

          if ( ix == grd%nx ) then
            do isp = 1, spl%nsp
              var%ni_3d(1,iy,1:grd%nz,isp) = var%ni(1:grd%nz,isp)
            end do
          end if

          write(*,*) grd%iy, '/', grd%ny, ',', iday, '/', set%nday, ',', grd%ix, '/', grd%nx

        end do ! end of x : local time
      end do ! end of rotation

      call cpu_time(var%t4)

      call p__io_progress(spl, var, grd, set) ! in

    end do ! end of y : latitude

    !----------------------------------------------------------------------------------------------------------
    !
    !                                         Output 2D or 3D density
    !
    !----------------------------------------------------------------------------------------------------------
    
    call cpu_time(t3)

    call p__io_rotation__fin(set, spl, var, grd) ! in

  end if


  !----------------------------------------------------------------------------------------------------------
  !
  !                                            Finalize: deallocate
  !
  !----------------------------------------------------------------------------------------------------------

  call p__io_deallocate(var, spl, xct, flx) ! inout

  write(*,*) '----------------------------------------------------------------------------------------------------------'

  if (t2-t1 < 60.0_dp) then 
    write(*,'("  cpu time  = ",f6.2," [s] to read data.")') t2-t1
  else if (t2-t1 >= 60.0_dp) then 
    write(*,'("  cpu time  = ",f6.2," [min] to read data.")') (t2-t1)/60.0_dp
  end if

  if (t3-t2 < 60.0_dp) then 
    write(*,'("  cpu time  = ",f6.2," [s] to calculate.")') t3-t2
  else if (t3-t2 >= 60.0_dp .and. t3-t2 < 3600.0_dp) then 
    write(*,'("  cpu time  = ",f6.2," [min] to calculate.")') (t3-t2)/60.0_dp
  else if (t3-t2 >= 3660.0_dp) then 
    write(*,'("  cpu time  = ",f6.2," [h] to calculate.")') (t3-t2)/3600.0_dp
  end if

  write(*,*) 'Finish.'


end program e__main
