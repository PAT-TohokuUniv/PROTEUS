!==============================================================================================================
!
!                　Photochemical and RadiatiOn Transport model for Extensive USe (PROTEUS)
!
!==============================================================================================================
!
! 2 July 2020, Yuki Nakamura
!
!         This photochemical model is designed for application to many planets like Earth, Mars, Venus, Jupiter,
!       Titan, exoplanets, etc. For the flexibility to select a planet, to add and to remove chemical reactions,
!       Graphical User Interface run on Python tkinter is connected to this model. The usage of this model for all
!       users is as follows.
!
!        1 Run PROTEUS.py on Python 3.
!        2 On the Planet Selection window, please select a planet intended to run.
!        3 On the Chemical reaction lists window, please select chemical reactions intended to include the photo-
!            chemical model. You can search certain reactions by species, references and labels on the search window.
!        4 Please set input density paths, vertical grids, boundary conditions, and calculation settings.
!        5 Press 'Output f90 module & Run model' to run the model.
!


!
! NOT RECOMMENDED TO USE "INCLUDE", BUT JUST IN CASE CMAKE DOES NOT WORK
!

!include "v__tdec.f90"
!include "c__prm.f90"
!include "p__io.f90"
!include "p__search.f90"
!include "v__Venus.f90"
!include "v__Earth.f90"
!include "v__Mars.f90"
!include "v__Jupiter.f90"
!include "p__EUVAC.f90"
!include "p__UV.f90"
!include "p__eddy_diffusion.f90"
!include "p__molecular_diffusion.f90"
!include "p__photochem_opticaldepth.f90"
!include "p__photochem_rate.f90"
!include "p__photochem_transport.f90"
!include "p__photochem_scheme.f90"
!
!! Project directory and v__in.f90 therein
!include "./Mars/Neutral_atmosphere/v__in.f90"
!include "./Jupiter/Nakamura_22_Case1/v__in.f90"

program e__main

  use v__tdec,                   only : set_, grd_, var_, cst_, xct_, spl_, flx_
  use c__prm,                    only : c__prm__ini, c__prm__planet
  use p__search,                 only : p__search_reactant, p__search_product, ch_identify, sp_index
  use v__Earth,                  only : v__Earth__ini, v__Earth__exe
  use v__Venus,                  only : v__Venus__ini, v__Venus__exe
  use v__Mars,                   only : v__Mars__ini, v__Mars__exe
  use v__Jupiter,                only : v__Jupiter__ini, v__Jupiter__exe
  use v__in,                     only : v__in__ini, v__in__exe
  use p__EUVAC,                  only : p__EUVAC_flux, p__EUVAC_cross_section
  use p__UV,                     only : p__UV_flux, p__UV_cross_section_dat, p__UV_cross_section_exe, &
    &                                   p__UV_EUV_xct_convergence_exe
  use p__photochem_opticaldepth, only : p__photochem_opticaldepth__exe
  use p__photochem_rate,         only : p__photochem_rate__exe
  use p__photochem_transport,    only : p__photochem_transport__exe
  use p__photochem_scheme,       only : p__photochem_scheme__exe
  use p__io,                     only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
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
  integer i, j, k, ix, xs, iy, iz, ich, isp, jsp, is, iday
  integer iout, iout_tmp, nsp_out
  real(dp) dt_out, t_out(0:10000)
  real(dp) cln_tmp
  real(dp) tmp
  character(len=256) fname, command
  character(len=256), allocatable :: species_out(:)
  real(dp), allocatable :: cln_out(:,:), n_out(:,:,:)

  real(dp) t1, t2, t3, t4

  ! Call cpu time
  call cpu_time(t1)

  ! for output of timevariations -----------------------------
  dt_out  = 1.0e2_dp ! sec
  nsp_out = 1
  allocate(species_out(nsp_out))
  species_out(1) = 'O3'


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
  print *, '----------------------------------------------------------------------------------------------------------'
  print *, '|                                                                                                        |'
  print *, '|        * Welcome to PROTEUS (Photochemical and RadiatiOn Transport model for Extensive USe) *          |'
  print *, '|                                                                                                        |'
  print *, '----------------------------------------------------------------------------------------------------------'
  write(*,'("  Selected planet       :   ",a80)') spl%planet
  write(*,'("  Selected project dir  :   ",a80)') set%dir_name

  !-----------------------------------------------------
  ! Calculation settings for the selected planet: auto
  !-----------------------------------------------------
  if (      spl%planet == 'Venus' ) then
    call v__Venus__ini(spl, set)
  else if ( spl%planet == 'Earth' ) then
    call v__Earth__ini(spl, set)
  else if ( spl%planet == 'Mars' ) then
    call v__Mars__ini(spl, set)
  else if ( spl%planet == 'Jupiter' ) then
    call v__Jupiter__ini(spl, set)
  end if

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
  print *, '----------------------------------------------------------------------------------------------------------'

  !-----------------------------------------------------
  ! Special treatment for selected planet
  !-----------------------------------------------------
  if (      spl%planet == 'Venus' ) then
    call v__Venus__exe(spl, cst, grd, flx, & ! in
      &                var                 ) ! inout
  else if ( spl%planet == 'Earth' ) then
    call v__Earth__exe(spl, cst, grd, flx, & ! in
      &                var                 ) ! inout
  else if ( spl%planet == 'Mars' ) then
    call v__Mars__exe(spl, cst, grd, flx, & ! in
      &               var                 ) ! inout
  else if ( spl%planet == 'Jupiter' ) then
    call v__Jupiter__exe(spl, cst, grd, flx, set, & ! in
      &                  var                      ) ! inout
  end if

  !-----------------------------------------------------
  ! EUVAC solar EUV flux Model
  !-----------------------------------------------------
  call p__EUVAC_flux(spl, grd, set, & ! in
    &                var, xct, flx  ) ! inout

  !-----------------------------------------------------
  ! Woods et al., solar UV flux reference
  !-----------------------------------------------------
  call p__UV_flux(spl, grd, cst,& ! in
    &             var, xct, flx ) ! inout

  !-----------------------------------------------------
  ! UV cross section data 
  !-----------------------------------------------------
  call p__UV_cross_section_dat(spl,    & ! in
    &                          xct, flx) ! inout

  var%istep = 0
  ! absorption cross sections
  xct%type = 'absorption'
  call p__EUVAC_cross_section(spl, & ! in
    &                         xct  ) ! inout
  call p__UV_cross_section_exe(var%Tn, grd%nz, spl, flx, var, & ! in
    &                          xct                            ) ! inout

  ! photoionization cross sections
  xct%type = 'ionization'
  call p__EUVAC_cross_section(spl, & ! in
    &                         xct  ) ! inout

  ! photodissociation cross sections
  xct%type = 'dissociation'
  call p__UV_cross_section_exe(var%Tn, grd%nz, spl, flx, var, & ! in
    &                          xct                            ) ! inout

  ! Insert EUV absorption cross section into UV cross section when the data is absent below 105 nm.
  call p__UV_EUV_xct_convergence_exe(grd%nz, spl, var, flx, & ! in
    &                                xct                    ) ! inout


  allocate(cln_out(nsp_out,0:set%nstep),n_out(nsp_out,0:set%nstep,grd%nz))

  ! output ==================================
  t_out(0) = 0.0_dp
  do i = 1, nsp_out
    isp = sp_index(spl,trim(ADJUSTL(species_out(i))))
    if (isp >= 1 .and. isp <= spl%nsp) then 
      cln_tmp = 0.0_dp
      do iz = 1, grd%nz
        cln_tmp = cln_tmp + var%ni(isp,iz)*grd%dalt(iz)
        n_out(i,0,iz) = var%ni(isp,iz)
      end do 
      cln_out(i,0) = cln_tmp
    end if
  end do

  ! Call cpu time
  call cpu_time(t2)

  print *, '----------------------------------------------------------------------------------------------------------'
  print *, 'Start calculation!'


  !----------------------------------------------------------------------------------------------------------
  !
  !                          Start photochemical calculation: 1D and 2D stable solution
  !
  !----------------------------------------------------------------------------------------------------------
  if (set%mode == '1D' .or. set%mode == '2D Lat') set%calc_stable = 1
  if (set%calc_stable == 1) then

    call cpu_time(var%t1)

    !!! For calculation at a certain sza
    !set%test_loc = 1 !!!!!!!!!!!!
    !set%sza = 75.0_dp *cst%pi/180.0_dp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    flx%mode_factor = 0.5_dp ! average over LT
    if ( set%test_loc == 1 ) then
      flx%mode_factor = 1.0_dp
      grd%sza = set%sza
      grd%sza_xact = set%sza
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    grd%ix = (grd%nx-1)/2+1  ! 12:00 LT

    do iz = 1, grd%nz
      do iy = 1, grd%ny
        do isp = 1, spl%nsp
          var%ni_stable(isp,iy,iz) = var%ni(isp,iz)
        end do
      end do
    end do

    do iy = 1, grd%ny
      grd%iy = iy

      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          var%ni(isp,iz) = var%ni_stable(isp,iy,iz)
        end do
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
        call p__photochem_opticaldepth__exe(spl, cst, flx, grd, set, & ! in
          &                                 xct, var                 ) ! inout
        call p__photochem_rate__exe(spl, cst, flx, grd, set, & ! in
          &                         xct, var                 ) ! inout
        call p__photochem_transport__exe(spl, cst, grd, set, & ! in
          &                              var                 ) ! inout
        call p__photochem_scheme__exe(spl, grd, set, & ! in
          &                           var            ) ! inout

        if ( var%sum_time > set%fin_sec ) then
          exit loopt
        end if

        var%sum_time = var%sum_time + var%dtime

        write(*,'("  time step = ",i6,"  dt = ",e9.3," s  t = ",e9.3," s  : max(dN/N) = ",e9.3,"  @ zgrid = ",i3,"  by ",a10)') &
          &             var%istep, &
          &             var%dtime, var%sum_time, &
          &             var%max_dn_n(3), nint(var%max_dn_n(2)), spl%species(nint(var%max_dn_n(1)))


        iout_tmp = int(var%sum_time / dt_out)
        if (iout_tmp == iout + 1 .or. iout_tmp == iout + 2 .or. var%dtime > dt_out) then 
          iout = iout + 1
          ! output ==================================
          t_out(iout) = var%sum_time
          do i = 1, nsp_out
            isp = sp_index(spl,trim(ADJUSTL(species_out(i))))
            if (isp >= 1 .and. isp <= spl%nsp) then 
              cln_tmp = 0.0_dp
              do iz = 1, grd%nz
                cln_tmp = cln_tmp + var%ni(isp,iz)*grd%dalt(iz)
                n_out(i,iout,iz) = var%ni(isp,iz)
              end do 
              cln_out(i,iout) = cln_tmp
            end if
          end do
          
        end if


      end do loopt ! end of time step

      do isp = 1, spl%nsp
        do iz = 1, grd%nz
          var%ni_stable(isp,iy,iz) = var%ni(isp,iz)
        end do
      end do

      call cpu_time(var%t2)

      call p__io_progress(spl, var, grd, set) ! in
    end do ! end y

    !do ich = 1, spl%nch
    !  if (var%Pij(i,25,ich)/=0.0_dp) then 
    !    write(*,'(i3, e9.2)'), ich, var%Pij(1,25,ich)
    !  end if
    !end do 


    !----------------------------------------------------------------------------------------------------------
    !
    !                                           Output stable solution
    !
    !----------------------------------------------------------------------------------------------------------

    call p__io_stable__fin(set, spl, var, grd) ! in

    fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries'
    write(command,*) 'if [ ! -d ', trim(fname), ' ]; then mkdir -p ', trim(fname), '; fi'
    call system(command)

    fname = './'//trim(ADJUSTL(set%dir_name))//'/output/timeseries/alt_km.dat'
    open(11, file = fname, status = 'unknown' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz) / 1.0e3_dp
      end do
    close(11)

    fname = './'//trim(ADJUSTL(set%dir_name))//'/output/timeseries/timestamp.dat'
    open(11, file = fname, status = 'unknown' )
      do i = 0, iout
        write(11, *) t_out(i)
      end do
    close(11)

    do isp = 1, nsp_out
      jsp = sp_index(spl,trim(ADJUSTL(species_out(isp))))
      if (jsp >= 1 .and. jsp <= spl%nsp) then 
        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/timeseries/n_'//trim(ADJUSTL(species_out(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do i = 0, iout
            write(11, *) (n_out(isp,i,iz), iz = 1, grd%nz)
          end do
        close(11)
        fname = './'//trim(ADJUSTL(set%dir_name))//'/output/timeseries/cln_'//trim(ADJUSTL(species_out(isp)))//'.dat'
        open(11, file = fname, status = 'unknown' )
          do i = 0, iout
            write(11, *) t_out(i), cln_out(isp,i)
          end do
        close(11)
      end if
    end do

  end if



  !----------------------------------------------------------------------------------------------------------
  !
  !            Start photochemical calculation: 2D and 3D rotation (not exactly the global model)
  !
  !                  !!! YOU MUST CALCULATE OR INPUT STABLE DENSITY BEFORE THIS MODE !!!
  !
  !----------------------------------------------------------------------------------------------------------
  set%start_rot   = 1

  if (set%mode == '2D Rot' .or. set%mode == '3D Rot') then

    if (set%mode == '2D Rot') then
      if (set%calc_stable == 0 .and. set%read_stable == 1) then
        xs = (grd%nx-1)/2+1
        do isp = 1, spl%nsp
          fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/num/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'unknown' )
            do iz = 1, grd%nz
              read(11, *) tmp, var%ni_stable(isp,1,iz) ! tmp: altitude in km
            end do
          close(11)
        end do
      end if
    else if (set%mode == '3D Rot') then
      if (set%calc_stable == 0 .and. set%read_stable == 1) then
        xs = (grd%nx-1)/2+1
        do isp = 1, spl%nsp
          fname = './'//trim(ADJUSTL(set%dir_name))//'/output/density/2Dstable/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'unknown' )
            do iz = 1, grd%nz
              read(11, *) tmp, (var%ni_stable(isp,iy,iz), iy = 1, grd%ny) ! tmp: altitude in km
            end do
          close(11)
        end do
      end if
    end if

    var%ni_3d = 0.0_dp
    do iz = 1, grd%nz
      do iy = 1, grd%ny
        ix = (grd%nx-1)/2+1
        do isp = 1, spl%nsp
          var%ni_3d(isp,ix,iy,iz) = var%ni_stable(isp,iy,iz)
        end do
      end do
    end do

    call cpu_time(var%t3)

    flx%mode_factor = 1.0_dp

    var%dtime = cst%daysec / grd%nx
    set%dtime_limit = cst%daysec / grd%nx

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
            do iz = 1, grd%nz
              do isp = 1, spl%nsp
                var%ni(isp,iz) = var%ni_3d(isp,ix,iy,iz)
              end do
            end do
          else if ( ix >= xs+1 ) then
            do iz = 1, grd%nz
              do isp = 1, spl%nsp
                var%ni(isp,iz) = var%ni_3d(isp,ix-1,iy,iz)
              end do
            end do
          end if

          ! dt upper limit
          if ( var%dtime  >= set%dtime_limit ) then
            var%dtime = set%dtime_limit ! dt DO NOT excess set%dtime_limit 
          end if

          !-----------------------------------------------------
          !            Photochemical calculation
          !-----------------------------------------------------
          call p__photochem_opticaldepth__exe(spl, cst, flx, grd, set, & ! in
            &                                 xct, var                 ) ! inout
          call p__photochem_rate__exe(spl, cst, flx, grd, set, & ! in
            &                         xct, var                 ) ! inout
          call p__photochem_transport__exe(spl, cst, grd, set, & ! in
            &                              var                 ) ! inout
          call p__photochem_scheme__exe(spl, grd, set, & ! in
            &                           var            ) ! inout

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              var%ni_3d(isp,ix,iy,iz) = var%ni(isp,iz)
            end do
          end do

          if ( ix == grd%nx ) then
            do iz = 1, grd%nz
              do isp = 1, spl%nsp
                var%ni_3d(isp,1,iy,iz) = var%ni(isp,iz)
              end do
            end do
          end if

          write(*,*) iday, grd%iy, grd%ix

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

    call p__io_rotation__fin(set, spl, var, grd) ! in

  end if


  !----------------------------------------------------------------------------------------------------------
  !
  !                                            Finalize: deallocate
  !
  !----------------------------------------------------------------------------------------------------------

  call p__io_deallocate(var, spl, xct, flx) ! inout

  ! Call cpu time
  call cpu_time(t3)

  print *, '----------------------------------------------------------------------------------------------------------'

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

  print *, 'Finish.'
  print *, ''


end program e__main
