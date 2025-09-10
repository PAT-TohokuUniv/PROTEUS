module p__io

  use v__tdec, only : set_, var_, spl_, xct_, flx_, grd_
  use p__search,  only : p__search_reactant, p__search_product, ch_identify, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private

  ! Private variables -----------------------------------------------------------------------------
  character(len=256), allocatable, private :: reaction_out(:,:)
  integer,            allocatable, private :: ch_out(:)
  real(dp),           allocatable, private :: t_out(:), cln_out(:,:), n_out(:,:,:), vmr_out(:,:,:)
  real(dp),           allocatable, private :: rate_out(:,:,:), Pi_out(:,:,:), Li_out(:,:,:)
  real(dp),           allocatable, private :: flux_down_out(:,:,:), flux_up_out(:,:,:)

  ! Public interface ------------------------------------------------------------------------------
  public :: p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &       p__io_progress, p__io_timeseries__ini, p__io_timeseries__exe, p__io_timeseries__fin


contains


  subroutine p__io_stable__fin(set, spl, var, grd)
    type(set_), intent(in) :: set
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    character(len=256) fname, num
    integer iy, iz, isp, jsp, ich
    real(dp) alt, ni
    character(len=256) outdir, command

    outdir = trim(ADJUSTL(set%dir_name))//'/output'
    write(command,*) 'mkdir ', trim(outdir)
    call system(command)

    outdir = trim(ADJUSTL(set%dir_name))//'/output/density'
    write(command,*) 'mkdir ', trim(outdir)
    call system(command)
    

    if (set%mode == '1D') then 

      outdir = trim(ADJUSTL(set%dir_name))//'/output/density/num'
      write(command,*) 'mkdir ', trim(outdir)
      call system(command)
      
      outdir = trim(ADJUSTL(set%dir_name))//'/output/density/vmr'
      write(command,*) 'mkdir ', trim(outdir)
      call system(command)

      outdir = trim(ADJUSTL(set%dir_name))//'/output/flux'
      write(command,*) 'mkdir ', trim(outdir)
      call system(command)

      outdir = trim(ADJUSTL(set%dir_name))//'/output/rate'
      write(command,*) 'mkdir ', trim(outdir)
      call system(command)

      fname = trim(ADJUSTL(set%dir_name))//'/output/density/num/*.dat'
      write(command,*) 'rm ', trim(fname)
      call system(command)

      fname = trim(ADJUSTL(set%dir_name))//'/output/density/vmr/*.dat'
      write(command,*) 'rm ', trim(fname)
      call system(command)

      fname = trim(ADJUSTL(set%dir_name))//'/output/flux/*.dat'
      write(command,*) 'rm ', trim(fname)
      call system(command)

      fname = trim(ADJUSTL(set%dir_name))//'/output/rate/*.dat'
      write(command,*) 'rm ', trim(fname)
      call system(command)

      do isp = 1, spl%nsp
        fname = trim(ADJUSTL(set%dir_name))//'/output/density/num/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%ni(iz,isp)
          end do
        close(11)

        fname = trim(ADJUSTL(set%dir_name))//'/output/density/vmr/vmr_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%ni(iz,isp)/var%n_tot(iz)
          end do
        close(11)

        if (spl%label_fix(isp) == 0) then

          fname = trim(ADJUSTL(set%dir_name))//'/output/flux/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'replace' )
            do iz = 1, grd%nz+1
              if (iz==1) then 
                alt = (grd%alt(1) - grd%dalt(1)/2.0d0) /1.0e3_dp
              end if
              if (iz>=2 .and. iz<=grd%nz) then 
                alt = (grd%alt(iz) + grd%alt(iz-1)) /2.0d0 /1.0e3_dp
              end if
              if (iz==grd%nz+1) then 
                alt = (grd%alt(grd%nz) + grd%dalt(grd%nz)/2.0d0) /1.0e3_dp
              end if
              write(11, fmt='(f10.3, E20.10E4, E20.10E4)') alt, &
                & var%Fluxup(iz,spl%all_to_var(isp)), &
                & var%Fluxdwn(iz,spl%all_to_var(isp))
            end do
          close(11)

          fname = trim(ADJUSTL(set%dir_name))//'/output/flux/v_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'replace' )
            do iz = 1, grd%nz+1
              if (iz==1) then 
                alt = (grd%alt(1) - grd%dalt(1)/2.0d0) /1.0e3_dp
                ni = var%ni(isp,1)
              end if
              if (iz>=2 .and. iz<=grd%nz) then 
                alt = (grd%alt(iz) + grd%alt(iz-1)) /2.0d0 /1.0e3_dp
                ni = var%ni(iz,isp)
              end if
              if (iz==grd%nz+1) then 
                alt = (grd%alt(grd%nz) + grd%dalt(grd%nz)/2.0d0) /1.0e3_dp
                ni = var%ni(grd%nz,isp)
              end if
              write(11, fmt='(f10.3, E20.10E4, E20.10E4)') alt, &
                & var%Fluxup(iz,spl%all_to_var(isp))/ni, &
                & var%Fluxdwn(iz,spl%all_to_var(isp))/ni
            end do
          close(11)

          fname = trim(ADJUSTL(set%dir_name))//'/output/flux/D_'//trim(ADJUSTL(spl%species(isp)))//'.dat'
          open(11, file = fname, status = 'replace' )
            do iz = 1, grd%nz
              write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%D_mol(iz,isp)
            end do
          close(11)

        end if
      end do

      do ich = 1, spl%nch
        write(num,'(I3)') ich
        fname = trim(ADJUSTL(set%dir_name))//'/output/rate/k'//trim(ADJUSTL(num))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%ki(iz,ich)
          end do
        close(11)
      end do

      do ich = 1, spl%nch
        write(num,'(I3)') ich
        fname = trim(ADJUSTL(set%dir_name))//'/output/rate/rate'//trim(ADJUSTL(num))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%rate(iz,ich)
          end do
        close(11)
      end do

      do isp = 1, spl%nsp_i
        jsp = spl%var_to_all(isp)
        fname = trim(ADJUSTL(set%dir_name))//'/output/rate/P_'//trim(ADJUSTL(spl%species(jsp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%Pi(iz,isp) * 1.0e6_dp
          end do
        close(11)

        fname = trim(ADJUSTL(set%dir_name))//'/output/rate/L_'//trim(ADJUSTL(spl%species(jsp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%Li(iz,isp) * 1.0e6_dp
          end do
        close(11)
      end do

      fname = trim(ADJUSTL(set%dir_name))//'/output/flux/K_eddy.dat'
      open(11, file = fname, status = 'replace' )
        do iz = 1, grd%nz
          write(11, fmt='(f10.3, E15.5E4)') grd%alt(iz)/1.0e3_dp, var%K_eddy(iz)
        end do
      close(11)



    else if (set%mode == '2DLAT') then 

      outdir = trim(ADJUSTL(set%dir_name))//'/output/density/2DLAT'
      write(command,*) 'mkdir ', trim(outdir)
      call system(command)

      fname = trim(ADJUSTL(set%dir_name))//'/output/density/2DLAT/*.dat'
      write(command,*) 'rm ', trim(fname)
      call system(command)

      do isp = 1, spl%nsp
        fname = trim(ADJUSTL(set%dir_name))//'/output/density/2DLAT/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do iz = 1, grd%nz
            write(11, fmt='(f10.3)', advance='no') grd%alt(iz)/1.0e3_dp
            write(11, fmt='(E15.5E4)') (var%ni_stable(iy,iz,isp), iy = 1, grd%ny)
          end do
        close(11)
      end do

      fname = trim(ADJUSTL(set%dir_name))//'/output/density/2DLAT/grid.dat'
      open(11, file = fname, status = 'replace' )
        write(11, *) grd%nx, grd%ny, grd%nz
      close(11)

    end if


  end subroutine p__io_stable__fin


  subroutine p__io_timeseries__ini(set, grd, spl, var)
    implicit none
    type(set_), intent(in) :: set
    type(grd_), intent(in) :: grd
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    integer i, iz, isp, ich, nch_out
    real(dp) cln_tmp

    nch_out = 0!spl%nch
    allocate(reaction_out(nch_out,4),ch_out(nch_out))
    do ich = 1, nch_out
      ch_out(ich) = ich
    end do 
    allocate(t_out(0:set%nstep))
    allocate(cln_out(set%nsp_tout,0:set%nstep),n_out(grd%nz,set%nsp_tout,0:set%nstep),vmr_out(grd%nz,set%nsp_tout,0:set%nstep))
    allocate(rate_out(grd%nz,nch_out,0:set%nstep))
    allocate(Pi_out(grd%nz,set%nsp_tout,0:set%nstep),Li_out(grd%nz,set%nsp_tout,0:set%nstep))
    allocate(flux_down_out(grd%nz+1,set%nsp_tout,0:set%nstep),flux_up_out(grd%nz+1,set%nsp_tout,0:set%nstep))

    t_out(0) = 0.0_dp
    do i = 1, set%nsp_tout
      isp = sp_index(spl,trim(ADJUSTL(set%species_tout(i))))
      if (isp >= 1 .and. isp <= spl%nsp) then 
        cln_tmp = 0.0_dp
        do iz = 1, grd%nz
          cln_tmp = cln_tmp + var%ni(iz,isp)*grd%dalt(iz)
          n_out(iz,i,0) = var%ni(iz,isp)
          vmr_out(iz,i,0) = var%ni(iz,isp) / var%n_tot(iz)
        end do 
        cln_out(i,0) = cln_tmp
      end if
    end do

  end subroutine p__io_timeseries__ini


  subroutine p__io_timeseries__exe(set, grd, spl, var, iout)
    implicit none
    type(set_), intent(inout) :: set
    type(grd_), intent(in) :: grd
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    integer, intent(inout) :: iout
    integer i, isp, jsp

    if (set%dt_out < var%dtime) set%dt_out = var%dtime
      if (var%sum_time - t_out(iout) >= set%dt_out) then 
        iout = iout + 1
        ! output ==================================
        t_out(iout) = var%sum_time
        do i = 1, set%nsp_tout
          isp = sp_index(spl,trim(ADJUSTL(set%species_tout(i))))
          jsp = spl%all_to_var(isp)
          if (isp >= 1 .and. isp <= spl%nsp) then 
            cln_out(i,iout) = var%clm_ni(1,isp)
            n_out(1:grd%nz,i,iout) = var%ni(1:grd%nz,isp)
            vmr_out(1:grd%nz,i,iout) = var%ni(1:grd%nz,isp) / var%n_tot(1:grd%nz)
            Pi_out(1:grd%nz,i,iout) = var%Pi(1:grd%nz,jsp)
            Li_out(1:grd%nz,i,iout) = var%Li(1:grd%nz,jsp)
            flux_up_out(1:grd%nz+1,i,iout) = var%fluxup(1:grd%nz+1,jsp)
            flux_down_out(1:grd%nz+1,i,iout) = var%fluxdwn(1:grd%nz+1,jsp)
          end if
        end do

        !! Reaction rate timeseries output is not supported
        !do i = 1, nch_out
        !  do iz = 1, grd%nz
        !    ich = ch_out(i)
        !    rate_out(iz,i,iout) = var%rate(iz,ich)
        !  end do 
        !end do 
        
      end if
    

  end subroutine p__io_timeseries__exe


  subroutine p__io_timeseries__fin(set, grd, spl, iout) ! in
    implicit none
    type(set_), intent(in) :: set
    type(grd_), intent(in) :: grd
    type(spl_), intent(in) :: spl
    integer,    intent(in) :: iout
    character(len=256) fname, outdir, command, ci
    integer iz, i, isp, jsp, ich
    real(dp) tmp, n_tmp

    outdir = trim(ADJUSTL(set%dir_name))//'/output'
    write(command,*) 'mkdir ', trim(outdir)
    call system(command)

    outdir = trim(ADJUSTL(set%dir_name))//'/output/timeseries'
    write(command,*) 'mkdir ', trim(outdir)
    call system(command)

    fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/*.dat'
    write(command,*) 'rm ', trim(fname)
    call system(command)

    fname = trim(ADJUSTL(set%dir_name))//'/output/alt_density.dat'
    open(11, file = fname, status = 'replace' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz) / 1.0e3_dp
      end do
    close(11)

    fname = trim(ADJUSTL(set%dir_name))//'/output/alt_flux.dat'
    open(11, file = fname, status = 'replace' )
      do iz = 1, grd%nz+1
        if (iz == 1) tmp = grd%alt(1) - grd%dalt(1)/2.0_dp
        if (iz >= 2) tmp = grd%alt(iz-1) + grd%dalt(iz-1)/2.0_dp
        write(11, *) tmp / 1.0e3_dp
      end do
    close(11)

    fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/timestamp.dat'
    open(11, file = fname, status = 'replace' )
      do i = 0, iout
        write(11, *) t_out(i)
      end do
    close(11)

    do isp = 1, set%nsp_tout
      jsp = sp_index(spl,trim(ADJUSTL(set%species_tout(isp))))
      if (jsp >= 1 .and. jsp <= spl%nsp) then 
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/n_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            write(11, *) t_out(i), (n_out(iz,isp,i), iz = 1, grd%nz)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/cln_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            write(11, *) t_out(i), cln_out(isp,i)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/vmr_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            write(11, *) t_out(i), (vmr_out(iz,isp,i), iz = 1, grd%nz)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/Pi_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) write(11, *) t_out(i), (Pi_out(iz,isp,i)*1.0e6_dp, iz = 1, grd%nz)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/Li_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) write(11, *) t_out(i), (Li_out(iz,isp,i)*1.0e6_dp, iz = 1, grd%nz)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/flux_down_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) write(11, *) t_out(i), (flux_down_out(iz,isp,i), iz = 1, grd%nz+1)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/flux_up_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) write(11, *) t_out(i), (flux_up_out(iz,isp,i), iz = 1, grd%nz+1)
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/v_down_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) then 
              write(11, fmt='(E18.10E3)', advance='no') t_out(i)
              do iz = 1, grd%nz+1
                if (iz == 1) n_tmp = n_out(1,isp,i)
                if (iz >= 2 .and. iz <= grd%nz) n_tmp = ( n_out(iz,isp,i) + n_out(iz-1,isp,i) ) / 2.0_dp
                if (iz == grd%nz+1) n_tmp = n_out(grd%nz,isp,i)
                write(11, fmt='(E18.10E3)', advance='no') flux_down_out(iz,isp,i)/n_tmp
              end do
              write(11,*)
            end if
          end do
        close(11)
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/v_up_'//trim(ADJUSTL(set%species_tout(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            if (spl%all_to_var(jsp)/=0) then 
              write(11, fmt='(E18.10E3)', advance='no') t_out(i)
              do iz = 1, grd%nz+1
                if (iz == 1) n_tmp = n_out(1,isp,i)
                if (iz >= 2 .and. iz <= grd%nz) n_tmp = ( n_out(iz,isp,i) + n_out(iz-1,isp,i) ) / 2.0_dp
                if (iz == grd%nz+1) n_tmp = n_out(grd%nz,isp,i)
                write(11, fmt='(E18.10E3)', advance='no') flux_up_out(iz,isp,i)/n_tmp
              end do
              write(11,*)
            end if
          end do
        close(11)
      end if
    end do

    if (set%rate_tout == 1) then 
      do ich= 1, spl%nch
        write(ci,*) ich
        fname = trim(ADJUSTL(set%dir_name))//'/output/timeseries/rate_'//trim(ADJUSTL(ci))//'.dat'
        open(11, file = fname, status = 'replace' )
          do i = 0, iout
            write(11, *) t_out(i), (rate_out(iz,ich,i)*1.0e6_dp, iz = 1, grd%nz)
          end do
        close(11)
      end do
    end if

  end subroutine p__io_timeseries__fin

  subroutine p__io_rotation__fin(set, spl, var, grd)
    type(set_), intent(in) :: set
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    character(len=256) fname, command, outdir
    character(len=256) species(spl%nsp)
    real(dp) lat, lt, m(spl%nsp), q(spl%nsp)
    integer ix, iy, iz, isp, nx, ny, nz, nsp
    real(dp) alt(grd%nz), dalt(grd%nz)
    namelist /grid/ nx, ny, nz, nsp
    namelist /altitude/ alt, dalt
    namelist /species_list/ species, m, q

    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/altitude.dat'
    open(11, file = fname, status = 'replace' )
      do iz = 1, grd%nz
        write(11, *) grd%alt(iz) / 1.0e3_dp
      end do
    close(11)

    nsp = spl%nsp
    species(1:nsp) = spl%species(1:nsp)
    m(1:nsp) = var%m(1:nsp)
    q(1:nsp) = var%q(1:nsp)
    nx = grd%nx
    ny = grd%ny
    nz = grd%nz
    alt = grd%alt 
    dalt = grd%dalt 
    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/grid.out'
    open(11, file = fname, status = 'replace' )
    write(11,nml=grid)
    close(11)
    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/altitude.out'
    open(11, file = fname, status = 'replace' )
    write(11,nml=altitude)
    close(11)
    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/species.out'
    open(11, file = fname, status = 'replace' )
    write(11,nml=species_list)
    close(11)

    outdir = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))
    write(command,*) 'mkdir ', trim(outdir)
    call system(command)

    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/*.dat'
    write(command,*) 'rm ', trim(fname)
    call system(command)

    if (set%mode == '2DROT') then 
      do isp = 1, spl%nsp
        fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/'//trim(ADJUSTL(spl%species(isp)))//'.dat'
        open(11, file = fname, status = 'replace' )
          do ix = 1, grd%nx
            lt = 24.0_dp*dble(ix-1)/dble(grd%nx-1)
            do iy = 1, grd%ny
              if (grd%ny == 1) lat = set%latitude 
              if (grd%ny /= 1) lat = 180.0_dp*dble(iy-((grd%ny+1)/2))/dble(grd%ny-1)
              do iz = 1, grd%nz
                write(11,*) lt, lat, grd%alt(iz)/ 1.0e3_dp, var%ni_3d(ix,iy,iz,isp)
              end do
            end do
          end do
        close(11)
      end do 
    end if

    if (set%mode == '3DROT') then 
      fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/density.dat'
      open(11, file = fname, status = 'replace', action = 'write', form = 'unformatted')
        write(11) var%ni_3d(1:grd%nx,1:grd%ny,1:grd%nz,1:spl%nsp)
      close(11)
    end if

    fname = trim(ADJUSTL(set%dir_name))//'/output/density/'//trim(ADJUSTL(set%mode))//'/grid.dat'
    open(11, file = fname, status = 'replace' )
      write(11, *) grd%nx, grd%ny, grd%nz
    close(11)


  end subroutine p__io_rotation__fin


  subroutine p__io_progress(spl, var, grd, set)
    type(spl_), intent(in) :: spl
    type(var_), intent(in) :: var
    type(grd_), intent(in) :: grd
    type(set_), intent(in) :: set
    character(len=256) fname, char

    fname = trim(ADJUSTL(set%dir_name))//'/progress.dat'
    open(11, file = fname, status = 'replace' )

      write(11, *) '########################################'
      write(11, *)
      write(11, *) '            ',trim(ADJUSTL(spl%planet))
      write(11, *)
      write(11, *) '########################################'
      write(11, *)
      write(11, *)

      if ( set%mode == '1D' ) then

        write(11, *) '----------------------------------------'
        write(11, *) '         1D stable calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculate stable solution...'
          if ( grd%iy < grd%ny ) then
            write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
            write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
        end if

      else if ( set%mode == '2D Lat' ) then

        write(11, *) '----------------------------------------'
        write(11, *) '         2D stable calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculate stable solution...'
          if ( grd%iy < grd%ny ) then
            write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
            write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
        end if

      else if (set%mode == '2D Rot') then
        write(11, *) '----------------------------------------'
        write(11, *) '       2D rotational calculation     '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          if ( set%start_rot == 0 ) then
            write(11, *) 'calculating stable solution...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished stable calculation!'
            end if
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
          end if
          if ( set%start_rot == 1 ) then
            write(11, *) 'finished stable calculation!'
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
            if ( grd%iy < grd%ny ) then
              write(11, *) 'calculating rotational mode...'
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        end if

      else if (set%mode == '3D Rot') then
        write(11, *) '----------------------------------------'
        write(11, *) '       3D rotational calculation     '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          if ( set%start_rot == 0 ) then
            write(11, *) 'calculating stable solution...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished stable calculation!'
            end if
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
          end if
          if ( set%start_rot == 1 ) then
            write(11, *) 'finished stable calculation!'
            write(11,*) 'dt     = ',var%dtime, ' [s]'
            write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
            write(11,*)
            if ( grd%iy < grd%ny ) then
              write(11, *) 'calculating rotational mode...'
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t4-var%t3, ' [s]'
          end if
        end if

      else if ( set%mode == '3D Global' ) then
        write(11, *) '----------------------------------------'
        write(11, *) '         3D global calculation          '
        write(11, *) '----------------------------------------'
        write(11, *)

        if ( set%calc_stable == 0 ) then
          write(11, *) 'skipped stable solution...'
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(1,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          end if
        else if ( set%calc_stable == 1 ) then
          write(11, *) 'calculating stable solution...'
          if ( grd%iy < grd%ny ) then
            write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
            write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
          else if ( grd%iy == grd%ny ) then
            write(11, *) 'finished stable calculation!'
          end if
          write(11,*) 'dt     = ',var%dtime, ' [s]'
          write(11,*) 'sum_t  = ',var%sum_time, ' [s]'
          write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          write(11,*)
          if ( set%start_rot == 1 ) then
            write(11, *) 'calculating rotational mode...'
            if ( grd%iy < grd%ny ) then
              write(char,'(f5.1)') 100.0_dp * dble(grd%iy)/dble(grd%ny)    
              write(11, *) 'current progress = '//trim(ADJUSTL(char))//' %'
            else if ( grd%iy == grd%ny ) then
              write(11, *) 'finished rotational calculation!'
            end if
            write(11,*) 'cpu_t  = ',var%t2-var%t1, ' [s]'
          end if
        end if

      end if


    close(11)


  end subroutine p__io_progress


  subroutine p__io_deallocate(var, spl, xct, flx) ! inout
    type(var_),           intent(inout) :: var
    type(spl_),           intent(inout) :: spl
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx

    ! deallocate
    deallocate(var%ni)
    deallocate(var%ni_stable)
    deallocate(var%ni_3d)
    deallocate(var%ni_new)
    deallocate(var%clm_ni)
    deallocate(var%Ti,var%Te,var%Tn)
    deallocate(var%Ti_3d,var%Te_3d,var%Tn_3d)
    deallocate(var%m, var%q)
    deallocate(spl%reactant_list)
    deallocate(spl%product_list)
    deallocate(spl%species)
    deallocate(spl%label_fix)
    deallocate(spl%all_to_var)
    deallocate(spl%var_to_all)
    deallocate(spl%reaction_type_list)
    deallocate(var%ki)
    deallocate(var%Pi, var%Pij)
    deallocate(var%Li)
    deallocate(var%K_eddy,var%D_mol)
    deallocate(var%Fluxup,var%Fluxdwn)
    deallocate(var%Jmtx, var%tAmtx, var%tLmtx, var%Umtx)
    deallocate(spl%rate_cases)
    deallocate(spl%T_range)
    deallocate(spl%major_species)
    deallocate(var%n_tot,var%m_mean)

    deallocate(flx%solar_EUV)
    deallocate(var%tau_EUV,var%I_EUV)
    deallocate(xct%sigma_a_EUV, xct%sigma_i_EUV)


  end subroutine p__io_deallocate


end module p__io
