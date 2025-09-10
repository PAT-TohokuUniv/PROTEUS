module p__photochem_scheme
  
  use v__tdec,            only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__search,          only : p__search_reactant, p__search_product, sp_index
  use p__io,              only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                            p__io_progress
  use p__chemical_scheme, only : p__chemical_scheme__exe, p__LU_decomposition, p__LU_solver

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private

  ! Module-level variables ------------------------------------------------------------------------
  real(dp),           allocatable, private :: n_var(:,:), n_fix(:,:), delta_ni(:,:)
  real(dp),           allocatable, private :: dni_ni(:,:), n_lower_boundary(:,:)
  character(len=256), allocatable, private :: species_var(:), species_fix(:)
  character(len=256), allocatable, private :: var_species_list(:), fix_species_list(:)

  ! Public interface ------------------------------------------------------------------------------
  public :: p__photochem_scheme__ini, p__photochem_scheme__exe

contains


  subroutine p__photochem_scheme__ini(spl, grd) ! in
    implicit none
    type(spl_), intent(in) :: spl
    type(grd_), intent(in) :: grd
    integer nsp, nsp_i, nsp_f, nch, nz

    nsp = spl%nsp
    nsp_i = spl%nsp_i
    nsp_f = spl%nsp - nsp_i
    nch   = spl%nch
    nz    = grd%nz

    allocate(n_var(nz, nsp_i), n_fix(nz, nsp_f), delta_ni(nz, nsp_i))
    allocate(dni_ni(nz, nsp), n_lower_boundary(nsp_i, 2))
    allocate(species_var(nsp_i), species_fix(nsp_f))
    allocate(var_species_list(nsp_i), fix_species_list(nsp_f))

  end subroutine p__photochem_scheme__ini


  subroutine p__photochem_scheme__exe(spl, grd, set, & ! in
    &                                 var            ) ! inout
    implicit none
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var

    integer isp, jsp, ich, iz
    integer i
    real(dp) tmp
    real(dp) lambda ! implicit factor
    real(dp) eps, dt_rate
    real(dp) max_eps, dt_min
    character(len=256) model

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    model = set%inversion

    ! for conversion 
    eps     = set%dt_inc_eps
    dt_rate = 1.0_dp + set%dt_rate / 100.0_dp
    max_eps = set%max_eps ! go to next timestep without applying the changes and re-calculate with a smaller dt if max(dn/n) > max_eps
    dt_min  = 1.0_dp ! ignore max_eps limitation if dt < dt_min

    ! advance time step scheme ------------------------------------------------------------------

    lambda   = 1.0_dp
    var%tAmtx = 0.0_dp
    var%Jmtx = 0.0_dp
    var%barr = 0.0_dp
    var%xarr = 0.0_dp
    
    i = 0
    species_var = ''
    species_fix = ''
    n_var = 0.0_dp
    n_fix = 0.0_dp
    do isp = 1, spl%nsp
      jsp = spl%all_to_var(isp)
      if (jsp /= 0) then 
        n_var(1:grd%nz,jsp) = var%ni(1:grd%nz,isp) * 1.0e-6_dp
        species_var(jsp) = spl%species(isp)
      end if
      
      if (jsp == 0) then 
        i = i + 1
        n_fix(1:grd%nz,i) = var%ni(1:grd%nz,isp) * 1.0e-6_dp
        species_fix(i) = spl%species(isp)
      end if
    end do 

    do isp = 1, spl%nsp_i
      i = spl%var_to_all(isp)
      n_lower_boundary(isp,1) = var%Lower_n(i,1)
      n_lower_boundary(isp,2) = var%Lower_n(i,2)
    end do 

    var%Pi = var%Pi * 1.0e-6_dp
    var%Li = var%Li * 1.0e-6_dp
    var%dphi_dz = var%dphi_dz * 1.0e-6_dp

    call  p__chemical_scheme__exe(set%scheme,                                                 & ! in:  
      &                           species_var, species_fix,                                   & ! in:    name of variable species and fixed species
      &                           spl%nsp_i, spl%nsp-spl%nsp_i, spl%nch, grd%nz,              & ! in:    
      &                           n_var, n_fix,                                               & ! in:    number density of variable and fixed species
      &                           n_lower_boundary,                                           & ! in:    reaction rate coefficients
      &                           var%Pi, var%Li, var%ki, var%dphi_dz,                        & ! in:  
      &                           var%d_dniu_dphi_dz, var%d_dni0_dphi_dz, var%d_dnil_dphi_dz, & ! in:  
      &                           var%tAmtx, var%Umtx, var%tLmtx,                             & ! inout:  
      &                           var%dtime,                                                  & ! inout: timestep
      &                           delta_ni                                                    ) ! out:   change of number density of variable species

    do iz = 1, grd%nz
      do isp = 1, spl%nsp
        i = spl%all_to_var(isp)
        var%ni_new(iz,isp) = 0.0_dp
        if (spl%label_fix(isp) == 0) then
          var%ni_new(iz,isp) = var%ni(iz,isp) + delta_ni(iz,i) * 1.0e6_dp
        else if (spl%label_fix(isp) == 1) then
          var%ni_new(iz,isp) = var%ni(iz,isp)
        end if
      end do
    end do

    do iz = 1, grd%nz
      do isp = 1, spl%nsp
        if (var%ni_new(iz,isp) < 0.0_dp) then
          var%ni_new(iz,isp) = 0.0_dp
        end if
        if (var%ni_new(iz,isp) < 1.0e-20_dp) then
          var%ni_new(iz,isp) = 1.0e-20_dp
        end if
      end do
    end do
        

    ! electron density : charge neutrality ---------------------------------------------------------
    do iz = 1, grd%nz
      do isp = 1, spl%nsp
        if (spl%species(isp) == 'e-') then
          var%ni_new(iz,isp) = 0.0_dp
          do jsp = 1, spl%nsp
            if (var%q(jsp) > 0.0_dp .and. jsp /= isp) then
              var%ni_new(iz,isp) &
                &  = var%ni_new(iz,isp) + var%ni_new(iz,jsp)
            end if
            if (var%q(jsp) < 0.0_dp .and. jsp /= isp) then
              var%ni_new(iz,isp) &
                &  = var%ni_new(iz,isp) - var%ni_new(iz,jsp)
            end if
          end do
        end if
      end do
    end do
  

    ! advance time step ------------------------------------------------------------------
    ! lower boundary condition: density
    do isp = 1, spl%nsp
      if (nint(var%Lower_n(isp,1)) == 1) then 
        var%ni_new(1,isp) = var%Lower_n(isp,2)
      end if
    end do

   ! convergence check and increase time step
    if (set%calc_stable == 1 .and. set%start_rot == 0) then
      if ( var%istep >= 2 ) then

        var%max_dn_n = 0.0_dp
        do iz = 1, grd%nz
          do isp = 1, spl%nsp
            ! delta_N / N
            if ( var%ni(iz,isp) /= 0.0d0 .and. spl%label_fix(isp)==0) then
              if(var%ni(iz,isp) < var%ni_new(iz,isp)) then 
                dni_ni(iz,isp) = dabs(var%ni_new(iz,isp) - var%ni(iz,isp)) / var%ni(iz,isp)
              else if (var%ni(iz,isp) >= var%ni_new(iz,isp)) then
                dni_ni(iz,isp) = dabs(var%ni_new(iz,isp) - var%ni(iz,isp)) / var%ni_new(iz,isp)
              end if
            else if ( var%ni(iz,isp) == 0.0d0 .and. spl%label_fix(isp)==0) then
              dni_ni(iz,isp) = 0.0d0
            end if
            if ( dni_ni(iz,isp) > var%max_dn_n(3) .and. spl%label_fix(isp)==0) then
              var%max_dn_n(1) = isp
              var%max_dn_n(2) = iz
              var%max_dn_n(3) = dni_ni(iz,isp)
              tmp = dble(isp)
            end if
          end do
        end do

        ! if max(dN/N) < 10 %, then goes to larger dt (*10)
        ! if dt becomes 100000 [sec], unstable at HC layer
        if ( set%mode == '1D' .or. set%mode == '2DLAT' ) then 
          if ( var%max_dn_n(3) > max_eps .and. var%istep >= 2 .and. var%dtime >= dt_min) then
            var%dtime = var%dtime / (1.0_dp + max_eps)
          end if
          if ( var%dtime < set%dtime_limit * 0.99_dp ) then
            if ( var%max_dn_n(3) < eps ) then
              if (set%scheme == 'implicit') then
                var%dtime = var%dtime * dt_rate
              else if (set%scheme == 'semi-implicit') then
                var%dtime = var%dtime * dt_rate
              else if (set%scheme == 'explicit') then
                var%dtime = var%dtime * dt_rate
              end if
              var%iter = 0
              do isp = 1, spl%nsp
                var%ni_0(1:grd%nz,isp) = var%ni_new(1:grd%nz,isp)
              end do
            else if ( var%max_dn_n(3) >= eps ) then
              var%iter = var%iter + 1
              if ( var%iter > 5000 ) then
                var%dtime = 1.0e-8_dp
                var%iter = 0
                do isp = 1, spl%nsp
                  var%ni_0(1:grd%nz,isp) = var%ni_new(1:grd%nz,isp)
                end do
              end if
            end if
          end if
          if ( var%dtime  >= set%dtime_limit * 0.99_dp ) then
            var%dtime = set%dtime_limit ! dt DO NOT excess set%dtime_limit <- set at v__'planet'__ini
          end if
        end if

      end if
    end if



    ! advance time step ------------------------------------------------------------------
    do isp = 1, spl%nsp
      do iz = 1, grd%nz
        if (var%ni_new(iz,isp) /= var%ni_new(iz,isp)) then
          var%ni(iz,isp) = var%ni_0(iz,isp)
          var%ni_new(iz,isp) = var%ni_0(iz,isp)
          var%dtime = 1.0e-8_dp
          var%iter = 0
          jsp = spl%all_to_var(isp)
          print *, iz, trim(spl%species(isp)), var%Pi(iz,jsp), var%Li(iz,jsp) 
          do ich = 1, spl%nch
            print *, ich, var%ki(iz,ich)
          end do
          call p__io_stable__fin(set, spl, var, grd) ! in
          stop
        end if
        !var%ni_0(iz,isp) = var%ni(iz,isp)
      end do
    end do

    if ( var%max_dn_n(3) <= max_eps .or. var%istep == 1 .or. var%dtime < dt_min) then
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          var%ni(iz,isp) = var%ni_new(iz,isp)
        end do
      end do
      var%sum_time = var%sum_time + var%dtime
    end if


  end subroutine p__photochem_scheme__exe

end module p__photochem_scheme