! p__photochem_scheme.f90
!
!
!
!

module p__photochem_scheme
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__search,      only : p__search_reactant, p__search_product, sp_index
  use p__io,          only : p__io_stable__fin, p__io_rotation__fin, p__io_deallocate, &
    &                        p__io_progress

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__photochem_scheme__exe

contains


  subroutine p__photochem_scheme__exe(spl, grd, set, & ! in
    &                                 var            ) ! inout
    implicit none
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var

    integer isp, jsp, ich, jch, iz
    integer i, j, k, ijch, ilist, ii, jj
    integer bij
    real(dp) a0, nj, nk, kij, tmp
    real(dp) n0, nl, nu
    real(dp) dn(spl%nsp_i, grd%nz), dni_ni(spl%nsp,grd%nz)
    real(dp) lambda ! implicit factor
    real(dp) eps, dt_rate
    character(len=256) model

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    model = set%inversion

    ! for conversion 
    if (set%scheme == 'implicit') then
      eps     = 1.0e-1_dp 
      dt_rate = 1.0e1_dp
    else if (set%scheme == 'semi-implicit') then 
      eps     = 1.0e-2_dp 
      dt_rate = 1.01_dp
    else if (set%scheme == 'explicit') then 
      eps     = 1.0e-2_dp 
      dt_rate = 1.01_dp
    end if

    ! advance time step scheme ------------------------------------------------------------------

      !---------------------------------------------------------
      ! explicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'explicit') then

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              jsp = spl%all_to_var(isp)
              if (spl%label_fix(isp) == 0) then
                dn(jsp,iz) = ( var%Pi(jsp,iz) - var%Li(jsp,iz) - var%dPhi_dz(jsp,iz) ) * var%dtime
                var%ni_new(isp,iz) = var%ni(isp,iz) + dn(jsp,iz)
              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz)  = var%ni(isp,iz)
              end if
            end do
          end do

        end if

      !---------------------------------------------------------
      ! semi-implicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'semi-implicit') then

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              jsp = spl%all_to_var(isp)
              var%ni_new(isp,iz) = 0.0_dp
              if (spl%label_fix(isp) == 0) then
                if (var%ni(isp,iz) /= 0.0_dp) then
                  var%ni_new(isp,iz) = ( var%ni(isp,iz) + (var%Pi(jsp,iz) - var%dPhi_dz(jsp,iz)) * var%dtime ) &
                    &                  / (1.0_dp + (var%Li(jsp,iz) / var%ni(isp,iz)) * var%dtime )
                end if

                !n0 = var%ni(isp,iz)
                !if (iz>=2) nl = var%ni(isp,iz-1); if (iz==1) nl = 0.0_dp
                !if (iz<=grd%nz-1) nu = var%ni(isp,iz+1); if (iz==grd%nz) nu = 0.0_dp
                !if (var%ni(isp,iz) /= 0.0_dp) then
                !  var%ni_new(isp,iz) = ( var%ni(isp,iz) + ( var%Pi(jsp,iz) &
                !    &                                       - var%d_dnil_dPhi_dz(jsp,iz) * nl &
                !    &                                       - var%d_dniu_dPhi_dz(jsp,iz) * nu ) * var%dtime ) &
                !    &                  / (1.0_dp + (var%Li(jsp,iz) / var%ni(isp,iz) + var%d_dni0_dPhi_dz(jsp,iz)) * var%dtime )
                !end if

                if (var%ni(isp,iz) == 0.0_dp) then
                  var%ni_new(isp,iz) = ( var%ni(isp,iz) + (var%Pi(jsp,iz) - var%dPhi_dz(jsp,iz)) * var%dtime ) &
                    &                  / (1.0_dp)
                end if

                !if (var%ni(isp,iz) == 0.0_dp) then
                !  var%ni_new(isp,iz) = ( var%ni(isp,iz) + ( var%Pi(jsp,iz) &
                !    &                                       - var%d_dnil_dPhi_dz(jsp,iz) * nl &
                !    &                                       - var%d_dniu_dPhi_dz(jsp,iz) * nu ) * var%dtime ) &
                !    &                  / (1.0_dp)
                !end if

              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz)  = var%ni(isp,iz)
              end if

            end do
          end do

        end if

      !---------------------------------------------------------
      ! implicit scheme
      !---------------------------------------------------------
        if (set%scheme == 'implicit') then

          lambda   = 1.0_dp
          var%tAmtx = 0.0_dp
          var%Jmtx = 0.0_dp
          var%barr = 0.0_dp
          var%xarr = 0.0_dp

          do iz = 1, grd%nz

            ! generate Chemical Jacobian matrix
            !   Jmtx(i,j) = dFi/dnj, Fi = Pi - Li - dPhi_i/dz
            var%Jmtx = 0.0_dp
            do ilist = 1, spl%n_Jlist
              i = spl%Jmtx_list(ilist,1)
              j = spl%Jmtx_list(ilist,2)
              ii = spl%var_to_all(i)
              jj = spl%var_to_all(j)
              do ich = 1, spl%Jmtx_list(ilist,3)
                ijch = spl%Jmtx_list(ilist,2+2*ich)
                a0 = dble(spl%Jmtx_list(ilist,2+2*ich+1))
                kij = var%ki(ijch,iz)
                tmp = a0 * kij

                bij = 0
                do isp = 2, spl%reactant_list(ijch,1)+1
                  k = spl%reactant_list(ijch,isp)
                  if (jj /= k) then
                    tmp = tmp * var%ni(k,iz)
                  else if (jj == k) then
                    bij = bij + 1
                  end if
                end do
                nj = var%ni(jj,iz)
                tmp = tmp * dble(bij) * nj**dble(bij-1)

                var%Jmtx(i,j) = var%Jmtx(i,j) + tmp

              end do

              ! Lower boundary: density fix
              if (nint(var%LowerBC(ii,1)) == 1.0_dp .and. iz == 1) then 
                var%Jmtx(i,j) = 0.0_dp
              end if

            end do

            do isp = 1, spl%nsp_i
              do jsp = 1, spl%nsp_i
                i = (iz-1)*spl%nsp_i+isp
                j = (iz-1)*spl%nsp_i+jsp + spl%nsp_i + 1 - i
                var%tAmtx(j,i) = - var%Jmtx(isp,jsp)
              end do
            end do

            ! vector b = P - L - dPhi_dz 
            do isp = 1, spl%nsp_i
              ii = spl%var_to_all(isp)
              i = (iz-1)*spl%nsp_i+isp
              var%barr(i) = ( var%Pi(isp,iz) - var%Li(isp,iz) ) - var%dPhi_dz(isp,iz)

              ! Lower boundary: density fix
              if (nint(var%LowerBC(ii,1)) == 1.0_dp .and. iz == 1) then 
                var%barr(i) = 0.0_dp
              end if

            end do

          end do ! z

          !print *, 'maxJ = ',maxJ
          !stop

          ! matrix A = ( I/dt - λ*J )
          !   Combining Chemical Jacobian and Transport Jacobian term
          do iz = 1, grd%nz
            do isp = 1, spl%nsp_i
              ii = spl%var_to_all(isp)
              i = (iz-1)*spl%nsp_i+isp
              j = spl%nsp_i + 1
              var%tAmtx(j,i) = var%tAmtx(j,i) + var%d_dni0_dPhi_dz(isp,iz)
              if (nint(var%LowerBC(ii,1)) == 1.0_dp .and. iz == 1) then 
                var%tAmtx(j,i) = 0.0_dp
              end if
            end do
          end do
          do iz = 2, grd%nz
            do isp = 1, spl%nsp_i
              ii = spl%var_to_all(isp)
              i = (iz-1)*spl%nsp_i+isp
              j = (iz-2)*spl%nsp_i+isp + spl%nsp_i + 1 - i
              var%tAmtx(j,i) = var%d_dnil_dPhi_dz(isp,iz)
              i = (iz-2)*spl%nsp_i+isp
              j = (iz-1)*spl%nsp_i+isp + spl%nsp_i + 1 - i
              var%tAmtx(j,i) = var%d_dniu_dPhi_dz(isp,iz-1)
              if (nint(var%LowerBC(ii,1)) == 1.0_dp .and. iz == 2) then 
                var%tAmtx(j,i) = 0.0_dp
              end if
            end do
          end do

          do iz = 1, grd%nz
            do isp = 1, spl%nsp_i
              i = (iz-1)*spl%nsp_i+isp
              j = spl%nsp_i + 1
              var%tAmtx(j,i) = var%tAmtx(j,i) + 1.0_dp / var%dtime
            end do
          end do

          ! Solve 'A x = b' using LU decomposition
          CALL p__LU_CMP_CHEM_DIFF( spl, grd, var )
          CALL p__LU_SLV_CHEM_DIFF( spl, grd, var )
          !CALL p__LU_IMPROVE(spl, grd, var)

          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              i = (iz-1)*spl%nsp_i+spl%all_to_var(isp)
              var%ni_new(isp,iz) = 0.0_dp
              if (spl%label_fix(isp) == 0) then
                var%ni_new(isp,iz) = var%ni(isp,iz) + var%xarr(i)
              else if (spl%label_fix(isp) == 1) then
                var%ni_new(isp,iz) = var%ni(isp,iz)
              end if
            end do
          end do

        end if

        do iz = 1, grd%nz
          do isp = 1, spl%nsp
            if (var%ni_new(isp,iz) < 0.0_dp) then
              var%ni_new(isp,iz) = 0.0_dp
            end if
            if (var%ni_new(isp,iz) < 1.0e-20_dp) then
              var%ni_new(isp,iz) = 1.0e-20_dp
            end if
          end do
        end do
        

    ! electron density : charge neutrality ---------------------------------------------------------
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (spl%species(isp) == 'e-') then
            var%ni_new(isp,iz) = 0.0_dp
            do jsp = 1, spl%nsp
              if (var%q(jsp) > 0.0_dp .and. jsp /= isp) then
                var%ni_new(isp,iz) &
                  &  = var%ni_new(isp,iz) + var%ni_new(jsp,iz)
              end if
              if (var%q(jsp) < 0.0_dp .and. jsp /= isp) then
                var%ni_new(isp,iz) &
                  &  = var%ni_new(isp,iz) - var%ni_new(jsp,iz)
              end if
            end do
          end if
        end do
      end do
    

    ! advance time step ------------------------------------------------------------------
      ! lower boundary condition: density
      do isp = 1, spl%nsp
        if (nint(var%LowerBC(isp,1)) == 1.0_dp) then 
          var%ni_new(isp,1) = var%LowerBC(isp,2)
        end if
      end do

    ! convergence check and iuncrease time step
      if (set%calc_stable == 1 .and. set%start_rot == 0) then
        if ( var%istep >= 2 ) then

          var%max_dn_n = 0.0_dp
          do iz = 1, grd%nz
            do isp = 1, spl%nsp
              ! delta_N / N
              if ( var%ni(isp,iz) /= 0.0d0 ) then
                dni_ni(isp,iz) = dabs(var%ni_new(isp,iz) - var%ni(isp,iz)) / var%ni(isp,iz)
              else if ( var%ni(isp,iz) == 0.0d0 ) then
                dni_ni(isp,iz) = 0.0d0
              end if
              if ( dni_ni(isp,iz) > var%max_dn_n(3) ) then
                var%max_dn_n(1) = isp
                var%max_dn_n(2) = iz
                var%max_dn_n(3) = dni_ni(isp,iz)
                tmp = dble(isp)
              end if
            end do
          end do

          ! if max(dN/N) < 10 %, then goes to larger dt (*10)
          ! if dt becomes 100000 [sec], unstable at HC layer
          if ( set%mode == '1D' .or. set%mode == '2D Lat' ) then 
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
                do iz = 1, grd%nz
                  do isp = 1, spl%nsp
                    var%ni_0(isp,iz) = var%ni_new(isp,iz)
                  end do
                end do
              else if ( var%max_dn_n(3) >= eps ) then
                var%iter = var%iter + 1
                if ( var%iter > 5000 ) then
                  var%dtime = 1.0e-8_dp
                  var%iter = 0
                  do iz = 1, grd%nz
                    do isp = 1, spl%nsp
                      var%ni_0(isp,iz) = var%ni_new(isp,iz)
                    end do
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
      do iz = 1, grd%nz
        do isp = 1, spl%nsp
          if (var%ni_new(isp,iz) /= var%ni_new(isp,iz)) then
            var%ni(isp,iz) = var%ni_0(isp,iz)
            var%ni_new(isp,iz) = var%ni_0(isp,iz)
            var%dtime = 1.0e-8_dp
            var%iter = 0
            jsp = spl%all_to_var(isp)
            print *, iz, trim(spl%species(isp)), var%Pi(jsp,iz), var%Li(jsp,iz) 
            do ich = 1, spl%nch
              print *, ich, var%ki(ich,iz)
            end do

            do i = 1, grd%nz
              print *, i, var%tau_EUV(5,i)
            end do

            call p__io_stable__fin(set, spl, var, grd) ! in
            stop
          end if
          !var%ni_0(isp,iz) = var%ni(isp,iz)
          var%ni(isp,iz) = var%ni_new(isp,iz)
        end do
      end do


  end subroutine p__photochem_scheme__exe


  !--------------------------------------------------------------------
  ! LU decomposition : chemical & flux implicit
  !--------------------------------------------------------------------
  subroutine p__LU_CMP_CHEM_DIFF( spl, grd, var )
    implicit none
    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var

    integer i, j, k, is, ie, ks, ke, nsp1, nsp1_i, nsp1_j
    integer N, nsp, nz
    real(dp) rUmtx, sum0

    nsp = spl%nsp_i; nsp1 = nsp+1
    nz  = grd%nz
    N   = nsp * nz

    var%tLmtx = 0.0_dp
    var%Umtx  = 0.0_dp
    var%tLmtx(nsp1,:) = 1.0_dp

    do j = 1, N

      nsp1_j = nsp1-j

      is = j-nsp
      if (is<=1) var%Umtx(1+nsp1_j,j) = var%tAmtx(j+nsp,1)
      if (is<2) is = 2
      ie = j
      do i = is, ie
        nsp1_i = nsp1-i
        ks = j-nsp
        ke = i-1
        sum0 = DOT_PRODUCT(var%tLmtx(ks+nsp1_i:ke+nsp1_i,i),var%Umtx(ks+nsp1_j:ke+nsp1_j,j))
        var%Umtx(i+nsp1_j,j) = var%tAmtx(j+nsp1_i,i)-sum0
      end do

      is = j+1
      ie = j+nsp; if (ie>N) ie = N
      rUmtx =  1.0_dp / var%Umtx(nsp1,j)
      do i = is, ie
        nsp1_i = nsp1-i
        ks = i-nsp
        ke = j-1
        sum0 = DOT_PRODUCT(var%tLmtx(ks+nsp1_i:ke+nsp1_i,i),var%Umtx(ks+nsp1_j:ke+nsp1_j,j))
        var%tLmtx(j+nsp1_i,i) = (var%tAmtx(j+nsp1_i,i)-sum0)*rUmtx  
      end do

    end do

  end subroutine p__LU_CMP_CHEM_DIFF


  ! Solve simultaneous N eqations for N unknowns using LU decomposition
  subroutine p__LU_SLV_CHEM_DIFF( spl, grd, var )
    implicit none
    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var
    integer i, j, js, je
    integer N, nsp, nz, nsp1, nsp1_i
    real(dp) sum0, rUmtx, tmp

    nsp = spl%nsp_i; nsp1 = nsp+1
    nz  = grd%nz
    N   = nsp * nz

    var%yarr(1) = var%barr(1)
    do i = 2, N
      nsp1_i = nsp1-i
      js = i-nsp; if (js<1) js = 1
      je = i-1
      sum0 = DOT_PRODUCT(var%tLmtx(js+nsp1_i:je+nsp1_i,i),var%yarr(js:je))
      var%yarr(i) = var%barr(i)-sum0
    end do

    var%xarr(N) = var%yarr(N) / var%Umtx(nsp1,N)
    do i = N-1, 1, -1
      sum0 = 0.0_dp
      js = i+1
      je = i+nsp; if (je>N) je = N
      rUmtx = 1.0_dp / var%Umtx(nsp1,i)
      do j = js, je
        sum0 = sum0 + var%Umtx(i+nsp1-j,j) * var%xarr(j)
      end do
      var%xarr(i) = (var%yarr(i)-sum0) * rUmtx
    end do

  end subroutine p__LU_SLV_CHEM_DIFF


  ! Iteration using LU decomposition
  subroutine p__LU_IMPROVE(spl, grd, var)

    type(spl_), intent(in)    :: spl
    type(grd_), intent(in)    :: grd
    type(var_), intent(inout) :: var

    integer i, j, N, nsp, nz
    double precision sum

    nsp = spl%nsp_i
    nz  = grd%nz
    N   = nsp * nz

    do i = 1, N
      sum = 0.0_dp
      do j = i-nsp, i+nsp
        if (j >= 1 .and. j <= N) then 
          sum = sum + var%tAmtx(j+nsp+1-i,i) * var%xarr(j)
        end if 
      end do
      var%rarr(i) = sum - var%barr(i)
    end do

    !CALL p__LU_SLV_CHEM_DIFF( spl, grd, var )
    do i = 1, N
      var%xarr(i) = var%xarr(i) - var%dxarr(i)
    end do

  end subroutine p__LU_IMPROVE




end module p__photochem_scheme
