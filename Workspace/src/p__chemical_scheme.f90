module p__chemical_scheme
  
  use p__PROTEUS,       only : p__PROTEUS_Jacobian

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__chemical_scheme__exe, p__LU_decomposition, p__LU_solver

contains

  subroutine p__chemical_scheme__exe(scheme,                                         & ! in:    'explicit' or 'implicit'
    &                                var_species_list, fix_species_list,             & ! in:    name of variable species and fixed species
    &                                nsp_var, nsp_fix, nch, nz,                      & ! in:    number of variable species, fixed species, chemical reactions, and vertical grids
    &                                n_var, n_fix,                                   & ! in:    number density of variable and fixed species
    &                                n_lower_boundary,                               & ! in:    fixed density at lower boundary
    &                                prod, loss, k_coef, dphi_dz,                    & ! in:    production and loss rate of variable species, reaction rate coefficients, and vertical gradient of vertical flux
    &                                d_dniu_dphi_dz, d_dni0_dphi_dz, d_dnil_dphi_dz, & ! in:    Jacobian terms of vertical diffusion
    &                                tAmtx, Umtx, tLmtx,                             & ! inout: transposed chemical Jacobian matrix, LU decomposition matrices
    &                                dt,                                             & ! inout: timestep
    &                                delta_n_var                                     ) ! out:   change of number density of variable species
    ! < Input > ----------------------------------------------------------------------------------------------------------
    ! - scheme
    !   * String of 'explicit' or 'implicit' for chemical scheme.
    !
    ! - species(nsp)
    !   * Strings of chemical species and have the number of elements   
    !     equal to the number of chemical species to be treated as variables. 
    !   * Please note that all arrays to be given in this subroutine for variable and fixed species 
    !     should be aligned with the order of species_list string arrays.
    !
    ! - nsp_var, nsp_fix
    !   * The number of chemical species to be treated as variables and fixed species. 
    !
    ! - nch
    !   * The number of chemical reactions.
    !
    ! - nz
    !   * The number of grids in vertical direction. 
    !
    ! - n_var(nz,nsp_var), n_fix(nz,nsp_fix)
    !   * Number density of variable and fixed species in [cm^-3] in the simulation vertical grids.
    !
    ! - n_lower_boundary(nsp,2)
    !   * Lower boundary condition of number density of variable species in [cm^-3].
    !     Column 1: label 0.0_dp or 1.0_dp for fixed or free boundary condition for number density.
    !     Column 2: number density of variable species in [cm^-3] fixed at the lower boundary.
    !
    ! - prod(nz,nsp_var), loss(nz,nsp_var)
    !   * Production and loss rate of variable species in [cm^-3 s^-1] in the simulation vertical grids.
    !
    ! - k_coef(nz,nch)
    !   * Reaction rate coefficients in [cm^3 s^-1] in the simulation vertical grids.
    !
    ! - dphi_dz(nz,nsp), d_dniu_dphi_dz(nz,nsp), d_dni0_dphi_dz(nz,nsp), d_dnil_dphi_dz(nz,nsp)
    !   * Vertical gradient of vertical flux in [cm^-3 s^-1] and Jacobian terms of vertical diffusion
    !     in [s^-1] in the simulation vertical grids.
    !
    ! - tAmtx(2*nsp_var+1,nsp_var*nz), Umtx(nsp_var+1,nsp_var*nz), tLmtx(nsp_var+1,nsp
    !   * Jacobian matrix and LU decomposition matrices.
    !
    ! - dt
    !   * timestep in [s]
    !
    ! < Output > ---------------------------------------------------------------------------------------------------------
    ! - delta_n_var(nz,nsp_var)
    !   * Change of number density of variable species in [cm^-3] in the simulation vertical grids.
    !     
    !---------------------------------------------------------------------------------------------------------------------
    implicit none
    integer,          intent(in)    :: nsp_var, nsp_fix, nch, nz 
    character(len=*), intent(in)    :: scheme, var_species_list(nsp_var), fix_species_list(nsp_fix)
    real(dp),         intent(in)    :: n_var(nz,nsp_var), n_fix(nz,nsp_fix)
    real(dp),         intent(in)    :: n_lower_boundary(nsp_var,2)
    real(dp),         intent(in)    :: prod(nz,nsp_var), loss(nz,nsp_var), k_coef(nz,nch)
    real(dp),         intent(in)    :: dphi_dz(nz,nsp_var)
    real(dp),         intent(in)    :: d_dniu_dphi_dz(nz,nsp_var), d_dni0_dphi_dz(nz,nsp_var), d_dnil_dphi_dz(nz,nsp_var)
    real(dp),         intent(inout) :: tAmtx(2*nsp_var+1,nsp_var*nz), Umtx(nsp_var+1,nsp_var*nz), tLmtx(nsp_var+1,nsp_var*nz), dt
    real(dp),         intent(out)   :: delta_n_var(nz,nsp_var)

    integer isp, iz, i, j
    real(dp) Fvec(nsp_var*nz), xvec(nsp_var*nz)

    !
    ! Explicit scheme
    !
    if ( trim(ADJUSTL(scheme)) == 'explicit' ) then 

      do isp = 1, nsp_var
        delta_n_var(1:nz,isp) = ( prod(1:nz,isp) - loss(1:nz,isp)  - dphi_dz(1:nz,isp) ) * dt
      end do 

    end if

    !
    ! Semi-implicit scheme
    !
    if ( trim(ADJUSTL(scheme)) == 'semi-implicit' ) then 

      do isp = 1, nsp_var
        where (n_var(1:nz,isp) > 0.0_dp) 
          delta_n_var(1:nz,isp) = ( n_var(1:nz,isp) + (prod(1:nz,isp) - dPhi_dz(1:nz,isp)) * dt ) &
            &                 / (1.0_dp + (loss(1:nz,isp) / n_var(1:nz,isp)) * dt ) &
            &                 - n_var(1:nz,isp)
        end where
        where (n_var(1:nz,isp) == 0.0_dp) 
          delta_n_var(1:nz,isp) = ( n_var(1:nz,isp) + (prod(1:nz,isp) - dPhi_dz(1:nz,isp)) * dt ) &
            &                 - n_var(1:nz,isp)
        end where
      end do 

    end if

    !
    ! Implicit scheme
    !
    if ( trim(ADJUSTL(scheme)) == 'implicit' ) then 

      tAmtx = 0.0_dp
      ! Calculating chemical Jacobian matrix
      call p__PROTEUS_Jacobian(var_species_list, fix_species_list, & ! in:  name of variable species and fixed species
        &                      nsp_var, nsp_fix,                   & ! in:  number of variable species and fixed species
        &                      nz,                                 & ! in:  number of vertical grids
        &                      n_var, n_fix,                       & ! in:  number density of variable and fixed species
        &                      k_coef,                             & ! in:  reaction rate coefficients
        &                      tAmtx                               ) ! out: transposed chemical Jacobian matrix

      ! matrix A = ( I/dt - J )
      !   Combining chemical and transport terms of Jacobian matrix
      do isp = 1, nsp_var
        do iz = 1, nz
          i = (iz-1)*nsp_var+isp
          j = nsp_var + 1

          Fvec(i) = ( prod(iz,isp) - loss(iz,isp) ) - dphi_dz(iz,isp)
          if (nint(n_lower_boundary(isp,1)) == 1.0_dp .and. iz == 1) then 
            Fvec(i) = 0.0_dp
          end if

          tAmtx(j,i) = tAmtx(j,i) + d_dni0_dphi_dz(iz,isp)
          if (nint(n_lower_boundary(isp,1)) == 1.0_dp .and. iz == 1) then 
            tAmtx(j,i) = 0.0_dp
          end if
        end do
      end do

      do isp = 1, nsp_var
        do iz = 2, nz
          i = (iz-1)*nsp_var+isp
          j = (iz-2)*nsp_var+isp + nsp_var + 1 - i
          tAmtx(j,i) = d_dnil_dphi_dz(iz,isp)
          i = (iz-2)*nsp_var+isp
          j = (iz-1)*nsp_var+isp + nsp_var + 1 - i
          tAmtx(j,i) = d_dniu_dphi_dz(iz-1,isp)
          if (nint(n_lower_boundary(isp,1)) == 1.0_dp .and. iz == 2) then 
            tAmtx(j,i) = 0.0_dp
          end if
        end do
      end do

      do isp = 1, nsp_var
        do iz = 1, nz
          i = (iz-1)*nsp_var+isp
          j = nsp_var + 1
          tAmtx(j,i) = tAmtx(j,i) + 1.0_dp / dt
        end do
      end do

      ! Solve 'A x = b' using LU decomposition
      call p__LU_decomposition(nsp_var, nz,       & ! in:
        &                      tAmtx, Umtx, tLmtx ) ! inout:
        
      call p__LU_solver(nsp_var, nz,       & ! in:
        &               Umtx, tLmtx, Fvec, & ! in:
        &               xvec               ) ! out:
    
      do isp = 1, nsp_var
        do iz = 1, nz
          i = (iz-1)*nsp_var+isp
          delta_n_var(iz,isp) = xvec(i)
        end do
      end do

    end if


  end subroutine p__chemical_scheme__exe




  subroutine p__LU_decomposition(nsp, nz,           & ! in:
    &                            tAmtx, Umtx, tLmtx ) ! inout:
    implicit none
    integer,          intent(in)    :: nsp, nz 
    real(dp),         intent(inout) :: tAmtx(2*nsp+1,nsp*nz), Umtx(nsp+1,nsp*nz), tLmtx(nsp+1,nsp*nz)
  
    integer i, j, is, ie, ks, ke, nsp1, nsp1_i, nsp1_j
    integer N
    real(dp) rUmtx, sum0
  
    nsp1 = nsp+1
    N   = nsp * nz
  
    tLmtx = 0.0_dp
    Umtx  = 0.0_dp
    tLmtx(nsp1,:) = 1.0_dp
  
    do j = 1, N
  
      nsp1_j = nsp1-j
  
      is = j-nsp
      if (is<=1) Umtx(1+nsp1_j,j) = tAmtx(j+nsp,1)
      if (is<2) is = 2
      ie = j
      do i = is, ie
        nsp1_i = nsp1-i
        ks = j-nsp
        ke = i-1
        sum0 = DOT_PRODUCT(tLmtx(ks+nsp1_i:ke+nsp1_i,i),Umtx(ks+nsp1_j:ke+nsp1_j,j))
        Umtx(i+nsp1_j,j) = tAmtx(j+nsp1_i,i)-sum0
      end do
  
      is = j+1
      ie = j+nsp; if (ie>N) ie = N
      rUmtx =  1.0_dp / Umtx(nsp1,j)
      do i = is, ie
        nsp1_i = nsp1-i
        ks = i-nsp
        ke = j-1
        sum0 = DOT_PRODUCT(tLmtx(ks+nsp1_i:ke+nsp1_i,i),Umtx(ks+nsp1_j:ke+nsp1_j,j))
        tLmtx(j+nsp1_i,i) = (tAmtx(j+nsp1_i,i)-sum0)*rUmtx  
      end do
  
    end do
  
  end subroutine p__LU_decomposition




  subroutine p__LU_solver(nsp, nz,           & ! in:
    &                     Umtx, tLmtx, bvec, & ! in:
    &                     xvec               ) ! out:
    implicit none
    integer,          intent(in)    :: nsp, nz 
    real(dp),         intent(in)    :: Umtx(nsp+1,nsp*nz), tLmtx(nsp+1,nsp*nz), bvec(nsp*nz)
    real(dp),         intent(out)   :: xvec(nsp*nz)
  
    integer i, j, js, je
    integer N, nsp1, nsp1_i
    real(dp) sum0, rUmtx, yvec(nsp*nz)

    xvec = 0.0_dp
    yvec = 0.0_dp
  
    nsp1 = nsp+1
    N   = nsp * nz
  
    yvec(1) = bvec(1)
    do i = 2, N
      nsp1_i = nsp1-i
      js = i-nsp; if (js<1) js = 1
      je = i-1
      sum0 = DOT_PRODUCT(tLmtx(js+nsp1_i:je+nsp1_i,i),yvec(js:je))
      yvec(i) = bvec(i)-sum0
    end do
  
    xvec(N) = yvec(N) / Umtx(nsp1,N)
    do i = N-1, 1, -1
      sum0 = 0.0_dp
      js = i+1
      je = i+nsp; if (je>N) je = N
      rUmtx = 1.0_dp / Umtx(nsp1,i)
      do j = js, je
        sum0 = sum0 + Umtx(i+nsp1-j,j) * xvec(j)
      end do
      xvec(i) = (yvec(i)-sum0) * rUmtx
    end do
  
  end subroutine p__LU_solver



end module p__chemical_scheme


