module p__photochem_rate
  
  use v__tdec,                 only : spl_, var_, grd_, cst_, xct_, flx_, set_
  use p__PROTEUS,              only : p__PROTEUS_source
  use p__photolysis_rate,      only : photolysis_rate
  use p__photoionization_rate, only : photoionization_rate
  use p__EUV,                  only : p__EUVAC_photoionization_rate
  use p__search,               only : sp_index

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private

  ! Module-level variables --------------------------------------------------------------
  real(dp), allocatable, private :: Tn(:,:,:), Ti(:,:,:,:), Te(:,:,:)
  real(dp), allocatable, private :: n_var(:,:,:,:), n_fix(:,:,:,:)
  real(dp), allocatable, private :: v_var(:,:,:,:,:), v_fix(:,:,:,:,:)
  real(dp), allocatable, private :: prod(:,:,:,:), loss(:,:,:,:), kout(:,:,:,:)
  real(dp), allocatable, private :: J_rate_in(:,:,:,:)
  character(len=256), allocatable, private :: species_var(:), species_fix(:)
  
  ! Public interface --------------------------------------------------------------------
  public :: p__photochem_rate__ini, p__photochem_rate__exe

contains


  subroutine p__photochem_rate__ini(spl, grd) ! inout
    implicit none
    type(spl_), intent(in) :: spl
    type(grd_), intent(in) :: grd
    integer nsp_i, nsp_f, nch, nz
    integer nsp

    nsp = spl%nsp
    nsp_i = spl%nsp_i
    nsp_f = nsp - nsp_i
    nch   = spl%nch
    nz    = grd%nz

    allocate(Tn(1,1,nz), Ti(1,1,nz,nsp_i), Te(1,1,nz))
    allocate(n_var(1,1,nz,nsp_i), n_fix(1,1,nz,nsp_f))
    allocate(v_var(3,1,1,nz,nsp_i), v_fix(3,1,1,nz,nsp_f))
    allocate(prod(1,1,nz,nsp_i), loss(1,1,nz,nsp_i), kout(1,1,nz,nch))
    allocate(J_rate_in(1,1,nz,nch))
    allocate(species_var(nsp_i), species_fix(nsp_f))

  end subroutine p__photochem_rate__ini


  subroutine p__photochem_rate__exe(spl, cst, grd, flx, set, & ! in
    &                               xct, var                 ) ! inout
    implicit none
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(flx_),           intent(in)    :: flx
    type(set_),           intent(in)    :: set
    type(xct_),           intent(inout) :: xct
    type(var_),           intent(inout) :: var

    integer isp, jsp, ksp, ich, jch, iz
    integer i
    real(dp) tmp, tmp1, tmp2

    !-------------------------------------------------------------------------------------------------
    !
    !                                    Photochemical calculation
    !
    !-------------------------------------------------------------------------------------------------

    ! column density -------------------------------------------------------------------------------
      var%clm_ni = 0.0_dp
      do isp = 1, spl%nsp
        var%clm_ni(grd%nz,isp) = var%ni(grd%nz,isp)*grd%dalt(grd%nz)
        do iz = grd%nz, 2, -1
          var%clm_ni(iz-1,isp) = var%clm_ni(iz,isp) &
            &                   + (var%ni(iz-1,isp) + var%ni(iz,isp)) &
            &                   / 2.0_dp * grd%dalt(iz-1)
        end do
      end do

    ! 'M' density : sum of all ---------------------------------------------------------
      isp = sp_index(spl, 'M')
      if (isp >= 1 .and. isp <= spl%nsp) then 
        var%ni(1:grd%nz,isp) = 0.0_dp
        do jsp = 1, spl%nsp
          if (isp /= jsp) then
            var%ni(1:grd%nz,isp) = var%ni(1:grd%nz,isp) + var%ni(1:grd%nz,jsp)
          end if
        end do
      end if

    ! total density, mean mass ---------------------------------------------------------
      do iz = 1, grd%nz
        tmp1 = 0.0_dp
        tmp2 = 0.0_dp
        do isp = 1, spl%nsp
          if (spl%species(isp) /= 'M' .and. spl%species(isp) /= 'products' .and. spl%species(isp) /= 'hv') then
            tmp1 = tmp1 + var%ni(iz,isp) * var%m(isp)
            tmp2 = tmp2 + var%ni(iz,isp)
          end if
        end do
        var%m_mean(iz) = tmp1 / tmp2
        var%n_tot(iz) = tmp2
      end do

    ! reaction rate -------------------------------------------------------------------------------
      var%ki = 0.0_dp

      Tn(1,1,:) = var%Tn(:)
      do isp = 1, spl%nsp_i
        Ti(1,1,:,isp) = var%Ti(:)
      end do 
      Te(1,1,:) = var%Te(:)

      i = 0
      species_var = ''
      species_fix = ''
      n_var = 0.0_dp
      n_fix = 0.0_dp
      do isp = 1, spl%nsp
        jsp = spl%all_to_var(isp)
        if (jsp /= 0) then 
          n_var(1,1,1:grd%nz,jsp) = var%ni(1:grd%nz,isp) * 1.0e-6_dp
          species_var(jsp) = spl%species(isp)
        end if
        
        if (jsp == 0) then 
          i = i + 1
          n_fix(1,1,1:grd%nz,i) = var%ni(1:grd%nz,isp) * 1.0e-6_dp
          species_fix(i) = spl%species(isp)
        end if
      end do 
      v_var = 0.0_dp
      v_fix = 0.0_dp

      !
      ! photolysis and photoionization rate
      !

      J_rate_in = 0.0_dp
      call photolysis_rate(grd%alt(1:), grd%dalt(1:), var%ni(1:,1:), & ! in
        &                  var%Tn(1:), grd%sza(grd%ix,grd%iy),       & ! in
        &                  cst%R, cst%Mplanet,                       & ! in
        &                  J_rate_in(1,1,1:grd%nz,1:spl%nch)         ) ! out
      if (set%euv_input /= 'EUVAC') then 
        call photoionization_rate(grd%alt(1:), grd%dalt(1:), var%ni(1:,1:), & ! in
          &                       var%Tn(1:), grd%sza(grd%ix,grd%iy),       & ! in
          &                       cst%R, cst%Mplanet,                       & ! in
          &                       J_rate_in(1,1,1:grd%nz,1:spl%nch)         ) ! out
      else if (set%euv_input == 'EUVAC') then 
        call p__EUVAC_photoionization_rate(cst, grd, xct, flx, spl,               & ! in
          &                                var, J_rate_in(1,1,1:grd%nz,1:spl%nch) ) ! inout
      end if
      J_rate_in = J_rate_in * flx%irradiance_factor

      !
      ! Special reaction rate case
      !
      ! if you input production or loss rate explicitly, please input as var%ki_special in p__planet__exe in p__planet.f90
      do jch = 1, var%nspecial
        ich = nint(var%ich_special(jch))
        do iz = 1, grd%nz
          var%ni(iz,0) = 1.0_dp ! if there are no reactant like meteoroid ablation, P = k * ni(0,iz) = k
          if (    spl%reaction_type_char(ich) == 'electron impact' &
          &  .or. spl%reaction_type_char(ich) == 'proton impact' &
          &  .or. spl%reaction_type_char(ich) == 'H impact'  &
          &  .or. spl%reaction_type_char(ich) == 'Meteoroid ablation' &
          &  .or. spl%reaction_type_char(ich) == 'Rainout' ) then
            if (spl%reactant_list(ich,0) == 1 .and. spl%reactant_list(ich,1) == 0) then
              J_rate_in(1,1,iz,ich) = var%ki_special(jch,grd%ix,grd%iy,iz) * 1.0e-6_dp
            end if
            if (spl%reactant_list(ich,0) == 1 .and. spl%reactant_list(ich,1) /= 0) then
              J_rate_in(1,1,iz,ich) = var%ki_special(jch,grd%ix,grd%iy,iz) 
            end if
            if (spl%reactant_list(ich,0) == 2) then
              J_rate_in(1,1,iz,ich) = var%ki_special(jch,grd%ix,grd%iy,iz) * 1.0e6_dp
            end if
            if (spl%reactant_list(ich,0) == 3) then
              J_rate_in(1,1,iz,ich) = var%ki_special(jch,grd%ix,grd%iy,iz) * 1.0e12_dp
            end if
          end if
        end do
      end do

      do ich = 1, spl%nch
        if (spl%reaction_type_char(ich) == 'datafile') then ! rate from external datafile
          J_rate_in(1,1,:,ich) = set%rate_from_datafile(:,ich) ! cgs unit
        end if
      end do

      !
      ! Production and loss rates
      !
      
      prod = 0.0_dp
      loss = 0.0_dp
      kout = 0.0_dp

      call p__PROTEUS_source(species_var, species_fix,      & ! in:  name of variable species and fixed species
        &                    spl%nsp_i, spl%nsp-spl%nsp_i,  & ! in:  number of variable species and fixed species
        &                    1, 1, grd%nz,                  & ! in:  number of grids in x, y, z direction
        &                    Tn, Ti, Te,                    & ! in:  temperature of neutrals, ions and electrons
        &                    v_var, v_fix,                  & ! in:  three dimensional velocity vector of variable and fixed species
        &                    n_var, n_fix,                  & ! in:  number density of variable and fixed species
        &                    J_rate_in,                     & ! in:  photoionization / dissociation reaction list
        &                    prod, loss, kout               ) ! out: production and loss rates for variable species

      do ich = 1, spl%nch
        var%ki(1:grd%nz,ich) = kout(1,1,1:grd%nz,ich)
      end do 
      
      do isp = 1, spl%nsp_i
        var%Pi(1:grd%nz,isp) = prod(1,1,1:grd%nz,isp) * 1.0e6_dp
        var%Li(1:grd%nz,isp) = loss(1,1,1:grd%nz,isp) * 1.0e6_dp
      end do 

      do iz = 1, grd%nz
        do ich = 1, spl%nch
          tmp = var%ki(iz,ich)
          do jsp = 1, spl%reactant_list(ich,0)
            ksp = spl%reactant_list(ich,jsp)
            if (ksp > 0) then
              tmp = tmp * var%ni(iz,ksp) * 1.0e-6_dp
            end if
          end do
          var%rate(iz,ich) = tmp * 1.0e6_dp
        end do
      end do ! z


  end subroutine p__photochem_rate__exe



end module p__photochem_rate
