 module p__UV
  
  use v__tdec,            only : var_, grd_, cst_, xct_, spl_, flx_, set_
  use c__prm,             only : c__prm__ini
  use p__photolysis_rate, only : p__photolysis_rate__ini, load_cross_section_dat, get_cross_section
  use p__PROTEUS,         only : get_solar_flux, get_number_of_wavelength_bins, get_wavelength_bin
  use p__search,          only : p__search_reactant, p__search_product, ch_identify, sp_index


  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__UV_flux, p__UV_cross_section

contains


  !-------------------------------------------------
  ! solar UV model: Woods et al., 2016
  !-------------------------------------------------
  subroutine p__UV_flux(spl, grd, set, & ! in
    &                   var, xct, flx  ) ! inout
    implicit none
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(set_),           intent(in)    :: set
    type(var_),           intent(inout) :: var
    type(xct_),           intent(inout) :: xct
    type(flx_),           intent(inout) :: flx
    integer i
    character(len=256) fname

    !---------------------------------------------------
    ! wavelength bin definition
    !---------------------------------------------------    
    flx%nwl_UV = get_number_of_wavelength_bins()
    allocate(flx%solar_UV(flx%nwl_UV), flx%lambda_UV(flx%nwl_UV), flx%dlambda_UV(flx%nwl_UV))
    call get_wavelength_bin(flx%lambda_UV, flx%dlambda_UV)

    !---------------------------------------------------
    ! solar flux 
    !---------------------------------------------------
    !units: photons m^-2 s^-1 nm^-1
    fname = set%solar_flux
    call get_solar_flux(flx%lambda_UV, flx%dlambda_UV, flx%solar_UV)

    flx%solar_UV = flx%solar_UV / (flx%orbit * flx%orbit) 
    
  end subroutine p__UV_flux


  subroutine p__UV_cross_section(spl, grd, flx, var) ! inout
    type(spl_),   intent(in)     :: spl
    type(grd_),   intent(inout)  :: grd
    type(flx_),   intent(inout)  :: flx
    type(var_),   intent(inout)  :: var

    call p__photolysis_rate__ini(flx%nwl_UV, grd%nz, spl%nsp, spl%nch, &
      &                          spl%species(1:), var%m(1:), spl%reaction_type_list, &
      &                          spl%reactant_list, spl%product_list, &
      &                          flx%lambda_UV(1:), flx%dlambda_UV(1:), flx%solar_UV(1:)) ! inout
    call load_cross_section_dat('./UV/xsect_data') ! in
    call get_cross_section(var%Tn(1:), 'absorption', 1) ! out
    call get_cross_section(var%Tn(1:), 'photolysis', 1) ! out

  end subroutine p__UV_cross_section


end module p__UV