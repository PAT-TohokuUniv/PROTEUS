module v__Earth

  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__Earth__ini

contains

  
  subroutine v__Earth__ini(spl, cst, grd, flx, & ! in
    &                      var                 ) ! inout
    type(spl_),   intent(in)     :: spl
    type(cst_),   intent(in)     :: cst
    type(grd_),   intent(in)     :: grd
    type(flx_),   intent(in)     :: flx
    type(var_),   intent(inout)  :: var
    integer nspecial

    if ( spl%planet == 'Earth' ) then 
      nspecial = var%nspecial
      if (var%nspecial == 0) nspecial = 1
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))
      var%ki_special = 0.0_dp

    end if

  end subroutine v__Earth__ini


end module v__Earth
