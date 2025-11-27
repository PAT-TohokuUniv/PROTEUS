module v__Jupiter

  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__Jupiter__ini

contains


  subroutine v__Jupiter__ini(spl, cst, grd, flx, set, & ! in
    &                        var                      ) ! inout
    implicit none
    type(spl_),   intent(in)     :: spl
    type(cst_),   intent(in)     :: cst
    type(grd_),   intent(in)     :: grd
    type(flx_),   intent(in)     :: flx
    type(set_),   intent(in)     :: set
    type(var_),   intent(inout)  :: var
    integer i, j, ix, iy, iz, isp, ich, jch, nspecial
    real(dp) tmp, tmp1, tmp2, tmpzarr(grd%nz)
    character(len=256) fname

    if ( spl%planet == 'Jupiter' ) then

      if (var%nspecial == 0) nspecial = 1
      if (var%nspecial /= 0) nspecial = var%nspecial
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))

      var%ki_special = 0.0_dp

    end if

  end subroutine v__Jupiter__ini


end module v__Jupiter
