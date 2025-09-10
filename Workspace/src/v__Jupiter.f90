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

      isp = sp_index(spl, 'H2')

      if (var%nspecial == 0) nspecial = 1
      if (var%nspecial /= 0) nspecial = var%nspecial
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))

      var%ki_special = 0.0_dp

      ! auroral electron precipitation
      if (var%nspecial >= 1) then
        do ich = 1, spl%nch
          jch = 1
          if ( spl%reaction_type_char(ich) == 'electron impact' ) then
            var%ich_special(1) = ich

            fname = './Jupiter/k_aurora.dat'
            open(11, file = fname, status = 'unknown' )
              do iz = 1, grd%nz
              do iy = 1, grd%ny
              do ix = 1, grd%nx
                read(11,*) tmp, tmp, tmp, tmp1
                var%ki_special(jch,ix,iy,iz) = tmp1 / var%ni(iz,isp)
              end do
              end do
              end do
            close(11)
          end if
          jch = jch + 1
        end do
      end if

    end if

  end subroutine v__Jupiter__ini


end module v__Jupiter
