module v__pref

  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__pref__ini

contains


  subroutine v__pref__ini(spl, cst, grd, flx, set, & ! in
    &                     var                      ) ! inout
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

    allocate(var%ki_special(grd%nx,grd%ny,grd%nz,spl%nch))

    select case (spl%planet)

    case ('Earth')

      var%ki_special = 0.0_dp

      do ich = 1, spl%nch
        ! Users can define special reaction treatment here
        ! If users specified "FILE=path_to_datafile" in the reaction_list.txt, 
        ! the datafile is loaded and registered in v__in.f90, so users do not need to load it again.
        if (spl%reaction_type_list(ich) == 1000) then
          ! Load datafile for user-defined special reaction
          do iz = 1, grd%nz
            var%ki_special(:,:,iz,ich) = 0.0_dp ! This is dummy example.
          end do
        end if
      end do
    
    case ('Venus')

      var%ki_special = 0.0_dp

    case ('Mars')

      var%ki_special = 0.0_dp

    case ('Jupiter')

      var%ki_special = 0.0_dp

    end select

  end subroutine v__pref__ini


end module v__pref
