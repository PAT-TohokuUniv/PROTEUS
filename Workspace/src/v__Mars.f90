module v__Mars
  
  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__Mars__ini

contains


  subroutine v__Mars__ini(spl, cst, grd, flx, & ! in
    &                     var                 ) ! inout
    type(spl_),   intent(in)     :: spl
    type(cst_),   intent(in)     :: cst
    type(grd_),   intent(in)     :: grd
    type(flx_),   intent(in)     :: flx
    type(var_),   intent(inout)  :: var
    integer i, j, ix, iy, iz, isp, ich, jch, nspecial
    real(dp) tmp, tmparr(grd%nz)
    character(len=256) fname

    if ( spl%planet == 'Mars' ) then 
      nspecial = var%nspecial
      if (var%nspecial == 0) nspecial = 1
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))
      var%ki_special = 0.0_dp

      jch = 0
      do ich = 1, spl%nch
        if ( spl%reaction_type_char(ich) == 'proton impact' ) then

          ! HOx production
          if ( p__search_product(spl, ich, 'OH') == 1 ) then

            fname = './Mars/SEP_input/P_HOx.dat'
            open(11, file = fname, status = 'old' )
              do iz = 1, grd%nz
                read(11,*) tmp, tmparr(iz)
              end do
            close(11)

            jch = jch + 1
            var%ich_special(jch) = ich

            isp = sp_index(spl,'CO2')

            do iz = 1, grd%nz
            do iy = 1, grd%ny
            do ix = 1, grd%nx
              var%ki_special(jch,ix,iy,iz) = tmparr(iz) 
            end do
            end do
            end do

          end if

        end if

      end do 

    end if

  end subroutine v__Mars__ini


end module v__Mars
