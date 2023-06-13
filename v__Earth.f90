!------------------------------------------------------------------------------
!
!                                       Earth
!
!------------------------------------------------------------------------------


module v__Earth
  use v__tdec,     only : set_, cst_, spl_, var_, grd_, flx_
  use p__search,   only : p__search_reactant, p__search_product, sp_index

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: v__Earth__ini, v__Earth__exe

contains

  !----------------------------------------------------------------------------------------------------------
  !
  !                                            Calculation settings
  !
  !----------------------------------------------------------------------------------------------------------
  subroutine v__Earth__ini(spl, & ! in
    &                      set  ) ! inout
    type(spl_),   intent(in)      :: spl
    type(set_),   intent(inout)   :: set

    real(dp), parameter :: Eday = 86400.0_dp

    character(len=256) fname

    if ( spl%planet == 'Earth' ) then

      !----------------------------------------------------------------
      !                    SETTING CALCULATION MODE
      !----------------------------------------------------------------
      ! mode selection :
      !      1 : 1D stable solution only (average over LT, latitude is fixed)
      !      2 : 2D stable solution only (average over LT)
      !
      !      3 : 2D rotation (longitude, latitude are fixed)
      !      4 : 3D rotation (longitude is fixed)
      !      5 : 3D global calculation # NOT YET #
      !
      !         mode 1, 3    => Set latitude below !
      !         mode 3, 4, 5 => Select CALCULATE or SKIP stable calculation below !
      !----------------------------------------------------------------
      !        SELECT 'CALCULATE' OR 'SKIP' STABLE CALCULATION
      !----------------------------------------------------------------
      ! calc_stable :
      !      1 : calclate stable solution
      !      0 : skip stable calclation
      !----------------------------------------------------------------
      !               SETTING INPUT STABLE DENSITY PATH
      !----------------------------------------------------------------
      ! if 'calc_stable =  0' (= skip stable calclation)
      !  -> please fill in the path of input density data
      !  -> grids of density data should be the same as the current grids
      !----------------------------------------------------------------
      !                   SETTING PLANETARY DAYS
      !----------------------------------------------------------------
      ! nday  : number of planetary days for rotational calculation
      !         (calculation days = nday - 0.5  [days]
      !----------------------------------------------------------------
      !                   SETTING CALCULATION SCHEME
      !----------------------------------------------------------------
      ! scheme = 1 : 'implicit'      : calculate chemical-flux Jacobian
      !                                to use implicit method
      !          2 : 'semi-implicit' : dt must be small
      !          3 : 'explicit'      : dt must be very small
      !

      ! calculation mode
      set%calc_stable =  0
      set%read_stable =  1
      set%test_loc    =  0
      set%sza         =  0.0_dp

      !----------------------------------------------------------------
      !      SETTING HORIZONTAL RESOLUTION FOR 2D, 3D CALCULATION
      !----------------------------------------------------------------
      ! latitude, local time grid number settings
      !  
      !      resx : resolution of local time, longitude [deg]
      !      resy : resolution of latitude [deg]
      !
      !      nx  : numer of local time grids : adapted automatically
      !      ny  : number of latitude grids  : adapted automatically
      !
      !  latitude, local time is calculated in 'p__prm.f90'
      !      local time :  00 - 24 LT (1: 00LT, nx: 24LT)
      !      latitude   : -90 ~ 90    (1: -90,  ny: +90)
      !--------------------------------
      !  nx, ny : MUST BE ODD NUMBER
      !--------------------------------

      set%resx = 3 ! [deg] localtime resolution
      set%resy = 3 ! [deg] latitude resolution

      !                                           END Calculation settings
      !----------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------
      !                  RESOLUTION ERROR MESSAGE
      !----------------------------------------------------------------
      if ( mod(360,set%resx) /= 0 .or. mod(180,set%resy) /= 0) then
        fname = './progress.dat'
        open(11, file = fname, status = 'unknown' )
          write(11,*) 'change resolution!'
          write(11,*)'!------------------------------------------'
          write(11,*)'!  360/resx, 180/resy : MUST BE INTEGER'
          write(11,*)'!------------------------------------------'
        close(11)
        stop
      else if ( mod(360/set%resx+1,2)==0 .or. mod(360/set%resy+1,2)==0 ) then
        fname = './progress.dat'
        open(11, file = fname, status = 'unknown' )
          write(11,*) 'change resolution!', 360/set%resx+1, 360/set%resy+1
          write(11,*)'!------------------------------------------'
          write(11,*)'!  nx, ny : MUST BE ODD NUMBER'
          write(11,*)'!------------------------------------------'
        close(11)
        stop
      end if

      ! number of latitude and local time grids are automatically adapted with the selected mode
      if ( set%mode == '1D' ) then
        set%nx     =  1
        set%ny     =  1
      else if ( set%mode == '2D Lat' ) then
        set%nx     =  1
        set%ny     =  180/set%resy + 1
      else if( set%mode == '2D Rot' ) then
        set%nx     =  360/set%resx + 1
        set%ny     =  1
      else if( set%mode == '3D Rot' .or. set%mode == '3D Global' ) then
        set%nx     =  360/set%resx + 1
        set%ny     =  180/set%resy + 1
      end if

    end if


  end subroutine v__Earth__ini



  !----------------------------------------------------------------------------------------------------------
  !
  !                                           Special treatment for Mars
  !
  !----------------------------------------------------------------------------------------------------------
  subroutine v__Earth__exe(spl, cst, grd, flx, & ! in
    &                      var                 ) ! inout
    type(spl_),   intent(in)     :: spl
    type(cst_),   intent(in)     :: cst
    type(grd_),   intent(in)     :: grd
    type(flx_),   intent(in)     :: flx
    type(var_),   intent(inout)  :: var
    integer i, j, ix, iy, iz, isp, ich, nspecial
    real(dp) tmp, tmp1, tmp2, tmpzarr1(grd%nz), tmpzarr2(grd%nz)
    character(len=256) fname

    if ( spl%planet == 'Earth' ) then 
      nspecial = var%nspecial
      if (var%nspecial == 0) nspecial = 1
      allocate(var%ich_special(nspecial), var%ki_special(nspecial,grd%nx,grd%ny,grd%nz))
      var%ki_special = 0.0_dp

    end if

  end subroutine v__Earth__exe


end module v__Earth
