module p__search
  
  use v__tdec, only : cst_, spl_, var_, grd_, flx_

  implicit none
  integer(4), parameter         :: sp = 4, dp = 8

  private
  public :: p__search_reactant, p__search_product, ch_identify, sp_index

contains

  function p__search_reactant(spl, ich, species)
    type(spl_),       intent(in) :: spl
    integer,          intent(in) :: ich
    character(len=*), intent(in) :: species
    integer p__search_reactant
    integer i

    p__search_reactant = 0
    do i = 1, spl%reactant_list(ich,0)
      if( spl%species(spl%reactant_list(ich,i)) == species ) then
        p__search_reactant = 1
      end if
    end do

  end function p__search_reactant

  function p__search_product(spl, ich, species)
    type(spl_),       intent(in) :: spl
    integer,          intent(in) :: ich
    character(len=*), intent(in) :: species
    integer p__search_product
    integer i

    p__search_product = 0
    do i = 1, spl%product_list(ich,0)
      if( spl%species(spl%product_list(ich,i)) == species ) then
        p__search_product = 1
      end if
    end do

  end function p__search_product

  function ch_identify(spl, ich, type, in_reactant, in_product)
    ! this function should not be used for searching reactions other than photoionization and photodissociation.
    implicit none
    type(spl_),       intent(in) :: spl
    character(len=*), intent(in) :: type, in_reactant(:), in_product(:)
    integer,          intent(in) :: ich
    integer ch_identify
    integer i
    character(len=256), allocatable :: reactants(:),  products(:)

    ch_identify = 0

    if (spl%reaction_type_char(ich) == type) then

      allocate(reactants(spl%reactant_list(ich,0)),  products(spl%product_list(ich,0)))

      do i = 1, spl%reactant_list(ich,0)
        reactants(i) = trim(spl%species(spl%reactant_list(ich,i)))
      end do

      do i = 1, spl%product_list(ich,0)
        products(i) = trim(spl%species(spl%product_list(ich,i)))
      end do

      ! if the number of products is 1: 1 case
      if (spl%product_list(ich,0) == 1 .and. size(in_product) == 1) then 
        if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1)))) then 
          ch_identify = 1
        end if
      ! if the number of products is 2: 2 cases
      else if (spl%product_list(ich,0) == 2 .and. size(in_product) == 2) then 
        if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1)))) then 
          ch_identify = 1
        end if
      ! if the number of products is 3: 6 cases
      else if (spl%product_list(ich,0) == 3 .and. size(in_product) == 3) then 
        if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(3)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(1))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(2))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(2))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(3))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(1)))) then 
          ch_identify = 1
        else if (    trim(ADJUSTL(reactants(1)))==trim(ADJUSTL(in_reactant(1))) &
        & .and. trim(ADJUSTL(products(1)))==trim(ADJUSTL(in_product(3))) &
        & .and. trim(ADJUSTL(products(2)))==trim(ADJUSTL(in_product(1))) &
        & .and. trim(ADJUSTL(products(3)))==trim(ADJUSTL(in_product(2)))) then 
          ch_identify = 1
        end if

      end if

    end if

  end function ch_identify
  

  function sp_index(spl, species)
    type(spl_),       intent(in) :: spl
    character(len=*), intent(in) :: species
    integer sp_index
    integer isp

    sp_index = 0
    loop: do isp = 1, spl%nsp
      if( spl%species(isp) == species ) then
        sp_index = isp
        exit loop
      end if
    end do loop

  end function sp_index


end module p__search

! count number of lines of headers
subroutine p__count_header(nh, fname)
  implicit none
  integer nh
  character(len=256) fname
  character(len=256) strm
  open(11, file = fname, status = 'unknown' )
    nh = 0
    do
      read(11,'(A)') strm
      if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
        nh = nh + 1
      else
        exit
      end if
    end do
  close(11)
end subroutine p__count_header

! count number of lines without headers
subroutine p__count_lines(nl, fname)
  implicit none
  integer nl, ios
  character(len=256) fname
  character(len=256) strm
  open(11, file = fname, status = 'unknown' )
    nl = 0
    do
      read(11,'(A)',iostat = ios) strm
      if (ios < 0) exit
      if (strm(1:1) == "#" .or. strm(1:1) == "!" .or. strm(1:1) == ";") then
        cycle
      else
        nl = nl + 1
      endif
    end do
  close(11)
end subroutine p__count_lines

! solving max value
subroutine p__max(arrin, N, imax, max)
  implicit none
  integer N, i, imax
  real(8) arrin(N), tmp, max

  tmp = arrin(1)
  imax = 1
  do i = 2, N
    if (tmp < arrin(i)) then
      tmp = arrin(i)
      imax = i
    end if
  end do
  max = tmp

end subroutine p__max
