program sigma_O3
  implicit none
  integer i, j
  real(8) sigma(88738,11), lambda(88738)
  character(len=256) fname

  fname = 'sigma_O3_tmp.dat'
  open(11, file = fname, status = 'unknown')
    do i = 1, 88738
      read(11,*) lambda(i), (sigma(i,j), j = 1, 11)
    end do 
  close(11)

  fname = 'sigma_a_O3_new.dat'
  open(11, file = fname, status = 'unknown')
    do i = 1, 88738
      write(11,'(f8.3, 11e14.6)') lambda(i), (sigma(i,12-j), j = 1, 11)
    end do 
  close(11)

end program sigma_O3 