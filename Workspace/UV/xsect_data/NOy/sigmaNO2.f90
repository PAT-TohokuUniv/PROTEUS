program sigma_NO2
  implicit none
  integer i, j
  real(8) sigma(27993,2), lambda(27993)
  character(len=256) fname

  fname = 'sigma_a_NO2_1.dat'
  open(11, file = fname, status = 'unknown')
    do i = 1, 27993
      read(11,*) lambda(i), sigma(i,1)
    end do 
  close(11)

  fname = 'sigma_a_NO2_2.dat'
  open(11, file = fname, status = 'unknown')
    do i = 1, 27993
      read(11,*) lambda(i), sigma(i,2)
    end do 
  close(11)

  fname = 'sigma_a_NO2.dat'
  open(11, file = fname, status = 'unknown')
    do i = 1, 27993
      write(11,'(f10.5, 2e14.6)') lambda(i), sigma(i,1), sigma(i,2)
    end do 
  close(11)

end program sigma_NO2