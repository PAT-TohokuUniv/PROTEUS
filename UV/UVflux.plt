reset
set term postscript enhanced color

set xlabel font "Arial,25"
set ylabel font "Arial,25"
set cblabel font "Arial,20"
set tics font "Arial,25"
set key font "Arial,20"
set title font"Arial,20"

set xlabel 'Wavelength [nm]'
set ylabel 'Solar flux'

set xlabel offset 0,-2
set ylabel offset -6,0

set log y
set format y "10^{%L}"

set xr [100:250]

set out "UVflux.eps"
plot "UVflux_1nm.dat"   u 1:2 w l, \
	 "UVflux_small.dat" u 1:2 w l