reset

color_list="grey10 gray20 grey20 gray30 grey30 gray40 grey40 gray50 \
grey50 gray60 grey60 gray70 grey70 gray80 grey80 gray90 grey90 gray100 grey100 \
gray grey light-gray light-grey dark-gray dark-grey red light-red dark-red \
yellow light-yellow dark-yellow green light-green dark-green spring-green \
forest-green sea-green blue light-blue dark-blue midnight-blue navy \
medium-blue royalblue skyblue cyan light-cyan dark-cyan magenta light-magenta \
dark-magenta turquoise light-turquoise dark-turquoise pink light-pink \
dark-pink coral light-coral orange-red salmon light-salmon dark-salmon \
aquamarine khaki dark-khaki goldenrod light-goldenrod dark-goldenrod gold \
beige brown orange dark-orange violet dark-violet plum purple"

color_list_16="#1a1a1a #333333 #333333 #4d4d4d #4d4d4d #666666 #666666 #7f7f7f \
#7f7f7f #999999 #999999 #b3b3b3 #b3b3b3 #cccccc #cccccc #e5e5e5 #e5e5e5 \
#ffffff #ffffff #bebebe #bebebe #d3d3d3 #d3d3d3 #a9a9a9 #a9a9a9 #ff0000 \
#f03232 #8b0000 #ffff00 #ffffe0 #c8c800 #00ff00 #90ee90 #006400 #00ff7f \
#228b22 #2e8b57 #0000ff #add8e6 #00008b #191970 #000080 #0000cd #4169e1 \
#87ceeb #00ffff #e0ffff #008b8b #ff00ff #f055f0 #8b008b #40e0d0 #afeeee \
#00ced1 #ffc0cb #ffb6c1 #ff1493 #ff7f50 #f08080 #ff4500 #fa8072 #ffa07a \
#e9967a #7fffd4 #f0e68c #bdb76b #daa520 #eedd82 #b8860b #ffd700 #f5f5dc \
#a52a2a #ffa500 #ff8c00 #ee82ee #9400d3 #dda0dd #a020f0"

################################ output density #############################

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'Mixing ratio [ppm]'
set ylabel 'altitude [km]'

#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_vmr_major.eps"
set xr [1E-7:5E-4]
set yr [0:48]
plot "./output/density/vmr/vmr_SO2.dat"   u 2:1 w linespoints lw 3 title "SO2", \
     "./output/density/vmr/vmr_CO.dat"    u 2:1 w linespoints lw 3 title "CO", \
     "./output/density/vmr/vmr_H2O.dat"    u 2:1 w linespoints lw 3 title "H2O", \
     "./output/density/vmr/vmr_OCS.dat"     u 2:1 w linespoints lw 3 title "OCS", \
     "./output/density/vmr/vmr_H2SO4.dat" u 2:1 w linespoints lw 3 title "H2SO4", \
     "./output/density/vmr/vmr_S2.dat"   u 2:1 w linespoints lw 3 title "S2", \
     "./output/density/vmr/vmr_SO3.dat"    u 2:1 w linespoints lw 3 title "SO3", \


##############################

reset

color_list="grey10 gray20 grey20 gray30 grey30 gray40 grey40 gray50 \
grey50 gray60 grey60 gray70 grey70 gray80 grey80 gray90 grey90 gray100 grey100 \
gray grey light-gray light-grey dark-gray dark-grey red light-red dark-red \
yellow light-yellow dark-yellow green light-green dark-green spring-green \
forest-green sea-green blue light-blue dark-blue midnight-blue navy \
medium-blue royalblue skyblue cyan light-cyan dark-cyan magenta light-magenta \
dark-magenta turquoise light-turquoise dark-turquoise pink light-pink \
dark-pink coral light-coral orange-red salmon light-salmon dark-salmon \
aquamarine khaki dark-khaki goldenrod light-goldenrod dark-goldenrod gold \
beige brown orange dark-orange violet dark-violet plum purple"

color_list_16="#1a1a1a #333333 #333333 #4d4d4d #4d4d4d #666666 #666666 #7f7f7f \
#7f7f7f #999999 #999999 #b3b3b3 #b3b3b3 #cccccc #cccccc #e5e5e5 #e5e5e5 \
#ffffff #ffffff #bebebe #bebebe #d3d3d3 #d3d3d3 #a9a9a9 #a9a9a9 #ff0000 \
#f03232 #8b0000 #ffff00 #ffffe0 #c8c800 #00ff00 #90ee90 #006400 #00ff7f \
#228b22 #2e8b57 #0000ff #add8e6 #00008b #191970 #000080 #0000cd #4169e1 \
#87ceeb #00ffff #e0ffff #008b8b #ff00ff #f055f0 #8b008b #40e0d0 #afeeee \
#00ced1 #ffc0cb #ffb6c1 #ff1493 #ff7f50 #f08080 #ff4500 #fa8072 #ffa07a \
#e9967a #7fffd4 #f0e68c #bdb76b #daa520 #eedd82 #b8860b #ffd700 #f5f5dc \
#a52a2a #ffa500 #ff8c00 #ee82ee #9400d3 #dda0dd #a020f0"

################################ output density #############################

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'Mixing ratio [ppm]'
set ylabel 'altitude [km]'

#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_vmr_minor.eps"
set xr [1E-13:5E-6]
set yr [0:48]
plot "./output/density/vmr/vmr_H2.dat"   u 2:1 w linespoints lw 3 title "H2", \
     "./output/density/vmr/vmr_H2S.dat"    u 2:1 w linespoints lw 3 title "H2S", \
     "./output/density/vmr/vmr_S3.dat"    u 2:1 w linespoints lw 3 title "S3", \
     "./output/density/vmr/vmr_SO.dat"     u 2:1 w linespoints lw 3 title "SO", \
     "./output/density/vmr/vmr_SO2Cl2.dat" u 2:1 w linespoints lw 3 title "SO2Cl2", \
     "./output/density/vmr/vmr_ClSO2.dat"   u 2:1 w linespoints lw 3 title "ClSO2", \
     "./output/density/vmr/vmr_SH.dat"    u 2:1 w linespoints lw 3 title "SH", \



##############################

reset

color_list="grey10 gray20 grey20 gray30 grey30 gray40 grey40 gray50 \
grey50 gray60 grey60 gray70 grey70 gray80 grey80 gray90 grey90 gray100 grey100 \
gray grey light-gray light-grey dark-gray dark-grey red light-red dark-red \
yellow light-yellow dark-yellow green light-green dark-green spring-green \
forest-green sea-green blue light-blue dark-blue midnight-blue navy \
medium-blue royalblue skyblue cyan light-cyan dark-cyan magenta light-magenta \
dark-magenta turquoise light-turquoise dark-turquoise pink light-pink \
dark-pink coral light-coral orange-red salmon light-salmon dark-salmon \
aquamarine khaki dark-khaki goldenrod light-goldenrod dark-goldenrod gold \
beige brown orange dark-orange violet dark-violet plum purple"

color_list_16="#1a1a1a #333333 #333333 #4d4d4d #4d4d4d #666666 #666666 #7f7f7f \
#7f7f7f #999999 #999999 #b3b3b3 #b3b3b3 #cccccc #cccccc #e5e5e5 #e5e5e5 \
#ffffff #ffffff #bebebe #bebebe #d3d3d3 #d3d3d3 #a9a9a9 #a9a9a9 #ff0000 \
#f03232 #8b0000 #ffff00 #ffffe0 #c8c800 #00ff00 #90ee90 #006400 #00ff7f \
#228b22 #2e8b57 #0000ff #add8e6 #00008b #191970 #000080 #0000cd #4169e1 \
#87ceeb #00ffff #e0ffff #008b8b #ff00ff #f055f0 #8b008b #40e0d0 #afeeee \
#00ced1 #ffc0cb #ffb6c1 #ff1493 #ff7f50 #f08080 #ff4500 #fa8072 #ffa07a \
#e9967a #7fffd4 #f0e68c #bdb76b #daa520 #eedd82 #b8860b #ffd700 #f5f5dc \
#a52a2a #ffa500 #ff8c00 #ee82ee #9400d3 #dda0dd #a020f0"

################################ output density #############################

set term postscript enhanced color
#set term aqua
set log x
set xlabel 'Number density [m^{-3}]'
set ylabel 'altitude [km]'

#set xr [0.01:5]

set xlabel font "Arial,30"
set ylabel font "Arial,30"
set tics font "Arial,25"
set key font "Arial,15"
set title font"Arial,20"

set key spacing 1.5

set format x "10^{%L}"
#unset key

set xlabel offset 0,-2
set xtics offset 0,-1
set ylabel offset -2,0

set out "output_num_radical.eps"
set xr [1e6:1e13]
set yr [0:48]
plot "./output/density/num/SCl.dat"   u 2:1 w linespoints lw 3 title "SCl", \
     "./output/density/num/HSCl.dat"    u 2:1 w linespoints lw 3 title "HSCl", \
     "./output/density/num/Cl2.dat"    u 2:1 w linespoints lw 3 title "Cl2", \
     "./output/density/num/Cl.dat"     u 2:1 w linespoints lw 3 title "Cl", \
     "./output/density/num/S.dat" u 2:1 w linespoints lw 3 title "S", \
     "./output/density/num/H.dat"   u 2:1 w linespoints lw 3 title "H", \
     "./output/density/num/OH.dat"    u 2:1 w linespoints lw 3 title "OH", \
