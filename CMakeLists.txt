cmake_minimum_required(VERSION 3.19)
enable_language(Fortran) 
project(PhotoChemistry)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)
add_executable(e__main
    e__main.f90
    v__tdec.f90
    c__prm.f90
    p__io.f90
    p__search.f90
    v__Venus.f90
    v__Earth.f90
    v__Mars.f90
    v__Jupiter.f90
    p__EUVAC.f90
    p__UV.f90
    p__eddy_diffusion.f90
    p__molecular_diffusion.f90
    p__photochem_opticaldepth.f90
    p__photochem_rate.f90
    p__photochem_transport.f90
    p__photochem_scheme.f90
    Mars/Neutral_atmosphere/v__in.f90
)

#ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
endif()
#gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
endif()
