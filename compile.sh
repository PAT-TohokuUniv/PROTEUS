#!/bin/csh -f
cd build
cmake ..
cmake --build .
cd ..
./e__main