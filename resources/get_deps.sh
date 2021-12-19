#!/bin/bash

cd external

if [ ! -d "yaml-cpp" ] ; then
    git clone https://github.com/jbeder/yaml-cpp
fi

if [ ! -d "docopt.cpp" ] ; then
    git clone https://github.com/docopt/docopt.cpp
fi

if [ ! -d "fmt" ] ; then
    git clone https://github.com/fmtlib/fmt
fi

if [ ! -d "sundials" ] ; then
    git clone https://github.com/LLNL/sundials
fi

echo "Building yaml-cpp"
cd yaml-cpp
cmake -B build -DCMAKE_INSTALL_PREFIX=../ -DYAML_CPP_BUILD_TESTS=NO
cd build
make
make install
cd ../../

echo "Building docopt"
cd "docopt.cpp"
cmake -B build -DCMAKE_INSTALL_PREFIX=../ -DWITH_TESTS=NO
cd build
make
make install
cd ../../

echo "Building fmt"
cd fmt
cmake -B build -DCMAKE_INSTALL_PREFIX=../ -DFMT_TEST=NO
cd build
make
make install
cd ../../

echo "Building sundials"
cd sundials
cmake -B build -DCMAKE_INSTALL_PREFIX=../ \
               -DEXAMPLES_INSTALL_PATH=../sundials_examples \
               -DENABLE_MPI=YES \
               -DENABLE_LAPACK=YES
cd build
make
make install
cd ../../

echo "Done"
