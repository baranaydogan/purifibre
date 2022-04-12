#!/bin/bash

cmakeExe=cmake
staticBuild=OFF

c_compiler=gcc
cxx_compiler=g++

rm -rf build

mkdir -p build
cd build

${cmakeExe} -DCMAKE_C_COMPILER=${c_compiler} -DCMAKE_CXX_COMPILER=${cxx_compiler} -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC_EXE=${staticBuild} ..
${cmakeExe} --build . --config Release --target install --parallel 16

cd ..

