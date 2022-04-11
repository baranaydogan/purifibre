#!/bin/bash

cmakeExe=cmake
buildType=Release #Release or Debug
staticBuild=OFF


c_compiler=gcc
cxx_compiler=g++

rm -rf build

mkdir -p build
cd build

${cmakeExe} -DCMAKE_C_COMPILER=${c_compiler} -DCMAKE_CXX_COMPILER=${cxx_compiler} -DCMAKE_BUILD_TYPE=${buildType} -DBUILD_STATIC_EXE=${staticBuild} ..
${cmakeExe} --build . --config ${buildType} --target install --parallel 16

cd ..

