@SET cmakeExe="C:/Program Files/CMake/bin/cmake"
@SET staticBuild=ON

#@rmdir build /s /q
@mkdir build
@cd build
@%cmakeExe% -DBUILD_STATIC_EXE=%staticBuild% ..
@%cmakeExe% --build . --config Release --target install --parallel 8
@cd ..
