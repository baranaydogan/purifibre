cmake_minimum_required(VERSION 3.15.0)

PROJECT(Purifibre)

# Using C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_STATIC_EXE "Build static executable" OFF)

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

if(UNIX AND NOT APPLE)

   # Building for Linux or Apple
   message(STATUS "Building for Linux")

   # Add zlib support
   find_package(ZLIB REQUIRED) 
   add_definitions(-D HAVE_ZLIB)

   # USE THE FOLLOWING FOR STATIC LINKING. 
   # THIS REQUIRES A SYSTEM WITH libpthread.a and libz.a. ARCH LINUX BASED SYSTEMS DO NOT HAVE THESE.
   if (BUILD_STATIC_EXE)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -static -static-libgcc -static-libstdc++ -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive")
   endif()

	
elseif(MSVC)

   # Building for Windows
   message(STATUS "Building for Windows")
   add_definitions(-D BUILD_FOR_WINDOWS)

   # Enable static compilation
   # set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   # set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Release>:Release>")
   set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>")

   # Disable some compiler warnings for cl
   add_compile_options(/wd4244) # conversion from 'double' to 'x', possible loss of data
   add_compile_options(/wd4267) # conversion from 'size_t' to 'x', possible loss of data
   add_compile_options(/wd4996) # 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
   add_compile_options(/wd4305) # truncation from 'double' to 'float'
   add_compile_options(/wd4101) # unreferenced local variable

   # No zlib support for Windows at the moment

else()
 message(FATAL_ERROR "This operating system is not supported")
 return()
endif()

# List all source files
file(GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c*")

# Remove purifibre executable
list(REMOVE_ITEM SRCS "src/main.cpp")

add_library(Purifibre STATIC ${SRCS})

set(CMDNAME purifibre)
add_executable(${CMDNAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)

if (BUILD_STATIC_EXE)
   if(UNIX AND NOT APPLE)
      target_link_libraries(${CMDNAME} Purifibre)
   else()
      target_link_libraries(${CMDNAME} Purifibre PRIVATE Threads::Threads)
   endif()
else()
   if(UNIX AND NOT APPLE)
      target_link_libraries(${CMDNAME} Purifibre Threads::Threads pthread z)
   else()
      target_link_libraries(${CMDNAME} Purifibre Threads::Threads)
   endif()
endif()

INSTALL(TARGETS ${CMDNAME} DESTINATION ${CMAKE_BINARY_DIR}/install/bin)
INSTALL(CODE "configure_file(${CMAKE_BINARY_DIR}/libPurifibre.a ${CMAKE_BINARY_DIR}/install/lib/libPurifibre.a COPYONLY)")
