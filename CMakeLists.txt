cmake_minimum_required(VERSION 3.15.0)

PROJECT(Purifibre_project)

set(CMDNAME purifibre)

# Using C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_STATIC_EXE "Build static executable" OFF)

find_package(Threads REQUIRED)

# List all source files
file(GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c*")

# Remove purifibre executable
list(REMOVE_ITEM SRCS "src/main.cpp")

# Set executable
add_executable(${CMDNAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
INSTALL(TARGETS ${CMDNAME} DESTINATION ${CMAKE_BINARY_DIR}/install/bin)

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

      add_library(Purifibre STATIC ${SRCS})
      target_link_libraries(${CMDNAME} Purifibre)
      INSTALL(CODE "configure_file(${CMAKE_BINARY_DIR}/libPurifibre.a ${CMAKE_BINARY_DIR}/install/lib/libPurifibre.a COPYONLY)")
   else()
      add_library(Purifibre SHARED ${SRCS})
      target_link_libraries(${CMDNAME} Purifibre Threads::Threads pthread z)
      INSTALL(CODE "configure_file(${CMAKE_BINARY_DIR}/libPurifibre.so ${CMAKE_BINARY_DIR}/install/lib/libPurifibre.so COPYONLY)")
   endif()


elseif(MSVC)

   # Building for Windows
   message(STATUS "Building for Windows")
   add_definitions(-D BUILD_FOR_WINDOWS)
   add_definitions(/bigobj)

   # Disable some compiler warnings for cl
   add_compile_options(/wd4244) # conversion from 'double' to 'x', possible loss of data
   add_compile_options(/wd4267) # conversion from 'size_t' to 'x', possible loss of data
   add_compile_options(/wd4996) # 'sprintf': This function or variable may be unsafe. Consider using sprintf_s instead.
   add_compile_options(/wd4305) # truncation from 'double' to 'float'
   add_compile_options(/wd4101) # unreferenced local variable
   add_compile_options(/wd4661) # no suitable definition provided
   add_compile_options(/wd4477) # sprintf format string for %lu requires argument of type 'unsigned long'

   # No zlib support for Windows at the moment

   if (BUILD_STATIC_EXE)
     add_library(libPurifibre STATIC ${SRCS})
     set_property(TARGET libPurifibre PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
     set_property(TARGET ${CMDNAME}   PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
     target_link_libraries(${CMDNAME} libPurifibre)
     INSTALL(CODE "configure_file(${CMAKE_BINARY_DIR}/Release/libPurifibre.lib ${CMAKE_BINARY_DIR}/install/lib/libPurifibre.lib COPYONLY)")
   else()
     # It is better to link against .lib instead of DLL so we will static compile libPurifibre.
     add_library(libPurifibre STATIC ${SRCS})
     set_property(TARGET libPurifibre PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
     set_property(TARGET ${CMDNAME}   PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
     target_link_libraries(${CMDNAME} libPurifibre Threads::Threads)
     INSTALL(CODE "configure_file(${CMAKE_BINARY_DIR}/Release/libPurifibre.lib ${CMAKE_BINARY_DIR}/install/lib/libPurifibre.lib COPYONLY)")
   endif()

else()
 message(FATAL_ERROR "This operating system is not supported")
 return()
endif()
