cmake_minimum_required(VERSION 3.16)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Plugin build type" FORCE)
endif()
project(dseams_lammps_plugin LANGUAGES CXX)

# Out-of-tree LAMMPS runtime plugin. Point LAMMPS_SRC at the src directory
# of the LAMMPS that will load the plugin. The version string in that
# tree's version.h is written into the plugin and checked by LAMMPS.

if(NOT LAMMPS_SRC)
  message(FATAL_ERROR "Set LAMMPS_SRC to the LAMMPS src directory")
endif()
if(NOT SEAMS_INCLUDE)
  message(FATAL_ERROR "Set SEAMS_INCLUDE to the directory that holds seams_insitu.h")
endif()
if(NOT SEAMS_LIBRARY)
  message(FATAL_ERROR "Set SEAMS_LIBRARY to libyodaLib.so")
endif()

find_file(LAMMPS_VERSION_H version.h PATHS "${LAMMPS_SRC}" NO_DEFAULT_PATH REQUIRED)
find_file(LAMMPS_FIX_H fix.h PATHS "${LAMMPS_SRC}" NO_DEFAULT_PATH REQUIRED)

option(DSEAMS_LAMMPS_MPI "Build for an MPI LAMMPS. Off uses the LAMMPS MPI stubs." ON)

add_library(dseams_lammps MODULE fix_dseams.cpp)
target_compile_definitions(dseams_lammps PRIVATE DSEAMS_LAMMPS_PLUGIN)
target_compile_features(dseams_lammps PRIVATE cxx_std_17)
if(DSEAMS_LAMMPS_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(dseams_lammps PRIVATE MPI::MPI_CXX)
else()
  # The stub header must win over any system mpi.h.
  target_include_directories(dseams_lammps BEFORE PRIVATE "${LAMMPS_SRC}/STUBS")
endif()
target_include_directories(dseams_lammps PRIVATE "${LAMMPS_SRC}" "${SEAMS_INCLUDE}")
target_link_libraries(dseams_lammps PRIVATE "${SEAMS_LIBRARY}")
set_target_properties(dseams_lammps PROPERTIES
  PREFIX "lib"
  OUTPUT_NAME "dseams_lammps"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

get_filename_component(SEAMS_LIBDIR "${SEAMS_LIBRARY}" DIRECTORY)
set_target_properties(dseams_lammps PROPERTIES
  BUILD_RPATH "${SEAMS_LIBDIR}"
  INSTALL_RPATH "${SEAMS_LIBDIR}")

message(STATUS "d-SEAMS LAMMPS plugin uses headers in ${LAMMPS_SRC}")
message(STATUS "Load it with: plugin load ${CMAKE_BINARY_DIR}/libdseams_lammps.so")
