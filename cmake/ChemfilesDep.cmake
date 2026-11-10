# Optional chemfiles dependency. On by default.
# find_package uses an installed copy. Otherwise FetchContent builds
# tag 0.10.4, which carries XYZ, PDB, GRO, LAMMPS, XTC, TRR, DCD,
# Amber NetCDF and mmCIF.

option(SEAMS_WITH_CHEMFILES "Read trajectories with chemfiles" ON)

if(SEAMS_WITH_CHEMFILES)
  find_package(chemfiles QUIET CONFIG)
  if(NOT chemfiles_FOUND)
    include(FetchContent)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
    FetchContent_Declare(
      chemfiles
      GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
      GIT_TAG 0.10.4
      GIT_SHALLOW TRUE)
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(CHFL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(CHFL_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(chemfiles)
  endif()
endif()
