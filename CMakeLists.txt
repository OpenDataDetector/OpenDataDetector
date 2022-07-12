# Open Data Dector project
#
# (c) 2021 CERN for the benefit of the ODD project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.14 )
project( OpenDataDetector VERSION 0.1 LANGUAGES CXX )

# Standard CMake include(s).
include( GNUInstallDirs )

# Set a better default for installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/install" CACHE PATH "default install path" FORCE)
endif()

# Find the DD4Hep package with corresponding library
find_package(DD4hep 1.11 REQUIRED CONFIG COMPONENTS DDCore DDDetectors)

include(FetchContent)

if(NOT TARGET ActsDD4hep::ActsDD4hep)
  # ensure glue library is available if not already part of the build
  FetchContent_Declare(actsdd4hep
    GIT_REPOSITORY https://github.com/acts-project/acts-dd4hep.git
    GIT_TAG v1.0.1)
  FetchContent_MakeAvailable(actsdd4hep)
endif()


# Add the main factory 
add_subdirectory(factory)


configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
  ${PROJECT_BINARY_DIR}/this_odd.sh
  @ONLY)
install(
  FILES ${PROJECT_BINARY_DIR}/this_odd.sh
  DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install the compact files
install(DIRECTORY xml DESTINATION share/OpenDataDetector)
