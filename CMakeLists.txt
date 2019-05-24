#####################################################################
## brODErs solvers for reacting flows.                             ##
##                                                                 ##
## 2016 - Georgios Bellas-Chatzigeordis (Shocking)                 ##
## 2017 - Stefano Boccelli (LARSEN and MultiLARSEN)                ##
#####################################################################

# Basic setup of project and minimum required CMake version
cmake_minimum_required(VERSION 2.6)
project(brODErs++)
enable_language(CXX)
#set (Tutorial_VERSION_MAJOR 0) # Major version number
#set (Tutorial_VERSION_MINOR 1) # Minor version number

# ==================================
# =======   L I B R A R I E S
# ==================================

# --- Mutation++ library ---

set(MPP_INSTALL_DIR "" CACHE PATH "Mutation++ installation directory")
include_directories(${MPP_INSTALL_DIR}/include/mutation++)
#include_directories(${MPP_INSTALL_DIR}/include/eigen3)
include_directories(${MPP_INSTALL_DIR}/../thirdparty/eigen)
include_directories(${MPP_INSTALL_DIR}/include/)
link_directories(${MPP_INSTALL_DIR}/lib)

# --- Boost Library ---
# 
# Asks to the user the explicit path.
# BOOST_ROOT is the suggested location.
# Boost_NO_SYSTEM_PATHS makes sure that boost is only searched in the specified dir.

set(BOOST_DIRECTORY "" CACHE STRING "Boost library directory")
if(BOOST_DIRECTORY)
  set(BOOST_ROOT ${BOOST_DIRECTORY})
  set(Boost_NO_SYSTEM_PATHS ON)
  find_package( Boost REQUIRED )
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# ==================================
# =======   S O U R C E S
# ==================================

# Include Directories
include_directories(src/general)
include_directories(src/larsen)
include_directories(src/shocking)

# Descend into the src directory to build all targets and libraries
add_subdirectory(src)

# Main executable
get_property(brODErs++_SRCS GLOBAL PROPERTY brODErs++_SRCS)
set_source_files_properties(${brODErs++_SRCS} PROPERTIES COMPILE_FLAGS -cpp)
add_executable(brODErs++ ${brODErs++_SRCS})

###target_link_libraries(brODErs++ -lmutation++ ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES} )
target_link_libraries(brODErs++ -lmutation++ ${VTK_LIBRARIES} ${Boost_LIBRARIES})

# Installation
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS brODErs++ DESTINATION bin)

