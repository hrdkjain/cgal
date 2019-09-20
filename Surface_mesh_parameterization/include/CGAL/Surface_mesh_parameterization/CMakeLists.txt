#cmake -G Eclipse\ CDT4\ -\ Unix\ Makefiles -DCGAL_DIR=~/tmp/cgal/ .
PROJECT(Surface_mesh_parameterization)
SET( CMAKE_BUILD_TYPE Debug)
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules/" ${CMAKE_MODULE_PATH})
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations -pg -std=c++0x")



# Boost and its components
find_package(Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()
# CGAL and its components
FIND_PACKAGE( OpenCV REQUIRED )
find_package( PCL REQUIRED )
# Find CGAL and CGAL Qt5
find_package( CGAL REQUIRED )
find_package( CGAL QUIET COMPONENTS ImageIO)# MPFI)
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()
find_package( GMP REQUIRED ./Find)
if(NOT GMP_FOUND)
  message(STATUS "This project requires GMP and is not found, and will not be compiled.")
  return()  
endif()
find_package(PythonLibs 2.7 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
include_directories(${GMP_INCLUDE_DIR})
include( ${CGAL_USE_FILE} )
include_directories(../../)
include_directories(${PYTHON_INCLUDE_DIRS})

#set(CMAKE_BUILD_TYPE Release)
file(GLOB SOURCE_FILES *.cpp *.h)
ADD_EXECUTABLE(Main ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(Main ${OpenCV_LIBS} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${PCL_LIBRARIES} ${CPPL_LIBS} ${SMI_LIBS} ${PYTHON_LIBRARIES} ${ND_LIBS} -lCGAL -lm -lmpfr -lboost_filesystem -lboost_system -lgmpxx -lgmp)
# add the install targets
install (TARGETS Main DESTINATION ~/bin)

