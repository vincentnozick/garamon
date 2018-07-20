cmake_minimum_required(VERSION 3.5)

project(garamon_generator)

# add a "output" directory in the project binary directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# include eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()

# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O2 -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")

# files to compile
file(GLOB_RECURSE source_files src/*.cpp src/*.hpp)
message(STATUS "  files " ${source_files})

# compile
add_executable(garamon_generator ${source_files})

# link libs
include_directories(${EIGEN3_INCLUDE_DIR})

