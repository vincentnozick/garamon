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

# files to compile
file(GLOB_RECURSE source_files src/*.cpp src/*.hpp)
message(STATUS "  files " ${source_files})

# compile
add_executable(garamon_generator ${source_files})

# link libs
include_directories(${EIGEN3_INCLUDE_DIR})

# compilation flags
target_compile_features(garamon_generator PRIVATE cxx_std_14) 
target_compile_options(garamon_generator PRIVATE -Wall -O2) 
