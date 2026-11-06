# 3.24: FetchContent_Declare(... FIND_PACKAGE_ARGS) is used in extern/*/CMakeLists.txt
cmake_minimum_required(VERSION 3.24)

project(garamon_generator)

set(CMAKE_CXX_STANDARD 20)
message(DEBUG "CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD}")

# CMAKE_CXX_STANDARD 20 is a requirement:
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(DEBUG "CMAKE_CXX_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED}")

# To avoid portability issues, language extensions (enabled by default)
# are disabled before any target is created. This allows to stick to the
# language standard such as "-std=c++20" instead of "-std=gnu++20" (which
# is c++20 with GNU extensions):
set(CMAKE_CXX_EXTENSIONS OFF)
message(DEBUG "CMAKE_CXX_EXTENSIONS ${CMAKE_CXX_EXTENSIONS}")


option(GARAMON_BUILD_TESTS "Build the unit tests of the generator and of some generated algebras (requires Catch2)" ON)

# include Eigen (C++ template library for linear algebra)
add_subdirectory( "extern/Eigen" )

# include CLI11 (command line parser for C++11)
add_subdirectory( "extern/CLI11" )

# add a "output" directory in the project binary directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# files to compile
file(GLOB_RECURSE source_files src/*.cpp src/*.hpp)
list(REMOVE_ITEM source_files ${PROJECT_SOURCE_DIR}/src/main.cpp)
message(STATUS "${CMAKE_PROJECT_NAME} source files:")
foreach(src_file ${source_files})
  message(STATUS "  ${src_file}")
endforeach()

# the generator itself, as a library (also used by the unit tests)
add_library(garamon_core STATIC ${source_files})
target_include_directories(garamon_core PUBLIC src)
target_link_libraries(garamon_core PUBLIC Eigen3::Eigen)
if(MSVC)
    target_compile_options(garamon_core PRIVATE /W3)
else()
    target_compile_options(garamon_core PRIVATE -Wall)
endif()

# the command line generator
add_executable(garamon_generator src/main.cpp)

# link libs
target_link_libraries(garamon_generator PRIVATE
  garamon_core
  CLI11::CLI11
)

# compilation flags
if(MSVC)
    target_compile_options(garamon_generator PRIVATE /W3)
else()
    target_compile_options(garamon_generator PRIVATE -Wall)
endif()

# unit tests: 'ctest' from the build directory
if(GARAMON_BUILD_TESTS)
  # include Catch2 (modern, C++-native, header-only testing framework)
  add_subdirectory( "extern/Catch2" )

  enable_testing()
  add_subdirectory(tests)
endif()

# install

include(GNUInstallDirs) # GNU standard directories for installation

install(TARGETS garamon_generator RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(DIRECTORY data DESTINATION .)
install(DIRECTORY conf DESTINATION .)