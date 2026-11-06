cmake_minimum_required(VERSION 3.14)

project(cmake_project_name_sample)

# default build type of the single-config generators (Ninja, Makefiles), overridable with -DCMAKE_BUILD_TYPE=...
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# no compiler extensions (-std=c++20, not -std=gnu++20;
# the C++ standard comes from the cmake_project_name_original_case::cmake_project_name_original_case target)
set(CMAKE_CXX_EXTENSIONS OFF)

# require lib cmake_project_name_original_case (header-only, brings Eigen along)
# > set CMAKE_PREFIX_PATH to the cmake_project_name_original_case install prefix (and to Eigen if not found otherwise)
find_package(cmake_project_name_original_case CONFIG REQUIRED)
message(STATUS "lib cmake_project_name_original_case ${cmake_project_name_original_case_VERSION} found: ${cmake_project_name_original_case_DIR}")


# compile
add_executable(cmake_project_name_sample src/main.cpp)

# warnings (the optimization level comes from the build type)
if(MSVC)
    target_compile_options(cmake_project_name_sample PRIVATE /W3)
else()
    target_compile_options(cmake_project_name_sample PRIVATE -Wall)
endif()

# link
target_link_libraries(cmake_project_name_sample PRIVATE cmake_project_name_original_case::cmake_project_name_original_case)
