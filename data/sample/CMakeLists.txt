cmake_minimum_required(VERSION 3.5)

project(cmake_project_name_sample)

# set the module path for "findLib.cmake"
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# require Eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()

# require lib cmake_project_name_original_case
find_package(cmake_project_name_upper_case)
if(${cmake_project_name_upper_case_FOUND})
    message(STATUS "lib cmake_project_name_upper_case found")
    message(STATUS "  include " ${cmake_project_name_upper_case_INCLUDE_DIRS})
    message(STATUS "  lib     " ${cmake_project_name_upper_case_LIBRARIES})
else()
    message(STATUS "lib cmake_project_name_upper_case should be set manually, please uncomment the commented lines in the CMakeLists.txt file")
    #set(cmake_project_name_upper_case_INCLUDE_DIRS "[path to]/garamon_cmake_project_name_original_case/src/")
    #set(cmake_project_name_upper_case_LIBRARIES "[path to]/garamon_cmake_project_name_original_case/build/libcmake_project_name_original_case.dylib")

    # for manual setting, check if valid
    if(EXISTS "${cmake_project_name_upper_case_INCLUDE_DIRS}")
        message(STATUS "include for lib cmake_project_name_original_case found: ${cmake_project_name_upper_case_INCLUDE_DIRS}" )
    else()
        message(FATAL_ERROR "include for lib cmake_project_name_original_case invalid or not found: ${cmake_project_name_upper_case_INCLUDE_DIRS}" )
    endif()

    # for manual setting, check if valid
    if(EXISTS "${cmake_project_name_upper_case_LIBRARIES}")
        message(STATUS "library cmake_project_name_original_case found: ${cmake_project_name_upper_case_LIBRARIES}" )
    else()
        message(FATAL_ERROR "library cmake_project_name_original_case not found: ${cmake_project_name_upper_case_LIBRARIES}" )
    endif()
endif()


# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++14")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# compile
add_executable(cmake_project_name_sample src/main.cpp)

# includes
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${cmake_project_name_upper_case_INCLUDE_DIRS})

# link
target_link_libraries(cmake_project_name_sample ${cmake_project_name_upper_case_LIBRARIES})




