cmake_minimum_required(VERSION 3.10)

project(cmake_project_name_original_case)
set(CMAKE_BUILD_TYPE Release)

# include eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
else()
    include_directories("/usr/include/eigen3") # manually specify the include location
endif()

# call the CMakeLists.txt to make the documentation (Doxygen)
# > 'make html' to generate the documentation
add_subdirectory(doc)

# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++14 -Wno-return-local-addr")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-report ") # comilation profiling:  -Q -ftime-report
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-elide-constructors") # force the compilor to use our move 'copy constructor'

# files to compile
file(GLOB_RECURSE source_files src/cmake_project_name_original_case/*.cpp)
file(GLOB_RECURSE header_files src/cmake_project_name_original_case/*.hpp)

# display info
message(STATUS "  sources")
foreach(src_file ${source_files})
    message("          " ${src_file})
endforeach()
message(STATUS "  headers")
foreach(src_file ${header_files})
    message("          " ${src_file})
endforeach()

# generate a library
add_library(cmake_project_name_original_case SHARED ${source_files} ${header_files})
#add_library(cmake_project_name_original_case STATIC ${source_files} ${header_files})

# include directory path
include_directories(src)
include_directories(${EIGEN3_INCLUDE_DIR})

# install lib
install(FILES ${header_files} ${source_files} DESTINATION /usr/local/include/cmake_project_name_original_case)
install(TARGETS cmake_project_name_original_case
        RUNTIME DESTINATION /usr/local/bin
        LIBRARY DESTINATION /usr/local/lib
        ARCHIVE DESTINATION /usr/local/lib)

