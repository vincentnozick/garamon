cmake_minimum_required(VERSION 3.14)

project(cmake_project_name_original_case VERSION 1.0.0 LANGUAGES CXX)

# default build type of the single-config generators (Ninja, Makefiles), overridable with -DCMAKE_BUILD_TYPE=...
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# no compiler extensions (-std=c++20, not -std=gnu++20;
# the C++ standard is set in the 'target_compile_features' command for the library target)
set(CMAKE_CXX_EXTENSIONS OFF)

# include eigen
find_package(Eigen3 REQUIRED)


# python binding
option(BUILD_PYTHON "Build the Python binding (requires Python 3 and pybind11)" OFF)


# call the CMakeLists.txt to make the documentation (Doxygen)
# > 'make html' to generate the documentation
find_package(Doxygen OPTIONAL_COMPONENTS QUIET)
if(DOXYGEN_FOUND)
    add_subdirectory(doc)
else()
    message(WARNING "Doxygen not found, skip dococumentation")
endif()


# GNU standard directories for installation
include(GNUInstallDirs)

# header files
file(GLOB_RECURSE header_files src/cmake_project_name_original_case/*.hpp)

# display info
message(STATUS "  headers")
foreach(src_file ${header_files})
    message("          " ${src_file})
endforeach()


# header-only library
add_library(cmake_project_name_original_case INTERFACE)
add_library(cmake_project_name_original_case::cmake_project_name_original_case ALIAS cmake_project_name_original_case)
target_include_directories(cmake_project_name_original_case INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(cmake_project_name_original_case INTERFACE cxx_std_20)

if(TARGET Eigen3::Eigen)
    message(STATUS "lib EIGEN3 ${Eigen3_VERSION} found")
    target_link_libraries(cmake_project_name_original_case INTERFACE Eigen3::Eigen)
else()
    target_include_directories(cmake_project_name_original_case INTERFACE ${EIGEN3_INCLUDE_DIR})
endif()

if (BUILD_PYTHON)
    if (CMAKE_VERSION VERSION_LESS 3.18)
        message(FATAL_ERROR "BUILD_PYTHON requires CMake 3.18 or more")
    endif()

    # the Python used to build the module is the one that will import it (an active venv is honoured)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    # locate pybind11 installed with 'pip install pybind11'
    execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
                    OUTPUT_VARIABLE pybind11_cmake_dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED HINTS "${pybind11_cmake_dir}")

    # cmake_project_name_original_case is header-only: the module has no runtime dependency
    pybind11_add_module(cmake_project_name_original_case_py src/cmake_project_name_original_case/PythonBindings.cpp)
    target_link_libraries(cmake_project_name_original_case_py PRIVATE cmake_project_name_original_case)

    # > 'cmake --install build --component python' to install the module into the site-packages
    # > -DPYTHON_USER_INSTALL=ON to install into the per-user site-packages (python -m site --user-site)
    option(PYTHON_USER_INSTALL "Install the Python module into the per-user site-packages" OFF)
    if (PYTHON_USER_INSTALL)
        execute_process(COMMAND "${Python_EXECUTABLE}" -c "import site; print(site.getusersitepackages())"
                        OUTPUT_VARIABLE python_default_install_dir
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        file(TO_CMAKE_PATH "${python_default_install_dir}" python_default_install_dir)
    else()
        set(python_default_install_dir "${Python_SITEARCH}")
    endif()
    set(PYTHON_INSTALL_DIR "${python_default_install_dir}" CACHE PATH "Installation directory of the Python module")
    install(TARGETS cmake_project_name_original_case_py
            LIBRARY DESTINATION "${PYTHON_INSTALL_DIR}" COMPONENT python
            RUNTIME DESTINATION "${PYTHON_INSTALL_DIR}" COMPONENT python)

    # > 'cmake --build build --target python_sample' to run sample/sample.py with the just-built module
    add_custom_target(python_sample
        COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=$<TARGET_FILE_DIR:cmake_project_name_original_case_py>"
                "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/sample/sample.py"
        DEPENDS cmake_project_name_original_case_py
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sample"
        VERBATIM)
endif()

# install lib
# > consumers use 'find_package(cmake_project_name_original_case CONFIG REQUIRED)' and 'target_link_libraries(... cmake_project_name_original_case::cmake_project_name_original_case)'
install(TARGETS cmake_project_name_original_case EXPORT cmake_project_name_original_caseTargets)

install(DIRECTORY src/cmake_project_name_original_case
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

set(cmake_project_name_original_case_config_install_dir ${CMAKE_INSTALL_DATADIR}/cmake/cmake_project_name_original_case)

install(EXPORT cmake_project_name_original_caseTargets
    NAMESPACE cmake_project_name_original_case::
    DESTINATION ${cmake_project_name_original_case_config_install_dir}
)

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/cmake_project_name_original_caseConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_project_name_original_caseConfig.cmake
    INSTALL_DESTINATION ${cmake_project_name_original_case_config_install_dir}
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cmake_project_name_original_caseConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_project_name_original_caseConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_project_name_original_caseConfigVersion.cmake
    DESTINATION ${cmake_project_name_original_case_config_install_dir}
)
