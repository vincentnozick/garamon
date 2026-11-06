# Try to find the Eigen package. If it’s not there,
# try downloading it and building it from source.
# This module exports target Eigen3::Eigen that you can link to with:
#
#   `target_link_libraries(<your_target> PRIVATE Eigen3::Eigen)`.
# 

set(_eigen_precis "Eigen (C++ template library for linear algebra)")

include(FetchContent)

FetchContent_Declare(
  eigen_project
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        bc3b39870ecb690a623a3f49149a358b95c5781d # 5.0.1, Nov 8, 2025
  GIT_SHALLOW    TRUE
  FIND_PACKAGE_ARGS 5.0 NAMES Eigen3 CONFIG GLOBAL
)

FetchContent_MakeAvailable(eigen_project)

if (eigen_project_FOUND)
  message(STATUS "${_eigen_precis}, version ${eigen_project_VERSION} found")
elseif(eigen_project_POPULATED)
  message(STATUS "${_eigen_precis}, version 5.0.1 populated.")
else()
  message(FATAL_ERROR "${_eigen_precis} not found")
endif()
