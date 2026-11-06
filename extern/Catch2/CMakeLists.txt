# Try to find the Catch2 package. If it’s not there,
# try downloading it and building it from source.
# This module exports:
# 
# - targets: Catch2::Catch2WithMain and Catch2::Catch2,
#   which you can refer to using the `target_link_libraries` command,
# 
# - the CATCH2_MODULE_PATH variable (set in the parent scope) for use in
#   the `tests` script, which needs the modules provided by Catch2.

set(_catch2_precis "Catch2 (unit testing framework for C++)")

include(FetchContent)

FetchContent_Declare(
  catch2_project
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        08092139210881f01a34e60c5f2fd6afa0ea4024   # v3.16.1, October 2026
  GIT_SHALLOW    TRUE
  FIND_PACKAGE_ARGS 3.16 NAMES Catch2 CONFIG GLOBAL
)

FetchContent_MakeAvailable(catch2_project)

if (catch2_project_FOUND)

  set(CATCH2_MODULE_PATH "${catch2_project_DIR}" PARENT_SCOPE)
  message(STATUS "${_catch2_precis}, version ${catch2_project_VERSION} found")

elseif(catch2_project_POPULATED)

  set(CATCH2_MODULE_PATH "${catch2_project_SOURCE_DIR}/extras" PARENT_SCOPE)

  # file(READ "${catch2_project_SOURCE_DIR}/src/catch2/catch_version_macros.hpp" _hdr)
  # string(REGEX MATCH "CATCH_VERSION_MAJOR ([0-9]+)" _ "${_hdr}")
  # set(_major "${CMAKE_MATCH_1}")
  # string(REGEX MATCH "CATCH_VERSION_MINOR ([0-9]+)" _ "${_hdr}")
  # set(_minor "${CMAKE_MATCH_1}")
  # string(REGEX MATCH "CATCH_VERSION_PATCH ([0-9]+)" _ "${_hdr}")
  # set(_patch "${CMAKE_MATCH_1}")
  # set(CATCH2_VERSION "${_major}.${_minor}.${_patch}")
  # message(STATUS "${_catch2_precis}, version ${CATCH2_VERSION} populated.")
  message(STATUS "${_catch2_precis}, version v3.16.1 populated.")

else()

  message(FATAL_ERROR "${_catch2_precis} not found")

endif()