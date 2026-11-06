# Try to find the CLI11 package. If it’s not there,
# try downloading it and building it from source.
# This module exports targets CLI11::CLI11 that you can link to with:
#
#   `target_link_libraries(<your_target> PRIVATE CLI11::CLI11)`.
#

set(_cli11_precis "CLI11 (command line parser for C++11)")

include(FetchContent)

FetchContent_Declare(
  cli11_project
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG        cbd58a3696887b34c70949aef21a71735a0c2ad5 # v2.7.2, Aug 2, 2026
  GIT_SHALLOW    TRUE
  FIND_PACKAGE_ARGS 2.7 NAMES CLI11 CONFIG GLOBAL
)

FetchContent_MakeAvailable(cli11_project)

if (cli11_project_FOUND)
  message(STATUS "${_cli11_precis}, version ${cli11_project_VERSION} found")
elseif(cli11_project_POPULATED)
  message(STATUS "${_cli11_precis}, version v2.7.2 populated.")
else()
  message(FATAL_ERROR "${_cli11_precis} not found")
endif()