# Unit tests of the Garamon generator and of some generated algebras.
#   cmake --build build && ctest --test-dir build --output-on-failure
#
# - generator_tests:  unit tests of the generator itself (src/), linked with garamon_core
# - <algebra>_tests:  the algebra is generated at build time by garamon_generator (from conf/ or tests/conf/)
#                     and its library is tested against an independent reference implementation
# - <algebra>_sample: the sample program generated with the algebra is compiled and run

# Catch.cmake (catch_discover_tests): shipped next to Catch2Config.cmake, or in 'extras' of the source tree

if(DEFINED CATCH2_MODULE_PATH)
  list(APPEND CMAKE_MODULE_PATH "${CATCH2_MODULE_PATH}")
endif()
include(Catch)

function(garamon_set_warnings target)
  if(MSVC)
      target_compile_options(${target} PRIVATE /W3)
  else()
      target_compile_options(${target} PRIVATE -Wall)
  endif()
endfunction()


# ------------------------------------------------------------------------------------------------ generator

add_executable(generator_tests
  generator/test_utilities.cpp
  generator/test_product_tools.cpp
  generator/test_metric_tools.cpp
  generator/test_config_parser.cpp
  generator/test_metadata.cpp
  generator/test_directory.cpp
  generator/test_product_to_string.cpp
  generator/test_generation.cpp
  generator/test_clifford_reference.cpp
)
target_link_libraries(generator_tests PRIVATE garamon_core Catch2::Catch2WithMain)
target_include_directories(generator_tests PRIVATE generator algebras/reference)
target_compile_definitions(generator_tests PRIVATE
  GARAMON_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
  GARAMON_TEST_OUTPUT_DIR="${CMAKE_CURRENT_BINARY_DIR}/generator_output"
)
garamon_set_warnings(generator_tests)
catch_discover_tests(generator_tests TEST_PREFIX "generator: " WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


# ------------------------------------------------------------------------------------------------ generated algebras

set(garamon_generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(GLOB_RECURSE garamon_template_files CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/data/*")

# garamon_add_algebra_tests(<name> <configuration file>)
# generates the library of the algebra at build time, then builds and registers its tests (algebras/<name>_tests.cpp
# and any other algebras/<name>_*.cpp) and its sample program.
function(garamon_add_algebra_tests name conf)
  set(library_dir "${garamon_generated_dir}/garamon_${name}")
  set(stamp "${garamon_generated_dir}/${name}.stamp")

  # the generator refuses to overwrite an existing library: remove it first
  add_custom_command(
    OUTPUT "${stamp}"
    COMMAND ${CMAKE_COMMAND} -E rm -rf "${library_dir}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${garamon_generated_dir}"
    COMMAND $<TARGET_FILE:garamon_generator> -c "${conf}" -t "${PROJECT_SOURCE_DIR}/data/" -o "${garamon_generated_dir}/"
    COMMAND ${CMAKE_COMMAND} -E touch "${stamp}"
    DEPENDS garamon_generator "${conf}" ${garamon_template_files}
    COMMENT "Generating the ${name} library with garamon_generator"
    VERBATIM
  )
  add_custom_target(generate_${name} DEPENDS "${stamp}")

  # tests of the generated library (C++ standard as set in the top CMake script)
  file(GLOB test_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/algebras/${name}_*.cpp")
  # (add_dependencies only orders the targets: without OBJECT_DEPENDS, the objects are not rebuilt in the build that regenerates the library)
  set_source_files_properties(${test_sources} PROPERTIES OBJECT_DEPENDS "${stamp}")
  add_executable(${name}_tests ${test_sources})
  add_dependencies(${name}_tests generate_${name})
  target_include_directories(${name}_tests PRIVATE "${library_dir}/src" algebras/common algebras/reference)
  target_link_libraries(${name}_tests PRIVATE Eigen3::Eigen Catch2::Catch2WithMain)
  if(MSVC)
    target_compile_options(${name}_tests PRIVATE /bigobj)
  endif()
  garamon_set_warnings(${name}_tests)
  catch_discover_tests(${name}_tests TEST_PREFIX "${name}: ")

  # the generated sample must compile and run
  set(sample_source "${library_dir}/sample/src/main.cpp")
  set_source_files_properties("${sample_source}" PROPERTIES GENERATED TRUE OBJECT_DEPENDS "${stamp}")
  add_executable(${name}_sample "${sample_source}")
  add_dependencies(${name}_sample generate_${name})
  target_include_directories(${name}_sample PRIVATE "${library_dir}/src")
  target_link_libraries(${name}_sample PRIVATE Eigen3::Eigen)
  add_test(NAME "${name}: sample" COMMAND ${name}_sample)
endfunction()

foreach(algebra c3ga c5ga e3ga p3ga)
  garamon_add_algebra_tests(${algebra} "${PROJECT_SOURCE_DIR}/conf/${algebra}.conf")
endforeach()

# test algebras covering other code paths of the generator (see the comments in the configuration files)
foreach(algebra m3ga n3ga rc3ga rm3ga rp3ga)
  garamon_add_algebra_tests(${algebra} "${CMAKE_CURRENT_SOURCE_DIR}/conf/${algebra}.conf")
endforeach()
