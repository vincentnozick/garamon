find_package (Doxygen)

message("cmake_project_name_original_case_SOURCE_DIR     " ${cmake_project_name_original_case_SOURCE_DIR})
message("cmake_project_name_original_case_BINARY_DIR     " ${cmake_project_name_original_case_BINARY_DIR})
message("DOXYGEN_EXECUTABLE  " ${DOXYGEN_EXECUTABLE})

if(DOXYGEN_FOUND)
  configure_file (${cmake_project_name_original_case_SOURCE_DIR}/doc/Doxyfile-html.cmake ${cmake_project_name_original_case_BINARY_DIR}/doc/Doxyfile-html)
  add_custom_target (html ${DOXYGEN_EXECUTABLE} ${cmake_project_name_original_case_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)
