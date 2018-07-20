# Find Libcmake_project_name_original_case
#
#  cmake_project_name_upper_case_FOUND        - True if Libcmake_project_name_original_case found.
#  cmake_project_name_upper_case_INCLUDE_DIRS - where to find Libcmake_project_name_original_case include files
#  cmake_project_name_upper_case_LIBRARIES    - where to find Libcmake_project_name_original_case binary libraries


# find include path
find_path (cmake_project_name_upper_case_INCLUDE_DIR 
		NAMES cmake_project_name_original_case
		HINTS /usr/local/include)

# find library file
find_library (cmake_project_name_upper_case_LIBRARY 
		NAMES cmake_project_name_original_case
		HINTS /usr/local/lib/ )


set(cmake_project_name_upper_case_LIBRARIES ${cmake_project_name_upper_case_LIBRARY})
set(cmake_project_name_upper_case_INCLUDE_DIRS ${cmake_project_name_upper_case_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set cmake_project_name_upper_case_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (cmake_project_name_upper_case DEFAULT_MSG 
                                   cmake_project_name_upper_case_LIBRARIES cmake_project_name_upper_case_INCLUDE_DIRS)

mark_as_advanced(cmake_project_name_upper_case_INCLUDE_DIR cmake_project_name_upper_case_LIBRARIE)





