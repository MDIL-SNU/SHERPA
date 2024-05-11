add_executable ( kmc kmc.c ../files.c )
target_include_directories ( kmc PRIVATE .. )
target_link_libraries ( kmc PRIVATE m )
add_custom_command ( TARGET kmc POST_BUILD COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:kmc>" "${CMAKE_SOURCE_DIR}/bin/$<TARGET_FILE_NAME:kmc>" )
