add_executable ( extractor extractor.c ../files.c )
target_include_directories ( extractor PRIVATE .. )
add_custom_command ( TARGET extractor POST_BUILD COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:extractor>" "${CMAKE_SOURCE_DIR}/bin/$<TARGET_FILE_NAME:extractor>" )
#install( TARGETS extractor DESTINATION ${CMAKE_SOURCE_DIR}/bin )
