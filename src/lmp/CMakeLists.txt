set ( SRC_FILES ../main.c ../art_nouveau.c ../dimer.c ../utils.c )
add_executable ( sherpa_lmp ${SRC_FILES} )
target_include_directories ( sherpa_lmp PRIVATE ../../include ../../lib/lmp )
target_link_libraries ( sherpa_lmp PRIVATE sherpa ${LAPACK_LIBRARIES} lmp )
add_custom_command ( TARGET sherpa_lmp POST_BUILD COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:sherpa_lmp>" "${CMAKE_SOURCE_DIR}/bin/$<TARGET_FILE_NAME:sherpa_lmp>" )
