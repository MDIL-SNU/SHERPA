set ( SRC_FILES ../main.c ../art_nouveau.c ../dimer.c ../utils.c )
add_executable ( sherpa_vasp ${SRC_FILES} )
target_include_directories ( sherpa_vasp PRIVATE ../../include ../../lib/vasp )
target_link_libraries ( sherpa_vasp PRIVATE sherpa ${LAPACK_LIBRARIES} vasp )
add_custom_command ( TARGET sherpa_vasp POST_BUILD COMMAND "${CMAKE_COMMAND}" -E copy "$<TARGET_FILE:sherpa_vasp>" "${CMAKE_SOURCE_DIR}/bin/$<TARGET_FILE_NAME:sherpa_vasp>" )
