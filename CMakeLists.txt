cmake_minimum_required ( VERSION 3.13 )
project ( SHERPA )
enable_language ( C )

find_package ( LAPACK REQUIRED )
set ( CMAKE_C_COMPILER "mpicc" )
set ( CMAKE_VERBOSE_MAKEFILE TRUE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
add_compile_options ( -O3 -g -Wall )
file ( MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

add_subdirectory ( utils/extractor )
add_subdirectory ( utils/kmc )
add_subdirectory ( lib/lmp )
add_subdirectory ( lib/vasp )
add_subdirectory ( lib/ase )
add_subdirectory ( src )

add_dependencies ( sherpa_lmp sherpa lmp )
add_dependencies ( sherpa_vasp sherpa vasp )
add_dependencies ( sherpa_ase sherpa ase )

install ( TARGETS extractor kmc sherpa_lmp sherpa_vasp sherpa_ase DESTINATION ${CMAKE_SOURCE_DIR}/bin )
