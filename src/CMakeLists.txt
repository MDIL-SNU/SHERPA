cmake_minimum_required ( VERSION 3.13 )
project ( Saddle_point_search )
enable_language ( C )

find_package ( MKL CONFIG REQUIRED )
set ( CMAKE_C_COMPILER "mpiicc" )
set ( CMAKE_VERBOSE_MAKEFILE TRUE )
add_compile_options ( -O3 -g -Wall )

add_executable ( EXTRACTOR extractor.c )
#add_executable ( KMC kmc.c )
add_library ( SPS config.c dataset.c linalg.c target.c )

set ( SRC_FILES main.c art_nouveau.c dimer.c kappa_dimer.c utils.c )

add_subdirectory ( LMP )
add_executable ( SPS_LMP ${SRC_FILES} )
target_include_directories ( SPS_LMP PRIVATE LMP )
target_link_libraries ( SPS_LMP LMP )

add_subdirectory ( VASP )
add_executable ( SPS_VASP ${SRC_FILES} )
target_include_directories ( SPS_VASP PRIVATE VASP )
target_link_libraries ( SPS_VASP VASP )
