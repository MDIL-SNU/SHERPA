cmake_minimum_required ( VERSION 3.13 )
project ( Saddle_point_search )
enable_language ( C )

find_package ( MKL CONFIG REQUIRED )
set ( CMAKE_C_COMPILER "mpiicc" )
set ( CMAKE_VERBOSE_MAKEFILE TRUE )
add_compile_options ( -O3 -g -Wall )

add_executable ( EXTRACTOR extractor.c files.c )
add_executable ( KMC kmc.c files.c )
add_library ( SHERPA config.c dataset.c linalg.c target.c )

set ( SRC_FILES main.c art_nouveau.c dimer.c kappa_dimer.c utils.c )

add_subdirectory ( LMP )
add_executable ( SHERPA_LMP ${SRC_FILES} )
target_include_directories ( SHERPA_LMP PRIVATE LMP )
target_link_libraries ( SHERPA_LMP LMP )

# PyTorch
#find_package ( Torch REQUIRED )
#set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}" )
#add_subdirectory ( LMP )
#add_executable ( SHERPA_LMP ${SRC_FILES} )
#target_include_directories ( SHERPA_LMP PRIVATE LMP )
#target_link_libraries ( SHERPA_LMP LMP "${TORCH_LIBRARIES}" )
#set_property ( TARGET SHERPA_LMP PROPERTY CXX_STANDARD 14 )

add_subdirectory ( VASP )
add_executable ( SHERPA_VASP ${SRC_FILES} )
target_include_directories ( SHERPA_VASP PRIVATE VASP )
target_link_libraries ( SHERPA_VASP VASP )

