add_library( sherpa config.c dataset.c linalg.c target.c )
target_include_directories ( sherpa PRIVATE ../include )

add_subdirectory ( lmp )
add_subdirectory ( vasp )
add_subdirectory ( ase )
