add_library ( vasp calculator.c input.c my_mpi.c )
target_include_directories ( vasp PRIVATE ../../include )
