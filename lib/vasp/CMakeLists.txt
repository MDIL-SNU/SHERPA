add_library ( vasp calculator.c my_mpi.c )
target_include_directories ( vasp PRIVATE ../../include )
#target_link_libraries ( vasp PUBLIC m )
