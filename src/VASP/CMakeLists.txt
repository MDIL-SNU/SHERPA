add_library ( VASP calculator.c input.c my_mpi.c )
target_include_directories ( VASP PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES> .. )
target_compile_options ( VASP PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS> )
target_link_libraries ( VASP PUBLIC $<LINK_ONLY:MKL::MKL> SPS m )
