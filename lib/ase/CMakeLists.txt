set ( Python_FIND_STRATEGY LOCATION )
find_package ( Python REQUIRED COMPONENTS Interpreter Development )

add_library ( ase calculator.c input.c)
target_include_directories ( ase PRIVATE ../../include ${Python_INCLUDE_DIRS} )
target_link_libraries ( ase PRIVATE ${Python_LIBRARIES} )
