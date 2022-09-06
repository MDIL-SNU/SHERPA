CMAKE_MINIMUM_REQUIRED ( VERSION 2.8.12 )

SET ( LMP_PATH /home/wjdwltn928/lammps-stable_29Oct2020 )
PROJECT ( "SPS" )
SET ( PROJECT_VERSION_MAJOR 0 )
SET ( PROJECT_VERSION_MINOR 1 )

SET ( CMAKE_VERBOSE_MAKEFILE true )
SET ( SRC_FILES 
        main.c 
        calculator.c
        config.c
        dataset.c
        dimer.c
        input.c
        target.c
        utils.c )
SET ( CMAKE_C_COMPILER "mpiicc" )

ADD_COMPILE_OPTIONS ( -O3 -g -Wall )
LINK_LIBRARIES( lammps m )
LINK_DIRECTORIES( ${LMP_PATH}/src )
INCLUDE_DIRECTORIES ( ${LMP_PATH}/src )

ADD_EXECUTABLE ( SPS ${SRC_FILES} )
