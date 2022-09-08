CMAKE_MINIMUM_REQUIRED ( VERSION 3.10 )
PROJECT ( Saddle_point_search )

SET ( LAMMPS_SRC_DIRECTORY "" CACHE PATH "Path for lammps source" )
IF ( NOT LAMMPS_SRC_DIRECTORY STREQUAL "" AND EXISTS ${LAMMPS_SRC_DIRECTORY}/cmake/CMakeLists.txt )
  OPTION( BUILD_LIB "Build LAMMPS library" ON )
  ADD_SUBDIRECTORY ( ${LAMMPS_SRC_DIRECTORY}/cmake lammps )
ELSE()
  FIND_PACKAGE( LAMMPS REQUIRED )
ENDIF()

SET ( SRC_FILES 
        main.c 
        calculator.c
        config.c
        dataset.c
        dimer.c
        input.c
        target.c
        utils.c )
SET ( CMAKE_C_COMPILER "mpiicc" CACHE PATH "" )
ADD_COMPILE_OPTIONS ( -O3 -g -Wall )
LINK_LIBRARIES ( m )

ENABLE_LANGUAGE ( C )
ADD_EXECUTABLE ( SPS ${SRC_FILES} )
TARGET_LINK_LIBRARIES ( SPS LAMMPS::lammps )
