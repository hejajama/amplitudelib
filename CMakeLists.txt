# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(amplitudelib CXX C Fortran)

# Optimization on (for debug flags, change to Debug)
set(CMAKE_BUILD_TYPE Debug)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)


if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(amplitudelib)
add_subdirectory(src)
add_subdirectory(fragmentation)
add_subdirectory(pdf)
#add_subdirectory(tools)


# specify where to place binaries and libraries
#set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY, ${CMAKE_BINARY_DIR}/lib)

set_target_properties(partondistribution PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
	ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
set_target_properties(fragmentation PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
	ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
set_target_properties(amplitude PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
	ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

set_target_properties(dipole_analyzer PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set_target_properties(cross_section PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set_target_properties(castor_spectra PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set_target_properties(tester PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
