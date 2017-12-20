
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	amplitude
	STATIC
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	gitsha1.h
	amplitudelib.cpp
	datafile.cpp
	wave_function.cpp
	virtual_photon.cpp
	ugd.cpp
	hybrid.cpp
	dis.cpp
	qcd.cpp
	../tools/tools.cpp
	../tools/interpolation.cpp
	../tools/interpolation2d.cpp
	../fourier/fourier.c
)

include_directories(amplitude .)
