add_executable(
	dipole_analyzer
	amplitude.cpp
	)

add_executable(
	cross_section
	cross_section.cpp
)
add_executable(
	dis
	dis.cpp
	data.cpp
)
add_executable(
	tester
	../tools/tester.cpp
)	

find_package(GSL REQUIRED)    #
target_link_libraries(
	dipole_analyzer
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	dis
	PRIVATE
	amplitude
	GSL::gsl
	GSL::gslcblas
)

target_link_libraries(
	cross_section
	PRIVATE
	amplitude
	partondistribution
	fragmentation
	GSL::gsl
	GSL::gslcblas
)


target_link_libraries(
        tester
        PRIVATE
        amplitude
        partondistribution
        fragmentation
        GSL::gsl
        GSL::gslcblas
)
