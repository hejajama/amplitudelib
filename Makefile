CXXFLAGS = `gsl-config --cflags` -O2 -Wall -pedantic -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lgfortran -lm 
FORTRANFLAGS = -O2 

include filelist.m

all: amplitude

amplitude: $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp src/amplitude.cpp src/cross_section.cpp
	#g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp -lgfortran -o amplitude
	ar cru libamplitude.a $(OBJECTS) $(FTOBJECTS) $(FOBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) src/amplitude.cpp libamplitude.a -o amplitude
	g++ $(CXXFLAGS) $(LDFLAGS) src/cross_section.cpp libamplitude.a -o cross_section

interpolator: tools/interpolator.o tools/interpolation.o
	g++ $(CXXFLAGS) $(LDFLAGS) tools/interpolator.o tools/interpolation.o -o interpolator

tester: tools/tester.cpp libamplitude.a
	g++ $(CXXFLAGS) $(LDFLAGS) tools/tester.cpp libamplitude.a -o tester 

fit: tools/f2fit.o tools/lhcfit.o
	g++ $(CXXFLAGS) $(LDFLAGS) tools/f2fit.o libamplitude.a -o f2fit
	g++ $(CXXFLAGS) $(LDFLAGS) tools/lhcfit.o libamplitude.a -o lhcfit

.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@ -I .
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
.f.o:
	gfortran $(FORTRANFLAGS) -c $< -c -o $@


clean:
	rm -f $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.o src/amplitude.o
	rm -f amplitude	
	rm libamplitude.a
