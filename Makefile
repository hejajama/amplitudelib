CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm 
FORTRANFLAGS = -O3 

include filelist.m

all: amplitude

amplitude: $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp -lgfortran -o amplitude
	ar cru libamplitude.a $(OBJECTS) $(FTOBJECTS) $(FOBJECTS)

interpolator: tools/interpolator.o tools/interpolation.o
	g++ $(CXXFLAGS) $(LDFLAGS) tools/interpolator.o tools/interpolation.o -o interpolator

.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
.f.o:
	gfortran $(FORTRANFLAGS) -c $< -c -o $@


clean:
	rm -f $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.o
	rm -f amplitude	
	rm libamplitude.a
