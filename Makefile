CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm 

include filelist.m

all: amplitude

amplitude: $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp libamplitude.a
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.cpp -lgfortran -o amplitude
	ar cru libamplitude.a $(OBJECTS) $(FTOBJECTS) $(FOBJECTS)
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
.f.o:
	gfortran -c $< -c -o $@


clean:
	rm -f $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) src/main.o
	rm -f amplitude	
	rm libamplitude.a
