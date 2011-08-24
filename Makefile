CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm 

SOURCES = src/main.cpp amplitudelib/amplitudelib.cpp \
	amplitudelib/datafile.cpp tools/interpolation.cpp tools/tools.cpp 
FTSOURCES = fourier/fourier.c
OBJECTS = $(SOURCES:.cpp=.o)
FTOBJECTS = $(FTSOURCES:.c=.o)

all: amplitudelib 

amplitudelib: $(OBJECTS) $(FTOBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) -o amplitude
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f amplitude	
