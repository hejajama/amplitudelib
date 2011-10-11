CXXFLAGS = `gsl-config --cflags` -O3 -Wall -pedantic -fopenmp# -I ./libbci-1.1.0/ 
LDFLAGS = `gsl-config --libs` -lm 

SOURCES = src/main.cpp amplitudelib/amplitudelib.cpp \
	amplitudelib/datafile.cpp tools/interpolation.cpp tools/tools.cpp \
	amplitudelib/wave_function.cpp amplitudelib/virtual_photon.cpp \
	amplitudelib/ugd.cpp \
	pdf/pdf.cpp pdf/cteq.cpp pdf/mrst.cpp pdf/mrst99.cpp \
	fragmentation/fragmentation.cpp fragmentation/kkp.cpp \
	amplitudelib/xs.cpp 
FTSOURCES = fourier/fourier.c
OBJECTS = $(SOURCES:.cpp=.o)
FTOBJECTS = $(FTSOURCES:.c=.o)
FSOURCES = fragmentation/fragmentation_kkp.f pdf/CT10Pdf.f
FOBJECTS = $(FSOURCES:.f=.o)

all: amplitudelib 

amplitudelib: $(OBJECTS) $(FTOBJECTS) $(FOBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(FTOBJECTS) $(FOBJECTS) -lgfortran -o amplitude
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
.f.o:
	gfortran -c $< -c -o $@


clean:
	rm -f $(OBJECTS) $(FTOBJECTS) $(FOBJECTS)
	rm -f amplitude	
