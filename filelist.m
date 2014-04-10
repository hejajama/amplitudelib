include pdf.m
include fragfun.m

AMPLITUDESOURCES = amplitudelib/amplitudelib.cpp \
	amplitudelib/datafile.cpp tools/interpolation.cpp tools/tools.cpp \
	amplitudelib/wave_function.cpp amplitudelib/virtual_photon.cpp \
	amplitudelib/ugd.cpp \
	tools/interpolation2d.cpp amplitudelib/qcd.cpp
SOURCES = $(AMPLITUDESOURCES) $(PDF_CPPSOURCES) $(FRAG_CPPSOURCES)
FTSOURCES = fourier/fourier.c
OBJECTS = $(SOURCES:.cpp=.o)
FTOBJECTS = $(FTSOURCES:.c=.o)
FSOURCES = $(PDF_FSOURCES) $(FRAG_FSOURCES)
FOBJECTS = $(FSOURCES:.f=.o)


