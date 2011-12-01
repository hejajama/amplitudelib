include pdf.m
include fragfun.m

AMPLITUDESOURCES = amplitudelib/amplitudelib.cpp \
	amplitudelib/datafile.cpp tools/interpolation.cpp tools/tools.cpp \
	amplitudelib/wave_function.cpp amplitudelib/virtual_photon.cpp \
	amplitudelib/ugd.cpp \
	amplitudelib/xs.cpp tools/interpolation2d.cpp 
SOURCES = $(AMPLITUDESOURCES) $(PDF_CPPSOURCES) $(FRAG_CPPSOURCES)
FTSOURCES = fourier/fourier.c
OBJECTS = $(SOURCES:.cpp=.o)
FTOBJECTS = $(FTSOURCES:.c=.o)
#FSOURCES = fragmentation/fragmentation_kkp.f pdf/CT10Pdf.f fragmentation/fragmentation_pkhff.f \
#	fragmentation/fragmentation_hknsff07.f fragmentation/dlib.f fragmentation/locate.f \
#	fragmentation/polint.f fragmentation/polin2.f fragmentation/grille_had_charged.f
FSOURCES = $(PDF_FSOURCES) $(FRAG_FSOURCES)
FOBJECTS = $(FSOURCES:.f=.o)


