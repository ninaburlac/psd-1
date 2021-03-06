#
# Makefile for PropNuclei package
#
# 


ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


CXXFLAGS      = -O -Wall -fPIC -D_REENTRANT $(CXXOPTS)
CXXFLAGS     += $(ROOTCFLAGS)  -I /software/Externals/fftw/3.3.4/include

ifeq ($(DEBUG),1)
CXXFLAGS      += -D_DEBUG
endif

LIBS          = $(ROOTLIBS) -lMinuit

CXX           = g++
LD            = g++
LDFLAGS       = 

CXXINCLU     +=  -I. $(ROOTCFLAGS) -I $(ROOTSYS)/include 
LIBS          = $(ROOTLIBS) 
GLIBS         = $(ROOTGLIBS)


LDFLAGS      +=  $(ROOTGLIBS) -lRooFit -lRooFitCore -lMinuit


FFTWLIBS      = -L/software/Externals/fftw/3.3.4/lib
FFTWFLAGS     =  -lfftw3


PROGRAM       = processPSD 

SRCS          = processPSD_v5.cpp

OBJ           = $(patsubst %.cc, %.o, $(SRCS)) 

src           = $(addprefix ./,$(SRCS))
obj           = $(addprefix ./,$(OBJ))
prg           = $(addprefix ./,$(PROGRAM))


PROGRAMII       = readTier 
SRCSII	      = readTier.cpp 
OBJII           = $(patsubst %.cc, %.o, $(SRCSII)) 
srcII           = $(addprefix ./,$(SRCSII))
objII           = $(addprefix ./,$(OBJII))
prgII           = $(addprefix ./,$(PROGRAMII))

PROGRAMIII       = plotSF_AoE
SRCSIII	      = plotSF_AoE.cpp 
OBJIII           = $(patsubst %.cc, %.o, $(SRCSIII)) 
srcIII           = $(addprefix ./,$(SRCSIII))
objIII           = $(addprefix ./,$(OBJIII))
prgIII           = $(addprefix ./,$(PROGRAMIII))

PROGRAMIV       = plotSF_AoE_weight
SRCSIV       = plotSF_AoE_weight.cpp
OBJIV           = $(patsubst %.cc, %.o, $(SRCSIV))
srcIV           = $(addprefix ./,$(SRCSIV))
objIV           = $(addprefix ./,$(OBJIV))
prgIV           = $(addprefix ./,$(PROGRAMIV))


# -----------------------------------------------------------------------------


all:: $(prg) $(prgII) $(prgIII) $(prgIV) $(con)

$(con): $(objc) 
	$(CXX) -o $@ $^ $(LDFLAGS)  $(FFTWFLAGS) $(FFTWLIBS)

$(prg):$(obj) 
#	@if ! [ -x bin ]; then mkdir bin ; fi
	$(CXX) -o $@ $^ $(LDFLAGS) $(FFTWFLAGS) $(FFTWLIBS) $(CXXINCLU)
#	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

#%.o : %.cc %.h
#$(CXX) $(CXXFLAGS) -c $<   -o $@ 
#%.o : %.cc 
#$(CXX) $(CXXFLAGS) -c $<   -o $@ 
#%.o : %.c %.h
#	$(CXX) $(CXXFLAGS) -c $<   -o $@ 

%.o : %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $< $(CXXINCLU)  -o $@
%.o : %.cc %.h
	$(CXX) $(CXXFLAGS) -c $< $(CXXINCLU)  -o $@
%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< $(CXXINCLU)  -o $@
%.o : %.c %.h
	$(CXX) $(CXXFLAGS) -c $< $(CXXINCLU)  -o $@


clean::
	rm -f $(con)
	rm -f $(prg)	
	rm -f $(prgII)	
	rm -f $(prgIII)	
	rm -f $(prgIV)

# -----------------------------------------------------------------------------

