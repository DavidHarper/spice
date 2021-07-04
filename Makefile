# Makefile for SPICE programs

SPICELIBDIR=/usr/local/spice/toolkit/lib

RM=/bin/rm -f

FC=gfortran
FFLAGS=-fPIC $(DEBUG)

LD=gfortran
LDOPTS=-L$(SPICELIBDIR) -Wl,-rpath=$(SPICELIBDIR)
LDLIBS=-lspice

all: analyseintegration cassini initialstate

ANALYSEINTEGRATION_OBJS=analyseintegration.o

analyseintegration: $(ANALYSEINTEGRATION_OBJS)
	$(LD) -o $(@) $(LDOPTS) $(ANALYSEINTEGRATION_OBJS) $(LDLIBS)

CASSINI_OBJS=cassini.o

cassini: $(CASSINI_OBJS)
	$(LD) -o $(@) $(LDOPTS) $(CASSINI_OBJS) $(LDLIBS)

INITIALSTATE_OBJS=initialstate.o

initialstate: $(INITIALSTATE_OBJS)
	$(LD) -o $(@) $(LDOPTS) $(INITIALSTATE_OBJS) $(LDLIBS)

clean:
	$(RM) *.o analyseintegration cassini initialstate
