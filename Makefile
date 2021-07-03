# Makefile for SPICE programs

SPICELIBDIR=/usr/local/spice/toolkit/lib

RM=/bin/rm -f

FC=gfortran
FFLAGS=-fPIC $(DEBUG)

LD=gfortran
LDOPTS=-L$(SPICELIBDIR) -Wl,-rpath=$(SPICELIBDIR)
LDLIBS=-lspice

CASSINI_OBJS=cassini.o

cassini: $(CASSINI_OBJS)
	$(LD) -o $(@) $(LDOPTS) $(CASSINI_OBJS) $(LDLIBS)

clean:
	$(RM) *.o cassini
