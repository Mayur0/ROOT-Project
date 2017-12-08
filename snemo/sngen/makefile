OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-g -O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) \
	-I../common

# Debug flags
# CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` `root-config --ldflags` -lGeom -lgfortran -lm
# LDFLAGS += -g

sngen: sngen.C 
	$(CXX) $(CFLAGS) sngen.C
	$(CXX) -o sngen sngen.o $(LDFLAGS)
	rm sngen.o

sngen_ideal: sngen_ideal.C 
	$(CXX) $(CFLAGS) sngen_ideal.C
	$(CXX) -o sngen_ideal sngen_ideal.o $(LDFLAGS)
	rm sngen_ideal.o

cassettegen: cassettegen.C 
	$(CXX) $(CFLAGS) cassettegen.C
	$(CXX) -o cassettegen cassettegen.o $(LDFLAGS)
	rm cassettegen.o

readOut: readOut.C 
	$(CXX) $(CFLAGS) readOut.C
	$(CXX) -o readOut readOut.o $(LDFLAGS)
	rm readOut.o

sngen_ideal_2: sngen_ideal_2.C 
	$(CXX) $(CFLAGS) sngen_ideal_2.C
	$(CXX) -o sngen_ideal_2 sngen_ideal_2.o $(LDFLAGS)
	rm sngen_ideal_2.o
