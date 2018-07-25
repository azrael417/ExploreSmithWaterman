# =========================================
# MOSAIK Banded Smith-Waterman Makefile
# (c) 2009 Michael Stromberg & Wan-Ping Lee
# =========================================

#kokkos path
KOKKOS_PATH=kokkos
KOKKOS_DEVICES=OpenMP
KOKKOS_ARCH=HSW

# ----------------
# compiler type
# ----------------
CXX = CC

#main target goes here
PROGRAM=SmithWaterman
all: $(PROGRAM)

#include kokkos
include kokkos/Makefile.kokkos

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= SWMain.cpp \
		BandedSmithWaterman.cpp \
		SmithWatermanGotoh.cpp \
		fasta_reader.cpp \
		fastq_reader.cpp \
		references.cpp \
		parameter_parser.cpp \
		$(KOKKOS_SRC)
OBJECTS= $(SOURCES:.cpp=.o)

# ----------------
# compiler options
# ----------------
export CXXFLAGS = -Wall -O3 $(KOKKOS_CXXFLAGS)
#PROGRAM=SmithWaterman
LIBS=$(KOKKOS_EXTRA_LIBS)

#all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean

clean: kokkos-clean
	@echo "Cleaning up."
	@rm -f $(OBJECTS) $(PROGRAM) *~
