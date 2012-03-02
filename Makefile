# =========================================
# MOSAIK Banded Smith-Waterman Makefile
# (c) 2009 Michael Stromberg & Wan-Ping Lee
# =========================================

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= SWMain.cpp \
		BandedSmithWaterman.cpp \
		SmithWatermanGotoh.cpp \
		fasta_reader.cpp \
		fastq_reader.cpp \
		references.cpp \
		parameter_parser.cpp
OBJECTS= $(SOURCES:.cpp=.o)

# ----------------
# compiler options
# ----------------

CFLAGS=-Wall -O3
LDFLAGS=-Wl,-s
#CFLAGS=-g
PROGRAM=BandedSmithWaterman
LIBS=

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(FLAGS) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	@echo "Cleaning up."
	@rm -f *.o $(PROGRAM) *~
