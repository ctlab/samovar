CXX ?= $(shell command -v g++ 2>/dev/null || command -v c++ 2>/dev/null || command -v clang++ 2>/dev/null || echo g++)
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
BIN := $(ROOT)/bin/samovar_combine_annotations
SRC := $(ROOT)/src/cpp/combine_annotations.cpp
SPEARMAN_BIN := $(ROOT)/bin/samovar_spearman
SPEARMAN_SRC := $(ROOT)/src/cpp/spearman.cpp
KMER_BIN := $(ROOT)/bin/samovar_kmer
KMER_SRC := $(ROOT)/src/cpp/kmer_profile.cpp

.PHONY: all spearman kmer clean

all: $(BIN) $(KMER_BIN)

$(BIN): $(SRC)
	mkdir -p $(ROOT)/bin
	$(CXX) $(CXXFLAGS) -o $@ $<

spearman: $(SPEARMAN_BIN)

$(SPEARMAN_BIN): $(SPEARMAN_SRC)
	mkdir -p $(ROOT)/bin
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $<

kmer: $(KMER_BIN)

$(KMER_BIN): $(KMER_SRC)
	mkdir -p $(ROOT)/bin
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $< -lz

clean:
	rm -f $(BIN) $(SPEARMAN_BIN) $(KMER_BIN)
