CXX ?= $(shell command -v g++ 2>/dev/null || command -v c++ 2>/dev/null || command -v clang++ 2>/dev/null || echo g++)
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
BIN := $(ROOT)/bin/samovar_combine_annotations
SRC := $(ROOT)/src/cpp/combine_annotations.cpp
SPEARMAN_BIN := $(ROOT)/bin/samovar_spearman
SPEARMAN_SRC := $(ROOT)/src/cpp/spearman.cpp

.PHONY: all spearman clean

all: $(BIN)

$(BIN): $(SRC)
	mkdir -p $(ROOT)/bin
	$(CXX) $(CXXFLAGS) -o $@ $<

spearman: $(SPEARMAN_BIN)

$(SPEARMAN_BIN): $(SPEARMAN_SRC)
	mkdir -p $(ROOT)/bin
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $<

clean:
	rm -f $(BIN) $(SPEARMAN_BIN)
