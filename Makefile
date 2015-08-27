#-------------------------------------------------------------------------------
# g++ and its options
#-------------------------------------------------------------------------------
GCC_VER := $(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.`)

CXXFLAGS = -O2 -Wall -Wno-unused-function -Wno-array-bounds -D__STDC_FORMAT_MACROS -funroll-loops -fprefetch-loop-arrays -fopenmp -I. -std=c++0x -static-libgcc
LIB = -lm -lz -lpthread

ifeq "4.5" "$(word 1, $(sort 4.5 $(GCC_VER)))"
	CXXFLAGS += -static-libstdc++
endif

ifneq ($(disablempopcnt), 1)
	CXXFLAGS += -mpopcnt
endif

#-------------------------------------------------------------------------------
# standalone headers
#-------------------------------------------------------------------------------
STANDALONE_H = a_star_node.h hash.h hash_map.h hash_map_single_thread.h \
			   hash_set.h hash_set_single_thread.h functional.h \
			   khash.h kseq.h pool.h hash_table.h hash_table_single_thread.h \
			   utils.h hmmer3b_parser.h kmer.h mem_file_checker-inl.h\
			   most_probable_path.h nucl_kmer_generator.h profile_hmm.h \
			   definitions.h prot_kmer_generator.h rank_and_select.h sdbg_multi_io.h

DEPS = Makefile $(STANDALONE_H)

#-------------------------------------------------------------------------------
# CPU version
#-------------------------------------------------------------------------------
all:  kingAssembler_find_seed kingAssembler_search

#-------------------------------------------------------------------------------
# Codon library
#-------------------------------------------------------------------------------
LIB_CODON_DIR = src/sequence
LIB_CODON = $(LIB_CODON_DIR)/AASequence.o
LIB_CODON += $(LIB_CODON_DIR)/AminoAcid.o
LIB_CODON += $(LIB_CODON_DIR)/CodingSequence.o
LIB_CODON += $(LIB_CODON_DIR)/Codon.o
LIB_CODON += $(LIB_CODON_DIR)/Mutation.o
LIB_CODON += $(LIB_CODON_DIR)/NTSequence.o
LIB_CODON += $(LIB_CODON_DIR)/Nucleotide.o

#-------------------------------------------------------------------------------
# CPU objectives
#-------------------------------------------------------------------------------
%.o: %.cpp %.h $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.C %.h $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# CPU Applications
#-------------------------------------------------------------------------------
kingAssembler_find_seed: fast_kmer_filter.cpp nucl_kmer.o prot_kmer.o city.o $(LIB_CODON) $(DEPS)
	$(CXX) $(CXXFLAGS) fast_kmer_filter.cpp nucl_kmer.o prot_kmer.o city.o $(LIB_CODON) $(LIB) -o kingAssembler_find_seed

kingAssembler_search: search.cpp succinct_dbg.o nucl_kmer.o node_enumerator.o codon.o hmm_graph_search.o city.o $(LIB_CODON) $(DEPS)
	$(CXX) $(CXXFLAGS) search.cpp succinct_dbg.o nucl_kmer.o node_enumerator.o codon.o hmm_graph_search.o city.o $(LIB) $(LIB_CODON) -o kingAssembler_search

# megahit_toolkit: $(TOOLKIT) $(DEPS)
# 	$(CXX) $(CXXFLAGS) $(TOOLKIT) $(LIB) -o megahit_toolkit

.PHONY:
clean:
	-rm -fr *.o \
		$(LIB_CODON) \
		kingAssembler_find_seed kingAssembler_search