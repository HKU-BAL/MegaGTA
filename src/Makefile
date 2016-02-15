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
STANDALONE_H = a_star_node.h hash.h hash_map_st.h \
			   pool_st.h hash_set_st.h functional.h \
			   khash.h kseq.h hash_table_st.h \
			   utils.h hmmer3b_parser.h kmer.h mem_file_checker-inl.h\
			   most_probable_path.h nucl_kmer_generator.h profile_hmm.h \
			   definitions.h prot_kmer_generator.h rank_and_select.h sdbg_multi_io.h \
			   node_enumerator.h

DEPS = Makefile $(STANDALONE_H)

#-------------------------------------------------------------------------------
# CPU version
#-------------------------------------------------------------------------------
# all:  kingAssembler_find_seed kingAssembler_search

#-------------------------------------------------------------------------------
# Codon library
#-------------------------------------------------------------------------------
LIB_CODON_DIR = sequence
LIB_CODON = $(LIB_CODON_DIR)/AASequence.o
LIB_CODON += $(LIB_CODON_DIR)/AminoAcid.o
LIB_CODON += $(LIB_CODON_DIR)/CodingSequence.o
LIB_CODON += $(LIB_CODON_DIR)/Codon.o
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

all: megagta

megagta: megagta.o assembler.o unitig_graph.o assembly_algorithms.o build_graph.o cx1_read2sdbg_s1.o cx1_read2sdbg.h cx1_read2sdbg_s2.o kthread.o \
            build_read_lib.o sequence_manager.o sequence_package.h \
			read_stat.o filter_by_len.o search.o fast_kmer_filter.o translate.o \
			succinct_dbg.o nucl_kmer.o codon.o hmm_graph_search.o city.o prot_kmer.o branch_group.o $(LIB_CODON) \
			options_description.o $(DEP)
	$(CXX) $(CXXFLAGS) megagta.o assembler.o unitig_graph.o assembly_algorithms.o build_graph.o cx1_read2sdbg_s1.o cx1_read2sdbg_s2.o kthread.o sequence_manager.o build_read_lib.o read_stat.o filter_by_len.o search.o fast_kmer_filter.o translate.o options_description.o succinct_dbg.o nucl_kmer.o codon.o hmm_graph_search.o city.o prot_kmer.o branch_group.o $(LIB) $(LIB_CODON) -o megagta

path_viewer: path_viewer.o succinct_dbg.o assembly_algorithms.o
	$(CXX) $(CXXFLAGS)  path_viewer.o succinct_dbg.o assembly_algorithms.o -o path_viewer $(LIB)

false_positive_checker: false_positive_checker.o succinct_dbg.o
	$(CXX) $(CXXFLAGS)  false_positive_checker.o succinct_dbg.o -o false_positive_checker $(LIB)

PHONY:
clean:
	-rm -fr *.o \
		$(LIB_CODON) \
		megagta
