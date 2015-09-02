#!/bin/bash
g++ -std=c++11 -O2 hmm_graph_search_tester.cpp succinct_dbg.cpp nucl_kmer.cpp codon.cpp hmm_graph_search.cpp city.o src/sequence/AASequence.o src/sequence/AminoAcid.o src/sequence/CodingSequence.o src/sequence/Codon.o src/sequence/Mutation.o src/sequence/NTSequence.o src/sequence/Nucleotide.o -o hmm_graph_search_test -fopenmp -static-libstdc++ -mpopcnt
