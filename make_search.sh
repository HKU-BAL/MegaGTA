#!/bin/bash
g++ -std=c++11 -lz -O2 search.cpp succinct_dbg.cpp nucl_kmer.cpp node_enumerator.cpp codon.cpp hmm_graph_search.cpp city.o AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o -o search  -static-libstdc++ -fopenmp

