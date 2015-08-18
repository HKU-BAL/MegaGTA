#!/bin/bash
g++ -std=c++11 hmm_graph_search_tester.cpp succinct_dbg.cpp nucl_kmer.cpp node_enumerator.cpp codon.cpp hmm_graph_search.cpp city.o AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o -o hmm_graph_search_test -fopenmp -static-libstdc++
