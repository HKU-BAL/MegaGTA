#!/bin/bash
g++ -std=c++11 -lz -O2 search.cpp succinct_dbg.cpp nucl_kmer.cpp node_enumerator.cpp codon.cpp hmm_graph_search.cpp city.o src/sequence/AASequence.o src/sequence/AminoAcid.o src/sequence/CodingSequence.o src/sequence/Codon.o src/sequence/Mutation.o src/sequence/NTSequence.o src/sequence/Nucleotide.o -o search  -static-libstdc++ -fopenmp

