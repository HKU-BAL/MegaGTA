#!/bin/bash
g++ -O2 -lz -std=c++11 -c fast_kmer_filter.cpp -o fast_kmer_filter.o -fopenmp
g++ -O2 -std=c++11 -c city.cpp -o city.o
g++ -std=c++11 -O2 -lz -o fast fast_kmer_filter.o nucl_kmer.cpp prot_kmer.cpp city.o AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o -fopenmp
