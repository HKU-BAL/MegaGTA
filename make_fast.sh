#!/bin/bash
g++ -O2 -lz -std=c++11 -c fast_kmer_filter.cpp -o fast_kmer_filter.o -fopenmp
g++ -O2 -c city.cpp -o city.o
g++ -O2 -lz -o fast fast_kmer_filter.o city.o AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o -fopenmp
