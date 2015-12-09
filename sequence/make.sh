#!/bin/bash
g++ -O2 -c translation_test.cpp -o translation.o
g++ -O2 AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o translation.o -o translation
