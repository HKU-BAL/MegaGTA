#!/bin/bash
g++ -c translation_test.cpp -o translation.o
g++ AASequence.o AminoAcid.o CodingSequence.o Codon.o Mutation.o NTSequence.o Nucleotide.o translation.o -o translation
