#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "kmer.h"
#include "city.h"

int main(int argc, char **argv) {


	Kmer<4, uint64_t> kmer1;
	int seq1[150];
	for (int i=0; i < 150; ++i){
		seq1[i] = rand() % 4 ;
	}

	for (int i=0; i < 150; ++i){
		std::cout << seq1[i] << " ";
	}

	std::cout << std::endl;

	for (int i=0; i < 45; ++i){
		kmer1.ShiftAppend(seq1[i],45);
	}

	for (int i=0; i < 45; ++i){
		std::cout << "0123"[kmer1.get_base(i)] << " ";
	}

	std::cout << std::endl;

	std::cout << kmer1.hash() <<std::endl;


	
	return 0;
}