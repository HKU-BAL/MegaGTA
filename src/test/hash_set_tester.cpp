#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "kmer.h"
#include "hash_set.h"

int main(int argc, char **argv) {


    Kmer<> kmer1, kmer2;
    int seq1[150];

    for (int i = 0; i < 150; ++i) {
        seq1[i] = rand() % 4 ;
    }

    for (int i = 0; i < 150; ++i) {
        std::cout << seq1[i] << " ";
    }

    std::cout << std::endl;

    for (int i = 0; i < 45; ++i) {
        kmer1.ShiftAppend(seq1[i], 45);
        kmer2.ShiftAppend(seq1[i + 1], 45);
    }

    for (int i = 0; i < 45; ++i) {
        std::cout << "ACGT"[kmer1.get_base(i)];
    }

    std::cout << std::endl;

    for (int i = 0; i < 45; ++i) {
        std::cout << "ACGT"[kmer2.get_base(i)];
    }

    std::cout << std::endl;

    // std::cout << kmer1.hash() <<std::endl;

    // std::cout << kmer2.hash() <<std::endl;

    HashSet<Kmer<> > kmerSet;
    std::pair<HashSet<Kmer<> >::iterator, bool> result = kmerSet.insert(kmer1);

    std::cout << "result: " << result.second << std::endl;

    HashSet<Kmer<> >::iterator iter1 = kmerSet.find(kmer1);
    HashSet<Kmer<> >::iterator iter2 = kmerSet.find(kmer2);

    if (*iter1 == kmer1) {
        std::cout << "test 1 passed" << std::endl;
    }

    if (iter2 == NULL) {
        std::cout << "test 2 passed" << std::endl;
    }


    return 0;
}