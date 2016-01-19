#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
// remove this part to refactor
#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif
// remove this part to refactor

using namespace std;

int main(int argc, char **argv) {
    vector<string> ref;

    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
        exit(1);
    }

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files

    while (kseq_read(seq) >= 0) {
        // printf("%s\n", seq->seq.s);
        ref.push_back(seq->seq.s);
    }

    // printf("FUCKFUCKFUCK\n");

    // cout << ref[ref.size()-1] << endl;

    kseq_destroy(seq);
    gzclose(fp);
}

