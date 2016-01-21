#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

int main(int argc, char **argv) {

    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
        exit(1);
    }

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files

    while (kseq_read(seq) >= 0) {
        // printf("name: %s\n", seq->name.s);
        // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        printf("%s\n", seq->seq.s);
        // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
    }

    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}