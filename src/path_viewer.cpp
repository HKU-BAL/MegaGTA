#include "succinct_dbg.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage %s <graph-prefix> <sequences-to-check>\n", argv[0]);
        exit(2);
    }

    SuccinctDBG dbg;
    // load graph
    dbg.LoadFromMultiFile(argv[1], true); // false means do not load the multiplicity to save memory
    printf("kmer_size: %d, graph_size: %lld\n", dbg.kmer_k, (long long)dbg.size);

    gzFile fp = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int kmer_size = dbg.kmer_k + 1;

    // ofstream myfile;
    // myfile.open (argv[3]);
    int dna_map[256];
    memset(dna_map, -1, sizeof(dna_map));

    for (int i = 0; i < 10; ++i) {
        dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
    }

    while (kseq_read(seq) >= 0) {
        for (int i = 0; i < seq->seq.l; ++i) {
            seq->seq.s[i] = dna_map[seq->seq.s[i]];
        }

        printf("%s:\n", seq->name.s);

        for (int i = 0; i + dbg.kmer_k + 1 < seq->seq.l; ++i) {
            printf("%d:", i);
            int64_t node_id = dbg.IndexBinarySearchEdge((uint8_t *)seq->seq.s + i);

            if (node_id == -1) {
                printf(" not found\n");
            }
            else {
                printf(" %d", dbg.EdgeMultiplicity(node_id));
                int64_t next[4];
                int outd = dbg.OutgoingEdges(node_id, next);
                printf(" %d children", outd);

                for (int j = 0; j < outd; ++j) {
                    printf(" %d", dbg.EdgeMultiplicity(next[j]));
                }

                printf("\n");
            }

        }
    }

    return 0;
}
