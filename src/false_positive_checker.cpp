#include "succinct_dbg.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
    fprintf(stderr, "kmer_size: %d, graph_size: %lld\n", dbg.kmer_k, (long long)dbg.size);

    gzFile fp = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files

    int dna_map[256];
    memset(dna_map, -1, sizeof(dna_map));

    for (int i = 0; i < 10; ++i) {
        dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
    }

    vector<vector<uint8_t> > seqs;
    vector<string> name;

    while (kseq_read(seq) >= 0) {
        name.push_back(seq->name.s);
        seqs.push_back(vector<uint8_t>());
        for (size_t i = 0; i < seq->seq.l; ++i) {
            seqs.back().push_back(dna_map[seq->seq.s[i]]);
        }
    }

#pragma omp parallel for
    for (size_t seq_id = 0; seq_id < seqs.size(); ++seq_id) {
        for (size_t i = 0; i + dbg.kmer_k + 1 < seqs[seq_id].size(); ++i) {
            int64_t node_id = dbg.IndexBinarySearchEdge(&seqs[seq_id][i]);

            if (node_id == -1) {
                printf("%s %zu %zu %zu\n", name[seq_id].c_str(), i, seqs[seq_id].size(), min(i, seqs[seq_id].size() - dbg.kmer_k - i));
            }

        }
    }

    return 0;
}
