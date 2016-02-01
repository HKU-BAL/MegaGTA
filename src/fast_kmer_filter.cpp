#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
#include "prot_kmer_generator.h"
#include "nucl_kmer.h"
#include "hash_set_st.h"
#include <string.h>
#include <string>
#include <vector>
#include <omp.h>
#include "sequence/NTSequence.h"
#include "sequence/AASequence.h"
#include "sequence_manager.h"

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

struct Seed {
    string nucl;
    string prot;
    int model_pos;

    Seed(const string &nucl = "", const string &prot = "", int model_pos = 0):
        nucl(nucl), prot(prot), model_pos(model_pos) {}

    bool operator == (const Seed &rhs) const {
        return rhs.nucl == nucl;
    }

    bool operator < (const Seed &rhs) const {
        return nucl < rhs.nucl;
    }
};

void ProcessSequenceMulti(const string &sequence, HashSetST<ProtKmer> &kmerSet, const int &kmer_size, vector<Seed> &candidates);

int find_start(int argc, char **argv) {
    ProtKmer::setUp();
    NuclKmer::setUp();

    if (argc == 1) {
        fprintf(stderr, "Usage: %s <ref_seq> <read.lib> <k_size> [num_threads=0]\n", argv[0]);
        exit(1);
    }

    int num_threads = 0;

    if (argc > 4) {
        num_threads = atoi(argv[4]);
    }

    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }

    omp_set_num_threads(num_threads);

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    int kmer_size = stoi(argv[3]);

    HashSetST<ProtKmer> kmerSet;

    //add kmers from reference to hash set
    while (kseq_read(seq) >= 0) {
        ProtKmerGenerator kmers = ProtKmerGenerator(seq->seq.s, kmer_size / 3, true);

        while (kmers.hasNext()) {
            ProtKmer temp = kmers.next();
            temp.model_position = kmers.getPosition();
            kmerSet.insert(temp);
        }
    }

    xlog("reference kmer set size: %lld\n", kmerSet.size());

    int count = 0;

    // read binary reads
    int64_t kMaxNumReads = 1 << 22;
    int64_t kMaxNumBases = 1 << 28;
    bool append = false;
    bool reverse = false;

    SequenceManager seq_manager;
    SequencePackage package;
    seq_manager.set_file_type(SequenceManager::kBinaryReads);

    seq_manager.set_file(argv[2]);
    seq_manager.set_readlib_type(SequenceManager::kSingle); // PE info not used
    seq_manager.set_package(&package);

    vector<vector<Seed> > seeds(num_threads);

    while ((count = seq_manager.ReadShortReads(kMaxNumReads, kMaxNumBases, append, reverse)) > 0) {
        xlog("Processing %d reads\n", count);
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < count; i++) {
                int len = package.length(i);
                if (len >= kmer_size) {
                    string s;
                    s.resize(len);

                    for (int j = 0; j < len; ++j) {
                        s[j] = "ACGT"[package.get_base(i, j)];
                    }
                    ProcessSequenceMulti(s, kmerSet, kmer_size, seeds[omp_get_thread_num()]);

                    for (int j = 0; j < len; ++j) {
                        s[j] = "ACGT"[3 - package.get_base(i, len - 1 - j)];
                    }

                    ProcessSequenceMulti(s, kmerSet, kmer_size, seeds[omp_get_thread_num()]);
                }
            }
    }

    if (argc > 5) {
        seq_manager.clear();
        seq_manager.set_file_type(SequenceManager::kFastxReads);
        seq_manager.set_file(argv[5]);
        seq_manager.set_readlib_type(SequenceManager::kSingle); // PE info not used
        seq_manager.set_package(&package);

        while ((count = seq_manager.ReadShortReads(kMaxNumReads, kMaxNumBases, append, reverse)) > 0) {
            xlog("Processing %d contigs\n", count);
    #pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < count; i++) {
                int len = package.length(i);
                if (len >= kmer_size) {
                    string s;
                    s.resize(len);

                    for (int j = 0; j < len; ++j) {
                        s[j] = "ACGT"[package.get_base(i, j)];
                    }
                    ProcessSequenceMulti(s, kmerSet, kmer_size, seeds[omp_get_thread_num()]);

                    for (int j = 0; j < len; ++j) {
                        s[j] = "ACGT"[3 - package.get_base(i, len - 1 - j)];
                    }

                    ProcessSequenceMulti(s, kmerSet, kmer_size, seeds[omp_get_thread_num()]);
                }
            }
        }
    }

    size_t total_size = 0;

    for (int i = 0; i < num_threads; ++i) {
        total_size += seeds[i].size();
    }

    seeds[0].reserve(total_size);

    for (int i = 1; i < num_threads; ++i) {
        seeds[0].insert(seeds[0].end(), seeds[i].begin(), seeds[i].end());
    }

    sort(seeds[0].begin(), seeds[0].end());
    total_size = unique(seeds[0].begin(), seeds[0].end()) - seeds[0].begin();
    random_shuffle(seeds[0].begin(), seeds[0].begin() + total_size);


    for (size_t i = 0; i < total_size; ++i) {
        printf("dump_gene_name\tdump_seq_name\tdump\t%s\ttrue\t%d\t%s\t%d\n", seeds[0][i].nucl.c_str(), 1, seeds[0][i].prot.c_str(), seeds[0][i].model_pos);
    }

    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

void ProcessSequenceMulti(const string &sequence, HashSetST<ProtKmer> &kmerSet, const int &kmer_size, vector<Seed> &candidates) {
    vector<ProtKmerGenerator> kmer_gens;
    seq::NTSequence nts = seq::NTSequence("", "", sequence);

    for (int i = 0; i < 3; i++) {
        kmer_gens.push_back(ProtKmerGenerator(seq::AASequence::translate(nts.begin() + i, nts.begin() + i + ((nts.size() - i) / 3) * 3).asString(), kmer_size / 3));
    }

    ProtKmer kmer;

    for (int gen = 0; gen < 3; gen++) {
        while (kmer_gens[gen].hasNext()) {
            kmer = kmer_gens[gen].next();
            HashSetST<ProtKmer>::iterator iter = kmerSet.find(kmer);

            if (iter != NULL) {
                int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;
                candidates.push_back(Seed(sequence.substr(nucl_pos, kmer_size), kmer.decodePacked(), iter->model_position));
                // printf("dump_gene_name\tdump_seq_name\tdump\t%s\ttrue\t%d\t%s\t%d\n", sequence.substr(nucl_pos, kmer_size).c_str(), 1, kmer.decodePacked().c_str(), iter->model_position);
            }
        }
    }
}
