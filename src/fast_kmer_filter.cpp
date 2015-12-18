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
#include "sequence/NTSequence.h"
#include "sequence/AASequence.h"

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

struct xtimer_t {
    struct timeval tv1, tv2;
    long long time_elapsed;

    xtimer_t() {
        reset();
    }
    void reset() {
        time_elapsed = 0;
    }
    void start() {
        gettimeofday(&tv1, NULL);
    }
    void stop() {
        gettimeofday(&tv2, NULL);
        time_elapsed += (long long)(tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    }
    double elapsed() {
        return time_elapsed / 1000000.0;
    }
};

struct Sequence {
	string name_;
	string comment_;
	string sequence_;

	Sequence(const string &name, const string &comment, const string &sequence)
		:name_(name), comment_(comment), sequence_(sequence) {
	}
};

void ProcessSequenceMulti(const string &sequence, const string &name, const string &comment, HashSetST<ProtKmer> &kmerSet, const int &kmer_size);
void ProcessSequence(const string &sequence, const string &name, const string &comment, HashSetST<ProtKmer> &kmerSet, const int &kmer_size);
char Comp(char c);
void RevComp(string &s);


int find_start(int argc, char **argv) {
	ProtKmer::setUp();
	NuclKmer::setUp();

	if (argc == 1) {
		fprintf(stderr, "Usage: %s <ref_seq> <read_seq> <k_size>\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzopen(argv[1], "r");
	gzFile fp2 = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    kseq_t *seq2 = kseq_init(fp2);

    // int kmer_size = 45;
    int kmer_size = stoi(argv[3]);
    int batch_size = 1 * 1024 * 1024;

    HashSetST<ProtKmer> kmerSet;
    //add kmers from reference to hash set
    while (kseq_read(seq) >= 0) { 
        ProtKmerGenerator kmers = ProtKmerGenerator(seq->seq.s, kmer_size/3, true);
        while (kmers.hasNext()) {
        	ProtKmer temp = kmers.next();
        	temp.model_position = kmers.getPosition();
        	kmerSet.insert(temp);
        }
    }

    int count = 0;
    vector<Sequence> sequence_storage;
	// xtimer_t timer;
	// timer.reset();
    while (kseq_read(seq2) >= 0) {
    	sequence_storage.push_back(Sequence(seq2->name.s, "", seq2->seq.s));
    	if (++count == batch_size) {
    		// timer.start();
    		count = 0;
    		#pragma omp parallel for schedule(dynamic, 1)
	    	for (int i = 0; i < batch_size; i++) {
	    	// vector<ProtKmerGenerator> kmer_gens;
			   	if ((int)sequence_storage[i].sequence_.size() >= kmer_size) {
			   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size);
			   		RevComp(sequence_storage[i].sequence_);
			   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size);
				}
		    }
	    	sequence_storage.clear();
	    	// timer.stop();
    	}
    }

    //do the remaining job
    if (count > 0) {
    	// timer.start();
    	#pragma omp parallel for schedule(dynamic, 1)
    	for (int i = 0; i < count; i++) {
		   	if ((int)sequence_storage[i].sequence_.size() >= kmer_size) {
		   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size);
		   		RevComp(sequence_storage[i].sequence_);
		   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size);
			}
	    }
	    // timer.stop();
    }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}
void ProcessSequenceMulti(const string &sequence, const string &name, const string &comment, HashSetST<ProtKmer> &kmerSet, const int &kmer_size) {
	vector<ProtKmerGenerator> kmer_gens;
	seq::NTSequence nts = seq::NTSequence("", "", sequence);
	for (int i = 0; i < 3; i++) {
	    kmer_gens.push_back(ProtKmerGenerator(seq::AASequence::translate(nts.begin() + i, nts.begin() + i + ((nts.size() - i) / 3) * 3).asString(), kmer_size/3));
	}
	ProtKmer kmer;
    for (int gen = 0; gen < 3; gen++) {
    	while (kmer_gens[gen].hasNext()) {
    		kmer = kmer_gens[gen].next();
    		HashSetST<ProtKmer>::iterator iter = kmerSet.find(kmer);
    		if (iter != NULL) {
    			int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;
    			printf("rplB\tSRR172903.7702200\t357259128\t%s\ttrue\t%d\t%s\t%d\n", sequence.substr(nucl_pos, kmer_size).c_str(), 1, kmer.decodePacked().c_str(), iter->model_position);
                // gene name, read name, 
                // printf("rplB\t%s\t357259128\t%s\ttrue\t%d\t%s\t%d\n", name.c_str(), sequence.substr(nucl_pos, kmer_size).c_str(), 1, kmer.decodePacked().c_str(), iter->model_position);
    		}
    	}
    }	
}

char Comp(char c) {
	switch (c) {
		case 'A':
		case 'a': return 'T';
		case 'C':
		case 'c': return 'G';
		case 'G':
		case 'g': return 'C';
		case 'T':
		case 't': return 'A';
		case 'N':
		case 'n': return 'N';
		default: assert(false);
	}
}

void RevComp(string &s) {
	reverse(s.begin(), s.end());
	for (unsigned i = 0; i < s.length(); ++i) {
		s[i] = Comp(s[i]);
	}
}
