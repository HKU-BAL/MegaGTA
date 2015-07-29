#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
#include "prot_kmer_generator.h"
#include "nucl_kmer.h"
#include "hash_set.h"
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include "src/sequence/NTSequence.h"
#include "src/sequence/AASequence.h"

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

//test
#include <bitset>

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

struct KmerHelper {
	ProtKmer kmer_;
	// string nucl_seq_;
	NuclKmer n_kmer_;
	int frame_;
	int position_;

	KmerHelper() {}

	KmerHelper(const ProtKmer &kmer, const string &nucl_seq, const int &frame, const int &position) {
		kmer_ = kmer;
		// nucl_seq_ = nucl_seq;
		n_kmer_ = NuclKmer(nucl_seq);
		frame_ = frame;
		position_ = position;
	}

	uint64_t hash() const {
		return n_kmer_.hash();
		// NuclKmer dna_kmer = NuclKmer(nucl_seq_);
		// return dna_kmer.hash();
		// return kmer_.hash();
	}

	bool operator ==(const KmerHelper &kmer_helper) const {
		if (kmer_.kmers[0] != kmer_helper.kmer_.kmers[0] || kmer_.kmers[1] != kmer_helper.kmer_.kmers[1]) {
			return false;
		}
		return true;
	}
};

void ProcessSequenceMulti(const string &sequence, const string &name, const string &comment, HashSet<ProtKmer> &kmerSet, const int &kmer_size, HashSet<KmerHelper> &starting_kmers);
void ProcessSequence(const string &sequence, const string &name, const string &comment, HashSet<ProtKmer> &kmerSet, const int &kmer_size);
char Comp(char c);
void RevComp(string &s);


int main(int argc, char **argv) {
	ProtKmer::setUp();
	NuclKmer::setUp();

	if (argc == 1) {
		fprintf(stderr, "Usage: %s <ref_seq> <read_seq>\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzopen(argv[1], "r");
	gzFile fp2 = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    kseq_t *seq2 = kseq_init(fp2);

    int kmer_size = 45;
    int batch_size = 1 * 1024 * 1024;

    HashSet<ProtKmer> kmerSet;

    while (kseq_read(seq) >= 0) { 
        // printf("%s\n", seq->seq.s);
        ProtKmerGenerator kmers = ProtKmerGenerator(seq->seq.s, kmer_size/3, true); // kmer = 45
        while (kmers.hasNext()) {
        	ProtKmer temp = kmers.next();
        	// cout << "seeds = " << temp.decodePacked();
        	kmerSet.insert(temp);
        	// cout << " result: " << result.second << endl;
        }
    }

  //   while (kseq_read(seq2) >= 0) {
  //   	// printf("%s\n", seq2->seq.s); 	
  //   	string string_seq(seq2->seq.s);
  //   	string string_name(seq2->name.s);
  //   	string string_comment;
  //   	if (seq2->comment.l) 
  //   		string_comment = string(seq2->comment.s); 
  //   	else 
  //   		string_comment = "";

  //   	string rc_string_seq = RevComp(string_seq);

  //   	if (string_seq.size() >= kmer_size) {
  //   		ProcessSequence(string_seq, string_name, string_comment, kmerSet, kmer_size);
  //   		ProcessSequence(rc_string_seq, string_name, string_comment, kmerSet, kmer_size);
		// }
  //   }
  //   return 0;

    	//multi-thread version
    int count = 0;
    vector<Sequence> sequence_storage;
    HashSet<KmerHelper> starting_kmers;
	// xtimer_t timer;
	// timer.reset();
    while (kseq_read(seq2) >= 0) {
    	sequence_storage.push_back(Sequence(seq2->name.s, seq2->comment.s, seq2->seq.s));
    	if (++count == batch_size) {
    		// timer.start();
    		count = 0;
    		#pragma omp parallel for schedule(dynamic, 1)
	    	for (int i = 0; i < batch_size; i++) {
	    	// vector<ProtKmerGenerator> kmer_gens;
			   	if (sequence_storage[i].sequence_.size() >= kmer_size) {
			   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size, starting_kmers);
			   		RevComp(sequence_storage[i].sequence_);
			   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size, starting_kmers);
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
	    	// vector<ProtKmerGenerator> kmer_gens;
		   	if (sequence_storage[i].sequence_.size() >= kmer_size) {
		   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size, starting_kmers);
		   		RevComp(sequence_storage[i].sequence_);
		   		ProcessSequenceMulti(sequence_storage[i].sequence_, sequence_storage[i].name_,sequence_storage[i].comment_, kmerSet, kmer_size, starting_kmers);
			}
	    }
	    // timer.stop();
    }

    // cerr << "multi-process time: " << timer.elapsed() << endl;

    for (HashSet<KmerHelper>::iterator i = starting_kmers.begin(); i != starting_kmers.end() ; i++) {
    	cout << "rplB\t" << "SRR172903.7702200\t" << "357259128\t";
    	printf("%s\ttrue\t%d\t%s\t%d\n", /*i->nucl_seq_.c_str()*/ "haha", i->frame_, i->kmer_.decodePacked().c_str(), i->position_);
    }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}
void ProcessSequenceMulti(const string &sequence, const string &name, const string &comment, HashSet<ProtKmer> &kmerSet, const int &kmer_size, HashSet<KmerHelper> &starting_kmers) {
	vector<ProtKmerGenerator> kmer_gens;
	seq::NTSequence nts = seq::NTSequence("", "", sequence);
	for (int i = 0; i < 3; i++) {
	    kmer_gens.push_back(ProtKmerGenerator(seq::AASequence::translate(nts.begin() + i, nts.begin() + i + ((nts.size() - i) / 3) * 3).asString(), kmer_size/3));
	}
	ProtKmer kmer;
    for (int gen = 0; gen < 3; gen++) {
    	while (kmer_gens[gen].hasNext()) {
    		kmer = kmer_gens[gen].next();
    		HashSet<ProtKmer>::iterator iter = kmerSet.find(kmer);
    		if (iter != NULL) {
    			// cout << kmer.decodePacked() << endl;
    			int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;

    			printf("rplB\tSRR172903.7702200\t357259128\t%s\ttrue\t%d\t%s\t%d\n", sequence.substr(nucl_pos, kmer_size).c_str(), gen + 1, kmer.decodePacked().c_str(), kmer_gens[gen].getPosition());
    			// starting_kmers.insert(KmerHelper(kmer, sequence.substr(nucl_pos, kmer_size), gen+1, kmer_gens[gen].getPosition()));
    		}
    	}
    }	
}

void ProcessSequence(const string &sequence, const string &name, const string &comment, HashSet<ProtKmer> &kmerSet, const int &kmer_size) {
	vector<ProtKmerGenerator> kmer_gens;
	for (int i = 0; i < 3; i++) {
	    string seq = sequence.substr(i);
	    seq::NTSequence nts = seq::NTSequence(name, comment, seq);
	    seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
	    // cout << ">" << string_name << endl;
	    // cout << aa.asString() << endl;

	    kmer_gens.push_back(ProtKmerGenerator(aa.asString(), kmer_size/3));
	}

	ProtKmer kmer;
	for (int gen = 0; gen < kmer_gens.size(); gen++) {
	  	while (kmer_gens[gen].hasNext()) {
	   		kmer = kmer_gens[gen].next();

	   		// cout << kmer.decodePacked() << endl;

	   		HashSet<ProtKmer>::iterator iter = kmerSet.find(kmer);
	   		if (iter != NULL) {
	   			int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;
	   			cout << "rplB\t" << "SRR172903.7702200\t" << "357259128\t";
	   			printf("%s\ttrue\t%d\t%s\t%d\n", sequence.substr(nucl_pos, kmer_size).c_str(), gen+1, kmer.decodePacked().c_str(), kmer_gens[gen].getPosition());
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
