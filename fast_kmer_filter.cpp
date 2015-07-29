#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
#include "prot_kmer_generator.h"
#include "hash_set.h"
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include "src/sequence/NTSequence.h"
#include "src/sequence/AASequence.h"

//test
#include <bitset>

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

struct Sequence {
	string name_;
	string comment_;
	string sequence_;

	Sequence(const string &name, const string &comment, const string &sequence) {
		name_ = name;
		comment_ = comment;
		sequence_ = sequence;
	}
};

struct KmerHelper {
	ProtKmer kmer_;
	string nucl_seq_;
	int frame_;
	int position_;

	KmerHelper() {}

	KmerHelper(const ProtKmer &kmer, const string &nucl_seq, const int &frame, const int &position) {
		kmer_ = kmer;
		nucl_seq_ = nucl_seq;
		frame_ = frame;
		position_ = position;
	}

	uint64_t hash() const {
		return kmer_.hash();
	}

	bool operator ==(const KmerHelper &kmer_helper) const {
		if (kmer_.kmers[0] != kmer_helper.kmer_.kmers[0] || kmer_.kmers[1] != kmer_helper.kmer_.kmers[1]) {
			return false;
		}
		return true;
	}
};

int main(int argc, char **argv) {

	if (argc == 1) {
		fprintf(stderr, "Usage: %s <ref_seq> <read_seq>\n", argv[0]);
		exit(1);
	}

	gzFile fp = gzopen(argv[1], "r");
	gzFile fp2 = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
    kseq_t *seq2 = kseq_init(fp2);

    int kmer_size = 45;
    int batch_size = 500000;

    HashSet<ProtKmer> kmerSet;

    while (kseq_read(seq) >= 0) { 
        printf("%s\n", seq->seq.s);
        string string_seq(seq->seq.s);
        ProtKmerGenerator kmers = ProtKmerGenerator(string_seq, kmer_size/3, true); // kmer = 45
        while (kmers.hasNext()) {
        	ProtKmer temp = kmers.next();
        	// cout << "seeds = " << temp.decodePacked();
        	std::pair<HashSet<ProtKmer>::iterator, bool> result = kmerSet.insert(temp);
        	// cout << " result: " << result.second << endl;
        }
    }

    while (kseq_read(seq2) >= 0) {
    	// printf("%s\n", seq2->seq.s); 	
    	string string_seq(seq2->seq.s);
    	string string_name(seq2->name.s);
    	string string_comment;
    	if (seq2->comment.l) 
    		string_comment = string(seq2->comment.s); 
    	else 
    		string_comment = "";

    	vector<ProtKmerGenerator> kmer_gens;
    	// ProtKmerGenerator kmer_gens[3];

    	if (string_seq.size() >= kmer_size) {

	    	for (int i = 0; i < 1; i++) {
	    		string seq = string_seq.substr(i);
	    		seq::NTSequence nts = seq::NTSequence(string_name, string_comment, seq);
	    		seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
	    		// cout << ">" << string_name << endl;
	    		// cout << aa.asString() << endl;
	    		kmer_gens.push_back(ProtKmerGenerator(aa.asString(), kmer_size/3));
	    	}

	    	ProtKmer kmer;
	    	// Kmer kmer;
		    for (int gen = 0; gen < kmer_gens.size(); gen++) {
		    	while (kmer_gens[gen].hasNext()) {
		    		kmer = kmer_gens[gen].next();
		    		HashSet<ProtKmer>::iterator iter = kmerSet.find(kmer);
		    		if (iter != NULL) {
		    			int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;
		    			// cout << "rplB\t" << "SRR172903.7702200\t" << "357259128\t";
		    			// printf("%s\ttrue\t%d\t%s\t%d\n", string_seq.substr(nucl_pos, kmer_size).c_str(), gen+1, kmer.decodePacked().c_str(), kmer_gens[gen].getPosition());
		    		}
		    	}
		    }
		}
    }


    	//multi-thread version
    // int count = 0;
    // vector<Sequence> sequence_storage;
    // HashSet<KmerHelper> starting_kmers;
    // while (int ret = kseq_read(seq2) >= 0) {
    // 	sequence_storage.push_back(Sequence(string(seq2->name.s), string(seq2->comment.s), string(seq2->seq.s)));
    // 	if (++count == batch_size || ret < 0) { //100000
    // 		cout << "count = "<< count << "\n" ;
    // 		if (ret >= 0)
    // 			count = 0;
    // 		else
    // 			batch_size = count;
    // 		#pragma omp parallel for
	   //  	for (int i = 0; i < batch_size; i++) {
	   //  		string string_seq = sequence_storage[i].sequence_;
	   //  		string string_name = sequence_storage[i].name_;
	   //  		string string_comment = sequence_storage[i].comment_;
	   //  		vector<ProtKmerGenerator> kmer_gens;
		  //   	if (string_seq.size() >= kmer_size) {
			 //    	for (int i = 0; i < 3; i++) {
			 //    		string seq = string_seq.substr(i);
			 //    		seq::NTSequence nts = seq::NTSequence(string_name, string_comment, seq);
			 //    		seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
			 //    		kmer_gens.push_back(ProtKmerGenerator(aa.asString(), kmer_size/3));
			 //    	}

			 //    	ProtKmer kmer;
				//     for (int gen = 0; gen < kmer_gens.size(); gen++) {
				//     	while (kmer_gens[gen].hasNext()) {
				//     		kmer = kmer_gens[gen].next();
				//     		HashSet<ProtKmer>::iterator iter = kmerSet.find(kmer);
				//     		if (iter != NULL) {
				//     			int nucl_pos = (kmer_gens[gen].getPosition() - 1) * 3 + gen;
				//     			starting_kmers.insert(KmerHelper(kmer, string_seq.substr(nucl_pos, kmer_size), gen+1, kmer_gens[gen].getPosition()));
				//     		}
				//     	}
				//     }
				// }
	   //  	}	    	
    // 	}    	
    // }

    // for (HashSet<KmerHelper>::iterator i = starting_kmers.begin(); i != starting_kmers.end() ; i++) {
    // 	cout << "rplB\t" << "SRR172903.7702200\t" << "357259128\t";
    // 	printf("%s\ttrue\t%d\t%s\t%d\n", i->nucl_seq_.c_str(), i->frame_, i->kmer_.decodePacked().c_str(), i->position_);
    // }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}
