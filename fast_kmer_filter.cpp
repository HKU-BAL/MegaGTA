#include <stdio.h>
#include <stdlib.h>
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

    HashSet<ProtKmer> kmerSet;

    while (kseq_read(seq) >= 0) { 
        // printf("%s\n", seq->seq.s);
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
		    			cout << "rplB\t" << "SRR172903.7702200\t" << "357259128\t";
		    			printf("%s\ttrue\t%d\t%s\t%d\n", string_seq.substr(nucl_pos, kmer_size).c_str(), gen+1, kmer.decodePacked().c_str(), kmer_gens[gen].getPosition());
		    		}
		    	}
		    }
		}
    }

    kseq_destroy(seq);
    gzclose(fp);
	return 0;
}
