#include "succinct_dbg.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include "kseq.h"
#include "nucl_kmer_generator.h"
#include "nucl_kmer.h"

#ifndef KSEQ_INITED
#define KSEQ_INITED
KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Usage %s <graph-prefix> <sequences-to-check> <split-to-kmer>\n", argv[0]);
		exit(2);
	}

	SuccinctDBG dbg;
	// load graph
	dbg.LoadFromMultiFile(argv[1], false); // false means do not load the multiplicity to save memory
	// printf("kmer_size: %d, graph_size: %lld\n", dbg.kmer_k + 1, (long long)dbg.size);

	NuclKmer::setUp();

	gzFile fp = gzopen(argv[2], "r");
    kseq_t *seq = kseq_init(fp); // kseq to read files
	int kmer_size = dbg.kmer_k + 1;

	ofstream myfile;
	myfile.open (argv[3]);

	while (kseq_read(seq) >= 0) {
        NuclKmerGenerator kmers = NuclKmerGenerator(seq->seq.s, kmer_size, false);

		myfile << ">" << seq->name.s << "\n";
		  
        while (kmers.hasNext()) {
        	NuclKmer temp = kmers.next();
        	myfile << temp.decodePacked() << "\n";
        }
    }

    myfile.close();

	ifstream sequence_file (argv[3]);

	// init dna_map
	int dna_map[256];
	memset(dna_map, -1, sizeof(dna_map));
	for (int i = 0; i < 10; ++i) {
		dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
	}

	if (sequence_file.is_open())
	{
		string line;
		getline (sequence_file, line);
		//the first line is the name

		//header
		while ( getline (sequence_file,line) ) { 

			bool ok = true;
			uint8_t seq[dbg.kmer_k + 1];
			for (int i = 0; i < dbg.kmer_k + 1; ++i) {
				if (dna_map[line[i]] == -1) {
					// printf("ACGTN only\n");
					ok = false;
					break;
				}
				seq[i] = dna_map[line[i]]; // $->0, A->1, C->2, G->3, T->4
			}

			if (!ok) { continue; }

		    int64_t node_id = dbg.IndexBinarySearchEdge(seq);
		    if (node_id == -1) {
		    	printf("No such seq: %s.\n", line.c_str());
		    } else {
		    	int64_t next[4];
		    	int outd = dbg.OutgoingEdges(node_id, next);

		    	printf("Kmer: %s\tOutdegree:\t%d\n", line.c_str(), outd);
		    	// printf("Kmer: %s\n", line.c_str());
		    	// for (int i = 0; i < outd; ++i) {
		    	// 	printf("Next #%d: label: %c\n", i, "$acgt"[dbg.GetEdgeOutLabel(next[i])]);
		    	// }
		    }
		}
	}


	

    return 0;
}