#include "succinct_dbg.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include "kseq.h"
#include "nucl_kmer_generator.h"
#include "nucl_kmer.h"
#include <vector>

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
    dbg.LoadFromMultiFile(argv[1], false); // false means do not load the multiplicity to save memory
    // printf("kmer_size: %d, graph_size: %lld\n", dbg.kmer_k + 1, (long long)dbg.size);

    NuclKmer::setUp();

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

    //container for kmers
    vector<string> query_sequence_container;
    vector<string> query_name_container;

    while (kseq_read(seq) >= 0) {
        query_sequence_container.push_back(seq->seq.s);
        query_name_container.push_back(seq->name.s);
    }

    #pragma omp parallel for schedule(dynamic, 1)

    for (int i = 0; i < query_sequence_container.size(); ++i) {
        NuclKmerGenerator kmers = NuclKmerGenerator(query_sequence_container[i], kmer_size, false);
        int kmer_count = 0, false_positive_kmer_count = 0;

        while (kmers.hasNext()) {
            kmer_count++;
            NuclKmer temp = kmers.next();
            string line = temp.decodePacked();
            // myfile << temp.decodePacked() << "\n";
            bool ok = true;
            uint8_t seq[dbg.kmer_k + 1];

            for (int j = 0; j < dbg.kmer_k + 1; ++j) {
                if (dna_map[line[j]] == -1) {
                    // printf("ACGTN only\n");
                    ok = false;
                    break;
                }

                seq[j] = dna_map[line[j]]; // $->0, A->1, C->2, G->3, T->4
            }

            if (!ok) {
                continue;
            }

            int64_t node_id = dbg.IndexBinarySearchEdge(seq);

            if (node_id == -1) {
                false_positive_kmer_count++;
                // printf("No such seq: %s.\n", line.c_str());
            }
            else {
                // int64_t next[4];
                // int outd = dbg.OutgoingEdges(node_id, next);

                // printf("Kmer: %s\tOutdegree:\t%d\n", line.c_str(), outd);
                // printf("Kmer: %s\n", line.c_str());
                // for (int i = 0; i < outd; ++i) {
                // 	printf("Next #%d: label: %c\n", i, "$acgt"[dbg.GetEdgeOutLabel(next[i])]);
                // }
            }
        }

        printf("%s\tfalse positive kmers: %d\ttotal number of kmers: %d\n", query_name_container[i].c_str(), false_positive_kmer_count, kmer_count);
    }

    // while (kseq_read(seq) >= 0) {
    //        NuclKmerGenerator kmers = NuclKmerGenerator(seq->seq.s, kmer_size, false);

    // 	myfile << ">" << seq->name.s << "\n";

    //        while (kmers.hasNext()) {
    //        	NuclKmer temp = kmers.next();
    //        	myfile << temp.decodePacked() << "\n";
    //        }
    //    }

    //refactor it to use multiple threads


    //    myfile.close();

    // ifstream sequence_file (argv[3]);

    // init dna_map
    // int dna_map[256];
    // memset(dna_map, -1, sizeof(dna_map));
    // for (int i = 0; i < 10; ++i) {
    // 	dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
    // }

    // if (sequence_file.is_open())
    // {
    // 	string line;
    // 	getline (sequence_file, line);
    // 	//the first line is the name

    // 	//header
    // 	while ( getline (sequence_file,line) ) {

    // 		bool ok = true;
    // 		uint8_t seq[dbg.kmer_k + 1];
    // 		for (int i = 0; i < dbg.kmer_k + 1; ++i) {
    // 			if (dna_map[line[i]] == -1) {
    // 				// printf("ACGTN only\n");
    // 				ok = false;
    // 				break;
    // 			}
    // 			seq[i] = dna_map[line[i]]; // $->0, A->1, C->2, G->3, T->4
    // 		}

    // 		if (!ok) { continue; }

    // 	    int64_t node_id = dbg.IndexBinarySearchEdge(seq);
    // 	    if (node_id == -1) {
    // 	    	printf("No such seq: %s.\n", line.c_str());
    // 	    } else {
    // 	    	int64_t next[4];
    // 	    	int outd = dbg.OutgoingEdges(node_id, next);

    // 	    	printf("Kmer: %s\tOutdegree:\t%d\n", line.c_str(), outd);
    // 	    	// printf("Kmer: %s\n", line.c_str());
    // 	    	// for (int i = 0; i < outd; ++i) {
    // 	    	// 	printf("Next #%d: label: %c\n", i, "$acgt"[dbg.GetEdgeOutLabel(next[i])]);
    // 	    	// }
    // 	    }
    // 	}
    // }




    return 0;
}