#include "succinct_dbg.h"

int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Usage %s <graph-prefix>\n", argv[0]);
		exit(2);
	}

	SuccinctDBG dbg;
	// load graph
	dbg.LoadFromMultiFile(argv[1], false); // false means do not load the multiplicity to save memory
	printf("kmer_size: %d, graph_size: %lld\n", dbg.kmer_k, (long long)dbg.size);

	// look look all nodes
	// uint8_t seq[1024];
	// for (int64_t i = 0; i < dbg.size; ++i) {
	// 	if (dbg.IsLast(i)) {
	// 		dbg.Label(i, seq);
	// 		for (int i = 0; i < dbg.kmer_k; ++i) {
	// 			printf("%c", "$ACGT"[seq[i]]);
	// 		}
	// 		puts("");
	// 	}
	// }

	// init dna_map
	int dna_map[256];
	memset(dna_map, -1, sizeof(dna_map));
	for (int i = 0; i < 10; ++i) {
		dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
	}

	// assume k is 21; find "AAAAAAAGAGTGTCTGATAGC"
	char buf[1024];
	while (gets(buf)) {
		bool ok = true;
		uint8_t seq[21];
		for (int i = 0; i < dbg.kmer_k; ++i) {
			if (dna_map[buf[i]] == -1) {
				printf("ACGTN only\n");
				ok = false;
				break;
			}
			seq[i] = dna_map[buf[i]]; // $->0, A->1, C->2, G->3, T->4
		}

		if (!ok) { continue; }

	    int64_t node_id = dbg.IndexBinarySearch(seq);
	    if (node_id == -1) {
	    	printf("No such seq: %s.\n", buf);
	    } else {
	    	int64_t next[4];
	    	int outd = dbg.NextNodes(node_id, next);

	    	printf("ID: %lld, Outdegree: %d\n", (long long)node_id, outd);
	    	for (int i = 0; i < outd; ++i) {
	    		printf("Next #%d: %lld, label: %c\n", i, (long long)next[i], "$ACGT"[dbg.GetNodeLastChar(next[i])]);
	    	}
	    }
	}

    return 0;
}