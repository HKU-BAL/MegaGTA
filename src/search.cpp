#include "hmm_graph_search.h"
#include "node_enumerator.h"
#include "profile_hmm.h"
#include "most_probable_path.h"
#include "hmmer3b_parser.h"
#include "succinct_dbg.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <omp.h>


using namespace std;

int search(int argc, char **argv) {

	if (argc < 5) {
		fprintf(stderr, "Usage: %s <succinct_dbg> <gene_list> <starting_kmers> <output_info> [num_threads=0]\n", argv[0]);
		exit(1);
	}

	int num_threads = 0;

	if (argc > 5) {
		num_threads = atoi(argv[5]);
	}
	if (num_threads == 0) {
		num_threads = omp_get_max_threads();
	}
	omp_set_num_threads(num_threads);

	int heuristic_pruning = 20; //this one should be able to adapt to user preference
	HMMGraphSearch::setUp();
	SuccinctDBG dbg;
	dbg.LoadFromMultiFile(argv[1], false);

	// -----refactor it to a list base read/write function
	ifstream gene_list_file (argv[2]);
	vector<string> gene_info;
	vector<vector<string>> gene_list;
	if (gene_list_file.is_open()) {
		string gene;
		while ( getline (gene_list_file, gene) ) {
			istringstream iss(gene);
			string gene_name, forward_hmm_path, reverse_hmm_path;
			iss >> gene_name >> forward_hmm_path >> reverse_hmm_path;
			gene_info.clear();
			gene_info.push_back(gene_name);
			gene_info.push_back(forward_hmm_path);
			gene_info.push_back(reverse_hmm_path);
			gene_list.push_back(gene_info);
		}
	}

	for (vector<string>& gene : gene_list) {
        fprintf(stderr, "START %s\n", gene[0].c_str());
		ifstream hmm_file (gene[1]);
		ProfileHMM forward_hmm = ProfileHMM(true);
		Parser::readHMM(hmm_file, forward_hmm);
		ifstream hmm_file_2 (gene[2]);
		ProfileHMM reverse_hmm = ProfileHMM(true);
		Parser::readHMM(hmm_file_2, reverse_hmm);
		MostProbablePath for_hcost = MostProbablePath(forward_hmm);
		MostProbablePath rev_hcost = MostProbablePath(reverse_hmm);

		pair<string, int> starting_kmer;
		vector<pair<string, int>> starting_kmer_storage;
		string sk = string(argv[3]) + "_" + gene[0] + "_starting_kmers.txt";
		ifstream starting_kmer_file (sk);

		if (starting_kmer_file.is_open()) {
			string line, line_array[8];
			while ( getline (starting_kmer_file, line) ) {
				istringstream iss(line);
				for (int i = 0; i < 8; ++i) {
					iss >> line_array[i];
				}
				transform(line_array[3].begin(), line_array[3].end(), line_array[3].begin(), ::tolower);
				starting_kmer = make_pair(line_array[3], stoi(line_array[7]) -1 );
				starting_kmer_storage.push_back(starting_kmer);
			}
		} else {
			exit(1);
		}
		FILE * out_file;
		string out_file_name = string(argv[4]) + "_raw_contigs_" + gene[0] + ".fasta";
		out_file = fopen(out_file_name.c_str(), "w");
		vector<HMMGraphSearch> search;
		vector<NodeEnumerator> for_node_enumerator, rev_node_enumerator;
		for (int i = 0; i < num_threads; ++i) {
			search.push_back(HMMGraphSearch(heuristic_pruning));
			for_node_enumerator.push_back(NodeEnumerator(forward_hmm, for_hcost));
			rev_node_enumerator.push_back(NodeEnumerator(reverse_hmm, rev_hcost));
		}

		for (int i = 0; i < num_threads; ++i) {
			search[i].constructPool();
		}

		HashMapST<AStarNode, AStarNode> term_nodes, term_nodes_rev;
		
		#pragma omp parallel for
		for (unsigned i = 0; i < starting_kmer_storage.size(); ++i) {
			search[omp_get_thread_num()].search(gene[0], starting_kmer_storage[i].first, forward_hmm, reverse_hmm, starting_kmer_storage[i].second, 
				for_node_enumerator[omp_get_thread_num()], rev_node_enumerator[omp_get_thread_num()], dbg, i, term_nodes, term_nodes_rev, out_file);
		}
		fclose(out_file);

        fprintf(stderr, "Done %s\n", gene[0].c_str());
	}

	return 0;
}
