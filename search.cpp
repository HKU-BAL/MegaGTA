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

int main(int argc, char **argv) {
	if (argc < 5) {
		fprintf(stderr, "Usage: %s <succinct_dbg> <forward_hmm> <reverse_hmm> <starting_kmers>\n", argv[0]);
		exit(1);
	}

	setvbuf ( stdout , NULL , _IOLBF , 1024 );

	int heuristic_pruning = 20;
	NuclKmer::setUp();
	HMMGraphSearch::setUp();
	SuccinctDBG dbg;
	dbg.LoadFromMultiFile(argv[1], false);
	ifstream hmm_file (argv[2]);
	ProfileHMM forward_hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file, forward_hmm);
	ifstream hmm_file_2 (argv[3]);
	ProfileHMM reverse_hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file_2, reverse_hmm);
	MostProbablePath for_hcost = MostProbablePath(forward_hmm);
	MostProbablePath rev_hcost = MostProbablePath(reverse_hmm);

	pair<string, int> starting_kmer;
	vector<pair<string, int>> starting_kmer_storage;
	ifstream starting_kmer_file (argv[4]);
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
	}
	vector<HMMGraphSearch> search;
	vector<NodeEnumerator> for_node_enumerator, rev_node_enumerator;
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		search.push_back(HMMGraphSearch(heuristic_pruning));
		for_node_enumerator.push_back(NodeEnumerator(forward_hmm, for_hcost));
		rev_node_enumerator.push_back(NodeEnumerator(reverse_hmm, rev_hcost));
	}
	
		#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < starting_kmer_storage.size(); ++i) {
		search[omp_get_thread_num()].search(starting_kmer_storage[i].first, forward_hmm, reverse_hmm, starting_kmer_storage[i].second, for_node_enumerator[omp_get_thread_num()], rev_node_enumerator[omp_get_thread_num()], dbg, i);
	}

	return 0;
}