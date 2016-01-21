#include "hmm_graph_search.h"
#include "node_enumerator.h"
#include "profile_hmm.h"
#include "most_probable_path.h"
#include "hmmer3b_parser.h"
#include "succinct_dbg.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    // NuclKmer::setUp();
    HMMGraphSearch::setUp();
    int pruning = 20;
    HMMGraphSearch search = HMMGraphSearch(pruning);
    SuccinctDBG dbg;
    dbg.LoadFromMultiFile(argv[1], false);
    ifstream hmm_file (argv[2]);
    ProfileHMM forward_hmm = ProfileHMM(true);
    Parser::readHMM(hmm_file, forward_hmm);
    ifstream hmm_file_2 (argv[3]);
    ProfileHMM reverse_hmm = ProfileHMM(true);
    Parser::readHMM(hmm_file_2, reverse_hmm);
    MostProbablePath for_hcost = MostProbablePath(forward_hmm);
    NodeEnumerator for_node_enumerator = NodeEnumerator(forward_hmm, for_hcost);
    MostProbablePath rev_hcost = MostProbablePath(reverse_hmm);
    NodeEnumerator rev_node_enumerator = NodeEnumerator(reverse_hmm, rev_hcost);

    HashMapSingleThread<AStarNode, AStarNode> term_nodes;

    int count = 0;

    // string starting_kmer = "atggctattaaaaagtataagccaataacaaatggtcgtcgtaat";
    // int starting_state = 0;
    // search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, for_node_enumerator, rev_node_enumerator, dbg, count);

    // string starting_kmer = "CACAAACGTCAATACCGTGTTATCGATTTT";
    // int starting_state = 57;
    // search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, for_node_enumerator, rev_node_enumerator, dbg, count, term_nodes);

    string starting_kmer = "AAGCGCCTCTACCGCGTCATCGACTTCAAG";
    int starting_state = 58;
    search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, for_node_enumerator, rev_node_enumerator, dbg, count, term_nodes);

    // starting_kmer = "cgtggtaaaaaatcatcagacaaacttatcgttcgtggacgtaag";
    // starting_state = 259;
    // search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, for_node_enumerator, rev_node_enumerator, dbg);

    return 0;
}