#include "hmm_graph_search.h"
#include "node_enumerator.h"
#include "profile_hmm.h"
#include "most_probable_path.h"
#include "hmmer3b_parser.h"
#include "succinct_dbg.h"

using namespace std;

int main(int argc, char **argv) {
	NuclKmer::setUp();
	HMMGraphSearch search = HMMGraphSearch(20);
	SuccinctDBG dbg;
	dbg.LoadFromMultiFile(argv[1], false);
	ifstream hmm_file (argv[2]);
	ProfileHMM forward_hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file, forward_hmm);
	ifstream hmm_file_2 (argv[3]);
	ProfileHMM reverse_hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file_2, reverse_hmm);
	MostProbablePath hcost = MostProbablePath(forward_hmm);
	NodeEnumerator node_enumerator = NodeEnumerator(forward_hmm, hcost);

	string starting_kmer = "gctggacgcagccgctggctcgggcgcaagccacaccagcgcggc";
	int starting_state = 207;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "ggccgctctcgttggatgcgtaaacgcccaactgttcgtggtagc";
	starting_state = 208;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "gctggccgctctcgttggatgcgtaaacgcccaactgttcgtggt";
	starting_state = 207;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "aaggctggacgcagccgctggctcgggcgcaagccacaccagcgc";
	starting_state = 206;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "ggcaaggctggacgcagccgctggctcgggcgcaagccacaccag";
	starting_state = 205;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "ggtaaagctggccgctctcgttggatgcgtaaacgcccaactgtt";
	starting_state = 205;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "aaagctggccgctctcgttggatgcgtaaacgcccaactgttcgt";
	starting_state = 206;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "ctgggcaaggctggacgcagccgctggctcgggcgcaagccacac";
	starting_state = 204;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "attggtaaagctggccgctctcgttggatgcgtaaacgcccaact";
	starting_state = 204;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "aacattggtaaagctggccgctctcgttggatgcgtaaacgccca";
	starting_state = 203;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "tctggcgaagttcgtatgatcttagcaacttgccgtgcaacaatc";
	starting_state = 177;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "ccttccggcgaactgcgccgcgttcacagcgagtgctacgccacc";
	starting_state = 176;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "aactctggcgaagttcgtatgatcttagcaacttgccgtgcaaca";
	starting_state = 176;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "attaaaaagtacaaacctaccacaaatggccgtcgtaacatgaca";
	starting_state = 2;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "gccgtcaaaaagtaccgtccctatacccccagccgtcgccagatg";
	starting_state = 1;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	starting_kmer = "atggccgtcaaaaagtaccgtccctatacccccagccgtcgccag";
	starting_state = 0;
	search.search(starting_kmer, forward_hmm, reverse_hmm, starting_state, node_enumerator, dbg);

	return 0;
}