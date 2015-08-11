#include "node_enumerator.h"
#include "profile_hmm.h"
#include "most_probable_path.h"
#include "hmmer3b_parser.h"
#include <iostream>
#include <fstream>
#include "succinct_dbg.h"
#include <vector>

using namespace std;

int main(int argc, char **argv) {
	NuclKmer::setUp();
	SuccinctDBG dbg;
	dbg.LoadFromMultiFile(argv[1], false);
	ifstream hmm_file (argv[2]);
	ProfileHMM hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file, hmm);
	MostProbablePath hcost = MostProbablePath(hmm);
	NodeEnumerator node_enumerator = NodeEnumerator(hmm, hcost);
	string kmer = "cggaagcgcaagacctcggaccgtttcatcgtcacccgtcgtaag";
	// string kmer = "cgtaataaaaaagctaaatcagacaaacttatcgttcgtcgtcgt";
	NuclKmer nucl_kmer = NuclKmer(kmer);
	AStarNode curr = AStarNode(NULL, nucl_kmer, 17, 'm');
	vector<AStarNode> result = node_enumerator.enumeratorNodes(curr, dbg);
	cout << "result size = " << result.size() << '\n';
	for (int i = 0; i < result.size(); i++) {
		cout << result[i].state << " "<<result[i].kmer.decodePacked() << '\n';
	}

	return 0;
}