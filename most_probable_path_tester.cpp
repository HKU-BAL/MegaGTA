#include "profile_hmm.h"
#include "most_probable_path.h"
#include "hmmer3b_parser.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
	ifstream hmm_file (argv[1]);
	ProfileHMM hmm = ProfileHMM(true);
	Parser::readHMM(hmm_file, hmm);
	MostProbablePath hcost = MostProbablePath(hmm);

	for (int i = 0; i <= hmm.modelLength(); i++) {
		cout << "hcost.computeHeuristicCost('m', " << i <<") = "<< hcost.computeHeuristicCost('m', i) << '\n';
		cout << "hcost.computeHeuristicCost('i', " << i <<") = "<< hcost.computeHeuristicCost('i', i) << '\n';
		cout << "hcost.computeHeuristicCost('d', " << i <<") = "<< hcost.computeHeuristicCost('d', i) << '\n';
	}
	
	return 0;
}