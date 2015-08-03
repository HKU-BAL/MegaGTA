#include <iostream>
#include <fstream>
#include "hmmer3b_parser.h"

using namespace std;

int main(int argc, char **argv) {
	ifstream hmm_file (argv[1]);
	ProfileHMM hmm;
	Parser::readHMM(hmm_file, hmm);

	for (int i = 0; i < hmm.model_length + 1; i++) {
		for (int j = 0; j < hmm.alphabet_length; j++) {
			cout <<hmm.emissions[i][j][0] << " ";
		}
		cout << "\n";
	}
	return 0;
}