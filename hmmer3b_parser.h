#ifndef PARSER_H__
#define PARSER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "profile_hmm.h"
#include <vector>
#include <string.h>
#include <ctype.h>

using namespace std;

class Parser {
	public:
		Parser() {};
		static void readHMM(ifstream &myfile, ProfileHMM &hmm) {
			if (myfile.is_open())
			{
				string line;
				//version
				getline (myfile, line);
				istringstream iss(line);
				string word, word2;
				iss >> word;
				hmm.version = word;

				//header
				while ( getline (myfile,line) ) {
					istringstream iss2(line);
					iss2 >> word >> word2;
					if (word == "NAME") {
						hmm.name = word2;
					}
					else if (word == "LENG") {
						hmm.model_length = stoi(word2);
					}
					else if (word == "ALPH") {
						if (word2 == "amino") {
							hmm.alphabet = ProfileHMM::protein;
						} else if (word2 == "RNA") {
							hmm.alphabet = ProfileHMM::nucleotide;
						} else if (word2 == "DNA") {
							hmm.alphabet = ProfileHMM::nucleotide;
						} else {

						}
					}
					else if (word == "HMM") {
						parseAlpha(line, hmm);
						break;
					}
				}

				//labels
				getline (myfile,line);

				//compo
				getline (myfile,line);
				istringstream iss3(line);
				string compo, probability;
				iss3 >> compo;
				if (compo == "COMPO") {
					for (int i = 0; i < hmm.alphabet_length; i++) {
						iss3 >> probability;
						double p = stod(probability);
						hmm.compo.push_back(p);
					}
				}

				for (int i = 0; i < hmm.NUM_TRANSITIONS + 1; i++) {
					for (int j = 0; j < hmm.model_length + 1; j++) {
						hmm.transitions[i].push_back(0.0);
					}
				}
				vector<vector<double>> vector_2d;
				vector<double> vector_1d;
				for (int i = 0; i < hmm.model_length + 1; i++) {
					hmm.emissions.push_back(vector_2d);
					for (int j = 0; j < hmm.alphabet_length; j++) {
						hmm.emissions[i].push_back(vector_1d);
						for (int k = 0; k < hmm.NUM_EMISSION_STATES; k++) {
							hmm.emissions[i][j].push_back(0.0);
						}
					}
				}

				//real
				string line_num;
				for (int i = 0; i <= hmm.model_length; i++) {			
					if (i > 0) {
						getline (myfile,line);
						istringstream iss4(line);
						iss4 >> line_num;
						for (int j = 0; j < hmm.alphabet_length; j++) {
							iss4 >> probability;
							double p;
							if (probability == "*") {
								p = 0.0;
							} else {
								p = stod(probability);
							}
							hmm.emissions[i][j][hmm.MSC] = p;
						}
					}
					getline (myfile,line);
					istringstream iss5(line);
					for (int j = 0; j < 20; j++) {
						iss5 >> probability;
						double p;
						if (probability == "*") {
							p = 0.0;
						} else {
							p = stod(probability);
						}
						hmm.emissions[i][j][hmm.ISC] = p;
					}

					getline (myfile,line);
					istringstream iss6(line);
					for (int j = 0; j < 7; j++) {
						iss6 >> probability;
						double p;
						if (probability == "*") {
							p = 0.0;
						} else {
							p = stod(probability);
						}
						hmm.transitions[j][i] = p;
					}
				}
				myfile.close();
			}
		}

		static void parseAlpha(string &line, ProfileHMM &hmm) {
			string alphabet;
			istringstream iss(line);
			iss >> alphabet;
			vector<char> alphabet_container;			
			fill_n(hmm.alpha_mapping, 127, -1);
			int count = 0;
			while (iss >> alphabet) {
				if (alphabet.size() != 1) {
					cout << "exception\n";
				}
				char alphabet_copy = alphabet[0];
				// strcpy(alphabet_copy, alphabet);
				alphabet_container.push_back(alphabet_copy);
				hmm.alpha_mapping[(int) toupper(alphabet_copy)] = count;
				hmm.alpha_mapping[(int) tolower(alphabet_copy)] = count;
				count++;
			}
			hmm.alphabet_length = alphabet_container.size();
		}
};

#endif