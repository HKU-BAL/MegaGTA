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
#include <math.h>

using namespace std;

class Parser {
  public:
    Parser() {};
    static void readHMM(ifstream &myfile, ProfileHMM &hmm) {
        if (myfile.is_open()) {
            string line;
            //version
            getline (myfile, line);
            istringstream iss(line);
            string word, word2;
            iss >> word;
            hmm.version = word;

            //header
            while ( getline (myfile, line) ) {
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
                    }
                    else if (word2 == "RNA") {
                        hmm.alphabet = ProfileHMM::nucleotide;
                    }
                    else if (word2 == "DNA") {
                        hmm.alphabet = ProfileHMM::nucleotide;
                    }
                    else {

                    }
                }
                else if (word == "HMM") {
                    parseAlpha(line, hmm);
                    break;
                }
            }

            //labels
            getline (myfile, line);

            //compo
            getline (myfile, line);
            istringstream iss3(line);
            string compo, probability;
            iss3 >> compo;

            if (compo == "COMPO") {
                for (int i = 0; i < hmm.alphabetLength(); i++) {
                    iss3 >> probability;
                    double p = exp(-1 * stod(probability));
                    hmm.compo.push_back(p);
                }
            }

            for (int i = 0; i < hmm.NUM_TRANSITIONS + 1; i++) {
                for (int j = 0; j < hmm.modelLength() + 1; j++) {
                    hmm.transitions[i].push_back(0.0);
                }
            }

            vector<vector<double>> vector_2d;
            vector<double> vector_1d;

            for (int i = 0; i < hmm.modelLength() + 1; i++) {
                hmm.emissions.push_back(vector_2d);

                for (int j = 0; j < hmm.alphabetLength(); j++) {
                    hmm.emissions[i].push_back(vector_1d);

                    for (int k = 0; k < hmm.NUM_EMISSION_STATES; k++) {
                        hmm.emissions[i][j].push_back(0.0);
                    }
                }
            }

            for (int i = 0; i < hmm.modelLength() + 1; i++) {
                hmm.max_match_emissions.push_back(- numeric_limits<double>::infinity());
            }

            //real
            string line_num;

            for (int i = 0; i <= hmm.modelLength(); i++) {
                if (i > 0) {
                    getline (myfile, line);
                    istringstream iss4(line);
                    iss4 >> line_num;

                    for (int j = 0; j < hmm.alphabetLength(); j++) {
                        iss4 >> probability;
                        double p;

                        if (probability == "*") {
                            p = 0.0;
                        }
                        else {
                            p = exp(-1 * stod(probability));
                        }

                        if (hmm.normalized == true) {
                            hmm.msc(i, j, log(p / hmm.compo[j]));
                        }
                        else {
                            hmm.msc(i, j, log(p));
                        }
                    }
                }

                getline (myfile, line);
                istringstream iss5(line);

                for (int j = 0; j < hmm.alphabetLength(); j++) {
                    iss5 >> probability;
                    double p;

                    if (probability == "*") {
                        p = 0.0;
                    }
                    else {
                        p = exp(-1 * stod(probability));
                    }

                    if (hmm.normalized == true) {
                        hmm.isc(i, j, 0);
                    }
                    else {
                        hmm.isc(i, j, log(p));
                    }
                }

                getline (myfile, line);
                istringstream iss6(line);

                for (int j = 0; j < hmm.NUM_TRANSITIONS; j++) {
                    iss6 >> probability;
                    double p;

                    if (probability == "*") {
                        p = 0.0;
                    }
                    else {
                        p = exp(-1 * stod(probability));
                    }

                    hmm.tsc(i, j, log(p));
                }

                for (int i = 0; i < hmm.alphabetLength(); i++) {
                    hmm.isc(hmm.modelLength(), i, - numeric_limits<double>::infinity());
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