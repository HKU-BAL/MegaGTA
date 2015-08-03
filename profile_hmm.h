#ifndef PROFILEHMM_H__
#define PROFILEHMM_H__

#include <string>
#include <vector>

using namespace std;

class ProfileHMM {
	public:
		static const int NUM_EMISSION_STATES = 2;
		static const int NUM_TRANSITIONS = 7;
		static const int NXSTATES = 5;
		static const int NXTRANS = 2;
		static const int MSC = 0;
		static const int ISC = 1;
		string version;
		string name;
		enum sequence_type {nucleotide, protein, unknown};
		sequence_type alphabet;
		vector<double> transitions[NUM_TRANSITIONS+1];
		vector<vector<vector<double>>> emissions;
		vector <double> compo;
		double xsc[NXSTATES][NXTRANS];
		int alpha_mapping[127];
		vector <double> max_match_emissions;

		int model_length;
		int alphabet_length;
		int l;

		bool normalized;


};

#endif