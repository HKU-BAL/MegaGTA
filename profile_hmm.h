#ifndef PROFILEHMM_H__
#define PROFILEHMM_H__

#include <string>
#include <vector>
#include <limits>

using namespace std;

class ProfileHMM {
	friend class Parser;
	public:
		ProfileHMM() {};
		ProfileHMM(const bool &normalized) : normalized(normalized) {};
		static const int NUM_EMISSION_STATES = 2;
		static const int NUM_TRANSITIONS = 7;
		static const int NXSTATES = 5;
		static const int NXTRANS = 2;
		static const int MSC = 0;
		static const int ISC = 1;
		string version;
		string name;
		enum sequence_type {nucleotide, protein, unknown};
		enum TSC {MM = 0, MI = 1, MD = 2, IM = 3, II = 4, DM = 5, DD = 6, BM = 7};
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

		bool normalized = true;

	public:
		int modelLength() {
			return model_length;
		}
		int alphabetLength() {
			return alphabet_length;
		}
		double getMaxMatchEmission(int i) {
			if (normalized) {
				return max_match_emissions[i];
			}
			else {
				return 0.0;
			}
		}
		sequence_type getAlphabet() {
			return alphabet;
		}
		double msc(int k, char b) {
			if (alpha_mapping[(int)b] == -1) {

			}
			return emissions[k][alpha_mapping[(int)b]][MSC];
		}
		double msc(int k, int b) {
			if (k == 0) {
				return - numeric_limits<double>::infinity();
			}
			return emissions[k][b][MSC];
		}
		void msc(int k, int i, double val) {
			emissions[k][i][MSC] = val;
			if (val > max_match_emissions[k]) {
				max_match_emissions[k] = val;
			}
		}
		double isc(int k, char b) {
			return emissions[k][alpha_mapping[(int)b]][ISC];
		}
		double isc(int k, int b) {
			return emissions[k][b][ISC];
		}
		void isc(int k, int i, double val) {
			emissions[k][i][ISC] = val;
		}
		double tsc(int k, TSC trans) {
			return transitions[(int)trans][k];
		}
		void tsc(int k, int trans, double val) {
			transitions[trans][k] = val;
		}
		void tsc(int k, TSC trans, double val) {
			transitions[(int)trans][k] = val;
		}



};

#endif