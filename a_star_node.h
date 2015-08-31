#ifndef A_STAR_NODE_H__
#define A_STAR_NODE_H__

#include "nucl_kmer.h"

using namespace std;

class AStarNode
{
public:
	int partial : 2;
	AStarNode *discovered_from;
	NuclKmer kmer; //128 bits
	int nucl_emission : 9;
	double score;
	char state;
	int16_t state_no;
	int fval;
	double this_node_score;
	int16_t indels;
	double real_score;
	int16_t length;
	char emission;
	int64_t node_id;

	int16_t negative_count = 0;
	double max_score = 0;

	AStarNode() : discovered_from(NULL) {nucl_emission = 0; partial = 0;};
	AStarNode(AStarNode *discovered_from, NuclKmer &kmer, int state_no, char state)
	: discovered_from(discovered_from), kmer(kmer), state_no(state_no), state(state) {};
	// ~AStarNode();
	bool operator< (const AStarNode &node2) const {
		if (fval < node2.fval) {
			return true;
		}
		else if (fval > node2.fval) {
			return false;
		}
		else {
			if (state_no > node2.state_no) {
				return true;
			}
			else if (state_no < node2.state_no) {
				return false;
			}
			else {
				int s1 = 0;
				int s2 = 0;
				switch (state) {
					case 'm': s1 = 3; break;
					case 'd': s1 = 2; break;
					case 'i': s1 = 1; break;
				}
				switch (node2.state) {
					case 'm': s2 = 3; break;
					case 'd': s2 = 2; break;
					case 'i': s2 = 1; break;
				}
				return s1 < s2;
			}
		}
	}
	bool operator== (const AStarNode &node2) const {
		if (!(kmer == node2.kmer)) {
			return false;
		}
		if (state != node2.state) {
			return false;
		}
		if (state_no != node2.state_no) {
			return false;
		}
		return true;
	}

	uint64_t hash() const {
		return CityHash64((const char*)kmer.kmers, sizeof(kmer.kmers[0]) * 2) + state + state_no;
	}
};

#endif