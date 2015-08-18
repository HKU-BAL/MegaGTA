#ifndef NODE_ENUMERATOR_H__
#define NODE_ENUMERATOR_H__

#include "profile_hmm.h"
#include "most_probable_path.h"
#include "nucl_kmer.h"
#include "succinct_dbg.h"
#include "a_star_node.h"
#include <iostream>
#include <vector>
#include "codon.h"
#include <set>

using namespace std;

class NodeEnumerator
{
private:
	static const int SCALE = 10000;
	static constexpr double hweight = 1.0;
	AStarNode next;
	uint8_t next_nucl;
	char emission;
	double match_trans;
	double ins_trans;
	double del_trans;
	NuclKmer next_kmer;
	int next_state;
	ProfileHMM *hmm;
	bool prot_search;
	MostProbablePath *hcost;
	static int dna_map[128];

public:
	NodeEnumerator(ProfileHMM &_hmm, MostProbablePath &_hcost) {
		hmm = &_hmm;
		prot_search = _hmm.getAlphabet() == ProfileHMM::protein;
		hcost = &_hcost;
		memset(dna_map, -1, sizeof(dna_map));
		for (int i = 0; i < 10; ++i) {
			dna_map["ACGTNacgtn"[i]] = "1234312343"[i] - '0';
		}
	};
	~NodeEnumerator() {};
	vector<AStarNode> enumeratorNodes(AStarNode &curr, bool forward, SuccinctDBG &dbg) {
		return enumeratorNodes(curr, forward, dbg, NULL);
	}
	vector<AStarNode> enumeratorNodes(AStarNode &curr, bool forward, SuccinctDBG &dbg, AStarNode *child_node) {
		//test
		// if (!forward) {
		// 	cout << "RC = "<<curr.kmer.decodePacked() << '\n';
		// }
		vector<AStarNode> ret;
		next_state = curr.state_no + 1;
		switch (curr.state) {
			case 'm':
				match_trans = hmm->tsc(curr.state_no, ProfileHMM::MM);
				ins_trans = hmm->tsc(curr.state_no, ProfileHMM::MI);
				del_trans = hmm->tsc(curr.state_no, ProfileHMM::MD);
				break;
			case 'd':
				match_trans = hmm->tsc(curr.state_no, ProfileHMM::DM);
				ins_trans = - numeric_limits<double>::infinity();
				del_trans = hmm->tsc(curr.state_no, ProfileHMM::DD);
				break;
			case 'i':
				match_trans = hmm->tsc(curr.state_no, ProfileHMM::IM);
				ins_trans = hmm->tsc(curr.state_no, ProfileHMM::II);
				del_trans = - numeric_limits<double>::infinity();
				break;
			default:
				assert(false);
		}
		double max_match_emission = hmm->getMaxMatchEmission(next_state);
		uint8_t seq[dbg.kmer_k];
		string buf = curr.kmer.decodePacked();
		for (int i = 0; i < dbg.kmer_k; ++i) {
			seq[i] = dna_map[buf[i]]; // $->0, A->1, C->2, G->3, T->4
		}
		int64_t node_id = dbg.IndexBinarySearch(seq);
		if (node_id == -1) {
	    	cout << "no such seq \n";
	    	return ret;
	    } else {
	    	int64_t next_node[4], next_node_2[4], next_node_3[4];
	    	vector<uint8_t> codons[64];
	    	int outd = dbg.NextNodes(node_id, next_node);
	    	if (outd == 0) {
	    		return ret;
	    	} else {
	    		for (int i = 0; i < outd; ++i) {
	    			for (int j = 0; j < 16; ++j) {
	    				codons[i * 16 + j].push_back(dbg.GetNodeLastChar(next_node[i])-1);
	    			}
	    			int outd_2 = dbg.NextNodes(next_node[i], next_node_2);
	    			if (outd_2 == 0) {
	    				continue;
	    			}
    				for (int j = 0; j < outd_2; ++j) {
    					for (int k = 0; k < 4; ++k) {
    						codons[i * 16 + j * 4 + k].push_back(dbg.GetNodeLastChar(next_node_2[j])-1);
    					}
    					int outd_3 = dbg.NextNodes(next_node_2[j], next_node_3);
    					if (outd_3 == 0) {
    						continue;
    					}
    					for (int k = 0; k < outd_3; ++k) {
    						codons[i * 16 + j * 4 + k].push_back(dbg.GetNodeLastChar(next_node_3[k])-1);
    					}
    				}
		    	}
	    	}

	    	//translate to aa
	    	set<char> aa;
	    	for (int i = 0; i < 64; ++i) {
	    		if (codons[i].size() == 3) {
	    			// cout << curr.kmer.decodePacked() << " " << curr.state_no <<'\n';
	    			next_kmer = curr.kmer.shiftLeftCopy(codons[i][0], codons[i][1], codons[i][2]);
	    			// cout << next_kmer.decodePacked() << " " << next_state << '\n';
	    			if (!forward) {
	    				emission = Codon::rc_codonTable[codons[i][0]][codons[i][1]][codons[i][2]];
	    			} else {
	    				emission = Codon::codonTable[codons[i][0]][codons[i][1]][codons[i][2]];
	    			}
	    			if (emission == '*') {
	    				continue;
	    			}
	    			set<char>::iterator it = aa.find(emission);
	    			aa.insert(emission);
	    			// cout << "emission = " << emission << '\n';
    				if (it == aa.end()) {
    					next = AStarNode(&curr, next_kmer, next_state, 'm');
			    		next.real_score = curr.real_score + match_trans + hmm->msc(next_state, emission);
			    		if (next.real_score >= curr.max_score) {
			    			next.max_score = next.real_score;
			    			next.negative_count = 0;
			    		} else {
			    			next.max_score = curr.max_score;
			    			next.negative_count = curr.negative_count + 1;
			    		}
			    		for (int j = 0; j < 3; j++) {
			    			next.nucl_emission[j] = "acgt"[codons[i][j]];
			    		}
			    		next.emission = emission;
			    		next.this_node_score = match_trans + hmm->msc(next_state, emission) - max_match_emission;
			    		next.length = curr.length + 1;
			    		next.score = curr.score + next.this_node_score;
			    		next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('m', next_state)));
			    		next.indels = curr.indels;

			    		ret.push_back(next);

			    		//insert node
			    		if (curr.state != 'd') {
			    			next = AStarNode(&curr, next_kmer, curr.state_no, 'i');
				    		next.real_score = curr.real_score + ins_trans + hmm->isc(next_state, emission);
				    		next.max_score = curr.max_score;
				    		next.negative_count = curr.negative_count + 1;
				    		for (int j = 0; j < 3; j++) {
				    			next.nucl_emission[j] = "acgt"[codons[i][j]];
				    		}
				    		next.emission = emission;
				    		next.this_node_score = ins_trans + hmm->isc(next_state, emission);
				    		next.length = curr.length + 1;
				    		next.score = curr.score + next.this_node_score;
				    		next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('i', curr.state_no)));
				    		next.indels = curr.indels + 1;

				    		ret.push_back(next);
			    		}
    				}
	    		}
	    	}

	    	//delete node
	    	if (curr.state != 'i') {
	    		next = AStarNode(&curr, curr.kmer, next_state, 'd');
	    		next.real_score = curr.real_score + del_trans;
	    		next.max_score = curr.max_score;
	    		next.negative_count = curr.negative_count + 1;
	    		for (int j = 0; j < 3; j++) {
	    			next.nucl_emission[j] = '-';
	    		}
	    		next.emission = '-';
	    		next.this_node_score = del_trans - max_match_emission;
	    		next.length = curr.length;
	    		next.score = curr.score + next.this_node_score;
	    		next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('d', next_state)));
	    		next.indels = curr.indels + 1;

	    		ret.push_back(next);
	    	}
	    }
	    return ret;
	}
	
};

#endif