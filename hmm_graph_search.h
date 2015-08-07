#include "profile_hmm.h"
#include "succinct_dbg.h"
#include "nucl_kmer.h"
#include <vector>
#include "hash_set.h"
#include <math.h>
#include <unordered_map>
#include "src/sequence/NTSequence.h"
#include "src/sequence/AASequence.h"

using namespace std;

class HMMGraphSearch
{
private:
	int heuristic_pruning = 20;
	static double exit_probabilities[3000];
	unordered_map<AStarNode, AStarNode> term_nodes;

public:
	HMMGraphSearch(arguments) {
		for (int i = 0; i < 3000; i++) {
            exit_probabilities[i] = log(2.0 / (i + 2)) * 2;
        }
	};
	~HMMGraphSearch();

	vector<>
	double scoreStart(ProfileHMM &hmm, string &starting_kmer, int starting_state) {
		double ret = 0;
		for (int i = 1; i <= starting_kmer.size(); i++) {
			ret += hmm.msc(starting_state + i, starting_kmer[i-1]) + hmm.tsc(starting_state + i - 1, ProfileHMM::MM) - hmm.getMaxMatchEmission(starting_state + i);
		}
		return ret;
	}
	
	double realScoreStart(ProfileHMM &hmm, string &starting_kmer, int starting_state) {
		double ret = 0;
		for (int i = 1; i <= starting_kmer.size(); i++) {
			ret += hmm.msc(starting_state + i, starting_kmer[i-1]) + hmm.tsc(starting_state + i - 1, ProfileHMM::MM);
		}
		return ret;
	}

	bool astarSearch(const ProfileHMM &hmm, AStarNode &starting_node, string &framed_word, SuccinctDBG &dbg, bool forward, NodeEnumerator &nodeEnumerator, AStarNode &goal_node) {
		string scoring_word;
		if (!forward) {
			if (hmm.getAlphabet() == ProfileHMM::protein) {
				seq::NTSequence nts = seq::NTSequence("", "", framed_word);
	    		seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
	    		scoring_word = aa.asString();
	    		reverse(scoring_word.begin(), scoring_word.end());
			}
		} else if (hmm.getAlphabet() == ProfileHMM::protein) {
			seq::NTSequence nts = seq::NTSequence("", "", framed_word);
	    	seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
			scoring_word = aa.asString();
		}
		NuclKmer kmer = NuclKmer(framed_word);
		AStarNode starting_node;
		if (hmm.getAlphabet() == ProfileHMM::protein) {
			starting_node = AStarNode(NULL, kmer, starting_state + (framed_word.size() / 3), 'm');
			starting_node.length = framed_word.size() / 3;
		} else {
			starting_node = AStarNode(NULL, kmer, starting_state, 'm');
			starting_node.length = framed_word.size();
		}
		starting_node.fval = 0;
		starting_node.score = scoreStart(hmm, scoring_word, starting_state);
		starting_node.real_score = realScoreStart(hmm, scoring_word, starting_state);

		return astarSearch(hmm, starting_node, dbg, nodeEnumerator, goal_node);
	}

	bool astarSearch(const ProfileHMM &hmm, AStarNode &starting_node, SuccinctDBG &dbg, NodeEnumerator &nodeEnumerator, AStarNode &goal_node) {
		if (starting_node.state_no >= hmm.modelLength()) {
			goal_node = starting_node;
			return true;
		}
		priority_queue<AStarNode, vector<AStarNode>> open;
		HashSet<AStarNode> closed;
		AStarNode curr;
		int opened_nodes = 1;

		unordered_map<AStarNode, AStarNode> open_hash;

		int repeated_nodes = 0;
		int replaced_nodes = 0;
		int pruned_nodes = 0;

		//need to add a cache here
		unordered_map<AStarNode, AStarNode>::const_iterator got = term_nodes.find(starting_node);
		if (got == term_nodes.end()) {
			for (AStarNode next : nodeEnumerator.enumeratorNodes(starting_node, dbg)) {
				open.push(next);
			}
		} else {
			for (AStarNode next : nodeEnumerator.enumeratorNodes(starting_node, dbg, got->second)) {
				open.push(next);
			}
		}
		

		if (open.empty()) {
			return false;
		}
		AStarNode inter_goal = starting_node;
		while (!open.empty()) {
			curr = open.top();
			open.pop();
			HashSet<AStarNode>::iterator iter = closed.find(kmer);
			if (iter != NULL) {
				continue;
			}
			if (curr.state_no >= hmm.modelLength()) {
				curr.partial = false;
				if ((curr.real_score + exit_probabilities[curr.length]) / log(2) 
						> (inter_goal.real_score + exit_probabilities[inter_goal.length]) / log(2)) {
					inter_goal = curr;
				}
				return getHighestScoreNode(inter_goal, goal_node);
			}

			closed.insert(curr);
			if ((curr.real_score + exit_probabilities[curr.length]) / log(2) 
					> (inter_goal.real_score + exit_probabilities[inter_goal.length]) / log(2)) {
				inter_goal = curr;
			}
			vector<AStarNode> temp_nodes_to_open;
			got = term_nodes.find(curr);
			if (got == term_nodes.end()) {
				temp_nodes_to_open = nodeEnumerator.enumeratorNodes(curr, dbg);
			} else {
				temp_nodes_to_open = nodeEnumerator.enumeratorNodes(curr, dbg, got->second);
			}
			for (AStarNode next : temp_nodes_to_open) {
				bool open_node = false;
				if (heuristic_pruning > 0) {
					if ((next.length < 5 || next.negative_count <= heuristic_pruning) && next.real_score > 0.0) {
						got = open_hash.find(next);
						if (got != open_hash.end()) {
							repeated_nodes++;
							if (got->second < next) {
								replaced_nodes++;
								open_node = true;
							}
						} else {
							open_node = true;
						}
					} else {
						pruned_nodes++;
					}
				} else {
					got = open_hash.find(next);
					if (got != open_hash.end()) {
						repeated_nodes++;
						if (got->second < next) {
							replaced_nodes++;
							open_node = true;
						}
					} else {
						open_node = true;
					}
				}

				if (open_node) {
					open_hash.insert(next, next);
					open_node++;
					open.push(next);
				}
			}
		}

		inter_goal.partial = true;
		return getHighestScoreNode(inter_goal, goal_node);
	}

	bool getHighestScoreNode(AStarNode &inter_goal, AStarNode &goal_node) {
		AStarNode temp_goal = inter_goal;
		goal_node = node;
		while (temp_goal.discovered_from != NULL) {
			temp_goal = *temp_goal.discovered_from;
			if (temp_goal.real_score > goal_node.real_score) {
				goal_node = temp_goal;
			}
		}
		return true;
	}
};