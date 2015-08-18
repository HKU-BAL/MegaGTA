#ifndef HMM_GRAPH_SEARCH_H__
#define HMM_GRAPH_SEARCH_H__

#include "profile_hmm.h"
#include "succinct_dbg.h"
#include "nucl_kmer.h"
#include <vector>
#include "hash_set.h"
#include "hash_map.h"
#include <math.h>
#include <queue>
#include "src/sequence/NTSequence.h"
#include "src/sequence/AASequence.h"
#include "node_enumerator.h"
#include <iostream>

using namespace std;

class HMMGraphSearch
{
private:
	int heuristic_pruning = 20;
	static double exit_probabilities[3000];
	HashMap<AStarNode, AStarNode> term_nodes;
	deque<AStarNode*> created_nodes;

public:
	HMMGraphSearch(int pruning) : heuristic_pruning(pruning) {
		for (int i = 0; i < 3000; i++) {
            exit_probabilities[i] = log(2.0 / (i + 2)) * 2;
        }
	};
	~HMMGraphSearch() {};
	void search(string &starting_kmer, ProfileHMM &forward_hmm, ProfileHMM &reverse_hmm, int &start_state, NodeEnumerator &node_enumerator, SuccinctDBG &dbg) {
		printf("left \n");
		int l_starting_state = reverse_hmm.modelLength() - start_state - starting_kmer.size() / (reverse_hmm.getAlphabet() == ProfileHMM::protein ? 3 : 1);
		AStarNode goal_node;
		astarSearch(reverse_hmm, l_starting_state, starting_kmer, dbg, false, node_enumerator, goal_node);
		string max_seq, nucl_seq;
		partialResultFromGoal(goal_node, false, max_seq);
		nucl_seq = max_seq + starting_kmer;
		cout << "left nucl_seq = "<<nucl_seq <<'\n';

		printf("right \n");
		astarSearch(forward_hmm, start_state, starting_kmer, dbg, true, node_enumerator, goal_node);
		max_seq = "";
		nucl_seq = "";

		partialResultFromGoal(goal_node, true, max_seq);
		delectAStarNodes();
		nucl_seq = starting_kmer + max_seq;
		cout << "right nucl_seq = "<<nucl_seq <<'\n';
	}
	void partialResultFromGoal(AStarNode &goal, bool forward, string &max_seq) {
		while (goal.discovered_from != NULL) {
			if (goal.state != 'd') {
				if (forward) {
					max_seq = goal.nucl_emission + max_seq;
				} else {
					max_seq = max_seq + goal.nucl_emission;
				}
			}
			goal = *goal.discovered_from;
		}		
	}


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

	bool astarSearch(ProfileHMM &hmm, int &starting_state, string &framed_word, SuccinctDBG &dbg, bool forward, NodeEnumerator &node_enumerator, AStarNode &goal_node) {
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

		return astarSearch(hmm, starting_node, dbg, forward, node_enumerator, goal_node);
	}

	bool astarSearch(ProfileHMM &hmm, AStarNode &starting_node, SuccinctDBG &dbg, bool forward, NodeEnumerator &node_enumerator, AStarNode &goal_node) {
		if (starting_node.state_no >= hmm.modelLength()) {
			goal_node = starting_node;
			return true;
		}
		priority_queue<AStarNode, vector<AStarNode>> open;
		HashSet<AStarNode> closed;
		AStarNode curr;
		int opened_nodes = 1;

		HashMap<AStarNode, AStarNode> open_hash;

		int repeated_nodes = 0;
		int replaced_nodes = 0;
		int pruned_nodes = 0;

		//need to add a cache here
		HashMap<AStarNode, AStarNode>::iterator got = term_nodes.find(starting_node);
		if (got == term_nodes.end()) {
			for (AStarNode next : node_enumerator.enumeratorNodes(starting_node, dbg)) {
				// cout << "enumeratorNodes = " << next.state_no<< '\n';
				// cout << "inter_goal discovered_from " << next.discovered_from->state_no <<'\n';
				open.push(next);
			}
		} else {
			for (AStarNode next : node_enumerator.enumeratorNodes(starting_node, dbg, &got->second)) {
				// cout << "enumeratorNodes = " << next.state_no<< '\n';
				// cout << "inter_goal discovered_from " << next.discovered_from->state_no <<'\n';
				open.push(next);
			}
		}
		

		if (open.empty()) {
			return false;
		}
		AStarNode inter_goal = starting_node;
		while (!open.empty()) {
			// cout << "open size = " <<open.size() << '\n';
			auto curr_ptr = new AStarNode();
			auto &curr = *curr_ptr;
			curr = open.top();

			created_nodes.push_back(curr_ptr);

			// cout << curr.kmer.decodePacked() << " " << curr.state_no <<" state = "<< curr.state <<" kmer = " << curr.kmer.decodePacked() <<'\n';
			// cout << "curr discovered_from " << curr.discovered_from->state_no <<'\n';
			open.pop();
			HashSet<AStarNode>::iterator iter = closed.find(curr);
			if (iter != NULL) {
				continue;
			}
			if (curr.state_no >= hmm.modelLength()) {
				// cout << curr.state_no << " " << hmm.modelLength() << " \n";
				curr.partial = false;
				if ((curr.real_score + exit_probabilities[curr.length]) / log(2) 
						> (inter_goal.real_score + exit_probabilities[inter_goal.length]) / log(2)) {
					inter_goal = curr;
				}
				getHighestScoreNode(inter_goal, goal_node);
				return true;
			}

			closed.insert(curr);
			if ((curr.real_score + exit_probabilities[curr.length]) / log(2) 
					> (inter_goal.real_score + exit_probabilities[inter_goal.length]) / log(2)) {
				inter_goal = curr;
			}
			vector<AStarNode> temp_nodes_to_open;
			got = term_nodes.find(curr);
			if (got == term_nodes.end()) {
				temp_nodes_to_open = node_enumerator.enumeratorNodes(curr, dbg);
			} else {
				temp_nodes_to_open = node_enumerator.enumeratorNodes(curr, dbg, &got->second);
			}
			for (AStarNode &next : temp_nodes_to_open) {
				// cout << "enumeratorNodes = " << next.state_no <<" state = "<< next.state <<" kmer = " << next.kmer.decodePacked() << endl;
				// cout << "next discovered_from " << next.discovered_from->state_no <<'\n';
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
					pair<AStarNode, AStarNode> pair_insert (next, next);
					open_hash.insert(pair_insert);
					open_node++;
					open.push(next);
				}
			}
		}

		inter_goal.partial = true;
		getHighestScoreNode(inter_goal, goal_node);
		return true;
	}

	void getHighestScoreNode(AStarNode &inter_goal, AStarNode &goal_node) {
		// cout << "inter_goal discovered_from " << inter_goal.discovered_from->state_no <<'\n';
		AStarNode temp_goal = inter_goal;
		// cout << "temp_goal discovered_from " << temp_goal.discovered_from->state_no <<'\n';
		goal_node = inter_goal;
		while (temp_goal.discovered_from != NULL) {
			// cout << "curr " << temp_goal.state_no<< '\n';
			// cout << "discovered_from " << temp_goal.discovered_from->state_no << " state = " << temp_goal.discovered_from->state<< '\n';
			temp_goal = *temp_goal.discovered_from;
			if (temp_goal.real_score > goal_node.real_score) {
				goal_node = temp_goal;
			}
		}
	}

	void delectAStarNodes() {
		for (auto ptr : created_nodes) {
			delete ptr;
		}
		created_nodes.clear();
	}
};

#endif