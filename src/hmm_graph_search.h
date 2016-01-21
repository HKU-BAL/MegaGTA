#ifndef HMM_GRAPH_SEARCH_H__
#define HMM_GRAPH_SEARCH_H__

#include "profile_hmm.h"
#include "succinct_dbg.h"
#include <vector>
#include "hash_set_st.h"
#include "hash_map_st.h"
#include "pool_st.h"
#include <math.h>
#include <queue>
#include <algorithm>
#include "sequence/NTSequence.h"
#include "sequence/AASequence.h"
#include "node_enumerator.h"
#include <iostream>
#include <stdio.h>

using namespace std;

class HMMGraphSearch {
  private:
    int heuristic_pruning = 20;
    static double exit_probabilities[3000];
    static int dna_map[128];

    // HashMapST<AStarNode, AStarNode> term_nodes;
    HashSetST<AStarNodePtr> closed;
    HashMapST<AStarNodePtr, AStarNodePtr> open_hash;

    PoolST<AStarNode> *pool_;

  public:
    HMMGraphSearch(int &pruning) : heuristic_pruning(pruning), pool_(NULL) {}
    void constructPool() {
        assert(!pool_);
        pool_ = new PoolST<AStarNode>;
        assert(pool_);
    }

    ~HMMGraphSearch() {
        if (pool_) {
            delete pool_;
        }
    }

    static void setUp() {
        for (int i = 0; i < 3000; i++) {
            exit_probabilities[i] = log(2.0 / (i + 2)) * 2;
        }

        memset(dna_map, -1, sizeof(dna_map));

        for (int i = 0; i < 10; ++i) {
            dna_map[int("ACGTNacgtn"[i])] = "1234312343"[i] - '0';
        }
    }

    void search(string &gene_name, string &starting_kmer, ProfileHMM &forward_hmm, ProfileHMM &reverse_hmm, int &start_state, NodeEnumerator &forward_enumerator,
                NodeEnumerator &reverse_enumerator, SuccinctDBG &dbg, int count, HashMapST<AStarNode, AStarNode> &term_nodes, HashMapST<AStarNode, AStarNode> &term_nodes_rev, FILE *out_file) {

        // if (start_state + starting_kmer.size() <= forward_hmm.modelLength() + 1) {
        //right, forward search
        AStarNode *goal_node = pool_->construct(), *goal_node2 = pool_->construct();
        string right_max_seq = "", left_max_seq = "";
        astarSearch(forward_hmm, start_state, starting_kmer, dbg, true, forward_enumerator, *goal_node, term_nodes);
        partialResultFromGoal(*goal_node, true, right_max_seq, term_nodes);

        // cout << "right start_state = " << start_state << endl;

        //left, reverse search
        int l_starting_state = reverse_hmm.modelLength() - start_state - starting_kmer.size() / (reverse_hmm.getAlphabet() == ProfileHMM::protein ? 3 : 1);
        astarSearch(reverse_hmm, l_starting_state, starting_kmer, dbg, false, reverse_enumerator, *goal_node2, term_nodes_rev);
        partialResultFromGoal(*goal_node2, false, left_max_seq, term_nodes_rev);
        deleteAStarNodes();
        RevComp(left_max_seq);

        fprintf(out_file, ">%s_contig_%d_contig_%d\n%s%s%s\n", gene_name.c_str(), count * 2, count * 2 + 1, left_max_seq.c_str(), starting_kmer.c_str(), right_max_seq.c_str());
        // }
    }

    void partialResultFromGoal(AStarNode &goal, bool forward, string &max_seq, HashMapST<AStarNode, AStarNode> &term_nodes) {
        int nucl_emission_mask = 0x7;
        auto ptr = &goal;
        max_seq.clear();

        static volatile int lock_ = 0;

        while (ptr->discovered_from != NULL) {
            // printf("in while %p\n", ptr->discovered_from);
            if (ptr->state != 'd') {
                for (int i = 0; i < 3; i++) {
                    max_seq.push_back("acgt-"[(ptr->nucl_emission >> 3 * i) & nucl_emission_mask]);
                }
            }

            pair<AStarNode, AStarNode> pair_to_cache (*(ptr->discovered_from), *ptr);

            while (__sync_lock_test_and_set(&lock_, 1)) while (lock_);

            term_nodes.insert(pair_to_cache);
            __sync_lock_release(&lock_);

            ptr = ptr->discovered_from;
        }

        reverse(max_seq.begin(), max_seq.end());
    }

    //bugs
    double scoreStart(ProfileHMM &hmm, string &starting_kmer, int starting_state) {
        double ret = 0;

        for (int i = 1; i <= (int)starting_kmer.size(); i++) {
            ret += hmm.msc(starting_state + i, starting_kmer[i - 1]) + hmm.tsc(starting_state + i - 1, ProfileHMM::MM) - hmm.getMaxMatchEmission(starting_state + i);
        }

        return ret;
    }

    double realScoreStart(ProfileHMM &hmm, string &starting_kmer, int starting_state) {
        double ret = 0;

        for (int i = 1; i <= (int)starting_kmer.size(); i++) {
            ret += hmm.msc(starting_state + i, starting_kmer[i - 1]) + hmm.tsc(starting_state + i - 1, ProfileHMM::MM);
        }

        return ret;
    }

    bool astarSearch(ProfileHMM &hmm, int &starting_state, string &framed_word, SuccinctDBG &dbg, bool forward, NodeEnumerator &node_enumerator,
                     AStarNode &goal_node, HashMapST<AStarNode, AStarNode> &term_nodes) {
        string scoring_word;

        if (!forward) {
            if (hmm.getAlphabet() == ProfileHMM::protein) {
                seq::NTSequence nts = seq::NTSequence("", "", framed_word);
                seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
                scoring_word = aa.asString();
                reverse(scoring_word.begin(), scoring_word.end());
            }
        }
        else if (hmm.getAlphabet() == ProfileHMM::protein) {
            seq::NTSequence nts = seq::NTSequence("", "", framed_word);
            seq::AASequence aa = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3);
            scoring_word = aa.asString();
        }

        // NuclKmer kmer;
        uint8_t seq[dbg.kmer_k + 1];

        if (!forward) {
            string rc_frame_word = framed_word;
            RevComp(rc_frame_word);

            // kmer = NuclKmer(rc_frame_word);
            for (int i = 0; i < dbg.kmer_k + 1; ++i) {
                seq[i] = dna_map[(int)rc_frame_word[i]];
            }
        }
        else {
            // kmer = NuclKmer(framed_word);
            for (int i = 0; i < dbg.kmer_k + 1; ++i) {
                seq[i] = dna_map[(int)framed_word[i]];
            }
        }

        AStarNode *starting_node_ptr = pool_->construct();
        AStarNode &starting_node = *starting_node_ptr;

        if (hmm.getAlphabet() == ProfileHMM::protein) {
            starting_node = AStarNode(NULL, starting_state + (framed_word.size() / 3), 'm');
            starting_node.length = framed_word.size() / 3;
        }
        else {
            starting_node = AStarNode(NULL, starting_state, 'm');
            starting_node.length = framed_word.size();
        }

        starting_node.fval = 0;
        starting_node.score = scoreStart(hmm, scoring_word, starting_state);
        starting_node.real_score = realScoreStart(hmm, scoring_word, starting_state);

        int64_t node_id = dbg.IndexBinarySearchEdge(seq);
        starting_node.node_id = node_id;

        return astarSearch(hmm, starting_node, dbg, forward, node_enumerator, goal_node, term_nodes);
    }

    bool astarSearch(ProfileHMM &hmm, AStarNode &starting_node, SuccinctDBG &dbg, bool forward, NodeEnumerator &node_enumerator,
                     AStarNode &goal_node, HashMapST<AStarNode, AStarNode> &term_nodes) {
        if (starting_node.state_no >= hmm.modelLength()) {
            goal_node = starting_node;
            fprintf(stderr, "\t-\t-\t-\t-\t-\tfalse\n");
            return true;
        }

        static const double log2 = log(2);

        priority_queue<AStarNodePtr> open;
        closed.clear();
        open_hash.clear();

        int opened_nodes = 1;

        int repeated_nodes = 0;
        int replaced_nodes = 0;
        int pruned_nodes = 0;

        //need to add a cache here
        HashMapST<AStarNode, AStarNode>::iterator got_term_node = term_nodes.find(starting_node);
        HashMapST<AStarNodePtr, AStarNodePtr>::iterator got;
        vector<AStarNode> temp_nodes_to_open;

        // printf("curr state: %c\n", starting_node.state);
        if (got_term_node == term_nodes.end()) {
            node_enumerator.enumerateNodes(temp_nodes_to_open, starting_node, forward, dbg);

            for (auto &next : temp_nodes_to_open) {
                auto next_ptr = AStarNodePtr(pool_->construct());
                next_ptr.get() = next;
                open.push(next_ptr);
            }
        }
        else {
            node_enumerator.enumerateNodes(temp_nodes_to_open, starting_node, forward, dbg, &got_term_node->second);

            for (auto &next : temp_nodes_to_open) {
                auto next_ptr = AStarNodePtr(pool_->construct());
                next_ptr.get() = next;
                open.push(next_ptr);
            }
        }

        if (open.empty()) {
            return false;
        }

        auto inter_goal_ptr = &starting_node;
        int closed_nodes = 0;

        while (!open.empty()) {
            auto curr_ptr = open.top();
            auto &curr = curr_ptr.get();
            closed_nodes++;
            // printf("%d: Curr: %lld, %d, %c\n", closed_nodes, curr.node_id, curr.state_no, curr.state);
            // printf("curr state: %c\n", curr.state);


            open.pop();

            // cout << "queue size = "<<open.size() << '\n';
            if (closed.find(curr_ptr) != closed.end()) {
                // puts("Skip");
                continue;
            }

            if (curr.state_no >= hmm.modelLength()) {
                curr.partial = 0;

                if ((curr.real_score + exit_probabilities[curr.length]) / log2
                        > (inter_goal_ptr->real_score + exit_probabilities[inter_goal_ptr->length]) / log2) {
                    inter_goal_ptr = &curr;
                }

                getHighestScoreNode(*inter_goal_ptr, goal_node);
                fprintf(stderr, "%d\t%zu\t%zu\t%d\t%d\t%d\tfalse\n", opened_nodes, open.size(), closed.size(), repeated_nodes, replaced_nodes, pruned_nodes);
                return true;
            }

            closed.insert(curr_ptr);

            if ((curr.real_score + exit_probabilities[curr.length]) / log2
                    > (inter_goal_ptr->real_score + exit_probabilities[inter_goal_ptr->length]) / log2) {
                inter_goal_ptr = &curr;
            }

            got_term_node = term_nodes.find(curr);

            if (got_term_node == term_nodes.end()) {
                node_enumerator.enumerateNodes(temp_nodes_to_open, curr, forward, dbg);
            }
            else {
                node_enumerator.enumerateNodes(temp_nodes_to_open, curr, forward, dbg, &got_term_node->second);
            }

            for (auto &next : temp_nodes_to_open) {
                auto next_ptr = AStarNodePtr(&next);
                bool open_node = false;

                if (heuristic_pruning > 0) {
                    if ((next.length < 5 || next.negative_count <= heuristic_pruning) && next.real_score > 0.0) {
                        got = open_hash.find(next_ptr);

                        if (got != open_hash.end()) {
                            repeated_nodes++;

                            if (got->second < next_ptr) {
                                replaced_nodes++;
                                open_node = true;
                            }
                        }
                        else {
                            open_node = true;
                        }
                    }
                    else {
                        pruned_nodes++;
                    }
                }
                else {
                    got = open_hash.find(next_ptr);

                    if (got != open_hash.end()) {
                        repeated_nodes++;

                        if (got->second < next_ptr) {
                            replaced_nodes++;
                            open_node = true;
                        }
                    }
                    else {
                        open_node = true;
                    }
                }

                if (open_node) {
                    next_ptr = AStarNodePtr(pool_->construct());
                    next_ptr.get() = next;
                    open_hash[next_ptr] = next_ptr;
                    opened_nodes++;
                    // printf("Next's discovered_from: %p\n", next.discovered_from);
                    open.push(next_ptr);
                }
            }
        }

        inter_goal_ptr->partial = 1;
        getHighestScoreNode(*inter_goal_ptr, goal_node);
        fprintf(stderr, "%d\t%zu\t%zu\t%d\t%d\t%d\ttrue\n", opened_nodes, open.size(), closed.size(), repeated_nodes, replaced_nodes, pruned_nodes);
        return true;
    }

    void getHighestScoreNode(AStarNode &inter_goal, AStarNode &goal_node) {
        AStarNode temp_goal = inter_goal;
        goal_node = inter_goal;

        while (temp_goal.discovered_from != NULL) {
            temp_goal = *temp_goal.discovered_from;

            if (temp_goal.real_score > goal_node.real_score) {
                goal_node = temp_goal;
            }
        }
    }

    void deleteAStarNodes() {
        pool_->clear();
    }

    char Comp(char c) {
        switch (c) {
        case 'A':
        case 'a':
            return 't';

        case 'C':
        case 'c':
            return 'g';

        case 'G':
        case 'g':
            return 'c';

        case 'T':
        case 't':
            return 'a';

        case 'N':
        case 'n':
            return 'n';

        case '-':
            return '-';

        default:
            assert(false);
        }
    }

    void RevComp(string &s) {
        reverse(s.begin(), s.end());

        for (unsigned i = 0; i < s.length(); ++i) {
            s[i] = Comp(s[i]);
        }
    }
};

#endif
