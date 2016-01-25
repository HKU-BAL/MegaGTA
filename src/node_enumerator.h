#ifndef NODE_ENUMERATOR_H__
#define NODE_ENUMERATOR_H__

#include "profile_hmm.h"
#include "most_probable_path.h"
#include "succinct_dbg.h"
#include "a_star_node.h"
#include <iostream>
#include <vector>
#include "codon.h"
#include <set>
#include <cmath>
#include <algorithm>

using namespace std;

#define MIN3(x, y, z) std::min(x, std::min(y, z))

class NodeEnumerator {
  private:
    static const int SCALE = 10000;
    static constexpr double hweight = 2.0;
    static constexpr double kLowCovSibling = 0; // -log(0.5);
    uint8_t next_nucl;
    char emission;
    double match_trans;
    double ins_trans;
    double del_trans;
    // NuclKmer next_kmer;
    int next_state;
    ProfileHMM *hmm;
    bool prot_search;
    MostProbablePath *hcost;
    AStarNode next;
    double low_cov_penalty;

  public:
    NodeEnumerator(ProfileHMM &_hmm, MostProbablePath &_hcost, double low_cov_penalty) {
        hmm = &_hmm;
        prot_search = _hmm.getAlphabet() == ProfileHMM::protein;
        hcost = &_hcost;
        low_cov_penalty = -log(low_cov_penalty);
    };
    ~NodeEnumerator() {};
    void enumerateNodes(vector<AStarNode> &ret, AStarNode &curr, bool forward, SuccinctDBG &dbg) {
        enumerateNodes(ret, curr, forward, dbg, NULL);
    }

    int calLowCovSibling(SuccinctDBG &dbg, int64_t a[4], int n) {
        int isMulti1 = 0;
        int numHighCov = 0;

        for (int i = 0; i < n; ++i) {
            if (dbg.IsMulti1(a[i])) {
                isMulti1 |= 1 << i;
            }
            else {
                numHighCov++;
            }
        }

        return numHighCov > 0 ? isMulti1 : 0;
    }

    void enumerateNodes(vector<AStarNode> &ret, AStarNode &curr, bool forward, SuccinctDBG &dbg, AStarNode *child_node) {
        ret.clear();
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

        if (curr.node_id == -1) {
            return;
        }
        else {
            int64_t next_node[4], next_node_2[4], next_node_3[4];
            vector<int64_t> packed_codon;
            int outd = dbg.OutgoingEdges(curr.node_id, next_node);
            int lowCovSibling = calLowCovSibling(dbg, next_node, outd);

            for (int i = 0; i < outd; ++i) {
                int outd_2 = dbg.OutgoingEdges(next_node[i], next_node_2);
                int lowCovSibling2 = calLowCovSibling(dbg, next_node_2, outd_2);

                for (int j = 0; j < outd_2; ++j) {
                    int outd_3 = dbg.OutgoingEdges(next_node_2[j], next_node_3);
                    int lowCovSibling3 = calLowCovSibling(dbg, next_node_3, outd_3);

                    for (int k = 0; k < outd_3; ++k) {
                        int64_t packed = (next_node_3[k] << 16) |
                                         ((dbg.GetEdgeOutLabel(next_node[i]) - 1) << 6) |
                                         ((dbg.GetEdgeOutLabel(next_node_2[j]) - 1) << 3) |
                                         (dbg.GetEdgeOutLabel(next_node_3[k]) - 1);
                        packed |= (dbg.IsMulti1(next_node[i]) &&
                                   dbg.IsMulti1(next_node_2[j]) &&
                                   dbg.IsMulti1(next_node_3[k])) << 9;
                        packed_codon.push_back(packed);

                        if ((lowCovSibling & (1 << i)) ||
                                (lowCovSibling2 & (1 << j)) ||
                                (lowCovSibling3 & (1 << k))) {
                            packed |= 1 << 10;
                        }
                    }
                }
            }

            //translate to aa
            for (int i = 0; i < (int)packed_codon.size(); ++i) {
                int64_t packed = packed_codon[i];

                // next_kmer = curr.kmer.shiftLeftCopy(codons[i][0], codons[i][1], codons[i][2]);
                if (!forward) {
                    emission = Codon::rc_codonTable[packed >> 6 & 7][packed >> 3 & 7][packed & 7];
                }
                else {
                    emission = Codon::codonTable[packed >> 6 & 7][packed >> 3 & 7][packed & 7];
                }

                if (emission == '*') {
                    continue;
                }

                if (child_node != NULL && child_node->node_id != (packed >> 16)) {
                    continue;
                }

                double lowCovPenalty = (packed & (1 << 9)) ? low_cov_penalty : 0;

                if (packed & (1 << 10)) lowCovPenalty += kLowCovSibling;

                next = AStarNode(&curr, next_state, 'm');

                next.real_score = curr.real_score + (match_trans + hmm->msc(next_state, emission)) - lowCovPenalty;

                if (next.real_score >= curr.max_score) {
                    next.max_score = next.real_score;
                    next.negative_count = 0;
                }
                else {
                    next.max_score = curr.max_score;
                    next.negative_count = curr.negative_count + 1;
                }

                next.nucl_emission = packed & ((1 << 9) - 1);

                next.emission = emission;
                double this_node_score = (match_trans + hmm->msc(next_state, emission)) - lowCovPenalty - max_match_emission;
                next.length = curr.length + 1;
                next.score = curr.score + this_node_score;
                next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('m', next_state)));
                next.indels = curr.indels;

                next.node_id = packed >> 16;

                if (child_node != NULL && *child_node == next) {
                    ret.push_back(next);
                    return;
                }
                else {
                    ret.push_back(next);
                }

                // ret.push_back(next);

                //insert node
                if (curr.state != 'd') {
                    next = AStarNode(&curr, curr.state_no, 'i');

                    next.real_score = curr.real_score + (ins_trans + hmm->isc(next_state, emission)) - lowCovPenalty;
                    next.max_score = curr.max_score;
                    next.negative_count = curr.negative_count + 1;

                    next.nucl_emission = packed & ((1 << 9) - 1);

                    next.emission = emission;
                    this_node_score = (ins_trans + hmm->isc(next_state, emission)) - lowCovPenalty;
                    next.length = curr.length + 1;
                    next.score = curr.score + this_node_score;
                    next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('i', curr.state_no)));
                    next.indels = curr.indels + 1;

                    next.node_id = packed >> 16;

                    if (child_node != NULL && *child_node == next) {
                        ret.push_back(next);
                        return;
                    }
                    else {
                        ret.push_back(next);
                    }
                }
            }

            //delete node
            if (curr.state != 'i') {
                next = AStarNode(&curr, next_state, 'd');

                next.real_score = curr.real_score + del_trans;
                next.max_score = curr.max_score;
                next.negative_count = curr.negative_count + 1;

                next.nucl_emission = (4 << 6) | (4 << 3) | 4;
                next.emission = '-';
                double this_node_score = del_trans - max_match_emission;
                next.length = curr.length;
                next.score = curr.score + this_node_score;
                next.fval = (int) (SCALE * (next.score + hweight * hcost->computeHeuristicCost('d', next_state)));
                next.indels = curr.indels + 1;

                next.node_id = curr.node_id;

                if (child_node != NULL && *child_node == next) {
                    ret.push_back(next);
                    return;
                }
                else {
                    ret.push_back(next);
                }

                // ret.push_back(next);
            }
        }
    }

};

#endif
