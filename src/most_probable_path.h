#ifndef MOST_PROBABLE_PATH_H__
#define MOST_PROBABLE_PATH_H__

#include "profile_hmm.h"
#include <vector>
#include <limits>
#include <assert.h>

using namespace std;

class MostProbablePath {
  private:
    ProfileHMM hmm_;
    vector<double> most_prob_from_state[3];

  public:
    MostProbablePath() {};
    MostProbablePath(ProfileHMM &hmm) : hmm_(hmm) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j <= hmm_.modelLength(); j++) {
                most_prob_from_state[i].push_back(- numeric_limits<double>::infinity());
            }
        }

        for (int i = 0; i <= hmm_.modelLength(); i++) {
            most_prob_from_state[0][i] = computeCostInternal('m', i);
            most_prob_from_state[1][i] = computeCostInternal('i', i);
            most_prob_from_state[2][i] = computeCostInternal('d', i);
        }
    };
    // ~MostProbablePath();
    double computeHeuristicCost(char state, int state_no) {
        switch (state) {
        case 'm':
            return most_prob_from_state[0][state_no];

        case 'i':
            return most_prob_from_state[1][state_no];

        case 'd':
            return most_prob_from_state[2][state_no];

        default:
            assert(false);
        }
    }
  private:
    double computeCostInternal(char pre_state, int state_no) {

        double h = 0;
        double match_trans = - numeric_limits<double>::infinity();
        double ins_trans = - numeric_limits<double>::infinity();
        double del_trans = - numeric_limits<double>::infinity();
        double best_match_prob = - numeric_limits<double>::infinity();
        double best_ins_prob = - numeric_limits<double>::infinity();

        for (int i = state_no + 1; i <= hmm_.modelLength(); i++) {
            switch (pre_state) {
            case 'm':
                match_trans = hmm_.tsc(i - 1, ProfileHMM::MM);
                ins_trans = hmm_.tsc(i - 1, ProfileHMM::MI);
                del_trans = hmm_.tsc(i - 1, ProfileHMM::MD);
                break;

            case 'd':
                match_trans = hmm_.tsc(i - 1, ProfileHMM::DM);
                ins_trans = - numeric_limits<double>::infinity();
                del_trans = hmm_.tsc(i - 1, ProfileHMM::DD);
                break;

            case 'i':
                match_trans = hmm_.tsc(i - 1, ProfileHMM::IM);
                ins_trans = hmm_.tsc(i - 1, ProfileHMM::II);
                del_trans = - numeric_limits<double>::infinity();
                break;

            default:
                assert(false);
            }

            //two more strange lines
            best_match_prob = - numeric_limits<double>::infinity();
            best_ins_prob = - numeric_limits<double>::infinity();

            for (int j = 0; j < hmm_.alphabetLength(); j++) {
                if (hmm_.msc(i, j) > best_match_prob) {
                    best_match_prob = hmm_.msc(i, j);
                }

                if (hmm_.isc(i, j) > best_ins_prob) {
                    best_ins_prob = hmm_.isc(i, j);
                }
            }

            match_trans += best_match_prob - hmm_.getMaxMatchEmission(i);
            del_trans -= hmm_.getMaxMatchEmission(i);
            ins_trans += best_ins_prob;

            //the incredible infinity step comes
            ins_trans = - numeric_limits<double>::infinity();

            if (ins_trans > match_trans && ins_trans > del_trans) {
                h += ins_trans;
                pre_state = 'i';
                i--;
            }
            else if (del_trans > match_trans && del_trans > ins_trans) {
                h += del_trans;
                pre_state = 'd';
            }
            else {
                h += match_trans;
                pre_state = 'm';
            }
        }

        return h;
    }


};

#endif