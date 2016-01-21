/*
 *  MEGAHIT
 *  Copyright (C) 2014 The University of Hong Kong
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include "branch_group.h"

bool BranchGroup::Search() {
    if (!sdbg_->IsValidEdge(begin_node_)) return false;
    int outd = sdbg_->EdgeOutdegree(begin_node_);
    if (outd <= 1 || outd > max_branches_) {
        return false;
    }

    branches_.push_back(BranchRecord(1, begin_node_));
    multiplicities_.push_back(0);

    bool converged = false;

    for (int j = 1; j < max_length_; ++j) {
        int num_branches = branches_.size();
        for (int i = 0; i < num_branches; ++i) {
            int64_t current = branches_[i].back();
            int64_t outgoings[4];
            int out_degree = sdbg_->OutgoingEdges(current, outgoings);

            if (out_degree >= 1) {
                // append this node the last of current branch
                branches_[i].push_back(outgoings[0]);
                multiplicities_[i] += sdbg_->EdgeMultiplicity(outgoings[0]);

                if ((int)branches_.size() + out_degree - 1 > max_branches_) {
                    // too many branches
                    return false;
                } else {
                    BranchRecord curr_branch = branches_[i];
                    int curr_branch_multiplicity = multiplicities_[i] - sdbg_->EdgeMultiplicity(outgoings[0]);
                    for (int x = 1; x < out_degree; ++x) {
                        curr_branch.pop_back();
                        curr_branch.push_back(outgoings[x]);
                        branches_.push_back(curr_branch);
                        multiplicities_.push_back(curr_branch_multiplicity + sdbg_->EdgeMultiplicity(outgoings[x]));
                    }
                }
            }
        }

        // check whether all branches's last nodes are coming from this branch group
        for (unsigned i = 0; i < branches_.size(); ++i) {
            int64_t last_node = branches_[i].back();
            int64_t incomings[4];
            int in_degree = sdbg_->IncomingEdges(last_node, incomings);

            if (in_degree == 1) {
                continue;
            } else {
                for (int x = 0; x < in_degree; ++x) {
                    bool exist_in_group = false;
                    for (auto it = branches_.begin(); it != branches_.end(); ++it) {
                        if ((*it)[j - 1] == incomings[x]) {
                            exist_in_group = true;
                            break;
                        }
                    }
                    if (!exist_in_group) {
                        return false;
                    }
                }
            }
        }

        // check converge
        end_node_ = branches_[0].back();
        if (sdbg_->EdgeOutdegree(end_node_) == 1) {
            converged = true;
            for (unsigned i = 1; i < branches_.size(); ++i) {
                if (branches_[i].back() != end_node_) {
                    converged = false;
                    break;
                }
            }
            if (converged) {
                break;
            }
        }
    }

    return converged && begin_node_ != end_node_;
}

bool BranchGroup::Pop(AtomicBitVector &marked) {
    int best_multiplicity = multiplicities_[0];
    int best_path = 0;
    std::vector<int64_t> locked_nodes;

    for (unsigned i = 1; i < branches_.size(); ++i) {
        int curr_multiplicity = multiplicities_[i];
        if (curr_multiplicity >= best_multiplicity) {
            best_path = i;
            best_multiplicity = curr_multiplicity;
        }
    }

    for (unsigned i = 0; i < branches_.size(); ++i) {
        for (unsigned j = 1; j + 1 < branches_[i].size(); ++j) {
            if (!marked.try_lock(branches_[i][j])) {
                for (auto it = locked_nodes.begin(); it != locked_nodes.end(); ++it) {
                    marked.unset(*it);
                }
                return false;
            }
            locked_nodes.push_back(branches_[i][j]);
            sdbg_->SetInvalidEdge(branches_[i][j]);
        }
    }

    for (unsigned j = 1; j + 1 < branches_[best_path].size(); ++j) {
        sdbg_->SetValidEdge(branches_[best_path][j]);
    }

    for (unsigned j = 1; j + 1 < branches_[best_path].size(); ++j) {
        marked.unset(branches_[best_path][j]);
    }

    return true;
}
