#include "a_star_node.h"
#include "nucl_kmer.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
    NuclKmer::setUp();
    NuclKmer kmer;
    long hash = 1000;
    int model = 1, model_2 = 2, model_3 = 3;
    char state = 'm';
    AStarNode node = AStarNode(NULL, kmer, model, state);
    AStarNode node_2 = AStarNode(&node, kmer, model_2, state);
    AStarNode node_3 = AStarNode(&node_2, kmer, model_3, state);
    AStarNode node_copy = node_2;
    cout << "node_2 discovered_from " << node_2.discovered_from->state_no << '\n';
    cout << "node_copy discovered_from " << node_copy.discovered_from->state_no << '\n';

    while (node_3.discovered_from != NULL) {
        cout << "state = " << node_3.state_no << '\n';
        node_3 = *node_3.discovered_from;
        cout << "not null " << endl;
    }
}