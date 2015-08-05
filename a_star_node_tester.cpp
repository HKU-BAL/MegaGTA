#include "a_star_node.h"
#include "nucl_kmer.h"

using namespace std;

int main(int argc, char **argv) {
	NuclKmer::setUp();
	NuclKmer kmer;
	long hash = 1000;
	int model = 1;
	char state = 'm';
	AStarNode node = AStarNode(NULL, kmer, hash, hash, model, state);
}