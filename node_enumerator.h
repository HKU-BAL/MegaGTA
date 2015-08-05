#include "profile_hmm.h"
#include "most_probable_path.h"
#include "nucl_kmer.h"
#include <vector>

using namespace std;

class NodeEnumerator
{
private:
	const double hweight = 1.0;
	AStarNode *next;
	char emission;
	double match_trans;
	double ins_trans;
	double del_trans;
	NuclKmer *next_kmer;
	ProfileHMM *hmm;
	bool prot_search;
	MostProbablePath *hcost;

public:
	NodeEnumerator(ProfileHMM &_hmm, MostProbablePath &_hcost, double _hweight) {
		hmm = &_hmm;
		prot_search = _hmm.getAlphabet() == ProfileHMM::protein;
		hcost = &_hcost;
		hweight = _hweight;
	};
	~NodeEnumerator();
	vector<AStarNode> enumeratorNodes(AStarNode &curr) {
		
	}


	
};