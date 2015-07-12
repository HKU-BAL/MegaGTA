#include <iostream>
#include <string>
#include <bitset>
#include "prot_kmer.h"

using namespace std;

int main(){
	char kmer_str[] = "ARNDCQEGHI";
	ProtKmer kmer1 = ProtKmer(kmer_str);
	cout << kmer1.hashCode() << endl;
	bitset<64> x(kmer1.kmers[0]);
	bitset<64> y(kmer1.kmers[1]);
	// cout << kmer1.kmers[0] << endl;
	cout << x << endl;
	cout << y << endl;


	return 0;
}