#include <iostream>
#include <string>
#include <bitset>
#include "nucl_kmer_generator.h"

using namespace std;

int main() {
    string seq = "ATGGCCGTCAAAAAGTACCGTCCCTATACCCCCAATGGCCGTCAAAAAGTACCGTCCCTATACCCA";
    NuclKmerGenerator kmers = NuclKmerGenerator(seq, 45, false);

    while (kmers.hasNext()) {
        NuclKmer temp = kmers.next();
        cout << "Kmer = " << temp.decodePacked() << " position = " << kmers.getPosition() << endl;
    }

    char kmer_str[] = "ATGGCCGTCAAAAAGTACCGTCCCTATACCCC";
    NuclKmer kmer1 = NuclKmer(kmer_str);
    cout << kmer1.hash() << endl;
    bitset<64> x(kmer1.kmers[0]);
    bitset<64> y(kmer1.kmers[1]);
    cout << x << endl;
    cout << y << endl;
    cout << kmer1.decodePacked() << endl;

    char kmer_str2[] = "ATGGCCGTCAAAAAGTACCGTCCCTATACCCA";
    NuclKmer kmer2 = NuclKmer(kmer_str2);
    cout << kmer2.hash() << endl;
    bitset<64> a(kmer2.kmers[0]);
    bitset<64> b(kmer2.kmers[1]);
    cout << a << endl;
    cout << b << endl;
    cout << kmer2.decodePacked() << endl;

    return 0;
}
