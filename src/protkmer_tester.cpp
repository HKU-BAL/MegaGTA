#include <iostream>
#include <string>
#include <bitset>
#include "prot_kmer.h"
#include "prot_kmer_generator.h"

using namespace std;

int main(){
	string seq = "..................................................................................................................................................................................................................................MAIKKYKPTS.NGRRGMTV.L..DF.SE...ITTDQ...............PEKS......LLA.PL..K..K.K..AGRN.N.QGKITVRHQ.GGGHKRQYRIIDF...KR.....D.KD..........GI..P.....................................................G.....R...VATIEYDPNRSANIALI.N.Y....A.D....G....E...K..............................R..............................Y.........ILA..........PK.NLKVGMEI...M...SG.................................................................P.NA...D...I......KV..........GNALPLE.............NIPVGTLVHNI...ELKPG....R..G........G.QLV..RAAGT.SAQVLGK...........EG..........................................................................KYVIIRLASGEVRMILGKCRATVGEVGNEQHELV..NIGKAGRARWL.GIRPT...VRGSVM..NPVDHPHGG.GE.......GKA..P.I.GR..kSPMTP..WG.KP...TL.G.YKTRKK..KNKSDKFI..IRRRKK-................................................................................................................................................................................................................................................................";

	ProtKmerGenerator kmers = ProtKmerGenerator(seq, 10 , true);

	while (kmers.hasNext()) {
		ProtKmer temp = kmers.next();
		cout << "Kmer = " << temp.decodePacked() << " position = " <<kmers.getPosition() << endl;
	}

	char kmer_str[] = "ARNDCQEGHI";
	ProtKmer kmer1 = ProtKmer(kmer_str);
	cout << kmer1.hash() << endl;
	bitset<64> x(kmer1.kmers[0]);
	bitset<64> y(kmer1.kmers[1]);
	cout << x << endl;
	cout << y << endl;
	cout << kmer1.decodePacked() << endl;

	char kmer_str2[] = "ARNDCQEGHIARNDC";
	ProtKmer kmer2 = ProtKmer(kmer_str2);
	cout << kmer2.hash() << endl;
	bitset<64> a(kmer2.kmers[0]);
	bitset<64> b(kmer2.kmers[1]);
	cout << a << endl;
	cout << b << endl;
	cout << kmer2.decodePacked() << endl;

	return 0;
}