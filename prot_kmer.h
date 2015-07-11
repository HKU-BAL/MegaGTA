#include "kmer_1.h"
#include <math.h>

class ProtKmer : public Kmer {
	public:
		ProtKmer(char char_kmer[]) {
			ProtKmer(char_kmer.size());
			initialize(char_kmer);
		}

	private:
		ProtKmer(int k){
			if (k > MAX_PROT_KMER_SIZE) {
				throw std::invalid_argument("k-mer lenght must be <= " + MAX_PROT_KMER_SIZE);
			}

			setUp();
			Kmer::k = k;
			l = (int) ceil(k / ((double)items_per_bucket));
			last_fill = items_per_bucket - (items_per_bucket * l - k);
			uint64_t tmp = 0;
			for (int index = 0; index < last_fill; index++) {
				tmp = (tmp << bits_to_shift) | char_mask;
			}
			last_mask = tmp;
		}

	public:
		ProtKmer(ProtKmer k) {
			Kmer::k = k.k;
			Kmer::l = k.l;
			memcpy(Kmer::kmers, k.kmers, k.kmers.size());
			Kmer::last_fill = k.last_fill;
			Kmer::last_mask = k.last_mask;
			setUp();
		}

	private:
		uint8_t ascii_map[127];
		char int_to_char[127];

	private:
		void setUp() {
			bits_to_shift = 0x5;
			char_mask = 0x1F;
			items_per_bucket = 12;

			for (int i = 0; i < 20; ++i) {
				// ARNDCQEGHI
				// LKMFPSTWYV
				ascii_map[(int)("ARNDCQEGHIarndcqeghi"[i])] = (uint8_t) ("01234567890123456789"[i] - '0'); 
				ascii_map[(int)("LKMFPSTWYVlkmfpstwyv"[i])] = (uint8_t) ("01234567890123456789"[i] - '0' + 10);
				if (i <= 10){
					int_to_char[(int)("0123456789"[i] - '0')] = "arndcqeghi"[i];
					int_to_char[(int)("01234567890123456789"[i] - '0' + 10)] = "lkmfpstwyv"[i];
				}				
			}
			ascii_map['*'] = 20;
			int_to_char[20] = '*';
		}

	public:
		void shiftRight(char c) {
			shiftRight(charToByte(c));
		}
		void shiftLeft(char c) {
			shiftLeft(charToByte(c));
		}
		uint8_t charToByte(char c) {
			return ascii_map[(int)c];
		}
		char intToChar(int i) {
			return int_to_char[i];
		}
};