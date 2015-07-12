#include "kmer_1.h"
#include <math.h>
#include <string.h>
#include <iostream>

#include <bitset>


class ProtKmer : public Kmer {
	public:
		ProtKmer(char char_kmer[]) : Kmer(0x5, 0x1F, 12, strlen(char_kmer)) {
			setUp();
			initialize(char_kmer);
		}

	// public:
	// 	ProtKmer(const ProtKmer &k) {
	// 		Kmer::k = k.k;
	// 		Kmer::l = k.l;
	// 		memcpy(Kmer::kmers, k.kmers, sizeof(k.kmers));
	// 		Kmer::last_fill = k.last_fill;
	// 		Kmer::last_mask = k.last_mask;
	// 		setUp();
	// 	}

	private:
		uint8_t ascii_map[127];
		char int_to_char[127];

	private:
		void setUp() {

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
			ascii_map['*'] = 10;
			int_to_char[10] = '*';
		}

	public:
		virtual void shiftRight(char c) {
			shiftRight(charToByte(c));
		}
		virtual void shiftLeft(char c) {
			shiftLeft(charToByte(c));
		}
		virtual uint8_t charToByte(char c) {
			return ascii_map[(int)c];
		}
		virtual char intToChar(int i) {
			return int_to_char[i];
		}
};