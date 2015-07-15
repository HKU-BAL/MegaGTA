#ifndef PROT_KMER_H__
#define PROT_KMER_H__

#include "kmer_1.h"
#include <math.h>
#include <string.h>
#include <algorithm>
#include <bitset>


class ProtKmer : public Kmer {
	public:
		ProtKmer() : Kmer() {setUp();}
		ProtKmer(const string &char_kmer) : Kmer(0x5, 0x1F, 12, char_kmer.size()) {
			setUp();
			initialize(char_kmer);
		}

	public:
		uint8_t ascii_map[127];
	private:
		char int_to_char[127];

	private:
		void setUp() {
			fill_n(ascii_map, 127, 31);
			for (int i = 0; i < 20; ++i) {
				// ARNDCQEGHI
				// LKMFPSTWYV
				ascii_map[(int)("ARNDCQEGHIarndcqeghi"[i])] = (uint8_t) ("01234567890123456789"[i] - '0'); 
				ascii_map[(int)("LKMFPSTWYVlkmfpstwyv"[i])] = (uint8_t) ("01234567890123456789"[i] - '0' + 10);
				if (i <= 9){
					int_to_char[(int)("0123456789"[i] - '0')] = "arndcqeghi"[i];
					int_to_char[(int)("01234567890123456789"[i] - '0' + 10)] = "lkmfpstwyv"[i];
				}				
			}
			ascii_map['*'] = 20;
			int_to_char[20] = '*';
		}

	public:
		virtual void shiftRight(char c) {
			Kmer::shiftRight(charToByte(c));
		}
		virtual void shiftLeft(char c) {
			Kmer::shiftLeft(charToByte(c));
		}
		virtual uint8_t charToByte(char c) {
			return ascii_map[(int)c];
		}
		virtual char intToChar(int i) {
			return int_to_char[i];
		}
};

#endif