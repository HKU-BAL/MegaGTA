#ifndef NUCL_KMER_H__
#define NUCL_KMER_H__

#include "kmer.h"
#include <math.h>
#include <string.h>
#include <algorithm>


class NuclKmer : public Kmer {
	public:
		NuclKmer() : Kmer() {}
		NuclKmer(const string &char_kmer) : Kmer(0x2, 0x3, 32, char_kmer.size()) {
			// setUp();
			initialize(char_kmer);
		}

	public:
		static uint8_t ascii_map[127];
	private:
		static char int_to_char[127];

	public:
		static void setUp() {
			fill_n(ascii_map, 127, 31);
			for (int i = 0; i < 10; ++i) {
				ascii_map[(int)("ACGTNacgtn"[i])] = (uint8_t) ("0123201232"[i] - '0');
				if (i <= 3) {
					int_to_char[(int)("0123"[i] - '0')] = "acgt"[i];
				}				
			}
		}

	public:
		virtual void shiftRight(char c) {
			Kmer::shiftRight(charToByte(c));
		}
		virtual void shiftLeft(char c) {
			Kmer::shiftLeft(charToByte(c));
		}
		NuclKmer shiftLeftCopy(uint8_t b) {
			NuclKmer ret = *this;
			ret.Kmer::shiftLeft(b);
			return ret;
		}
		NuclKmer shiftLeftCopy(uint8_t b1, uint8_t b2, uint8_t b3) {
			NuclKmer ret = *this;
			ret.Kmer::shiftLeft(b1);
			ret.Kmer::shiftLeft(b2);
			ret.Kmer::shiftLeft(b3);
			return ret;
		}
		virtual uint8_t charToByte(char c) {
			return ascii_map[(int)c];
		}
		virtual char intToChar(int i) {
			return int_to_char[i];
		}
};

#endif