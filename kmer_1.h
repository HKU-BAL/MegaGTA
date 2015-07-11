#include <stdint.h>
#include <string>
#include <stdexcept>

using namespace std;

class Kmer {
	public:
		static const int MAX_NUCL_KMER_SIZE = 64;
		static const int MAX_PROT_KMER_SIZE = 24;

	protected:
		uint64_t kmers[2];
		int k;
		int l;
		int last_fill;
		uint64_t last_mask; 

		//The following three default are for nucleotides
		int bits_to_shift = 0x2; //2 bits
		int char_mask = 0x3; // 11
		int items_per_bucket = 32;

	public:
		virtual uint8_t charToByte(char c);
		virtual char intToChar(int i);
		virtual void shiftRight(char c);
		virtual void shiftLeft(char c);
		// virtual void shiftRight(uint8_t b);
		// virtual void shiftLeft(uint8_t b);

		int length() {
			return k;
		}
		int packedLength() {
			return l;
		}
		long getPart(int index) {
			return kmers[index];
		}
		long* getPackedKmers() {
			return kmers;
		}

		void initialize(char char_kmer[]) {
			uint8_t b;

			for (int index = 0; index < sizeof(char_kmer); ++index) {
				b = charToByte(char_kmer[index]);
				if (b == -1) {
					throw std::invalid_argument("Kmer contains one or more invalid bases");
				}

				if (index < items_per_bucket) {
					kmers[0] = (kmers[0] << bits_to_shift) | (b & char_mask);
				} else {
					kmers[1] = (kmers[1] << bits_to_shift) | (b & char_mask);
				}
			}
		}

	protected:
		void shiftRight(uint8_t b) {
			if (l == 2) {
				uint8_t overflow = (uint8_t) (kmers[0] & char_mask);
				kmers[0] = (kmers[0] >> bits_to_shift) | ((uint64_t) b << ( (items_per_bucket-1) * bits_to_shift));

				kmers[1] = (kmers[1] >> bits_to_shift) | ((uint64_t) overflow << ((last_fill - 1) * bits_to_shift));
				kmers[1] &= last_mask;
			} else {
				kmers[0] = (kmers[0] >> bits_to_shift) | ((uint64_t) b << ((last_fill - 1) * bits_to_shift));
				kmers[0] &= last_mask;
			}
		}

		void shiftLeft(uint8_t b) {
			if (l == 2) {
				uint8_t overflow = (uint8_t) ((kmers[1] >> ((last_fill - 1) * bits_to_shift)) & char_mask);
				kmers[1] = (kmers[1] << bits_to_shift) | ((uint64_t) b & char_mask);
				kmers[1] &= last_mask;

				kmers[0] = (kmers[0] << bits_to_shift) | (overflow);
			} else {
				kmers[0] = (kmers[0] << bits_to_shift) | ((uint64_t) b & char_mask);
				kmers[0] &= last_mask;
			}
		}

	public:
		string decodePacked(uint64_t in_kmer[]){
			string buf;
			for (int seg = l; seg > 0; seg--) {
				int index = items_per_bucket;
				if (seg == 1) {
					index = last_fill;
				}
				uint64_t this_k = (seg == 2)? in_kmer[1] : in_kmer[0];
				for (; index > 0; index--) {
					string.append( intToChar((int) (this_k & char_mask)) );
					this_k = this_k >> bits_to_shift;
				}
			}
			return buf.reverse();
		}

		bool equals(const Kmer &kmer) const {
			if (kmers[0] != kmer.kmers[0] || kmers[1] != kmer.kmers[1]) {
				return false;
			}
			return true;
		}

		uint64_t hashCode() {
			return CityHash64((const char*)kmers, sizeof(kmers[0]) * 2);
		}
};

