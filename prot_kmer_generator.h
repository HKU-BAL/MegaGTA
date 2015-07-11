#include <stdint.h>
#include <string>
#include <stdexcept>
#include <ctype.h>

class ProtKmerGenerator {
	private:
		string bases_;
		int k_;
		Kmer next_;
		bool has_next;

		int index_;
		int position_;
		int cur_model_position_;
		bool model_only_ = false;

	public:
		ProtKmerGenerator(string &seq, int k) {
			ProtKmerGenerator(seq, k, false);
		}

		ProtKmerGenerator(string &seq, int k, bool model_only) {
			if (k > Kmer.max_prot_kmer_size) {
				throw std::invalid_argument("K-mer size cannot be larger than 24");
			}

			if (seq.length() < k) {
				throw std::invalid_argument("Sequence length is less than the kmer length");
			}

			bases_ = seq;
			k_ = k;
			model_only_ = model_only;
			index_ = 0;
			position_ = 1;
			has_next = getFirstKmer(0);
		}

		bool hasNext() {
			return has_next;
		}

		Kmer next() {
			// Kmer ret = next_;
			cur_model_position_ = position_;
			findNextKmer(k-1);
			return next_;
		}

	private:
		bool getFirstKmer(int klength) {
			char kmer_str[k_];
			while (index_ < bases_.length()) {
				char base = bases_[index_++];

				if (model_only_ && (islower(base) || base == '-' || base == 'X' || base == 'x')) {
					if (base == '-' || base == 'X') {
						position_++;
					}
					klength = 0;
				} else {
					if (!model_only_ || (model_only_ && (base != '.' && proteinAlphabet.contains(base) && base != '*'))) {
						if (ProtBinMapping.asciiMap[base] == -1) {
							throw std::invalid_argument("Unknown prot base" + base);
						}
						kmer_str[klength] = base;
						position_++;
						klength++;
					}

					if (klength == k_) {
						cur_model_position_ = position_;
						next_ = ProtKmer(kmer_str);
						return true;
						// return new ProtKmer(kmer_str);
					}
				}
			}
			return false;
		}

		void findNextKmer(int klength) {
			if (has_next == false) {
				return;
			}

			while (index_ < bases_.length()) {
				char base = bases_[index_++];

				if (model_only_ && (islower(base) || base == '-' || base == 'X' || base == 'x')) {
					if (base == '-' || base == 'X') {
						position_++;
					}
					klength = 0;
				} else {
					if (!model_only_ || (model_only_ && (base != '.' && proteinAlphabet.contains(base) && base != '*'))) {
						if (ProtBinMapping.asciiMap[base] == -1) {
							throw std::invalid_argument("Unknown prot base" + base);
						}
						next_ = next_.shiftLeft(base);
						position_++;
						klength++;
					}

					if (klength == k_) {
						return;
					}
				}
			}

			if (klength != k) {
				has_next = false;
			}
		}
	public:
		int getPosition() {
			return cur_model_position_ - k_;
		}


}