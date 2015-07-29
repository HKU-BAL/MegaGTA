#ifndef PROT_KMER_GENERATOR_H__
#define PROT_KMER_GENERATOR_H__


#include <stdint.h>
#include <string>
#include <stdexcept>
#include <ctype.h>
#include "kmer_1.h"
#include "prot_kmer.h"
#include <string.h>


class ProtKmerGenerator {
	private:
		std::string bases_;
		int k_;
		ProtKmer next_;
		bool has_next;

		int index_;
		int position_;
		int cur_model_position_;
		bool model_only_ = false;

	public:
		ProtKmerGenerator() {}

		ProtKmerGenerator(const std::string &seq, int k) : ProtKmerGenerator(seq, k, false) {}

		ProtKmerGenerator(const std::string &seq, int k, bool model_only) {
			if (k > Kmer::MAX_PROT_KMER_SIZE) {
				throw std::invalid_argument("K-mer size cannot be larger than 24");
			}

			if (seq.length() < k) {
				// throw std::invalid_argument("Sequence length is less than the kmer length");
				has_next = false;
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

		ProtKmer next() {
			ProtKmer ret = next_;
			cur_model_position_ = position_;
			findNextKmer(k_ -1);
			return ret;
		}

	private:
		bool getFirstKmer(int klength) {
			std::string kmer_str(k_, '\0');
			while (index_ < bases_.length()) {
				char base = bases_[index_++];

				if (model_only_ && (islower(base) || base == '-' || base == 'X' || base == 'x')) {
					if (base == '-' || base == 'X') {
						position_++;
					}
					klength = 0;
				} else {
					if (!model_only_ || (model_only_ && (base != '.' && next_.ascii_map[base] != 31 && base != '*'))) {
						if (next_.ascii_map[base] == 31 && base != 'X') {
							throw std::domain_error("Unknown prot base");
						}
						kmer_str[klength] = base;
						position_++;
						klength++;
					}

					if (base == 'X') {
						klength = 0;
					}

					if (klength == k_) {
						cur_model_position_ = position_;
						next_ = ProtKmer(kmer_str);
						return true;
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
					if (!model_only_ || (model_only_ && (base != '.' && next_.ascii_map[base] != 31 && base != '*'))) {
						if (next_.ascii_map[base] == 31 && base != 'X') {
							throw std::invalid_argument("Unknown prot base");
						}
						if (base != 'X') {
							next_.shiftLeft(base);
							klength++;
						} else {
							klength = 0;
						}
						position_++;						
					}

					if (klength == k_) {
						return;
					}
				}
			}

			if (klength != k_) {
				has_next = false;
			}
		}
	public:
		int getPosition() {
			return cur_model_position_ - k_;
		}


};

#endif