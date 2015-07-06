#include <stdint.h>
#include <string>
#include <stdexcept>

class ProtKmerGenerator {
	private:
		string bases_;
		int k_;
		Kmer next_;

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
				throw std::logic_error("Sequence length is less than the kmer length");
			}

			bases_ = seq;
			k_ = k;
			model_only_ = model_only;
			index_ = 0;
			position_ = 1;
			next_ = getFirstKmer(0);
		}

		bool hasNext() {
			return next_ != NULL;
		}

		Kmer next() {
			Kmer *ret = &next_;
			cur_model_position_ = position_;
			findNextKmer(k-1);
			return ret;
		}

	private:
		Kmer getFirstKmer(int klength) {
			
		}
}