#include "NTSequence.h"
#include "AASequence.h"

#include <iterator>
#include <fstream>

using namespace seq;

int main(int argc, char **argv)
{
  std::ifstream s(argv[1]);

  /*
   * Iterate over all nucleotide sequences in the file.
   */
  try {
    for (std::istream_iterator<NTSequence> i(s);
	 i != std::istream_iterator<NTSequence>();
	 ++i) {
      const NTSequence& seq = *i;

      AASequence aa = AASequence::translate(seq.begin(),
					    seq.begin() + (seq.size() / 3) * 3);
      aa.setName(seq.name());
      aa.setDescription("translated " + seq.description());

      // write out a FASTA file entry
      std::cout << aa;
    }
  } catch (ParseException& e) {
    std::cerr << "Error reading " << argv[1] << ": "
	      << e.message() << std::endl;
  }
}