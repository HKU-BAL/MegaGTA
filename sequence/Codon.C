#include "Codon.h"

namespace seq {

AminoAcid Codon::translate(const NTSequence::const_iterator triplet)
{
  const AminoAcid codonTable[4][4][4] = {
  { { AminoAcid::K /* AAA */,
      AminoAcid::N /* AAC */,
      AminoAcid::K /* AAG */,
      AminoAcid::N /* AAT */
    },
    { AminoAcid::T /* ACA */,
      AminoAcid::T /* ACC */,
      AminoAcid::T /* ACG */,
      AminoAcid::T /* ACT */
    },
    { AminoAcid::R /* AGA */,
      AminoAcid::S /* AGC */,
      AminoAcid::R /* AGG */,
      AminoAcid::S /* AGT */
    },
    { AminoAcid::I /* ATA */,
      AminoAcid::I /* ATC */,
      AminoAcid::M /* ATG */,
      AminoAcid::I /* ATT */
    }
  },
  { { AminoAcid::Q /* CAA */,
      AminoAcid::H /* CAC */,
      AminoAcid::Q /* CAG */,
      AminoAcid::H /* CAT */
    },
    { AminoAcid::P /* CCA */,
      AminoAcid::P /* CCC */,
      AminoAcid::P /* CCG */,
      AminoAcid::P /* CCT */
    },
    { AminoAcid::R /* CGA */,
      AminoAcid::R /* CGC */,
      AminoAcid::R /* CGG */,
      AminoAcid::R /* CGT */
    },
    { AminoAcid::L /* CTA */,
      AminoAcid::L /* CTC */,
      AminoAcid::L /* CTG */,
      AminoAcid::L /* CTT */
    }
  },
  { { AminoAcid::E /* GAA */,
      AminoAcid::D /* GAC */,
      AminoAcid::E /* GAG */,
      AminoAcid::D /* GAT */
    },
    { AminoAcid::A /* GCA */,
      AminoAcid::A /* GCC */,
      AminoAcid::A /* GCG */,
      AminoAcid::A /* GCT */
    },
    { AminoAcid::G /* GGA */,
      AminoAcid::G /* GGC */,
      AminoAcid::G /* GGG */,
      AminoAcid::G /* GGT */
    },
    { AminoAcid::V /* GTA */,
      AminoAcid::V /* GTC */,
      AminoAcid::V /* GTG */,
      AminoAcid::V /* GTT */
    }
  },
  { { AminoAcid::STP /* TAA */,
      AminoAcid::Y /* TAC */,
      AminoAcid::STP /* TAG */,
      AminoAcid::Y /* TAT */
    },
    { AminoAcid::S /* TCA */,
      AminoAcid::S /* TCC */,
      AminoAcid::S /* TCG */,
      AminoAcid::S /* TCT */
    },
    { AminoAcid::STP /* TGA */,
      AminoAcid::C /* TGC */,
      AminoAcid::W /* TGG */,
      AminoAcid::C /* TGT */
    },
    { AminoAcid::L /* TTA */,
      AminoAcid::F /* TTC */,
      AminoAcid::L /* TTG */,
      AminoAcid::F /* TTT */
    }
  } };

  if (*triplet == Nucleotide::GAP
      && (*(triplet + 1) == Nucleotide::GAP)
      && (*(triplet + 2) == Nucleotide::GAP))
    return AminoAcid::GAP;

  if (triplet->isAmbiguity()
      || (triplet + 1)->isAmbiguity()
      || (triplet + 2)->isAmbiguity())
    return AminoAcid::X;

  return
    codonTable[triplet->intRep()]
              [(triplet + 1)->intRep()]
              [(triplet + 2)->intRep()];
}

std::set<AminoAcid>
Codon::translateAll(const NTSequence::const_iterator triplet)
{
  std::set<AminoAcid> result;

  NTSequence s(triplet, triplet + 3);

  std::vector<NTSequence> possibilities;
  s.nonAmbiguousSequences(possibilities);

  for (unsigned i = 0; i < possibilities.size(); ++i)
    result.insert(translate(possibilities[i].begin()));

  return result;
}

};
