class Codon
{
public:
	Codon() {};
	~Codon() {};
	constexpr static char codonTable[4][4][4] = {
  { { 'K' /* AAA */,
      'N' /* AAC */,
      'K' /* AAG */,
      'N' /* AAT */
    },
    { 'T' /* ACA */,
      'T' /* ACC */,
      'T' /* ACG */,
      'T' /* ACT */
    },
    { 'R' /* AGA */,
      'S' /* AGC */,
      'R' /* AGG */,
      'S' /* AGT */
    },
    { 'I' /* ATA */,
      'I' /* ATC */,
      'M' /* ATG */,
      'I' /* ATT */
    }
  },
  { { 'Q' /* CAA */,
      'H' /* CAC */,
      'Q' /* CAG */,
      'H' /* CAT */
    },
    { 'P' /* CCA */,
      'P' /* CCC */,
      'P' /* CCG */,
      'P' /* CCT */
    },
    { 'R' /* CGA */,
      'R' /* CGC */,
      'R' /* CGG */,
      'R' /* CGT */
    },
    { 'L' /* CTA */,
      'L' /* CTC */,
      'L' /* CTG */,
      'L' /* CTT */
    }
  },
  { { 'E' /* GAA */,
      'D' /* GAC */,
      'E' /* GAG */,
      'D' /* GAT */
    },
    { 'A' /* GCA */,
      'A' /* GCC */,
      'A' /* GCG */,
      'A' /* GCT */
    },
    { 'G' /* GGA */,
      'G' /* GGC */,
      'G' /* GGG */,
      'G' /* GGT */
    },
    { 'V' /* GTA */,
      'V' /* GTC */,
      'V' /* GTG */,
      'V' /* GTT */
    }
  },
  { { '*' /* TAA */,
      'Y' /* TAC */,
      '*' /* TAG */,
      'Y' /* TAT */
    },
    { 'S' /* TCA */,
      'S' /* TCC */,
      'S' /* TCG */,
      'S' /* TCT */
    },
    { '*' /* TGA */,
      'C' /* TGC */,
      'W' /* TGG */,
      'C' /* TGT */
    },
    { 'L' /* TTA */,
      'F' /* TTC */,
      'L' /* TTG */,
      'F' /* TTT */
    }
  } };

};