class Codon {
  public:
    Codon() {};
    ~Codon() {};
    constexpr static char codonTable[4][4][4] = {
        {   {
                'K' /* AAA */,
                'N' /* AAC */,
                'K' /* AAG */,
                'N' /* AAT */
            },
            {
                'T' /* ACA */,
                'T' /* ACC */,
                'T' /* ACG */,
                'T' /* ACT */
            },
            {
                'R' /* AGA */,
                'S' /* AGC */,
                'R' /* AGG */,
                'S' /* AGT */
            },
            {
                'I' /* ATA */,
                'I' /* ATC */,
                'M' /* ATG */,
                'I' /* ATT */
            }
        },
        {   {
                'Q' /* CAA */,
                'H' /* CAC */,
                'Q' /* CAG */,
                'H' /* CAT */
            },
            {
                'P' /* CCA */,
                'P' /* CCC */,
                'P' /* CCG */,
                'P' /* CCT */
            },
            {
                'R' /* CGA */,
                'R' /* CGC */,
                'R' /* CGG */,
                'R' /* CGT */
            },
            {
                'L' /* CTA */,
                'L' /* CTC */,
                'L' /* CTG */,
                'L' /* CTT */
            }
        },
        {   {
                'E' /* GAA */,
                'D' /* GAC */,
                'E' /* GAG */,
                'D' /* GAT */
            },
            {
                'A' /* GCA */,
                'A' /* GCC */,
                'A' /* GCG */,
                'A' /* GCT */
            },
            {
                'G' /* GGA */,
                'G' /* GGC */,
                'G' /* GGG */,
                'G' /* GGT */
            },
            {
                'V' /* GTA */,
                'V' /* GTC */,
                'V' /* GTG */,
                'V' /* GTT */
            }
        },
        {   {
                '*' /* TAA */,
                'Y' /* TAC */,
                '*' /* TAG */,
                'Y' /* TAT */
            },
            {
                'S' /* TCA */,
                'S' /* TCC */,
                'S' /* TCG */,
                'S' /* TCT */
            },
            {
                '*' /* TGA */,
                'C' /* TGC */,
                'W' /* TGG */,
                'C' /* TGT */
            },
            {
                'L' /* TTA */,
                'F' /* TTC */,
                'L' /* TTG */,
                'F' /* TTT */
            }
        }
    };

    constexpr static char rc_codonTable[4][4][4] = {
        {   {
                'F' /* AAA = TTT*/,
                'V' /* AAC = GTT*/,
                'L' /* AAG = CTT*/,
                'I' /* AAT = ATT*/
            },
            {
                'C' /* ACA = TGT*/,
                'G' /* ACC = GGT*/,
                'R' /* ACG = CGT*/,
                'S' /* ACT = AGT*/
            },
            {
                'S' /* AGA = TCT*/,
                'A' /* AGC = GCT*/,
                'P' /* AGG = CCT*/,
                'T' /* AGT = ACT*/
            },
            {
                'Y' /* ATA = TAT*/,
                'D' /* ATC = GAT*/,
                'H' /* ATG = CAT*/,
                'N' /* ATT = AAT*/
            }
        },
        {   {
                'L' /* CAA = TTG*/,
                'V' /* CAC = GTG*/,
                'L' /* CAG = CTG*/,
                'M' /* CAT = ATG*/
            },
            {
                'W' /* CCA = TGG*/,
                'G' /* CCC = GGG*/,
                'R' /* CCG = CGG*/,
                'R' /* CCT = AGG*/
            },
            {
                'S' /* CGA = TCG*/,
                'A' /* CGC = GCG*/,
                'P' /* CGG = CCG*/,
                'T' /* CGT = ACG*/
            },
            {
                '*' /* CTA = TAG*/,
                'E' /* CTC = GAG*/,
                'Q' /* CTG = CAG*/,
                'K' /* CTT = AAG*/
            }
        },
        {   {
                'F' /* GAA = TTC*/,
                'V' /* GAC = GTC*/,
                'L' /* GAG = CTC*/,
                'I' /* GAT = ATC*/
            },
            {
                'C' /* GCA = TGC*/,
                'G' /* GCC = GGC*/,
                'R' /* GCG = CGC*/,
                'S' /* GCT = AGC*/
            },
            {
                'S' /* GGA = TCC*/,
                'A' /* GGC = GCC*/,
                'P' /* GGG = CCC*/,
                'T' /* GGT = ACC*/
            },
            {
                'Y' /* GTA = TAC*/,
                'D' /* GTC = GAC*/,
                'H' /* GTG = CAC*/,
                'N' /* GTT = AAC*/
            }
        },
        {   {
                'L' /* TAA = TTA*/,
                'V' /* TAC = GTA*/,
                'L' /* TAG = CTA*/,
                'I' /* TAT = ATA*/
            },
            {
                '*' /* TCA = TGA*/,
                'G' /* TCC = GGA*/,
                'R' /* TCG = CGA*/,
                'R' /* TCT = AGA*/
            },
            {
                'S' /* TGA = TCA*/,
                'A' /* TGC = GCA*/,
                'P' /* TGG = CCA*/,
                'T' /* TGT = ACA*/
            },
            {
                '*' /* TTA = TAA*/,
                'E' /* TTC = GAA*/,
                'Q' /* TTG = CAA*/,
                'K' /* TTT = AAA*/
            }
        }
    };

};