#include <zlib.h>
#include "sequence/NTSequence.h"
#include "sequence/AASequence.h"
#include "kseq.h"
#include <algorithm>

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

using namespace std;

int translate(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <nucl_seq> \n", argv[0]);
        exit(1);
    }

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        cout << ">" << seq->name.s << '\n';
        seq::NTSequence nts = seq::NTSequence("", "", seq->seq.s);
        string trans = seq::AASequence::translate(nts.begin(), nts.begin() + (nts.size() / 3) * 3).asString();
        transform(trans.begin(), trans.end(), trans.begin(), ::tolower);
        cout << trans << '\n';
    }

    return 0;
}
