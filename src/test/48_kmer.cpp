#include <map>
#include <vector>
#include <string>
#include "nucl_kmer_generator.h"
#include <iostream>
#include <zlib.h>
#include "kseq.h"

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif


int main(int argc, char **argv) {

    gzFile fp = gzopen(argv[1], "r");
    kseq_t *seq = kseq_init(fp);

    std::map<std::string, std::vector<string>> kmer_table;
    typedef std::map<std::string, std::vector<string>>::iterator it_type;
    it_type it;

    while (kseq_read(seq) >= 0) {
        // printf("%s\n", seq->seq.s);
        std::string string_seq(seq->seq.s);
        std::string string_name(seq->name.s);
        NuclKmerGenerator kmers = NuclKmerGenerator(string_seq, 48, false);

        while (kmers.hasNext()) {
            NuclKmer temp = kmers.next();
            std::string kmer = temp.decodePacked();
            it = kmer_table.find(kmer);

            if (it != kmer_table.end()) {
                it->second.push_back(string_name);
            }
            else {
                std::vector<string> vec;
                vec.push_back(string_name);
                kmer_table.insert(std::pair<std::string, std::vector<string>>(kmer, vec));
            }
        }
    }

    for (it_type it2 = kmer_table.begin(); it2 != kmer_table.end(); it2++) {
        if (it2->second.size() > 100) {
            std::cout << it2->first << " count = " << it2->second.size()  << "\n";
            typedef std::vector<string>::iterator vector_it;

            for (vector_it it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
                std::cout << *it3 << "\n";
            }

            std::cout << "\n";
        }

    }
}
