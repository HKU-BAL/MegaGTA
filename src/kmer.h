#ifndef KMER_1_H__
#define KMER_1_H__

#include <stdint.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include "city.h"
#include <string.h>

using namespace std;

class Kmer {
  public:
    static const int MAX_NUCL_KMER_SIZE = 64;
    static const int MAX_PROT_KMER_SIZE = 24;

  public: //protected
    uint64_t kmers[2] = {};
    int k;
    int l;
    int last_fill;
    uint64_t last_mask;

    //The following three default are for nucleotides
    int bits_to_shift = 0x2; //2 bits
    int char_mask = 0x3; // 11
    int items_per_bucket = 32;

  protected:
    Kmer() {}
    Kmer(int bits_to_shift, int char_mask, int items_per_bucket, int k) : k(k), bits_to_shift(bits_to_shift), char_mask(char_mask), items_per_bucket(items_per_bucket) {
        l = (k + items_per_bucket - 1) / items_per_bucket;
        last_fill = items_per_bucket - (items_per_bucket * l - k);
        uint64_t tmp = 0;

        for (int index = 0; index < last_fill; index++) {
            tmp = (tmp << bits_to_shift) | char_mask;
        }

        last_mask = tmp;
    }


  public:
    virtual uint8_t charToByte(char c) = 0;
    virtual char intToChar(int i) = 0;
    virtual void shiftRight(char c) = 0;
    virtual void shiftLeft(char c) = 0;

    int length() {
        return k;
    }
    int packedLength() {
        return l;
    }
    uint64_t getPart(int index) {
        return kmers[index];
    }
    uint64_t *getPackedKmers() {
        return kmers;
    }

    void initialize(const string &char_kmer) {
        uint8_t b;

        for (int index = 0; index < (int)char_kmer.size(); ++index) {
            b = charToByte(char_kmer[index]);

            if (b == 31) {
                throw std::invalid_argument("Kmer contains one or more invalid bases");
            }

            if (index < items_per_bucket) {
                kmers[0] = (kmers[0] << bits_to_shift) | (b & char_mask);
            }
            else {
                kmers[1] = (kmers[1] << bits_to_shift) | (b & char_mask);
            }
        }
    }

  protected:
    void shiftRight(uint8_t b) {
        if (l == 2) {
            uint8_t overflow = (uint8_t) (kmers[0] & char_mask);
            kmers[0] = (kmers[0] >> bits_to_shift) | ((uint64_t) b << ( (items_per_bucket - 1) * bits_to_shift));

            kmers[1] = (kmers[1] >> bits_to_shift) | ((uint64_t) overflow << ((last_fill - 1) * bits_to_shift));
            kmers[1] &= last_mask;
        }
        else {
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
            kmers[0] = kmers[0] << (64 - items_per_bucket * bits_to_shift) >> (64 - items_per_bucket * bits_to_shift);
        }
        else {
            kmers[0] = (kmers[0] << bits_to_shift) | ((uint64_t) b & char_mask);
            kmers[0] &= last_mask;
        }
    }

  public:
    string decodePacked() {
        string buf;

        for (int seg = l; seg > 0; seg--) {
            int index = items_per_bucket;

            if (seg == l) {
                index = last_fill;
            }

            uint64_t this_k = (seg == 2) ? kmers[1] : kmers[0];

            for (; index > 0; index--) {
                buf.append(1, intToChar((int) (this_k & char_mask)) );
                this_k = this_k >> bits_to_shift;
            }
        }

        reverse(buf.begin(), buf.end());
        return buf;
    }

    bool equals(const Kmer &kmer) const {
        if (kmers[0] != kmer.kmers[0] || kmers[1] != kmer.kmers[1]) {
            return false;
        }

        return true;
    }

    bool operator ==(const Kmer &kmer) const {
        if (kmers[0] != kmer.kmers[0] || kmers[1] != kmer.kmers[1]) {
            return false;
        }

        return true;
    }

    uint64_t hash() const {
        return CityHash64((const char *)kmers, sizeof(kmers[0]) * 2);
    }
};

#endif

