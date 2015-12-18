// This may look like C code, but it's really -*- C++ -*-
#ifndef MUTATION_H_
#define MUTATION_H_

#include <iostream>
#include <set>

#include "Nucleotide.h"
#include "AminoAcid.h"

namespace seq {

/**
 * A mutation in a sequence of type Char.
 *
 * \sa AAMutation, NTMutation
 */
template<typename Char>
class Mutation
{
public:
  Mutation(int pos, Char from, Char to)
    : pos_(pos),
      from_(from),
      to_(to)
  { }

  Mutation(int pos, Char to)
    : pos_(pos),
      to_(to)
  { }

  Mutation()
    : pos_(-1)
  { }

  int pos()   const { return pos_; }
  Char from() const { return from_; }
  Char to()   const { return to_; }

  Mutation<Char> reverse() const {
    return Mutation<Char>(pos_, to_, from_);
  }

  bool isValid() const { return pos_ >= 0; }

  bool operator== (const Mutation<Char>& other) const {
    return (pos_ == other.pos_) && (to_ == other.to_);
  }

  bool operator< (const Mutation<Char>& other) const {
    return (pos_ < other.pos_)
      || ((pos_ == other.pos_)
	  && (to_.intRep() < other.to_.intRep()));
  }

 private:
  int  pos_;
  Char from_, to_;
};

/**
 * A typedef for nucleotide mutations.
 */
typedef class Mutation<Nucleotide> NTMutation;
/**
 * A typedef for amino acid mutations.
 */
typedef class Mutation<AminoAcid> AAMutation;

extern std::set<AAMutation> readMutations(std::istream& mutationFile,
					  std::string prefix)
  throw (ParseException);

};

#endif // MUTATION_H_
