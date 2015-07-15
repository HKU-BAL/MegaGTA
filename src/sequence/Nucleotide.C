#include <ctype.h>
#include <stdlib.h>

#include "ParseException.h"
#include "Nucleotide.h"

namespace {

int sampleUniform(int one, int two)
{
  double d = drand48();
  
  return (d < 0.5 ? one : two);
}

int sampleUniform(int one, int two, int three)
{
  double d = drand48() * 3.;

  return (d < 1. ? one : (d < 2. ? two : three));
}

int sampleUniform(int one, int two, int three, int four)
{
  double d = drand48() * 4.;

  return (d < 1. ? one : (d < 2. ? two : (d < 3. ? three : four)));
}
};

namespace seq {

const char Nucleotide::NT_CHAR[] = {'A', 'C', 'G', 'T',
				    'M', 'R', 'W', 'S',
				    'Y', 'K', 'V', 'H',
				    'D', 'B', 'N', '-' };

const Nucleotide Nucleotide::A(Nucleotide::NT_A);
const Nucleotide Nucleotide::C(Nucleotide::NT_C);
const Nucleotide Nucleotide::G(Nucleotide::NT_G);
const Nucleotide Nucleotide::T(Nucleotide::NT_T);
const Nucleotide Nucleotide::M(Nucleotide::NT_M);
const Nucleotide Nucleotide::R(Nucleotide::NT_R);
const Nucleotide Nucleotide::W(Nucleotide::NT_W);
const Nucleotide Nucleotide::S(Nucleotide::NT_S);
const Nucleotide Nucleotide::Y(Nucleotide::NT_Y);
const Nucleotide Nucleotide::K(Nucleotide::NT_K);
const Nucleotide Nucleotide::V(Nucleotide::NT_V);
const Nucleotide Nucleotide::H(Nucleotide::NT_H);
const Nucleotide Nucleotide::D(Nucleotide::NT_D);
const Nucleotide Nucleotide::B(Nucleotide::NT_B);
const Nucleotide Nucleotide::N(Nucleotide::NT_N);
const Nucleotide Nucleotide::GAP(Nucleotide::NT_GAP);

Nucleotide::Nucleotide()
  : rep_(NT_N)
{ }

Nucleotide::Nucleotide(char c)
  throw (ParseException)
{
  switch (toupper(c)) {
  case 'A': rep_ = NT_A; break;
  case 'C': rep_ = NT_C; break;
  case 'G': rep_ = NT_G; break;
  case 'T': case 'U': rep_ = NT_T; break;
  case 'M': rep_ = NT_M; break;
  case 'R': rep_ = NT_R; break;
  case 'W': rep_ = NT_W; break;
  case 'S': rep_ = NT_S; break;
  case 'Y': rep_ = NT_Y; break;
  case 'K': rep_ = NT_K; break;
  case 'V': rep_ = NT_V; break;
  case 'H': rep_ = NT_H; break;
  case 'D': rep_ = NT_D; break;
  case 'B': rep_ = NT_B; break;
  case 'N': rep_ = NT_N; break;
  case '-': rep_ = NT_GAP; break;
  default:
    throw ParseException
      (std::string("Invalid nucleotide character: '") + c + "'");
  }
}

void Nucleotide::sampleAmbiguity()
{
  switch (rep_) {
  case NT_A:
  case NT_C:
  case NT_G:
  case NT_T:
  case NT_GAP:
    break;
  case NT_M:
    rep_ = sampleUniform(NT_A, NT_C); break;
  case NT_R:
    rep_ = sampleUniform(NT_A, NT_G); break;  
  case NT_W:
    rep_ = sampleUniform(NT_A, NT_T); break;
  case NT_S:
    rep_ = sampleUniform(NT_C, NT_G); break;
  case NT_Y:
    rep_ = sampleUniform(NT_C, NT_T); break;
  case NT_K:
    rep_ = sampleUniform(NT_G, NT_T); break;
  case NT_V:
    rep_ = sampleUniform(NT_A, NT_C, NT_G); break;
  case NT_H:
    rep_ = sampleUniform(NT_A, NT_C, NT_T); break;
  case NT_D:
    rep_ = sampleUniform(NT_A, NT_G, NT_T); break;
  case NT_B:
    rep_ = sampleUniform(NT_C, NT_G, NT_T); break;
  case NT_N:
    rep_ = sampleUniform(NT_A, NT_C, NT_G, NT_T); break;
  default:
    std::cerr << rep_ << std::endl;
    assert(false);
  }
}

/**
 * Get all non ambiguous nucleotides represented by this nucleotide.
 */ 
void Nucleotide::nonAmbiguousNucleotides(std::vector<Nucleotide>& result) const
{
  switch (rep_) {
  case NT_A:
  case NT_C:
  case NT_G:
  case NT_T:
  case NT_GAP:
    result.push_back(*this);
    break;
  case NT_M:
    result.push_back(A);
    result.push_back(C);
    break;
  case NT_R:
    result.push_back(A);
    result.push_back(G);
    break;
  case NT_W:
    result.push_back(A);
    result.push_back(T);
    break;
  case NT_S:
    result.push_back(C);
    result.push_back(G);
    break;
  case NT_Y:
    result.push_back(C);
    result.push_back(T);
    break;
  case NT_K:
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_V:
    result.push_back(A);
    result.push_back(C);
    result.push_back(G);
    break;
  case NT_H:
    result.push_back(A);
    result.push_back(C);
    result.push_back(T);
    break;
  case NT_D:
    result.push_back(A);
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_B:
    result.push_back(C);
    result.push_back(G);
    result.push_back(T);
    break;
  case NT_N:
    result.push_back(A);
    result.push_back(C);
    result.push_back(G);
    result.push_back(T);
    break;
  default:
    std::cerr << rep_ << std::endl;
    assert(false);
  }
}


std::ostream& operator<< (std::ostream& s, const Nucleotide nt)
{
  return s << nt.toChar();
}

};
