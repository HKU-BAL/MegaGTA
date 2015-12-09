// This may look like C code, but it's really -*- C++ -*-
#ifndef NUCLEOTIDESUBSTITUTIONMODEL_H_
#define NUCLEOTIDESUBSTITUTIONMODEL_H_

#include "Nucleotide.h"

namespace seq {

/**
 * Describes the mutation behaviour of (the whole of) replication
 * enzymes for one replication cycle. In particular, it models the
 * nucleotide-dependent error rate behaviour of the enzyme.
 */
class NucleotideSubstitutionModel
{
public:
  /**
   * Construct a stationary substitution model
   * based on:
   * <ul>
   *   <li> the 4 stationary nucleotide frequencies
   *   <li> the (relative) 6 symmetrical substitution rates.
   *   <li> the average error rate per generation per site (mu)
   * </ul>
   */
  NucleotideSubstitutionModel(double piA, double piC, double piG, double piT,
			      double rAC, double rAG, double rAT,
			      double rCG, double rCT, double rGT,
			      double errorRate);

  /**
   * Construct a (non-stationary) substitution model based on:
   * <ul>
   *   <li> the 12 relative substitution rates
   *   <LI> the 4 stationary nucleotide frequencies used to calibrate
   *        the average error rate per generation per site
   *   <li> the average error rate per generation per site
   * </ul>
   */
  NucleotideSubstitutionModel(double rAC, double rAG, double rAT,
			      double rCA, double rCG, double rCT,
			      double rGA, double rGC, double rGT,
			      double rTA, double rTC, double rTG,
			      double piA, double piC, double piG, double piT,
			      double errorRate);

  /**
   * Retrieve average rate of copying a nucleotide into a particular
   * other nucleotide.
   *
   * \return the error rate in nucleotide changes per replication cycle
   */
  double getMu(Nucleotide fromNT, Nucleotide toNT) const;

  void print(std::ostream& s) const;

private:
  double matrix_[4][4];
};

};

#endif // NUCLEOTIDESUBSTITUTIONMODEL
