// This may look like C code, but it's really -*- C++ -*-
#ifndef PARSE_EXCEPTION_H_
#define PARSE_EXCEPTION_H_

#include <string>

namespace seq {

/**
 * Exception thrown when an error was encountered while parsing the
 * string representation of an nucleotide, nucleotide sequence, amino
 * acid, amino acid sequence, or a FASTA file.
 *
 * \sa Nucleotide::Nucleotide(char), AminoAcid::AminoAcid(char),
 * NTSequence::NTSequence(const std::string, const std::string, const
 * std::string, bool), AASequence::AASequence(const std::string, const
 * std::string, const std::string), operator>> (std::istream&,
 * NTSequence&), operator>> (std::istream&, AASequence&)
 */
class ParseException
{
public:
  ParseException(const std::string message)
    : message_(message) { }

  /**
   * The message describing the error.
   */
  std::string message() const { return message_; }

private:
  std::string message_;
};

};

#endif // PARSE_EXCEPTION_H_
