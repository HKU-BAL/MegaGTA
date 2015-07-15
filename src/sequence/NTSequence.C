#include <fstream>

#include "NTSequence.h"
#include "ParseException.h"

namespace seq {

NTSequence::NTSequence()
  : std::vector<Nucleotide>()
{ }

NTSequence::NTSequence(unsigned size)
  : std::vector<Nucleotide>(size)
{ }

NTSequence::NTSequence(const std::string name, const std::string description,
		       const std::string aSeqString,
		       bool sampleAmbiguities)
  throw (ParseException)
  : std::vector<Nucleotide>(aSeqString.length()),
    name_(name),
    description_(description)
{
  for (unsigned i = 0; i < aSeqString.length(); ++i) {
    Nucleotide nt(aSeqString[i]);
    if (sampleAmbiguities)
      nt.sampleAmbiguity();

    (*this)[i] = nt;
  }
}

NTSequence::NTSequence(const const_iterator first,
		       const const_iterator last)
  : std::vector<Nucleotide>(first, last)
{ }

void NTSequence::sampleAmbiguities()
{
  for (unsigned i = 0; i < size(); ++i) {
    (*this)[i].sampleAmbiguity();
  }
}

void NTSequence::nonAmbiguousSequences(std::vector<NTSequence>& result) const
{
  iterateNonAmbiguous(NTSequence(), result);
}

void NTSequence::iterateNonAmbiguous(const NTSequence& head,
				     std::vector<NTSequence>& result) const
{
  /*
   * find the next ambigous codon (if any)
   */
  NTSequence s = head;
  unsigned i = head.size();
  for (; i < size(); ++i)
    if ((*this)[i].isAmbiguity())
      break;
    else
      s.push_back((*this)[i]);

  if (i == size())
    result.push_back(s);
  else {
    std::vector<Nucleotide> ambiguities;
    (*this)[i].nonAmbiguousNucleotides(ambiguities);
    for (unsigned i = 0; i < ambiguities.size(); ++i) {
      s.push_back(ambiguities[i]);
      iterateNonAmbiguous(s, result);
      s.pop_back();
    }
  }
}

std::string NTSequence::asString() const
{
  std::string result(size(), '-');

  for (unsigned i = 0; i < size(); ++i) {
    result[i] = (*this)[i].toChar();
  }

  return result;
}

/// \cond

void readFastaEntry(std::istream& i,
		    std::string& name,
		    std::string& description,
		    std::string& sequence)
  throw (ParseException)
{
    char ch;
    char c[512];

    i.getline(c, 511);
    if (i) {
      if (c[0] != '>') {
	throw ParseException(std::string("FASTA file expected '>', got: '")
			     + c[0] + "'");
      }

      std::string nameDesc = c + 1;
      std::string::size_type spacepos = nameDesc.find(" ");
      name = nameDesc.substr(0, spacepos);
      description = (spacepos == std::string::npos
		     ? ""
		     : nameDesc.substr(spacepos));

      for (ch = i.get(); (ch != EOF) && (ch != '>'); ch = i.get()) {
	if ((ch != '\n') && (ch != '\r') && (ch != ' ')) {
	  if (((ch >= 'a') && (ch <= 'z'))
	      || ((ch >= 'A') && (ch <= 'Z'))
	      || (ch == '-') || (ch == '*')) {
	    sequence += ch;
	  } else {
	    throw ParseException
	      (std::string("Illegal character in FASTA file: '")
	       + (char)ch + "'");
	  }
	}

	if (i.peek() == EOF)
	  break;
      }

      if (ch == '>')
	i.putback(ch);
    }
}

void writeFastaEntry(std::ostream& o,
		     const std::string& name,
		     const std::string& description,
		     const std::string& sequence)
{
  o << ">" << name << " " << description << std::endl;
  if (sequence.size() == 0)
    o << std::endl;
  else {
    for (unsigned i = 0; i <= (sequence.size() - 1) / 60; ++i) {
      int s = i * 60;
      o << sequence.substr(s, 60) << std::endl;
    }
  }
}

/// \endcond

std::istream& operator>>(std::istream& i, NTSequence& sequence)
  throw (ParseException)
{
  std::string name, description, seqString;

  readFastaEntry(i, name, description, seqString);
  sequence = NTSequence(name, description, seqString);

  return i;
}

std::ostream& operator<<(std::ostream& o, const NTSequence& sequence)
{
  writeFastaEntry(o, sequence.name(), sequence.description(),
		  sequence.asString());
  return o;
}

};
