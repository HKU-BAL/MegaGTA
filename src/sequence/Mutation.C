#include <boost/tokenizer.hpp>

#include "Mutation.h"

namespace seq {

std::set<AAMutation> readMutations(std::istream& mutationFile,
				   std::string prefix)
  throw (ParseException)
{
  std::set<AAMutation> result;

  std::string line;
  typedef boost::tokenizer<boost::escaped_list_separator<char> > csv_tok;
  getline(mutationFile, line);
  csv_tok tok(line);

  for (csv_tok::iterator i = tok.begin(); i != tok.end(); ++i) {
    std::string mutation = *i;

    if (mutation.length() < prefix.length() + 2)
      throw ParseException("Error while parsing mutation '"
			   + mutation + "': too short for mutation with "
			   "prefix '" + prefix + "'");

    if (mutation.substr(0, prefix.length()) != prefix)
      throw ParseException("Error while parsing mutation '"
			   + mutation + "': expected to start with '"
			   + prefix + "'");

    try {
      AminoAcid aa(mutation[mutation.length() - 1]);

      char *endptr;
      int pos
	= strtol(mutation.substr(prefix.length(),
				 mutation.length()
				 - prefix.length() - 1).c_str(),
		 &endptr, 10);

      if (*endptr != 0)
	throw ParseException("could not parse position");

      result.insert(AAMutation(pos, aa));
    } catch (ParseException& e) {
      throw ParseException("Error while parsing mutation '"
			   + mutation + "': " + e.message());
    }

  }

  return result;
}

};
