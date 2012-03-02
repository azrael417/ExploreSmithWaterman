
#ifndef UTIL_ParameterParser_H_
#define UTIL_ParameterParser_H_

#include <string>

using std::string;

struct Parameters {
  // i/o parameters
  string fasta;    // -f  --fasta
  string fastq;    // -q  --fastq

  // operation parameters
	
  // command line
  string command_line;

  // default values
  Parameters()
      : fasta()
      , fastq()
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);

#endif // UTIL_ParameterParser_H_
