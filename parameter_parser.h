
#ifndef UTIL_ParameterParser_H_
#define UTIL_ParameterParser_H_

#include <string>

using std::string;

struct Parameters {
  // i/o parameters
  string fasta;    // -f  --fasta
  string fastq;    // -q  --fastq

  // operation parameters
  float match;
  float mismatch;
  float open_gap;
  float extend_gap;
	
  // command line
  string command_line;

  // default values
  Parameters()
      : fasta()
      , fastq()
      , match(10.0)
      , mismatch(-9.0)
      , open_gap(15.0)
      , extend_gap(6.66)
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);

#endif // UTIL_ParameterParser_H_
