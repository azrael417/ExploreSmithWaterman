#include "parameter_parser.h"

#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

bool CheckParameters(const Parameters& param);
void PrintHelp(const string& program);

template<typename T>
bool convert_from_string(const std::string& s, T& r) {
  std::istringstream iss(s);
  iss >> r;
  return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}

void ParseArgumentsOrDie(const int argc, 
                         char* const * argv, 
                         Parameters* param) {

  if (argc == 1) { // no argument
    PrintHelp(argv[0]);
    exit(1);
  }

  // record command line
  param->command_line = argv[0];
  for ( int i = 1; i < argc; ++i ) {
    param->command_line += " ";
    param->command_line += argv[i];
  }

  const char *short_option = "hf:q:m:x:o:e:";

  const struct option long_option[] = {
    {"help", no_argument, NULL, 'h'},
    {"fasta", required_argument, NULL, 'f'},
    {"fastq", required_argument, NULL, 'q'},

    {"match", required_argument, NULL, 'm'},
    {"mismatch", required_argument, NULL, 'x'},
    {"open-gap",required_argument, NULL, 'o'},
    {"extend-gap",required_argument, NULL, 'e'},

    { 0, 0, 0, 0 }
  };

  int c = 0;
  bool help  = false;
  bool error = false;
  while (true) {
    int optionIndex = 0;
    c = getopt_long( argc, argv, short_option, long_option, &optionIndex );

    if ( c == -1 ) // end of options
      break;
		
    switch ( c ) {
      // help
      case 'h':
        help = true; break;
      // i/o parameters
      case 'f':
        param->fasta = optarg; break;
      case 'q':
        param->fastq = optarg; break;

      case 'm':
        if(!convert_from_string(optarg, param->match)) {
	  cerr << "WARNING: Cannot parse -m --match." << endl
               << "         Set it to default 10." << endl;
	}
	break;
	
      case 'x':
        if(!convert_from_string(optarg, param->mismatch)) {
	  cerr << "WARNING: Cannot parse -M --mismatch." << endl
               << "         Set it to default 9." << endl;
	}
	break;

      case 'o':
        if(!convert_from_string(optarg, param->open_gap)) {
	  cerr << "WARNING: Cannot parse -o --open-gap." << endl
               << "         Set it to default 15." << endl;
	}
	break;

      case 'e':
        if(!convert_from_string(optarg, param->extend_gap)) {
	  cerr << "WARNING: Cannot parse -e --extend-gap." << endl
               << "         Set it to default 6.66." << endl;
	}
	break;

      default: break;
    }

  }

  if (error) exit(1);

  if (help) {
    PrintHelp(argv[0]);
    exit(1);
  }

  if (!CheckParameters(*param))
    exit(1);
}

// true: passing the checker
bool CheckParameters(const Parameters& param) {
	
	bool errorFound = false;
	// necessary parameters
	if (param.fasta.empty()) {
		cerr << "ERROR: Please specific a FASTA file, -f." << endl;
		errorFound = true;
	}

	if (param.fastq.empty()) {
		cerr << "ERROR: Please specific a FASTQ file, -q." << endl;
		errorFound = true;
	}

	
	return !errorFound;

}

void PrintHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " -f <FASTA> -q <FASTQ>"
		<< endl
		<< endl
		<< "Help:" << endl
		<< endl
		<< "   -h --help             Print this help dialog." << endl
		<< endl
		
		<< "Input & Output:" << endl
		<< endl
		<< "   -f --fasta <FASTA>    Input reference FASTA file." << endl
		<< "   -q --fastq <FASTQ>    Input read FASTQ file." << endl
		<< endl

		<< "Operation parameters:" << endl
		<< endl
		<< "   -m --match <float>    Match score [10]." << endl
		<< "   -x --mismatch <float> Mismatch score [9]." << endl
		<< "   -o --open-gap <float> Gap open penalty [15]." << endl
		<< "   -e --extend-gap <float> Gap extention penalty [6.66]." << endl
		<< endl;
}

