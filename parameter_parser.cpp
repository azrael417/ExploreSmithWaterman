#include "parameter_parser.h"

#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <iostream>
#include <string>


using std::cout;
using std::endl;
using std::string;

bool CheckParameters(const Parameters& param);
void PrintHelp(const string& program);

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

  const char *short_option = "hi:";

  const struct option long_option[] = {
    { "help", no_argument, NULL, 'h' },
    { "fasta", required_argument, NULL, 'f' },
    { "fastq", required_argument, NULL, 'q' },

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
		cout << "ERROR: Please specific a FASTA file, -f." << endl;
		errorFound = true;
	}

	if (param.fastq.empty()) {
		cout << "ERROR: Please specific a FASTQ file, -q." << endl;
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
		<< "   -h --help           Print this help dialog." << endl
		<< endl
		<< "Input & Output:" << endl
		<< endl
		<< "   -f --fasta <FASTA>   Input reference FASTA file." << endl
		<< "   -q --fastq <FASTQ>   Input read FASTQ file." << endl
		
		<< endl;
}

