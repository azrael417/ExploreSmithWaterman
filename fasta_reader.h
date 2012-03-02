#ifndef UTIL_FASTA_READER_H_
#define UTIL_FASTA_READER_H_

#include <fstream>
#include <string>

using std::ifstream;
using std::string;

class FastaReader{
 public:
  FastaReader();
  ~FastaReader();
  bool Open(const char* filename);
  bool Close();
  bool LoadNextRead(string* readname, string* sequence);
 private:
  ifstream file_;
  string readname_;
  int line_;
  bool error_;
};

#endif // UTIL_FASTA_READER_H_
