#ifndef UTIL_FASTQ_READER_H_
#define UTIL_FASTQ_READER_H_

#include <fstream>
#include <string>

using std::ifstream;
using std::string;

class FastqReader{
 public:
  FastqReader();
  ~FastqReader();
  bool Open(const char* filename);
  bool Close();
  bool LoadNextRead(string* readname, string* sequence, string* qual);
  bool LoadNextBatch(string** readnames, string** sequence, string** qual, int* readsize, const int& BatchSize);
 private:
  ifstream file_;
  string readname_;
  int line_;
  bool error_;
};

#endif // UTIL_FASTQ_READER_H_
