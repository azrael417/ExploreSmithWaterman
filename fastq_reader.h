#ifndef UTIL_FASTQ_READER_H_
#define UTIL_FASTQ_READER_H_

#include <fstream>
#include <string>
#include "types.h"

using std::ifstream;
using std::string;

class FastqReader{
 public:
  FastqReader();
  ~FastqReader();
  bool Open(const char* filename);
  bool Close();
  bool LoadNextRead(string* readname, string* sequence, string* qual);
  bool LoadNextBatch(const int& BatchSize);
  View1D<char> GetSequence(const int& id);
  View1D<char> GetRead(const int& id);
  int GetSequenceLength(const int& id) const;
 private:
  ifstream file_;
  string readname_;
  int line_;
  bool error_;
  int readsize;
  View2D<char> readnames, sequences, quals;
  View1D<int> sequences_end;
};

#endif // UTIL_FASTQ_READER_H_
