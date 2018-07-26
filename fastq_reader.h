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
  KOKKOS_INLINE_FUNCTION FastqReader (const FastqReader& ) = default;
  KOKKOS_INLINE_FUNCTION FastqReader & operator= (const FastqReader& ) = default;
  bool Open(const char* filename);
  bool Close();
  bool LoadNextRead(string* readname, string* sequence, string* qual);
  bool LoadNextBatch(const int& BatchSize);
  KOKKOS_INLINE_FUNCTION View1D<char> GetSequence(const int& id) const;
  KOKKOS_INLINE_FUNCTION View2D<char> GetSequences() const;
  View1D<char> GetRead(const int& id) const;
  KOKKOS_INLINE_FUNCTION int GetSequenceLength(const int& id) const;
  KOKKOS_INLINE_FUNCTION View1D<int> GetSequenceLengths() const;
 private:
  ifstream file_;
  string readname_;
  int line_;
  bool error_;
  int readsize;
  View2D<char> readnames, sequences, quals;
  View1D<int> sequences_end;
};

KOKKOS_INLINE_FUNCTION int FastqReader::GetSequenceLength(const int& id) const{ 
  if(id >= readsize) return 0;
  else return sequences_end(id); 
}

KOKKOS_INLINE_FUNCTION View1D<int> FastqReader::GetSequenceLengths() const{ 
  return sequences_end; 
}

KOKKOS_INLINE_FUNCTION View1D<char> FastqReader::GetSequence(const int& id) const {
  if(id >= readsize) return View1D<char>();
  else return Kokkos::subview(sequences, id, Kokkos::ALL);
}

KOKKOS_INLINE_FUNCTION View2D<char> FastqReader::GetSequences() const {
  return sequences;
}


#endif // UTIL_FASTQ_READER_H_
