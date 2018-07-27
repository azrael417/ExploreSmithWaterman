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
  //KOKKOS_INLINE_FUNCTION FastqReader (const FastqReader& ) = default;
  //KOKKOS_INLINE_FUNCTION FastqReader & operator= (const FastqReader& ) = default;
  bool Open(const char* filename);
  bool Close();
  bool LoadNextRead(string* readname, string* sequence, string* qual);
  bool LoadNextBatch(const int& BatchSize);
  KOKKOS_INLINE_FUNCTION View1D<char> GetSequence(const int& id) const;
  KOKKOS_INLINE_FUNCTION View2D<char> GetSequences() const;
  View1D<char> GetRead(const int& id) const;
  KOKKOS_INLINE_FUNCTION int GetSequenceLength(const int& id) const;
  KOKKOS_INLINE_FUNCTION View1D<int> GetSequenceLengths() const;
  KOKKOS_INLINE_FUNCTION int GetReadsize() const;
 private:
  ifstream file_;
  string readname_;
  int line_;
  bool error_;
  int readsize;
  View2D<char, Kokkos::LayoutRight, Kokkos::HostSpace> readnames, sequences, quals;
  View1D<int, Kokkos::HostSpace> sequences_end;
  View2D<char, Kokkos::LayoutRight> d_sequences;
  View1D<int> d_sequences_end;
};

KOKKOS_INLINE_FUNCTION int FastqReader::GetReadsize() const{ 
  return readsize;
}

KOKKOS_INLINE_FUNCTION int FastqReader::GetSequenceLength(const int& id) const{ 
  if(id >= readsize) return 0;
  else return d_sequences_end(id); 
}

KOKKOS_INLINE_FUNCTION View1D<int> FastqReader::GetSequenceLengths() const{ 
  return d_sequences_end; 
}

KOKKOS_INLINE_FUNCTION View1D<char> FastqReader::GetSequence(const int& id) const {
  if(id >= readsize) return View1D<char>();
  else return Kokkos::subview(d_sequences, id, Kokkos::ALL);
}

KOKKOS_INLINE_FUNCTION View2D<char> FastqReader::GetSequences() const {
  return d_sequences;
}


#endif // UTIL_FASTQ_READER_H_
