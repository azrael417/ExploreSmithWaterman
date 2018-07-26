#ifndef UTIL_REFERENCES_H_
#define UTIL_REFERENCES_H_

#include <iostream>
#include <string>
#include <vector>

#include "types.h"

using std::string;

class References {
 public:
  References();
  void LoadReferences(const char* filename);
  //bool GetSequence(const int& chr_id, const int& pos, const int& length, 
  //                 string* seq) const;
  inline int GetReferenceCount() const;
  inline View1D<char> GetReferenceSequence(const int& id) const;
  inline View1D<char> GetReferenceSequence(const int& id, int* length) const;
  inline int GetMaxSequenceLength() const;

 private:
  int max_sequence_length_, max_name_length_;
  View2D<char, Kokkos::LayoutRight> names_, sequences_;
  View1D<int> sequences_end_;
};

inline int References::GetMaxSequenceLength() const {
  return max_sequence_length_;
}

int References::GetReferenceCount() const {
  return static_cast<int>(names_.extent(0));
}


inline View1D<char> References::GetReferenceSequence(const int& id) const {
  if (id >= static_cast<int>(sequences_.extent(0))) {
    std::cerr << "ERROR: The requested id(" << id << ") is invalid."<< std::endl;
    return View1D<char>();
  } else {
    return Kokkos::subview(sequences_, id, Kokkos::ALL);
  }
}


inline View1D<char> References::GetReferenceSequence(const int& id, int* length) const {
  View1D<char> result = GetReferenceSequence(id);
  if (result.size() != 0) *length = sequences_end_(id);
  else length = 0;

  return result;
}
#endif //UTIL_REFERENCES_H_
