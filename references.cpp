#include "references.h"

#include <iostream>
#include "fasta_reader.h"

using std::cerr;
using std::endl;

References::References()
    : names_()
    , sequences_() 
{}

void References::LoadReferences(const char* filename) {
  FastaReader reader;
  reader.Open(filename);
  names_.reserve(100);
  sequences_.reserve(100);

  int chr_id = 0;
  names_.resize(chr_id + 1);
  sequences_.resize(chr_id + 1);
  while (reader.LoadNextRead(&names_[chr_id], &sequences_[chr_id])){
    ++chr_id;
    names_.resize(chr_id + 1);
    sequences_.resize(chr_id + 1);
  }
  names_.erase(names_.begin() + chr_id);
  sequences_.erase(sequences_.begin() + chr_id);
  reader.Close();
}

bool References::GetSequence(
    const int& chr_id, 
    const int& pos, 
    const int& length,
    string* seq) const {
  
  if (chr_id >= static_cast<int>(names_.size())) {
    cerr << "ERROR: The reference id is invalid." 
         << endl
         << "       The total references is " << names_.size() 
	 << "; the request is" << chr_id << "."  
	 << endl;
    return false;
  }

  if (pos >= static_cast<int>(sequences_[chr_id].size())) {
    cerr << "ERROR: The requested position is invalid." 
         << endl
         << "       The length of reference is " << sequences_[chr_id].size() 
	 << "; the request is" << pos << "."  
	 << endl;
    return false;
  }

  int remain = sequences_[chr_id].size() - pos;
  int len = (remain < length) ? remain : length;

  *seq = sequences_[chr_id].substr(pos, len);
  return true;
}

