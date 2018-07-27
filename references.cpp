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
  
  //temporary variables
  std::vector<string> tmpnames_, tmpsequences_;
  tmpnames_.reserve(100);
  tmpsequences_.reserve(100);

  int chr_id = 0;
  max_sequence_length_ = 0;
  tmpnames_.resize(chr_id + 1);
  tmpsequences_.resize(chr_id + 1);
  while (reader.LoadNextRead(&tmpnames_[chr_id], &tmpsequences_[chr_id])){
    //compute max name length:
    max_name_length_ = tmpnames_[chr_id].size() > max_name_length_ ? tmpnames_[chr_id].size() : max_name_length_;
    //compute max sequence length:
    max_sequence_length_ = tmpsequences_[chr_id].size() > max_sequence_length_ ? tmpsequences_[chr_id].size() : max_sequence_length_;
    ++chr_id;
    tmpnames_.resize(chr_id + 1);
    tmpsequences_.resize(chr_id + 1);
  }
  tmpnames_.erase(tmpnames_.begin() + chr_id);
  tmpsequences_.erase(tmpsequences_.begin() + chr_id);
  
  //create view:
  names_ = View2D<char, Kokkos::LayoutRight, Kokkos::HostSpace>("names_", tmpnames_.size(), max_name_length_);
  sequences_end_ = View1D<int, Kokkos::HostSpace>("sequences_end_", tmpsequences_.size());
  sequences_ = View2D<char, Kokkos::LayoutRight, Kokkos::HostSpace>("sequences_", tmpsequences_.size(), max_sequence_length_);
  
  //fill
  for(unsigned int i=0; i<tmpnames_.size(); i++){
    for(unsigned int j=0; j<tmpnames_[i].size(); j++){
      names_(i,j) = tmpnames_[i][j];
    }
    if(tmpnames_[i].size() < max_name_length_) names_(i,tmpnames_[i].size()) = 0;
  }
  for(unsigned int i=0; i<tmpsequences_.size(); i++){
    for(unsigned int j=0; j<tmpsequences_[i].size(); j++){
      sequences_(i,j) = tmpsequences_[i][j];
    }
    if(tmpsequences_[i].size() < max_sequence_length_){
      sequences_(i,tmpsequences_[i].size()) = 0;
      sequences_end_(i) = tmpsequences_[i].size();
    }
    else sequences_end_(i) = max_sequence_length_;
  }

  //create mirrors and upload
  d_sequences_end_ = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), sequences_end_);
  d_sequences_ = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), sequences_);
  
  reader.Close();
}

//bool References::GetSequence(
//    const int& chr_id, 
//    const int& pos, 
//    const int& length,
//    string* seq) const {
//  
//  if (chr_id >= static_cast<int>(names_.extent(0))) {
//    cerr << "ERROR: The reference id is invalid." 
//         << endl
//         << "       The total references is " << names_.extent(0)
//	 << "; the request is" << chr_id << "."  
//	 << endl;
//    return false;
//  }
//    
//  if (pos >= static_cast<int>(sequences_end_(chr_id))) {
//    cerr << "ERROR: The requested position is invalid." 
//         << endl
//         << "       The length of reference is " << sequences_end_(chr_id)
//	       << "; the request is" << pos << "."  
//        << endl;
//    return false;
//  }
//
//  int remain = sequences_end_(chr_id) - pos;
//  int len = (remain < length) ? remain : length;
//
//  ViewToString(*seq, Kokkos::subview(sequences_, chr_id, Kokkos::ALL));
//  *seq = (*seq).substr(pos, len);
//  return true;
//}

