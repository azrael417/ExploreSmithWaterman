#include "fastq_reader.h"

#include <iostream>
#include <vector>

using std::cerr;
using std::endl;

namespace {
void TrimString(const char delimiter, string* readname){
  size_t pos = readname->find_first_of(delimiter);
  if (pos != string::npos)
    readname->erase(pos, readname->size() - pos);
  else
    ;// nothing
}
}

FastqReader::FastqReader()
  : file_()
  , readname_()
  , line_(0)
  , error_(false)
{}

FastqReader::~FastqReader() {
  this->Close();
}

bool FastqReader::Open(const char* filename) {
  file_.open(filename, ifstream::in);
  if (file_.good()) return true;
  else return false;
}

bool FastqReader::Close() {
  file_.close();
  return true;
}

bool FastqReader::LoadNextBatch(string** readnames,
                                string** sequences,
                                string** quals,
                                int* readsize,
                                const int& batchsize){
  *readnames = new string[batchsize];
  *sequences = new string[batchsize];
  *quals = new string[batchsize];
  
  bool IsOK = true;
  int count = 0;
  while(IsOK & (count<batchsize) ){
    IsOK = LoadNextRead(&(*readnames)[count], &(*sequences)[count], &(*quals)[count]);
    count++;
  }
  *readsize=count-1;
  
  return IsOK;
}

bool FastqReader::LoadNextRead(
    string* readname, 
    string* sequence, 
    string* qual) {
  readname->clear();
  sequence->clear();
  qual->clear();

  if (error_) return false;
  if (file_.eof()) return false;

  if (!readname_.empty()) { // readname is already get and stored in readname_
    *readname = readname_;
  } else {
    if (file_.eof()) return false;
    getline(file_, *readname);
    readname->erase(0, 1); // erase '@'
    ++line_;
  }

  if (readname->empty()) {
    cerr << "ERROR fastq: The readname is empty." << endl;
    error_ = true;
    return false;
  }

  TrimString(' ', readname);

  string temp;
  int fastq_step = 2;
  while(!file_.eof()) {
    temp.clear();
    getline(file_, temp);
    ++line_;
    
    if (!temp.empty()) { // get a line in temp
      if (fastq_step == 1) { // get the next read
        if (temp.size() == 1) {
	  cerr << "ERROR: A readname is not seen after '@'. Line: " << line_ << endl;
	  error_ = true;
	  return false;
	}
  readname_.clear();
	readname_ = temp.substr(1, temp.size() - 1);
	fastq_step = 2;
	break; // loading the current read done
      } else if (fastq_step == 2) {
        *sequence += temp;
	fastq_step = 3;
      } else if (fastq_step == 3) {
	fastq_step = 4;
	continue;
      } else if (fastq_step == 4){
	*qual += temp;
	fastq_step = 1;
      } else {
        cerr << "ERROR: Unknown fastq step. Line: " << line_ << endl;
	error_ = true;
	break;
      }
    } // end if-else
  } // end while

  if (sequence->empty()) {
    cerr << "ERROR: The sequence of " << *readname << " is empty." << endl;
    error_ = true;
    return false;
  } else if (qual->empty()) {
    cerr << "ERROR: The qual of " << *readname << " is empty." << endl;
    error_ = true;
    return false;
  } else {
    return true;
  }
}
