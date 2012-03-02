#include "fasta_reader.h"

#include <iostream>

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

FastaReader::FastaReader()
  : file_()
  , readname_()
  , line_(0)
  , error_(false)
{}

FastaReader::~FastaReader() {
  this->Close();
}

bool FastaReader::Open(const char* filename) {
  file_.open(filename, ifstream::in);
  if (file_.good()) return true;
  else return false;
}

bool FastaReader::Close() {
  file_.close();
  return true;
}

bool FastaReader::LoadNextRead(string* readname, string* sequence) {
  readname->clear();
  sequence->clear();

  if (error_) return false;
  if (file_.eof()) return false;

  if (!readname_.empty()) { // readname is already get and stored in readname_
    *readname = readname_;
  } else {
    if (file_.eof()) return false;
    getline(file_, *readname);
    readname->erase(0, 1); // erase '>'
    ++line_;
  }

  if (readname->empty()) {
    cerr << "ERROR: The readname is empty." << endl;
    error_ = true;
    return false;
  }

  TrimString(' ', readname);

  string temp;
  while(!file_.eof()) {
    temp.clear();
    getline(file_, temp);
    ++line_;
    
    if (!temp.empty()) { // get a line in temp
      if (temp[0] == '>') { // stop
        if (temp.size() == 1) {
	  cerr << "ERROR: A readname is not seen after '>'. Line: " << line_ << endl;
	  error_ = true;
	  return false;
	}

        readname_.clear();
	readname_ = temp.substr(1, temp.size() - 1);
	break;
      }

      *sequence += temp;
    } // end if-else
  } // end while

  if (sequence->empty()) {
    cerr << "ERROR: The sequence of " << *readname << " is empty." << endl;
    error_ = true;
    return false;
  } else {
    return true;
  }
}
