#ifndef UTIL_REFERENCES_H_
#define UTIL_REFERENCES_H_

#include <string>
#include <vector>

using std::string;

class References {
 public:
  References();
  void LoadReferences(const char* filename);
  bool GetSequence(const int& chr_id, const int& pos, const int& length, 
                   string* seq) const;
 private:
  std::vector<string> names_;
  std::vector<string> sequences_;
};

#endif //UTIL_REFERENCES_H_
