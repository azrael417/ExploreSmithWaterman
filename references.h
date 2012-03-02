#ifndef UTIL_REFERENCES_H_
#define UTIL_REFERENCES_H_

#include <iostream>
#include <string>
#include <vector>

using std::string;

class References {
 public:
  References();
  void LoadReferences(const char* filename);
  bool GetSequence(const int& chr_id, const int& pos, const int& length, 
                   string* seq) const;
  inline int GetReferenceCount() const;
  inline const char* GetReferenceSequence(const int& id) const;
  inline const char* GetReferenceSequence(const int& id, int* length) const;

 private:
  std::vector<string> names_;
  std::vector<string> sequences_;
};

int References::GetReferenceCount() const {
  return static_cast<int>(names_.size());
}

inline const char* References::GetReferenceSequence(const int& id) const {
  if (id >= static_cast<int>(sequences_.size())) {
    std::cerr << "ERROR: The requested id(" << id << ") is invalid."<< std::endl;
    return NULL;
  } else {
    return sequences_[id].c_str();
  }
}

inline const char* References::GetReferenceSequence(const int& id, int* length) const {
  const char* p_ref = GetReferenceSequence(id);
  if (p_ref != NULL) *length = sequences_[id].size();
  else length = 0;

  return p_ref;
}
#endif //UTIL_REFERENCES_H_
