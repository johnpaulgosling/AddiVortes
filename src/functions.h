#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <algorithm>

inline bool in_vector(int value, const std::vector<int>& vec) {
  return std::find(vec.begin(), vec.end(), value) != vec.end();
}

#endif