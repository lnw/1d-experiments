#ifndef AUXILIARY_HH
#define AUXILIARY_HH

#include <cstdlib>
#include <cassert>
#include <iomanip> // required for setfill()
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <set>

using namespace std;

template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p);

#define container_output(container) \
  template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
  { \
  s << "{"; \
  for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
    s << *x; \
    if(++x!=v.end()) s << ","; \
  } \
  s << "}"; \
  return s; \
  }
container_output(vector);
container_output(set);

template <typename S, typename T>
ostream& operator<<(ostream& s, const pair<S,T>& p) {
  s << "{" << p.first << "," << p.second << "}";
  return s;
}

#define container_product(container) \
  template <typename T> T product(const container<T>& v) { \
  T p = 1; \
  for(auto x: v) p *= x;  \
  return p; \
  }
container_product(initializer_list);
container_product(vector);

template <typename T>
std::string to_string_with_precision(const T val, const int p=4) {
  std::stringstream ss;
  ss << std::setprecision(p) << val;
  return ss.str();
}

template <typename T>
std::string to_string_fixedwidth(const T val, const int n=3) {
  std::stringstream ss;
  ss << std::setw(n) << std::setfill('0') << val;
  return ss.str();
}

template <typename T>
double to_double(const T& s) {
  std::stringstream ss;
  double result;
  ss << s;
  ss >> result;
  return result;
}

template <typename T>
int to_int(const T& s) {
  std::stringstream ss;
  int result;
  ss << s.raw();
  ss >> result;
  return result;
}

template <typename T>
size_t to_st(const T& s) {
  std::stringstream ss;
  size_t result;
  ss << s.raw();
  ss >> result;
  return result;
}

#endif // auxiliary
