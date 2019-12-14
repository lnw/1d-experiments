#ifndef GRID_1D_HH
#define GRID_1D_HH

#include <vector>
#include <iostream>
#include <algorithm>

typedef enum {EQUIDISTANT, GAUSSIAN_LIN, GAUSSIAN_LOG, TGRID, GAUSS_LOBATTO} grid_t;


class grid1d{

public:
  vector<double> p;
  vector<double> w;

  grid1d(){}
  // vector of bounds, vector of point counts, grid type
  grid1d(const vector<double> bounds, const vector<size_t> n, const grid_t g);
  size_t size() const {return p.size();}

  grid1d  operator+(const grid1d& g) const { return grid1d(*this) += g; }
  grid1d& operator+=(const grid1d& g) {
    p.insert(p.end(), g.p.begin(), g.p.end());
    w.insert(w.end(), g.w.begin(), g.w.end());
    return *this;
  }

  friend ostream& operator<<(ostream& S, const grid1d& G) {
    S << "{" << G.p << ", " << G.w << "}";
    return S;
  }

};

#endif


