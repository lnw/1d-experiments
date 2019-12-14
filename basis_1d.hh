#ifndef BASIS_1D_HH
#define BASIS_1D_HH

#include <vector>
#include <tuple>
#include <iostream>

using namespace std;

class system1d;

class basis1d{
public: 
  size_t n_cells;
  size_t n_funcs_per_cell;
  size_t n_funcs; // total number of functions
  bool overlapping;
  // vector<tensor_square> basis; // can't use vector<tensor_square> because the cells might be overlapping.  fuck it, we grow the basis and ignore linear dependencies
  // bool has_overlap; // did we store the overlap yet?
  // tensor_square_sym overlap; // n_funcs x n_funcs
  // bool has_S_invsqrt; // did we store the overlap yet?
  // tensor_square_sym S_invsqrt; // n_funcs x n_funcs

  basis1d(){};
  basis1d(const system1d& sys);
  basis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim);

  // void calc_and_set_norms_sqrt(const system1d& sys);
  // void calc_and_set_overlap(const system1d& sys);
  // tuple<tensor_square_sym,tensor_square_sym> get_overlap_invsqrt_sqrt() const;

  bool is_orthogonal() const;
  bool is_normal() const;

  friend ostream& operator<<(ostream& S, const basis1d& B) {
    S << "{";
    S << B.n_cells << ", ";
    S << B.n_funcs_per_cell << ", ";
    S << B.n_funcs << ", ";
    // S << B.weights << ", ";
    // S << B.basis << ", ";
    S << B.overlapping;
    S << "}";
    return S;
  }

};

#endif

