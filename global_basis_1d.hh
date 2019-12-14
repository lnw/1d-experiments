#ifndef GLOBAL_BASIS_1D_HH
#define GLOBAL_BASIS_1D_HH

#include <vector>
#include <tuple>

#include "timing.hh"
#include "tensor.hh"
#include "integrate.hh"
#include "basis_1d.hh"

using namespace std;

class globalbasis1d: public basis1d {
public: 
  tensor_square basis; // can't use vector<tensor_square> because the cells might be overlapping.  fuck it, we grow the basis and ignore linear dependencies
  bool has_overlap; // did we store the overlap yet?
  tensor_square_sym overlap; // n_funcs x n_funcs
  bool has_S_invsqrt_sqrt; // did we store the overlap yet?
  tensor_square_sym S_invsqrt; // n_funcs x n_funcs
  tensor_square_sym S_sqrt; // n_funcs x n_funcs

  globalbasis1d(){};
  globalbasis1d(const system1d& sys);
  globalbasis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim);

  void calc_and_set_overlap(const system1d& sys);
  void calc_and_set_overlap_invsqrt_sqrt(timing& times);

  bool is_orthogonal() const;
  bool is_normal() const;

  friend ostream& operator<<(ostream& S, const globalbasis1d& LO) {
    S << "{";
    S << LO.n_cells << ", ";
    S << LO.n_funcs_per_cell << ", ";
    S << LO.n_funcs << ", ";
    S << LO.basis << ", ";
    //S << LO.overlap;
    S << "}";
    return S;
  }

};


#endif

