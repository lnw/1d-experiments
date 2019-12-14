#ifndef LOCAL_BASIS_1D_HH
#define LOCAL_BASIS_1D_HH

#include <vector>
#include <tuple>

#include "tensor.hh"
// #include "grid_1d.hh"
#include "integrate.hh"
#include "basis_1d.hh"

using namespace std;

class localbasis1d: public basis1d {
public: 
  // can't use vector<tensor_square> because the cells might be overlapping.  fuck it, we grow the basis and ignore linear dependencies
  // meaning: the first row of basis[0] contains the indices for obtaining the first local orbital as a linear combination from LIP. etc.
  // i.e.: transforms local funcs to their LIP expansion: loc*trafo = LIP
  vector<tensor_square> basis;
  bool has_basis_inv; // did we store the basis back trafo yet?
  // i.e.: transforms LIP coeffs to local funcs: lip*trafo = loc
  vector<tensor_square> basis_inv;
  bool has_overlap; // did we store the overlap yet?
  vector<tensor_square_sym> overlap; // N x n_funcs_per_cell x n_funcs_per_cell
  bool has_S_invsqrt; // did we store the overlap yet?
  vector<tensor_square_sym> S_invsqrt; // N x n_funcs_per_cell x n_funcs_per_cell

  localbasis1d(const system1d& sys);
  localbasis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim);

  // void calc_and_set_norms_sqrt(const system1d& sys);
  void calc_and_set_overlap(const system1d& sys);
  void calc_and_set_overlap_invsqrt();
  void calc_and_set_basis_inv();

  bool is_orthogonal() const;
  bool is_normal() const;

  friend ostream& operator<<(ostream& S, const localbasis1d& LO) {
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

