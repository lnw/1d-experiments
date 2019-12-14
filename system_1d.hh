#ifndef SYSTEM_1D_HH
#define SYSTEM_1D_HH

#include <cmath>
#include <vector>

#include "auxiliary.hh"
#include "polynomial.hh"
#include "grid_1d.hh"
#include "tensor.hh"
#include "global_basis_1d.hh"

using namespace std;

// class globalbasis1d;

class system1d {
public:

  size_t n_cells;
  double cell_width;
  double grid_begin;
  size_t n_funcs_per_cell; // eg 7
  size_t n_gridpoints; // = n_cells * (size(grid) - overlapping) + overlapping
  size_t n_gridpoints_fine; // = n_cells * (size(grid_fine) - overlapping) + overlapping
  size_t n_gp_per_cell;
  size_t n_gp_fine_per_cell;
  bool overlapping;
  grid1d grid; // per cell (if applicable)
  grid1d grid_fine; // two points more, used when order 2n+1 is not enough
  // grid1d tgrid;
  vector<polynomial> legendre_polynomials;
  vector<polynomial> d_legendre_polynomials;
  // vector<polynomial> legendre_polynomials_normed;
  // vector<polynomial> d_legendre_polynomials_normed;
  globalbasis1d globbas1d_plain;
  // globalbasis1d globbas1d_normalised;
  tensor_1 pot_coeffs; // n_cells * n_funcs_per_cell  or  n_cells * (funcs_per_cel-1) + 1

  // call the four functions in this order
  void calc_and_set_polynomials();
  void calc_and_set_derivatives();
  // void calc_and_set_normalised_polynomials();
  // void calc_and_set_normalised_derivatives();

  void calc_and_set_plain_basis(timing& times);
  // void calc_and_set_normalised_basis();

  friend ostream& operator<<(ostream& S, const system1d& s) {
    S << "n cells: " << s.n_cells << endl;
    S << "cell width: " << s.cell_width << endl;
    S << "lower boundary: " << s.grid_begin << endl;
    S << "functions per cell: " << s.n_funcs_per_cell << endl;
    S << "overlapping: " << s.overlapping << endl;
    S << "grid: " << setprecision(12) << s.grid << endl;
    S << "grid fine: " << s.grid_fine << endl;
    // S << "tgrid: " << s.tgrid << endl;
    S << "pot coeffs: " << s.pot_coeffs << endl;
    S << "pol= " << s.legendre_polynomials << endl;
    S << "dpol= " << s.d_legendre_polynomials;
    // S << "poln= " << s.legendre_polynomials_normed << endl;
    // S << "dpoln= " << s.d_legendre_polynomials_normed;
    // S << "polynomials= " << s.legendre_polynomials << endl;
    // S << "derivatives of polynomials= " << s.d_legendre_polynomials;
    // S << "polynomials= " << s.legendre_polynomials << endl;
    // S << "derivatives of polynomials= " << s.d_legendre_polynomials;
    return S;
  }
};


#endif
