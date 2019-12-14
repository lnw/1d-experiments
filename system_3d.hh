#ifndef SYSTEM_3D_HH
#define SYSTEM_3D_HH

#include <cmath>
#include <vector>

#include "auxiliary.hh"
#include "polynomial.hh"
#include "grid_1d.hh"
#include "tensor.hh"

using namespace std;


class system3d {
public:

  size_t n_cells1d;
  size_t n_cells3d;
  double cell_width;
  double grid_begin;
  size_t n_funcs_per_cell1d; // eg 7
  size_t n_funcs_per_cell3d;
  size_t n_funcs1d;
  size_t n_funcs3d;
  bool overlapping;
  grid1d grid; // lets assume cubes to begin with
  grid1d tgrid;
  vector<polynomial> legendre_polynomials;
  vector<polynomial> d_legendre_polynomials;
  // vector<double> pot_coeffs;
  tensor_6 pot_coeffs; // cell z, cell y, cell x, p z, p y , p x

  void calc_and_set_polynomials();
  void calc_and_set_derivatives();
  void calc_and_set_derivatives_the_hard_way();

  friend ostream& operator<<(ostream& S, const system3d& s) {
    S << "n cells: " << s.n_cells1d << "/" << s.n_cells3d <<endl;
    S << "cell width: " << s.cell_width << endl;
    S << "lower boundary: " << s.grid_begin << endl;
    S << "functions per cell: " << s.n_funcs_per_cell1d << "/" << s.n_funcs_per_cell3d << endl;
    S << "n funcs: " << s.n_funcs1d << "/" << s.n_funcs3d <<endl;
    S << "overlapping: " << s.overlapping << endl;
    S << "grid: " << s.grid << endl;
    S << "tgrid: " << s.tgrid << endl;
    S << "pot coeffs: " << s.pot_coeffs << endl;
    S << "polynomials: " << s.legendre_polynomials << endl;
    S << "derivatives of polynomials: " << s.d_legendre_polynomials;
    return S;
  }
};


#endif
