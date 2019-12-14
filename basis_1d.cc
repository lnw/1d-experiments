
#include <vector>
#include <tuple>

// #include "tensor.hh"
// #include "grid_1d.h"
// #include "integrate.hh"
#include "basis_1d.hh"
#include "system_1d.hh"

using namespace std;

basis1d::basis1d(const system1d& sys): basis1d(sys.overlapping, sys.n_cells, sys.n_funcs_per_cell) {};


basis1d::basis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim):
  n_cells(c_dim),
  n_funcs_per_cell(p_dim),
  n_funcs(overlapping_ ? c_dim*(p_dim-1)+1 : c_dim*p_dim),
  overlapping(overlapping_)//,
  // basis(n_cells),
  {};

