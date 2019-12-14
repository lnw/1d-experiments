#ifndef MO_3D_HH
#define MO_3D_HH

#include <vector>

#include "tensor.hh"

using namespace std;


class mo3d{
public: 
  size_t dim1d, dim3d; // number of functions
  double orbital_e;
  // vector<double> orbital_coeffs; // cells in the order z, y, x.  Then within each cell, z, y, x // one value per grid point
  tensor_6 orbital_coeffs; // cells in the order z, y, x.  Then within each cell, z, y, x // one value per grid point

  mo3d(size_t c_dim, size_t p_dim): dim1d(c_dim*p_dim), dim3d(dim1d*dim1d*dim1d), orbital_e(0.0), orbital_coeffs(tensor_6({c_dim, c_dim, c_dim, p_dim, p_dim, p_dim})){};

  void guess(const system3d& sys);
  void solve(const system3d& sys, const double thr, const size_t max_cycle);
  void do_helmholtz_scf_step_slow(const system3d& sys);
  void do_helmholtz_scf_step(const system3d& sys);
  tensor_6 poisson_slow(const system3d& sys) const;
  tensor_6 poisson(const system3d& sys) const;
  void update_e(const system3d& sys);
  void normalise(const system3d& sys);

  friend ostream& operator<<(ostream& S, const mo3d& T) {
    S << "{";
    S << T.dim1d << "/" << T.dim3d << ", ";
    S << T.orbital_e << ", ";
    S << T.orbital_coeffs;
    S << "}";
    return S;
  }

};

#endif

