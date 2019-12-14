#ifndef MO_1D_HH
#define MO_1D_HH

#include <vector>

#include "tensor.hh"

using namespace std;

class localbasis1d;

class mo1d{
public: 
  size_t dim; // number of functions
  double orbital_e;
  tensor_1 orbital_coeffs; // can't use tensor_2 because the cells might be overlapping

  mo1d(const bool overlapping, const size_t c_dim, const size_t p_dim): dim(overlapping ? c_dim*(p_dim-1)+1 : c_dim*p_dim), orbital_e(-1.0), orbital_coeffs({dim}){};

  void guess(const system1d& sys, const localbasis1d& lo1d);
  // void solve(const system1d& sys, const localbasis1d& lo1d, const double thr);
  // void do_helmholtz_scf_step(const system1d& sys);

  void update_e(const system1d& sys, const localbasis1d& lo1d);
  void normalise(const system1d& sys, const localbasis1d& lo1d);


  friend ostream& operator<<(ostream& S, const mo1d& T) {
    S << "{";
    S << T.dim << ", ";
    S << T.orbital_e << ", ";
    S << T.orbital_coeffs;
    S << "}";
    return S;
  }

};


#endif

