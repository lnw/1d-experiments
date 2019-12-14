
#include <cmath>
#include <vector>

#include "auxiliary.hh"
#include "aux_math.hh"
#include "polynomial.hh"
#include "grid_1d.hh"
#include "integrate.hh"
#include "system_1d.hh"

using namespace std;

// the polynomials in one cell, centred at zero (needs to be shifted, when used)
void system1d::calc_and_set_polynomials(){
  legendre_polynomials = vector<polynomial>(n_funcs_per_cell);
  for (size_t j=0; j<n_funcs_per_cell; j++){
    legendre_polynomials[j].resize(n_funcs_per_cell); // 0 -- n-1
    double den = 1;
    for (size_t m=0; m<n_funcs_per_cell; m++){
      if(m==j) continue;
      den *= grid.p[j] - grid.p[m];
    }
    // cout << "den: " << den << endl;
    for (size_t x=0; x<n_funcs_per_cell; x++){ // coeffs from 0 -- n-1
      legendre_polynomials[j][x] = sum_x_over_n_except_j(n_funcs_per_cell-x, vector<size_t>({j}), grid.p) / den;
    }
  }

  // for(int i=0; i<n_funcs_per_cell; i++) legendre_polynomials[i] /= sqrt(grid.w[i]);
}


void system1d::calc_and_set_derivatives(){
  d_legendre_polynomials.resize(n_funcs_per_cell);
  for(int i=0; i<n_funcs_per_cell; i++){
    d_legendre_polynomials[i] = legendre_polynomials[i].derive();
  }
}


// void system1d::calc_and_set_normalised_polynomials(){
//   legendre_polynomials_normed = legendre_polynomials; // set first, so the integration succeeds, and overwrite after that
//   const int c = 0; // irrelevant
//   for(int i=0; i<n_funcs_per_cell; i++){
//     const double norm = integrate_funcs(*this, c, i, i, NONE);
//     const double norm_sqrt = sqrt(norm);
//     legendre_polynomials_normed[i] = legendre_polynomials[i] / norm_sqrt;
//   }
// }


// void system1d::calc_and_set_normalised_derivatives(){
//   d_legendre_polynomials_normed.resize(n_funcs_per_cell);
//   const int c = 0; // irrelevant
//   for(int i=0; i<n_funcs_per_cell; i++){
//     const double norm = integrate_funcs(*this, c, i, i, NONE);
//     const double norm_sqrt = sqrt(norm);
//     d_legendre_polynomials_normed[i] = legendre_polynomials_normed[i].derive();
//   }
// }

void system1d::calc_and_set_plain_basis(timing& times){
  globbas1d_plain = globalbasis1d(*this);
  globbas1d_plain.calc_and_set_overlap(*this);
  // cout << globbas1d_plain.overlap << endl;
  globbas1d_plain.calc_and_set_overlap_invsqrt_sqrt(times);
  // cout << globbas1d_plain.S_invsqrt << endl;
}

