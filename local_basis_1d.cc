
#include <vector>
#include <tuple>

#include "tensor.hh"
#include "grid_1d.hh"
#include "local_basis_1d.hh"
#include "system_1d.hh"

using namespace std;

localbasis1d::localbasis1d(const system1d& sys): localbasis1d(sys.overlapping, sys.n_cells, sys.n_funcs_per_cell) {};

localbasis1d::localbasis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim):
    basis1d(overlapping_, c_dim, p_dim),
    basis(n_cells),
    has_basis_inv(false),
    basis_inv(n_cells),
    has_overlap(false),
    overlap(n_cells),
    has_S_invsqrt(false),
    S_invsqrt(n_cells) {
};

bool localbasis1d::is_orthogonal() const {
  const double eps = 1.e-12;
  for(int c=0; c<n_cells; c++){
    tensor_square prod = matmul(basis[c], basis[c], 'N', 'N');
    for (int i=0; i<n_funcs_per_cell; i++){
      for (int j=0; j<n_funcs_per_cell; j++){
        if(i==j) continue;
        if(abs(prod(i,j)) > eps) return false;
      }
    }
  }
  return true;
}


bool localbasis1d::is_normal() const {
  const double eps = 1.e-12;
  for(int c=0; c<n_cells; c++){
    tensor_square prod = matmul(basis[c],basis[c],'N','N');
    for(int i=0; i<n_funcs_per_cell; i++){
      if(abs(prod(i,i)-1.0) > eps) return false;
    }
  }
  return true;
}


void localbasis1d::calc_and_set_overlap(const system1d& sys){
  for(int c=0; c<n_cells; c++){
    for(int f1=0; f1<n_funcs_per_cell; f1++){
      for(int f2=0; f2<n_funcs_per_cell; f2++){
        const size_t gi = c*(n_funcs_per_cell-overlapping) + f1;
        const size_t gj = c*(n_funcs_per_cell-overlapping) + f2;
assert(false); // use local basis
        overlap[c](gi,gj) += integrate_funcs(sys, c, f1, f2, NONE);
      }
    }
  }
}


// void localbasis1d::calc_and_set_overlap_invsqrt() {
//   for(int c=0; c<n_cells; c++){
//     tensor_1 evals;
//     tensor_square evecs;
//     tie(evals, evecs) = overlap[c].eigensystem_sym(times);
// assert(false); // use local basis
// 
//     tensor_square_sym evals_mat({n_funcs, n_funcs});
//  
//     for(int i=0; i<n_funcs; i++) evals_mat(i,i) = 1.0/sqrt(evals[i]);
//     tensor_square aux1( matmul( evals_mat, evecs, 'N', 'T' ));
//     S_invsqrt[c] = matmul( evecs, aux1, 'N', 'N' );
// 
//   // for(int i=0; i<n_funcs; i++) evals_mat(i,i) = sqrt(evals[i]);
//   // tensor_square aux2( matmul(evals_mat, evecs, 'N', 'T'));
//   // tensor_square_sym S_sqrt( matmul( evecs, aux2, 'N', 'N'));
//   }
//   has_S_invsqrt = true;
// }


void localbasis1d::calc_and_set_basis_inv() {
  for(int c=0; c<n_cells; c++){
    basis_inv[c] = basis[c].inverse();
  }
  has_basis_inv = true;
}


