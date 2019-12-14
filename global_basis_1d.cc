
#include <vector>

#include "tensor.hh"
// #include "grid_1d.hh"
#include "global_basis_1d.hh"
#include "system_1d.hh"

using namespace std;


globalbasis1d::globalbasis1d(const system1d& sys): globalbasis1d(sys.overlapping, sys.n_cells, sys.n_funcs_per_cell) {};

globalbasis1d::globalbasis1d(const bool overlapping_, const size_t c_dim, const size_t p_dim):
    basis1d(overlapping_, c_dim, p_dim),
    basis({n_funcs,n_funcs}),
    has_overlap(false),
    overlap({n_funcs, n_funcs}),
    has_S_invsqrt_sqrt(false),
    S_invsqrt({n_funcs, n_funcs}),
    S_sqrt({n_funcs, n_funcs}) {
  for(int c=0; c<n_cells; c++){
    for(int f1=0; f1<n_funcs_per_cell; f1++){
      const size_t gi = c*(n_funcs_per_cell-overlapping) + f1;
      basis(gi,gi) = 1.0;
      // for(int f2=0; f2<n_funcs_per_cell; f2++){
        // const size_t gj = c*(n_funcs_per_cell-overlapping) + f2;
        // weights(gi) = grid.w[p];
        // overlap(gi,gj) = integrate_funcs( sys, sys.lo1d, c, f1, f2, NONE);
      // }
    }
  }
};

bool globalbasis1d::is_orthogonal() const {
  const double eps = 1.e-12;
  tensor_square prod = matmul(basis, basis, 'N', 'N');
  for (int i=0; i<n_funcs; i++){
    for (int j=0; j<n_funcs; j++){
      if(i==j) continue;
      if(abs(prod(i,j)) > eps) return false;
    }
  }
  return true;
}


bool globalbasis1d::is_normal() const {
  const double eps = 1.e-12;
  tensor_square prod = matmul(basis, basis, 'N', 'N');
  for(int i=0; i<n_funcs; i++){
    if(abs(prod(i,i)-1.0) > eps) return false;
  }
  return true;
}


void globalbasis1d::calc_and_set_overlap(const system1d& sys){
  for(int c=0; c<n_cells; c++){
    for(int f1=0; f1<n_funcs_per_cell; f1++){
      for(int f2=0; f2<n_funcs_per_cell; f2++){
        const size_t gi = c*(n_funcs_per_cell-overlapping) + f1;
        const size_t gj = c*(n_funcs_per_cell-overlapping) + f2;
        overlap(gi,gj) += integrate_funcs(sys, c, f1, f2, NONE);
      }
    }
  }
  has_overlap = true;
}


void globalbasis1d::calc_and_set_overlap_invsqrt_sqrt(timing& times) {
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = overlap.eigensystem_sym(times);

  tensor_square_sym evals_mat({n_funcs, n_funcs});
 
  for(int i=0; i<n_funcs; i++) evals_mat(i,i) = 1.0/sqrt(evals[i]);
  tensor_square aux1( matmul(evals_mat, evecs, 'N', 'T'));
  S_invsqrt = matmul( evecs, aux1, 'N', 'N');

  for(int i=0; i<n_funcs; i++) evals_mat(i,i) = sqrt(evals[i]);
  tensor_square aux2( matmul(evals_mat, evecs, 'N', 'T'));
  S_sqrt = matmul( evecs, aux2, 'N', 'N');

  has_S_invsqrt_sqrt = true;
}

