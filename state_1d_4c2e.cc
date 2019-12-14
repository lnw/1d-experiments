#include <vector>
#include <tuple>
#include <algorithm>

#include "timing.hh"
#include "tensor.hh"
#include "local_basis_1d.hh"
#include "system_1d.hh"
#include "state_1d.hh"

using namespace std;


// calc (ij|kl) where i,j,k,l are MOs which are expanded in (mu nu|lambda sigma) which are LIP
double get_4c2eintegral_mo_local_ortho_grid(const system1d& sys, const state1d& S, size_t mo_i, size_t mo_j, size_t mo_k, size_t mo_l, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);

  const double a = 0.15;
  double int_4c2e = 0;

  vector<tensor_1> D_mu(sys.n_cells); // expansion of i in LIP
  vector<tensor_1> D_nu(sys.n_cells); // expansion of j in LIP
  vector<tensor_1> D_lambda(sys.n_cells); // expansion of k in LIP
  vector<tensor_1> D_sigma(sys.n_cells); // expansion of l in LIP
  for(int c=0; c<sys.n_cells; c++){
    tensor_1 C_mu({sys.n_funcs_per_cell}); // expansion of i in local orbitals
    tensor_1 C_nu({sys.n_funcs_per_cell}); // expansion of j in local orbitals
    tensor_1 C_lambda({sys.n_funcs_per_cell}); // expansion of k in local orbitals
    tensor_1 C_sigma({sys.n_funcs_per_cell}); // expansion of l in local orbitals
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      C_mu[i] = S.local_mo_coeffs(mo_i, c*sys.n_funcs_per_cell+i);
      C_nu[i] = S.local_mo_coeffs(mo_j, c*sys.n_funcs_per_cell+i);
      C_lambda[i] = S.local_mo_coeffs(mo_k, c*sys.n_funcs_per_cell+i);
      C_sigma[i] = S.local_mo_coeffs(mo_l, c*sys.n_funcs_per_cell+i);
    }
    D_mu[c] = basis[c] * C_mu;
    D_nu[c] = basis[c] * C_nu;
    D_lambda[c] = basis[c] * C_lambda;
    D_sigma[c] = basis[c] * C_sigma;
  }

  for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
    // cout << "cell mu nu: " << cell_munu << endl;

    tensor_1 int_lambda_sigma({sys.n_gp_fine_per_cell});
    // for every point in int_lambda_mu (which uses the fine grid)
    for(int k=0; k<sys.n_gp_fine_per_cell; k++){ // grid point index in cell
      // const double r1 = sys.grid_begin + (i_c+0.5)*sys.cell_width + sys.grid_fine.p[i_p];
      const double r1 = cell_munu * sys.cell_width + sys.grid_fine.p[k];

      // for(int cell_lasi=cell_munu; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
      for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
        tensor_square int_lambda_sigma_tmp({sys.n_funcs_per_cell, sys.n_funcs_per_cell});

        for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
          for(int sigma=lambda; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
            for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
              // const double r2 = sys.grid_begin + (cell_lasi+0.5)*sys.cell_width + sys.grid_fine.p[k_p];
              const double r2 = cell_lasi * sys.cell_width + sys.grid_fine.p[k_p];
              const double r12 = r1-r2;
              // cout << k << ", " << lambda << ", " << sigma << ", " << r12 << endl;
              // cout << "  " << sys.grid_fine.w[k_p] * exp(-a*r12*r12) * polynomials_fine_grid[lambda][k_p] * polynomials_fine_grid[sigma][k_p] << endl;
              const double int_lambda_sigma_tmp_tmp = sys.grid_fine.w[k_p] * exp(-a*r12*r12) *
                                                      polynomials_fine_grid[lambda][k_p] *
                                                      polynomials_fine_grid[sigma][k_p];
              int_lambda_sigma_tmp(lambda, sigma) += int_lambda_sigma_tmp_tmp;
              if(lambda!=sigma) int_lambda_sigma_tmp(sigma, lambda) += int_lambda_sigma_tmp_tmp;
            }
          } // sigma
        } // lambda
        // cout << int_lambda_sigma_tmp << endl;
        const int fac = 1;//cell_munu==cell_lasi ? 1 : 2;
        int_lambda_sigma[k] += D_lambda[cell_lasi] * ( int_lambda_sigma_tmp * D_sigma[cell_lasi] ) * fac;
      } // cell lambda/sigma
    } // k
    // cout << int_lambda_sigma << endl;

    tensor_square int_two({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
      for(int nu=mu; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
        for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
          // cout <<  "  " << sys.grid_fine.w[k_p] * int_lambda_sigma[k_g] * polynomials_fine_grid[mu][k_p] * polynomials_fine_grid[nu][k_p] << endl;
          const double int_two_tmp = sys.grid_fine.w[k_p] * int_lambda_sigma[k_p] *
                                     polynomials_fine_grid[mu][k_p] *
                                     polynomials_fine_grid[nu][k_p];
          int_two(mu, nu) += int_two_tmp;
          if(mu!=nu) int_two(nu, mu) += int_two_tmp;
        }
      } // nu
    } // mu

    int_4c2e += D_mu[cell_munu] * ( int_two * D_nu[cell_munu] );
  }

  auto t1 = Clock::now();
  times.append_item("4c2e integral (local, ortho, grid)", t1 - t0);
  times.ind_down();

  return int_4c2e;
}


// calc (ij|kl) where i,j,k,l are MOs which are expanded in (mu nu|lambda sigma) which are LIP
double get_4c2eintegral_mo_precalc(const system1d& sys, const state1d& S, size_t mo_i, size_t mo_j, size_t mo_k, size_t mo_l, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_lip_4c2e);

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const tensor_6& lip_ints = S.lip_4c2e;

  vector<tensor_1> D_mu(sys.n_cells); // expansion of i in LIP
  vector<tensor_1> D_nu(sys.n_cells); // expansion of j in LIP
  vector<tensor_1> D_lambda(sys.n_cells); // expansion of k in LIP
  vector<tensor_1> D_sigma(sys.n_cells); // expansion of l in LIP
  for(int c=0; c<sys.n_cells; c++){
    tensor_1 C_mu({sys.n_funcs_per_cell}); // expansion of i in local orbitals
    tensor_1 C_nu({sys.n_funcs_per_cell}); // expansion of j in local orbitals
    tensor_1 C_lambda({sys.n_funcs_per_cell}); // expansion of k in local orbitals
    tensor_1 C_sigma({sys.n_funcs_per_cell}); // expansion of l in local orbitals
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      C_mu[i] = S.local_mo_coeffs(mo_i, c*sys.n_funcs_per_cell+i);
      C_nu[i] = S.local_mo_coeffs(mo_j, c*sys.n_funcs_per_cell+i);
      C_lambda[i] = S.local_mo_coeffs(mo_k, c*sys.n_funcs_per_cell+i);
      C_sigma[i] = S.local_mo_coeffs(mo_l, c*sys.n_funcs_per_cell+i);
    }
    D_mu[c] = basis[c] * C_mu;
    D_nu[c] = basis[c] * C_nu;
    D_lambda[c] = basis[c] * C_lambda;
    D_sigma[c] = basis[c] * C_sigma;
  }

  double int_4c2e = 0;
  for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
    //for(int cell_lasi=cell_munu; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
    for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
      const int fac_cell = 1;//cell_munu==cell_lasi ? 1 : 2;
      for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
        double tmp_nu=0;
        for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
          const int fac_nu = 1;//nu==mu ? 1 : 2;
#if 1
          double tmp_lambda=0;
          const size_t indx0 = ((((cell_munu*sys.n_cells+cell_lasi)*sys.n_funcs_per_cell+mu)*sys.n_funcs_per_cell+nu)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
          int i=0;
          for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
            double tmp_sigma=0;
            for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
              // tmp_sigma += lip_ints(cell_munu, cell_lasi, mu, nu, lambda, sigma) * D_sigma[cell_lasi][sigma];
              tmp_sigma += lip_ints[indx0+i] * D_sigma[cell_lasi][sigma];
              i++;
            }
            tmp_lambda += tmp_sigma * D_lambda[cell_lasi][lambda];
          }
#else
          const size_t indx = ((((cell_munu*sys.n_cells+cell_lasi)*sys.n_funcs_per_cell+mu)*sys.n_funcs_per_cell+nu)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
          // cout << "indx " << indx << endl;
          tensor_square_sym tmp({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
          int i=0;
          for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){
            for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){
              tmp[i] = lip_ints[indx+i];
              i++;
            }
          }
          double tmp_lambda = D_lambda[cell_lasi] * (tmp * D_sigma[cell_lasi]);
#endif
          tmp_nu += tmp_lambda * D_nu[cell_munu][nu] * fac_nu;
        }
        int_4c2e += tmp_nu * D_mu[cell_munu][mu] * fac_cell;
      }
    }
  }

  auto t1 = Clock::now();
  times.append_item("4c2e mo integral (precalc)", t1 - t0);
  times.ind_down();

  return int_4c2e;
}



// a given 4 centre 2 electron integral (ij|kl) that spans two spatial domains
// i and o can be decomposed into $(ij|kl) = (ij|kl)^(ii) + (ij|kl)^(io) +
// (ij|kl)^(oo)$, where (ij|kl)^(ii), (ij|kl)^(io), (ij|kl)^(oo) are the
// interactions within i, between i and io and within o respectively.  In order
// to optimise the orbitals an one cell in the field of their neighbouring
// cells we need (ij|kl)^(ii) + (ij|kl)^(io), but (ij|kl)^(io) cannot easily be
// calculated, so we calculate (ij|kl) - (ij|kl)^(oo)

// include that funny environment construction only if 'env', otherwise just the chosen determinents from two cells

double get_1cell_2electron_interaction(const system1d& sys, const state1d& S,
                                       const size_t cell1, const size_t cell2,
                                       const size_t loc_i, const size_t loc_j, const size_t loc_k, const size_t loc_l,
                                       const bool env, const bool normalised,
                                       timing& times){

  const int debug_level = 15;
  if(!env){
    const int mode = 0;
    const double one_cell_int = get_4c2eintegral_1cexc_precalc(sys, S, cell1, cell2, loc_i, loc_j, loc_k, loc_l, mode, normalised, debug_level, times);
    return one_cell_int;
  }
  else{ // include some 2e interaction with the HF_0 mo
    int mode = 1;
    const double tot = get_4c2eintegral_1cexc_precalc(sys, S, cell1, cell2, loc_i, loc_j, loc_k, loc_l, mode, normalised, debug_level, times);
    mode = 2;
    const double outer = get_4c2eintegral_1cexc_precalc(sys, S, cell1, cell2, -1, -1, -1, -1, mode, normalised, debug_level, times);
    if( (cell1==15 || cell2==16) && loc_i==0 && loc_j==0 && loc_k==0 && loc_l==1 )
      cout << "cell1, cell2, tot, outer: " << cell1 << ", " << cell2 << ": " << tot << ", " << outer << endl;
    return  tot - outer;
  }

}



// calc (ij|kl) where i,j,k,l are local functions which are expanded in (mu
// nu|lambda sigma) which are LIP here, the excitations are local to the cell
// 'cell' and loc_i... are the local orthogonal functions in that cell, so we
// exchange local orbitals in one cell and assume the rest of the mo to be the
// first mo.

// mode: 0 -- excitations within one cell
//       1 -- excitations within one, constant MO-0 in the remaining system
//       2 -- MO-0 in the remaining system, 0 in the cell
// normalised: assume an occupation of 1 in the cell, or take occupation from mo_0?
double get_4c2eintegral_1cexc_precalc(const system1d& sys, const state1d& S,
                                      const size_t cell1, const size_t cell2,
                                      const int loc_i, const int loc_j, const int loc_k, const int loc_l,
                                      const int mode,
                                      const bool normalised,
                                      int debug_level, // 0..20, printing
                                      timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  // setting i, j, k, l all to -1 is a magic value that sets all functions in
  // the selected cell to zero, useful for the complement
  bool zero_out;
  if(mode==2) zero_out=true;
  else if(loc_i<0 || loc_j<0 || loc_k<0 || loc_l<0) assert(false);
  else zero_out=false;

  assert(S.has_local_mo_coeffs);
  assert(S.has_lip_4c2e);

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const tensor_6& lip_ints = S.lip_4c2e;

  // first get coeffs of the lowest MO (= first local func), modify one cell later
  vector<tensor_1> D_mu(sys.n_cells); // expansion of local funcs in LIP
  const int mo_0 = 0; // the occupied local orbital
  if(mode==1 || mode==2)
  for(int c=0; c<sys.n_cells; c++){
    // cout << c << endl;
    tensor_1 C_mu({sys.n_funcs_per_cell}); // expansion of the first mo (or a local func) in local funcs
    C_mu[0] = S.local_mo_coeffs(mo_0, c*sys.n_funcs_per_cell+0); // and the rest are zero
    // cout << "t1 " << C_mu << endl;
    D_mu[c] = basis[c] * C_mu;
    // cout << "t2 " << D_mu[c] << endl;
  }
  vector<tensor_1> D_nu(D_mu);
  vector<tensor_1> D_lambda(D_mu);
  vector<tensor_1> D_sigma(D_mu);

  if(!zero_out){
    tensor_1 C_mu({sys.n_funcs_per_cell}),
             C_nu({sys.n_funcs_per_cell}),
             C_lambda({sys.n_funcs_per_cell}),
             C_sigma({sys.n_funcs_per_cell});
    double coeff1, coeff2;
    if(!normalised){
      coeff1 = S.local_mo_coeffs(mo_0, cell1*sys.n_funcs_per_cell+0); // coeff of mo0 in local funcs in cell1
      coeff2 = S.local_mo_coeffs(mo_0, cell2*sys.n_funcs_per_cell+0); // coeff of mo0 in local funcs in cell2
    }
    else{
      coeff1 = 1.0; coeff2 = 1.0;
    }
    C_mu[loc_i] = coeff1;
    C_nu[loc_j] = coeff1;
    C_lambda[loc_k] = coeff2;
    C_sigma[loc_l] = coeff2;
    D_mu[cell1] = basis[cell1] * C_mu;
    D_nu[cell1] = basis[cell1] * C_nu;
    D_lambda[cell2] = basis[cell2] * C_lambda;
    D_sigma[cell2] = basis[cell2] * C_sigma;
//    if( debug_level>=15 && (cell1==15 || cell2==16) && loc_i==0 && loc_j==0 && loc_k==0 && loc_l==1){
//      // cout << cell1*sys.n_funcs_per_cell << endl;
//      // cout << cell2*sys.n_funcs_per_cell << endl;
//      // cout << S.local_mo_coeffs << endl;
//      cout << "coeff1, coeff2: " << coeff1 << ", " << coeff2 << endl;
//      cout << "C (mo in loc), D(mo in lip) " << C_sigma << ", " << D_sigma[cell2] << endl;
//    }
  }
  else{ // zero out
    tensor_1 nada({sys.n_funcs_per_cell});
    D_mu[cell1] = nada;
    D_nu[cell1] = nada;
    D_lambda[cell2] = nada;
    D_sigma[cell2] = nada;
  }

  double int_4c2e = 0;
  if(mode==0){
    // for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
      { int cell_munu = cell1;
      // for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
        { int cell_lasi = cell2;
        for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
          double tmp_nu=0;
          for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
            // const int fac_nu = 1;//nu==mu ? 1 : 2;
            double tmp_lambda=0;
            const size_t indx0 = ((((cell_munu*sys.n_cells+cell_lasi)*sys.n_funcs_per_cell+mu)*sys.n_funcs_per_cell+nu)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
            int i=0;
            for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
              double tmp_sigma=0;
              for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
                // tmp_sigma += lip_ints(cell_munu, cell_lasi, mu, nu, lambda, sigma) * D_sigma[cell_lasi][sigma];
                tmp_sigma += lip_ints[indx0+i] * D_sigma[cell_lasi][sigma];
                i++;
              }
              tmp_lambda += tmp_sigma * D_lambda[cell_lasi][lambda];
            }
            tmp_nu += tmp_lambda * D_nu[cell_munu][nu]; // * fac_nu;
          }
          int_4c2e += tmp_nu * D_mu[cell_munu][mu];
        }
      }
    }
  }
  else if(mode==1 || mode==2){
    for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
      // { int cell_munu = cell1;
      for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
        // { int cell_lasi = cell2;
        for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
          double tmp_nu=0;
          for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
            // const int fac_nu = 1;//nu==mu ? 1 : 2;
            double tmp_lambda=0;
            const size_t indx0 = ((((cell_munu*sys.n_cells+cell_lasi)*sys.n_funcs_per_cell+mu)*sys.n_funcs_per_cell+nu)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
            int i=0;
            for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
              double tmp_sigma=0;
              for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
                // tmp_sigma += lip_ints(cell_munu, cell_lasi, mu, nu, lambda, sigma) * D_sigma[cell_lasi][sigma];
                tmp_sigma += lip_ints[indx0+i] * D_sigma[cell_lasi][sigma];
                i++;
              }
              tmp_lambda += tmp_sigma * D_lambda[cell_lasi][lambda];
            }
            tmp_nu += tmp_lambda * D_nu[cell_munu][nu]; // * fac_nu;
          }
          int_4c2e += tmp_nu * D_mu[cell_munu][mu];
        }
      }
    }
  }

  auto t1 = Clock::now();
  // times.append_item("4c2e 1cexc integral (precalc)", t1 - t0); // too fast, too many
  times.ind_down();

  return int_4c2e;
}


// calc all (mu nu|lambda sigma) which are LIP
// the order is c_munu, c_lasi, mu, nu, lambda, sigma
void state1d::calc_and_set_4c2e_lip(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);
  const double a = 0.15;

  double written_cells = 0;

  // grid, cells, lambda, sigma
  tensor_6 munulambdasigma({sys.n_cells, sys.n_cells, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell});
  for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
    // cout << "cell mu nu: " << cell_munu << endl;

    tensor_4 int_lambda_sigma({sys.n_gp_fine_per_cell, sys.n_cells, sys.n_funcs_per_cell, sys.n_funcs_per_cell}); // grid, cells, lambda, sigma
    // for every point in int_lambda_mu (which uses the fine grid)
    for(int k=0; k<sys.n_gp_fine_per_cell; k++){ // grid point index in cell
      // const double r1 = sys.grid_begin + (i_c+0.5)*sys.cell_width + sys.grid_fine.p[i_p];
      const double r1 = cell_munu * sys.cell_width + sys.grid_fine.p[k];

      // for(int cell_lasi=cell_munu; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
      for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma

        for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
          for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
            double tmp=0;
            for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
              // const double r2 = sys.grid_begin + (cell_lasi+0.5)*sys.cell_width + sys.grid_fine.p[k_p];
              const double r2 = cell_lasi * sys.cell_width + sys.grid_fine.p[k_p];
              const double r12 = r1-r2;
              // cout << k << ", " << lambda << ", " << sigma << ", " << r12 << endl;
              // cout << "  " << sys.grid_fine.w[k_p] * exp(-a*r12*r12) * polynomials_fine_grid[lambda][k_p] * polynomials_fine_grid[sigma][k_p] << endl;
              tmp += sys.grid_fine.w[k_p] * exp(-a*r12*r12) *
                     polynomials_fine_grid[lambda][k_p] *
                     polynomials_fine_grid[sigma][k_p];
            }
            int_lambda_sigma(k, cell_lasi, lambda, sigma) = tmp;
          } // sigma
        } // lambda
      } // cell lambda/sigma
    } // k

    for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
      for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
        for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
          for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, an LIP
            for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, an LIP // do only one triangle
              double tmp=0;
              for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
                tmp += sys.grid_fine.w[k_p] * int_lambda_sigma(k_p, cell_lasi, lambda, sigma) *
                       polynomials_fine_grid[mu][k_p] *
                       polynomials_fine_grid[nu][k_p];
              }
              munulambdasigma(cell_munu, cell_lasi, mu, nu, lambda, sigma) = tmp;
              written_cells++;
            } // sigma
          } // lambda
        } // nu
      } // mu
    }
  }

#if 0
  const size_t cell_munu_test = 10, cell_lasi_test = 12;
  const size_t a_test = 0, b_test = 1, c_test = 2, d_test = 3;
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_munu_test, cell_lasi_test, b_test, a_test, c_test, d_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, d_test, c_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_munu_test, cell_lasi_test, b_test, a_test, d_test, c_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_lasi_test, cell_munu_test, c_test, d_test, a_test, b_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_lasi_test, cell_munu_test, d_test, c_test, a_test, b_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_lasi_test, cell_munu_test, c_test, d_test, b_test, a_test)) < 1.e-10);
  assert( abs(munulambdasigma(cell_munu_test, cell_lasi_test, a_test, b_test, c_test, d_test) - munulambdasigma(cell_lasi_test, cell_munu_test, d_test, c_test, b_test, a_test)) < 1.e-10);
#endif

  lip_4c2e = munulambdasigma;
  has_lip_4c2e = true;

  cout << "cells written in (mu nu|lambda sigma): " << written_cells << " of " << pow(sys.n_funcs_per_cell, 4)*pow(sys.n_cells, 2) << endl;

  auto t1 = Clock::now();
  times.append_item("all (mu nu|lambda sigma) integrals (local, ortho, grid)", t1 - t0);
  times.ind_down();
}


