
#include <cassert>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <vector>

#include "timing.hh"
#include "lina.hh"
#include "aux_math.hh"
#include "integrate.hh"
#include "tensor.hh"
#include "mo_1d.hh"
#include "system_1d.hh"
#include "state_1d.hh"

using namespace std;


tensor_square_sym get_kin_lip(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  tensor_square_sym T_single_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
  for(int f1=0; f1<sys.n_funcs_per_cell; f1++){
    for(int f2=f1; f2<sys.n_funcs_per_cell; f2++){
      const int cell = 0; // does not matter
      // const double integral = integrate_funcs(sys, cell, f1, f2, NABLA2_REAL); // NABLA2_REAL because the polynomials are not zero at the cell boundaries
      const double integral = integrate_funcs(sys, cell, f1, f2, NABLA2_APPR); // fuck it?
      T_single_cell(f1,f2) = -0.5 * integral;
      T_single_cell(f2,f1) = -0.5 * integral;
    }
  }
  // cout << "T single cell: " << T_single_cell << endl;
  tensor_square_sym T({sys.n_gridpoints, sys.n_gridpoints}); // the kinetic energy contribution to H
  for(int c=0; c<sys.n_cells; c++){
    for(int f1=0; f1<sys.n_funcs_per_cell; f1++){
      for(int f2=0; f2<sys.n_funcs_per_cell; f2++){
        const int gi = c*(sys.n_funcs_per_cell - sys.overlapping) + f1;
        const int gj = c*(sys.n_funcs_per_cell - sys.overlapping) + f2;
        T(gi,gj) += T_single_cell(f1,f2);
      }
    }
  }
  // cout << "T= " << T << endl;
  auto t1 = Clock::now();
  times.append_item("kinetic (lip)", t1 - t0);
  times.ind_down();
  return T;
}



tensor_square_sym get_pot_lip(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  // interaction with the external potential
  tensor_square_sym V({sys.n_gridpoints, sys.n_gridpoints}); // the potential energy contribution to H
  for(int c=0; c<sys.n_cells; c++){
    for(int f1=0; f1<sys.n_funcs_per_cell; f1++){
      for(int f2=0; f2<sys.n_funcs_per_cell; f2++){
        const int gi = c*(sys.n_funcs_per_cell-sys.overlapping) + f1;
        const int gj = c*(sys.n_funcs_per_cell-sys.overlapping) + f2;
        const double integral =  integrate_funcs(sys, c, f1, f2, POTENTIAL_E);
// cout << "gigj: " << gi << ", " << gj << ", " << integral <<  endl;
        V(gi,gj) += integral;
        // cout << integral << endl;
      }
    }
  }
  // cout << "V= " << V << endl;
  auto t1 = Clock::now();
  times.append_item("potential (lip)", t1 - t0);
  times.ind_down();
  return V;
}



// kinetic energy matrix in local orbitals
tensor_square get_kin_local(const system1d& sys, const state1d& S, const operator_t op, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  tensor_square T_global({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell});
  for(int c=0; c<sys.n_cells; c++){
    tensor_square T_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});

    if(op==NABLA2_APPR){
      for(int i=0; i<sys.n_funcs_per_cell; i++){
        for(int j=i; j<sys.n_funcs_per_cell; j++){
          // const double integral = integrate_funcs(sys, c, i, j, NABLA2_REAL); // NABLA2_REAL because the polynomials are not zero at the cell boundaries
          // const double integral = integrate_funcs(sys, c, i, j, NABLA2_APPR); // fuck it?
          const double integral = integrate_funcs(sys, c, i, j, op);
          T_cell(i,j) = integral;
          T_cell(j,i) = integral;
          // cout << integral << endl;
        }
      }
    }
    else if(op==NABLA2_REAL){
      for(int i=0; i<sys.n_funcs_per_cell; i++){
        for(int j=0; j<sys.n_funcs_per_cell; j++){
          const double integral = integrate_funcs(sys, c, i, j, op);
          T_cell(i,j) = integral;
        }
      }
    }

    const tensor_square& trafo = S.locbas1d_HF_orthonorm.basis[c];
    tensor_square T_cell_local = matmul(trafo, matmul(T_cell, trafo, 'N', 'N'), 'T', 'N');

    for(int i=0; i<sys.n_funcs_per_cell; i++){
      const int gi = c*sys.n_funcs_per_cell + i;
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        const int gj = c*sys.n_funcs_per_cell + j;
        T_global(gi,gj) = -0.5 * T_cell_local(i,j);
      }
    }
  }
  auto t1 = Clock::now();
  times.append_item("kinetic (local)", t1 - t0);
  times.ind_down();
  return T_global;
}



// potential energy matrix in local orbitals
tensor_square_sym get_pot_local(const system1d& sys, const state1d &S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  tensor_square_sym V_global({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell});
  for(int c=0; c<sys.n_cells; c++){
    tensor_square_sym V_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        const double integral =  integrate_funcs(sys, c, i, j, POTENTIAL_E);
// cout << "gigj: " << gi << ", " << gj << ", " << integral <<  endl;
        V_cell(i,j) += integral;
        // cout << integral << endl;
      }
    }

    const tensor_square& trafo = S.locbas1d_HF_orthonorm.basis[c];
    tensor_square_sym V_cell_local = matmul(trafo, matmul(V_cell, trafo, 'N', 'N'), 'T', 'N');

    for(int i=0; i<sys.n_funcs_per_cell; i++){
      const int gi = c*sys.n_funcs_per_cell + i;
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        const int gj = c*sys.n_funcs_per_cell + j;
        V_global(gi,gj) = V_cell_local(i,j);
      }
    }
  }
  auto t1 = Clock::now();
  times.append_item("potential (local)", t1 - t0);
  times.ind_down();
  return V_global;
}


vector<vector<double>> tabulate_polynomial_coeffs(const vector<polynomial>& polys, const grid1d& grid){
  vector<vector<double>> polynomials_fine_grid(polys.size());
  for(int p=0; p<polys.size(); p++){
    polynomials_fine_grid[p] = vector<double>(grid.size());
    for(int k=0; k<grid.size(); k++) polynomials_fine_grid[p][k] = polys[p](grid.p[k]);
  }
  return polynomials_fine_grid;
}


// (mu nu | lambda sigma)
// coulomb, using LIP, assuming orthogonal functions that can be integrated on a grid
// applicable to Gauss-Lobatto (with the approximation that exp(-ar^2) is locally linear)
// not applicable to equidistant
tensor_square_sym get_coulomb_lip(const system1d& sys, const tensor_square_sym& dens, const double interaction_scaling, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);

  tensor_square_sym J_mat({sys.n_gridpoints, sys.n_gridpoints}); // coulomb contribution to H
  const double a = 0.15;

  for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell index
    for(int lambda_p=0; lambda_p<sys.n_funcs_per_cell; lambda_p++){ // function index in cell
      const int lambda = cell_lasi * (sys.n_funcs_per_cell - sys.overlapping) + lambda_p; // global index

      for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){ // function index in cell
        const int sigma = cell_lasi * (sys.n_funcs_per_cell - sys.overlapping) + sigma_p; // global index

        tensor_1 int_lambda_sigma({sys.n_gridpoints_fine});
        // for every point in int_lambda_sigma (which uses the fine grid)
        for(int i_c=0; i_c<sys.n_cells; i_c++){ // cell index
        for(int i_p=0; i_p<sys.n_gp_fine_per_cell; i_p++){ // grid point index in cell
          const int i = i_c * (sys.n_gp_fine_per_cell - sys.overlapping) + i_p; // global index
          const double r1 = i_c*sys.cell_width + sys.grid_fine.p[i_p];

          // integrate on the fine grid, ie, the product of lambda and sigma is non-zero at every grid point
          // calc tensor_1: \int e^-a r_12^2 * sigma * lambda d2
          double int_lambda_sigma_tmp = 0;
          for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){
            const double r2 = cell_lasi*sys.cell_width + sys.grid_fine.p[k_p];
            const double r12 = r1-r2;
            int_lambda_sigma_tmp += sys.grid_fine.w[k_p] * exp(-a*r12*r12) * interaction_scaling *
                                    polynomials_fine_grid[sigma_p][k_p] *
                                    polynomials_fine_grid[lambda_p][k_p];
            // cout <<  sys.grid.w[lambda_p] << ", " << sys.grid.w[sigma_p] << ", " << exp(-a*r12*r12) << ", " << r12 << endl;
          }
          int_lambda_sigma[i] = int_lambda_sigma_tmp;
        } }
        // cout << int_lambda_sigma << endl;

        for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){
          for(int mu_p=0; mu_p<sys.n_funcs_per_cell; mu_p++){
            const int mu = cell_munu * (sys.n_funcs_per_cell - sys.overlapping) + mu_p;

            for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){
              const int nu = cell_munu * (sys.n_funcs_per_cell - sys.overlapping) + nu_p;

              double int_two = 0;
              // integrate on a grid two basis functions and int_lambda_sigma
              for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){
                const int k = cell_munu * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; //  global index
                int_two += sys.grid_fine.w[k_p] * int_lambda_sigma[k] *
                           polynomials_fine_grid[mu_p][k_p] *
                           polynomials_fine_grid[nu_p][k_p];
              }
              J_mat(lambda, sigma) += dens(mu, nu) * int_two;
            }
          }
        }
      }
    }
  }
  // cout << "JmatLIP= " << J_mat << endl;

  auto t1 = Clock::now();
  times.append_item("coulomb (lip)", t1 - t0);
  times.ind_down();
  return J_mat;
}


// (mu nu | lambda sigma)
// coulomb, using LIP, assuming orthogonal functions that can be integrated on a grid
// applicable to Gauss-Lobatto (with the approximation that exp(-ar^2) is locally linear)
// not applicable to equidistant
tensor_square_sym get_coulomb_lip_precalc(const system1d& sys, const state1d& S, const double interaction_scaling, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_lip_4c2e);
  assert(S.has_density);
  const tensor_6& lip_4c2e = S.lip_4c2e;
  const tensor_square_sym& dens = S.density;

  tensor_square_sym J_mat({sys.n_gridpoints, sys.n_gridpoints}); // coulomb contribution to H

  for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell index
    for(int lambda_p=0; lambda_p<sys.n_funcs_per_cell; lambda_p++){ // function index in cell
      const int lambda = cell_lasi * (sys.n_funcs_per_cell - sys.overlapping) + lambda_p; // global index

      for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){ // function index in cell
        const int sigma = cell_lasi * (sys.n_funcs_per_cell - sys.overlapping) + sigma_p; // global index

        for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){
          const size_t ind0 = ((((cell_lasi*sys.n_cells+cell_munu)*sys.n_funcs_per_cell+lambda_p)*sys.n_funcs_per_cell+sigma_p)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
          size_t i = 0;
          for(int mu_p=0; mu_p<sys.n_funcs_per_cell; mu_p++){
            const int mu = cell_munu * (sys.n_funcs_per_cell - sys.overlapping) + mu_p;

            for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){
              const int nu = cell_munu * (sys.n_funcs_per_cell - sys.overlapping) + nu_p;

              // J_mat(lambda, sigma) += dens(mu, nu) * lip_4c2e(cell_munu, cell_lasi, lambda_p, sigma_p, mu_p, nu_p);
              J_mat(lambda, sigma) += dens(mu, nu) * lip_4c2e[ind0+i];
              i++;
            }
          }
        }
      }
    }
  }
  // cout << "JmatLIP= " << J_mat << endl;

  auto t1 = Clock::now();
  times.append_item("coulomb (lip, precalc)", t1 - t0);
  times.ind_down();
  return J_mat;
}


// coulomb, using local basis functions, that can be integrated on a grid (correct)
tensor_square_sym get_coulomb_local(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const tensor_square_sym& loc_dens = S.local_density;
  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);

  tensor_square_sym J_mat({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell}); // coulomb contribution to H
  const double a = 0.15;

  for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
    for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, a local function
      const size_t ind_lambda = cell_lasi*sys.n_funcs_per_cell+lambda;

      for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, a local function
        const size_t ind_sigma = cell_lasi*sys.n_funcs_per_cell+sigma;

        tensor_1 int_lambda_sigma({sys.n_gridpoints_fine});
        // for every point in int_lambda_sigma (which uses the fine grid)
        for(int i_c=0; i_c<sys.n_cells; i_c++){ // cell index
        for(int i_p=0; i_p<sys.n_gp_fine_per_cell; i_p++){ // grid point index in cell
          const int i = i_c * (sys.n_gp_fine_per_cell - sys.overlapping) + i_p; // global index
          const double r1 = i_c*sys.cell_width + sys.grid_fine.p[i_p];

          double int_lambda_sigma_tmp = 0;
          for(int l1_p=0; l1_p<sys.n_funcs_per_cell; l1_p++){ // local point index // expansion of mu/nu in LIP
          for(int l2_p=0; l2_p<sys.n_funcs_per_cell; l2_p++){ // local point index // expansion of mu/nu in LIP

            for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
              const double r2 = cell_lasi*sys.cell_width + sys.grid_fine.p[k_p];
              const double r12 = r1-r2;
              int_lambda_sigma_tmp += sys.grid_fine.w[k_p] * exp(-a*r12*r12) *
                                      polynomials_fine_grid[l1_p][k_p] * basis[cell_lasi](lambda, l1_p) *
                                      polynomials_fine_grid[l2_p][k_p] * basis[cell_lasi](sigma, l2_p);
            }
          } }
          int_lambda_sigma[i] = int_lambda_sigma_tmp;
        } }
        // cout << int_lambda_sigma << endl;

        for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // the cell_munu
#if 1
          tensor_square_sym lip_ints_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
          for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, a LIP
  
            for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, a LIP
  
              double int_two = 0;
              for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){ // integration grid
                const int k = cell_munu * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; // global index
                int_two += sys.grid_fine.w[k_p] * int_lambda_sigma[k] *
                           polynomials_fine_grid[mu][k_p] *
                           polynomials_fine_grid[nu][k_p];
              }
              lip_ints_cell(mu, nu) = int_two;
            }
          }

          tensor_square_sym local_ints_cell = matmul(basis[cell_munu], matmul(lip_ints_cell, basis[cell_munu], 'N', 'N'), 'T', 'N'); // LIP -> local
  
          for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // a local function
            const size_t ind_mu = cell_munu*sys.n_funcs_per_cell+mu;
            for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // a local function
              const size_t ind_nu = cell_munu*sys.n_funcs_per_cell+nu;
              J_mat(ind_lambda, ind_sigma) += loc_dens(ind_mu, ind_nu) * local_ints_cell(mu, nu);
            }
          }
#else
          for(int mu_p=0; mu_p<sys.n_funcs_per_cell; mu_p++){ // mu, a LIP
  
            for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){ // nu, a LIP
  
              double int_two = 0;
              for(int l1_p=0; l1_p<sys.n_funcs_per_cell; l1_p++){ // local point index // expansion of mu/nu in LIP
              for(int l2_p=0; l2_p<sys.n_funcs_per_cell; l2_p++){ // local point index // expansion of mu/nu in LIP
                for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){ // integration grid
                  const int k = cell_munu * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; // global index
                  int_two += sys.grid_fine.w[k_p] * int_lambda_sigma[k] *
                           sys.legendre_polynomials[l1_p](sys.grid_fine.p[k_p]) * basis[cell_munu](mu_p, l1_p) *
                           sys.legendre_polynomials[l2_p](sys.grid_fine.p[k_p]) * basis[cell_munu](nu_p, l2_p);
                }
              } }
              J_mat(cell_lasi*sys.n_funcs_per_cell+lambda_p, cell_lasi*sys.n_funcs_per_cell+sigma_p) += loc_dens(cell_munu*sys.n_funcs_per_cell+mu_p, cell_munu*sys.n_funcs_per_cell+nu_p) * int_two;
            }
          }
#endif
        }
      }
    }
  }

  auto t1 = Clock::now();
  times.append_item("coulomb (local)", t1 - t0);
  times.ind_down();

  return J_mat;
}


// coulomb, using local basis functions, that can be integrated on a grid (correct)
// using the precalculated (ab|cd)
tensor_square_sym get_coulomb_local_precalc(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const tensor_square_sym& loc_dens = S.local_density;
  const tensor_6& lip_4c2e = S.lip_4c2e;

  tensor_square_sym J_mat({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell}); // coulomb contribution to H

  for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
    tensor_square_sym J_mat_lip({sys.n_funcs_per_cell, sys.n_funcs_per_cell}); // coulomb contribution to H, for one cell_lasi
    for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, a LIP
      for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, a LIP
        for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu

          tensor_square_sym lip_ints_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
          const size_t ind0 = ((((cell_lasi*sys.n_cells+cell_munu)*sys.n_funcs_per_cell+lambda)*sys.n_funcs_per_cell+sigma)*sys.n_funcs_per_cell+0)*sys.n_funcs_per_cell+0;
          size_t i = 0;
          for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, a LIP
            for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, a LIP
              lip_ints_cell[i] = lip_4c2e[ind0+i];//  (cell_lasi, cell_munu, lambda, sigma, mu, nu);
              // lip_ints_cell(mu, nu) = lip_4c2e(cell_lasi, cell_munu, lambda, sigma, mu, nu);
              i++;
            }
          }

          tensor_square_sym local_ints_cell = matmul(basis[cell_munu], matmul(lip_ints_cell, basis[cell_munu], 'N', 'N'), 'T', 'N');

          for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // a local function
            const size_t ind_mu = cell_munu*sys.n_funcs_per_cell+mu;
            for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // a local function
              const size_t ind_nu = cell_munu*sys.n_funcs_per_cell+nu;
              J_mat_lip(lambda, sigma) += loc_dens(ind_mu, ind_nu) * local_ints_cell(mu, nu);
            }
          } 
        } // munu
      } // sigma
    } // lambda

    tensor_square_sym J_mat_local = matmul(basis[cell_lasi], matmul(J_mat_lip, basis[cell_lasi], 'N', 'N'), 'T', 'N');
    
    for(int lambda=0; lambda<sys.n_funcs_per_cell; lambda++){ // lambda, a local function
      const size_t ind_lambda = cell_lasi*sys.n_funcs_per_cell+lambda;
      for(int sigma=0; sigma<sys.n_funcs_per_cell; sigma++){ // sigma, a local function
        const size_t ind_sigma = cell_lasi*sys.n_funcs_per_cell+sigma;
        J_mat(ind_lambda, ind_sigma) = J_mat_local(lambda, sigma);
      }
    }
  } // lasi

  auto t1 = Clock::now();
  times.append_item("coulomb (local, precalc)", t1 - t0);
  times.ind_down();

  return J_mat;
}



// (sigma nu | lambda mu)
// exchange, using LIP, assuming orthogonal functions (wrong) that can be integrated on a grid (wrong)
tensor_square_sym get_exchange_lip(const system1d& sys, const tensor_square_sym& dens, const double interaction_scaling, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);

  tensor_square_sym K_mat({sys.n_gridpoints, sys.n_gridpoints}); // exchange contribution to H
  const double a = 0.15;

  for(int lambda_c=0; lambda_c<sys.n_cells; lambda_c++){
  for(int lambda_p=0; lambda_p<sys.n_funcs_per_cell; lambda_p++){
    const int lambda = lambda_c * (sys.n_funcs_per_cell - sys.overlapping) + lambda_p;

    const int mu_c = lambda_c;
    for(int mu_p=0; mu_p<sys.n_funcs_per_cell; mu_p++){
      const int mu = mu_c * (sys.n_funcs_per_cell - sys.overlapping) + mu_p;

      tensor_1 int_lambda_mu({sys.n_gridpoints_fine});
      for(int i_c=0; i_c<sys.n_cells; i_c++){ // cell index
      for(int i_p=0; i_p<sys.n_gp_fine_per_cell; i_p++){ // function index in cell
        const int i = i_c * (sys.n_gp_fine_per_cell - sys.overlapping) + i_p; // global index
        const double r1 = sys.grid_begin + (i_c+0.5)*sys.cell_width + sys.grid_fine.p[i_p];

        // integrate on the grid, ie, the product of lambda and sigma is non-zero at only one grid point
        // calc tensor_1 \int  e^-a r_12^2 * sigma * lambda d2
        double int_lambda_mu_tmp = 0;
        for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){
          const double r2 = sys.grid_begin + (lambda_c+0.5)*sys.cell_width + sys.grid_fine.p[k_p];
          const double r12 = r1-r2;
          int_lambda_mu_tmp += sys.grid_fine.w[k_p] * exp(-a*r12*r12) * interaction_scaling *
                               polynomials_fine_grid[mu_p][k_p] *
                               polynomials_fine_grid[lambda_p][k_p];
          // cout <<  sys.grid.w[lambda_p] << ", " << sys.grid.w[sigma_p] << ", " << exp(-a*r12*r12) << ", " << r12 << endl;
        }
        int_lambda_mu[i] = int_lambda_mu_tmp;
      } }
      // cout << int_lambda_mu << endl;

      for(int sigma_c=0; sigma_c<sys.n_cells; sigma_c++){
      for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){
        const int sigma = sigma_c * (sys.n_funcs_per_cell - sys.overlapping) + sigma_p;

        const int nu_c = sigma_c;
        for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){
          const int nu = nu_c * (sys.n_funcs_per_cell - sys.overlapping) + nu_p;

          double int_two = 0;
          // integrate on a grid two basis functions and int_lambda_sigma
          for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){
            const int k = sigma_c * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; //  global index
            int_two += sys.grid_fine.w[k_p] * int_lambda_mu[k] *
                       polynomials_fine_grid[sigma_p][k_p] *
                       polynomials_fine_grid[nu_p][k_p];
          }

      // cout << int_two << ", " <<  dens(lambda, sigma_point) << endl;
          const double dens_normed = dens(mu, sigma); // the LIP are not normalised, hence D is not normalised (and its trace is not the number of electrons)
          K_mat(lambda, nu) += int_two * dens_normed;
        }
      } }
    }
  } }
  // cout << "KmatLIP= "  << K_mat << endl;

  auto t1 = Clock::now();
  times.append_item("exchange (lip, ortho, grid)", t1 - t0);
  times.ind_down();

  return K_mat;
}


// coulomb, using local basis functions, assuming orthogonal functions (correct) that can be integrated on a grid (correct)
tensor_square_sym get_exchange_local(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  const tensor_square_sym& loc_dens = S.local_density;
  const vector<vector<double>> polynomials_fine_grid = tabulate_polynomial_coeffs(sys.legendre_polynomials, sys.grid_fine);

  tensor_square_sym K_mat({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell}); // coulomb contribution to H
  const double a = 0.15;

  for(int cell_lamu=0; cell_lamu<sys.n_cells; cell_lamu++){ // cell of lambda and mu
  for(int lambda_p=0; lambda_p<sys.n_funcs_per_cell; lambda_p++){ // lambda, a local function

    // const int sigma_c = cell_lamu;
    for(int mu_p=0; mu_p<sys.n_funcs_per_cell; mu_p++){ // sigma, a local function

      tensor_1 int_lambda_mu({sys.n_gridpoints_fine});
      // for every point in int_lambda_mu (which uses the fine grid)
      for(int i_c=0; i_c<sys.n_cells; i_c++){ // cell index
      for(int i_p=0; i_p<sys.n_gp_fine_per_cell; i_p++){ // grid point index in cell
        const int i = i_c * (sys.n_gp_fine_per_cell - sys.overlapping) + i_p; // global index
        const double r1 = sys.grid_begin + (i_c+0.5)*sys.cell_width + sys.grid_fine.p[i_p];

        double int_lambda_mu_tmp = 0;
        for(int l1_p=0; l1_p<sys.n_funcs_per_cell; l1_p++){ // local point index // expansion of mu/nu in LIP
        for(int l2_p=0; l2_p<sys.n_funcs_per_cell; l2_p++){ // local point index // expansion of mu/nu in LIP

          for(int k_p=0; k_p<sys.n_gp_fine_per_cell; k_p++){ // integration
            const double r2 = sys.grid_begin + (cell_lamu+0.5)*sys.cell_width + sys.grid_fine.p[k_p];
            const double r12 = r1-r2;
            int_lambda_mu_tmp += sys.grid_fine.w[k_p] * exp(-a*r12*r12) *
                                 polynomials_fine_grid[l1_p][k_p] * basis[cell_lamu](lambda_p, l1_p) *
                                 polynomials_fine_grid[l2_p][k_p] * basis[cell_lamu](mu_p, l2_p);
          }
        } }
        int_lambda_mu[i] = int_lambda_mu_tmp;
      } }
      // cout << int_lambda_mu << endl;

      for(int cell_sinu=0; cell_sinu<sys.n_cells; cell_sinu++){ // the cell_sinu
#if 1
        tensor_square_sym lip_ints_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
        for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){ // mu, a LIP

          for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){ // nu, a LIP

            double int_two = 0;
            for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){ // integration grid
              const int k = cell_sinu * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; // global index
              int_two += sys.grid_fine.w[k_p] * int_lambda_mu[k] *
                         polynomials_fine_grid[sigma_p][k_p] *
                         polynomials_fine_grid[nu_p][k_p];
            }
            lip_ints_cell(sigma_p, nu_p) = int_two;
          }
        }
        tensor_square_sym local_ints_cell = matmul(basis[cell_sinu], matmul(lip_ints_cell, basis[cell_sinu], 'N', 'N'), 'T', 'N');

        for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){ // a local function
          for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){ // a local function
            K_mat(cell_lamu*sys.n_funcs_per_cell+mu_p, cell_sinu*sys.n_funcs_per_cell+nu_p) += loc_dens(cell_lamu*sys.n_funcs_per_cell+lambda_p, cell_sinu*sys.n_funcs_per_cell+sigma_p) * local_ints_cell(sigma_p, nu_p);
          }
        }
#else
        for(int sigma_p=0; sigma_p<sys.n_funcs_per_cell; sigma_p++){ // mu, a LIP

          for(int nu_p=0; nu_p<sys.n_funcs_per_cell; nu_p++){ // nu, a LIP

            double int_two = 0;
            for(int l1_p=0; l1_p<sys.n_funcs_per_cell; l1_p++){ // local point index // expansion of mu/nu in LIP
            for(int l2_p=0; l2_p<sys.n_funcs_per_cell; l2_p++){ // local point index // expansion of mu/nu in LIP
              for(int k_p=0; k_p<sys.grid_fine.size(); k_p++){ // integration grid
                const int k = cell_sinu * (sys.n_gp_fine_per_cell - sys.overlapping) + k_p; // global index
                int_two += sys.grid_fine.w[k_p] * int_lambda_mu[k] *
                         sys.legendre_polynomials[l1_p](sys.grid_fine.p[k_p]) * basis[cell_sinu](sigma_p, l1_p) *
                         sys.legendre_polynomials[l2_p](sys.grid_fine.p[k_p]) * basis[cell_sinu](nu_p, l2_p);
              }
            } }
            K_mat(cell_lamu*sys.n_funcs_per_cell+mu_p, cell_sinu*sys.n_funcs_per_cell+nu_p) += loc_dens(cell_lamu*sys.n_funcs_per_cell+lambda_p, cell_sinu*sys.n_funcs_per_cell+sigma_p) * int_two;
          }
        }
#endif
      }
    }
  } }
  // cout << "Kmatlocal " << K_mat << endl;

  auto t1 = Clock::now();
  times.append_item("exchange (local, ortho, grid)", t1 - t0);
  times.ind_down();

  return K_mat;
}



// // guess occupations ... lets use sin(a*x) from 0 to pi, where a enumerates the mos
// void state1d::guess_occupied_mos(const system1d& sys){
//   for(int a=0; a<occupied_mos.size(); a++){
//     const size_t occ = occupied_mos[a];
//     for(int c=0; c<sys.n_cells; c++){
//       for(int p=0; p<sys.n_funcs_per_cell; p++){
//         const double x = sys.grid_begin + (c+0.5)*sys.cell_width + sys.grid.p[p];
//         const int gindex = c*(sys.n_funcs_per_cell - sys.overlapping) + p;
//         mo_coeffs(occ,gindex) = sin(((occ+1)*(x-sys.grid_begin))*M_PI/(-2*sys.grid_begin));
//       }
//     }
//   // cout << mo_coeffs.get_slice(a, XX) << endl;
//   }
// }


void state1d::set_mo_coeffs(const tensor_square& coeffs){
  mo_coeffs = coeffs;
}


void state1d::set_mo_energies(const tensor_1& en){
  mo_energies = en;
}


// update density from mo_coeffs and occupied_mo's
// assuming two electrons per orbital: P = 2*sum C*C
// assumes MO to be normalised
void state1d::update_density(const system1d& sys){
  tensor_square dens_new({n_funcs, n_funcs});
  // as the basis functions are not normalised, their norm needs to be taken into account
  // get norms
  vector<double> norms(sys.n_funcs_per_cell);
  for(int i=0; i<sys.n_funcs_per_cell; i++){
    // norms[i] = sqrt(sys.grid.w[i]);
    norms[i] = 1.0;
  }
  // get the actual density
  double trace=0;
  for(int a=0; a<occupied_mos.size(); a++){
    const size_t occ = occupied_mos[a];
    for(int c1=0; c1<sys.n_cells; c1++){
      for(int f1=0; f1<sys.n_funcs_per_cell; f1++){
        const size_t i = c1*(sys.n_funcs_per_cell-sys.overlapping) + f1;
        const double norm_i = norms[f1];
        trace += 2 * mo_coeffs(occ,i) * mo_coeffs(occ,i) * norm_i * norm_i;
        for(int c2=0; c2<sys.n_cells; c2++){
          for(int f2=0; f2<sys.n_funcs_per_cell; f2++){
            const size_t j = c2*(sys.n_funcs_per_cell-sys.overlapping) + f2;
            const double norm_j = norms[f2];
            // dens_new(i,j) += 2 * mo_coeffs(occ,i) * norm_i * mo_coeffs(occ,j) * norm_j;
            dens_new(i,j) = 2 * mo_coeffs(occ,i) * norm_i * mo_coeffs(occ,j) * norm_j;
          }
        }
      }
    }
  }

  // cout << "trace: " << trace << " (should be: " << 2*occupied_mos.size() << ")" << endl;
  // const double eps = 1.e-6;
  cout << "trace, LIP: " << trace << endl;
  //assert(abs(trace - 2*occupied_mos.size()) < eps); // trace should be the number of electrons
  density = dens_new;
  has_density = true;
  // cout << density << endl;
}


// normalise the nth orbital
void state1d::normalise_mo_lip(const system1d& sys, const size_t n){
  tensor_1 mo_work = mo_coeffs.get_slice(n, XX);
  const double norm = integrate_system(sys, mo_work, mo_work, NONE);
  const double norm_sqrt = sqrt(norm);
  for(double &d: mo_work) d /= norm_sqrt;
  // cout << n << ": " << norm << endl;
  mo_energies[n] /= norm;
  mo_coeffs.set_slice(n, XX, mo_work);
  // mo_work = mo_coeffs.get_slice(n, XX);
  // cout << integrate_system(sys, mo_work, mo_work, NONE) << endl;
}


void state1d::normalise_occupied_lip(const system1d& sys){
  for(const size_t mo: occupied_mos){
    normalise_mo_lip(sys, mo);
  }
}


void state1d::normalise_all_lip(const system1d& sys){
  for(size_t mo=0; mo<n_funcs; mo++){
    normalise_mo_lip(sys, mo);
  }
}


void state1d::do_scf(const system1d& sys, const tensor_square_sym& T, const tensor_square_sym& V, tensor_square_sym& J_mat, tensor_square_sym& K_mat, double e0, const double thr, const double interaction_scaling, const bool omit_coulomb, timing& times){

  times.ind_up();
  auto t0 = Clock::now();

  assert(has_lip_4c2e);
  assert(sys.globbas1d_plain.has_S_invsqrt_sqrt);
  const tensor_square& S_invsqrt = sys.globbas1d_plain.S_invsqrt;

  // tensor_square_sym J_mat({sys.n_gridpoints, sys.n_gridpoints}); // coumlomb contribution to H
  // tensor_square_sym K_mat({sys.n_gridpoints, sys.n_gridpoints}); // exchange contribution to H
  tensor_square_sym H({sys.n_gridpoints, sys.n_gridpoints});
  tensor_square_sym H_prime({sys.n_gridpoints, sys.n_gridpoints});
  double old_e;

  int cycle = 0;
  do{
    old_e = e0;
    normalise_occupied_lip(sys);
    update_density(sys);
    // const tensor_square_sym& density = S.density_ro();
    // J_mat = get_coulomb_lip(sys, density, lip_4c2e, interaction_scaling, times);
    J_mat = get_coulomb_lip_precalc(sys, *this, interaction_scaling, times);
    //if(cycle==0) cout << "J = " << J_mat << endl;
    // K_mat = get_exchange_lip(sys, density, lip_4c2e, interaction_scaling, times);
    //if(cycle==0) cout << "K = " << K_mat << endl;

    cout << "getting H ..." << endl;
    // H_prime = T + V + J_mat - K_mat * 0.5;
    if(omit_coulomb) 
      H_prime = T + V;
    else
      H_prime = T + V + J_mat * 0.5;
    H = matmul(S_invsqrt, matmul(H_prime, S_invsqrt, 'N', 'N'), 'N', 'N');
    // cout << "hp= " << H_prime << endl;
    // cout << "h= " << H << endl;

    cout << "diagonalising H ..." << endl;
    tensor_1 evals;
    tensor_square evecs;
    tie(evals, evecs) = H.eigensystem_sym(times);
    evecs = matmul(S_invsqrt, evecs, 'N', 'N');
    // cout << "eigenvalue 0, 1, 2: " << evals[0] << ", " << evals[1] << ", " << evals[2] << endl;
    // cout << "eigenvectors: " << evecs << endl;
    set_mo_coeffs(evecs);
    set_mo_energies(evals);

    e0=0;
    for(int i=0; i<occupied_mos.size(); i++){
      e0 += evals[occupied_mos[i]];
    }
    cout << "cycle " << cycle++ << " E=" << setprecision(13) << e0 << " dE=" << e0 - old_e <<endl;
//    // cout << S.density << endl;
  } while(fabs(e0 - old_e) > thr);

  update_density(sys);

  auto t1 = Clock::now();
  times.append_item("do scf", t1 - t0);
  times.ind_down();
}



tensor_square state1d::get_H_1el_local(const system1d& sys, const operator_t op, timing& times) const {
  // cout << "getting 1 electron Hamiltonian ..." << endl;
  const tensor_square T_local = get_kin_local(sys, *this, op, times);
  const tensor_square_sym V_local = get_pot_local(sys, *this, times);
  tensor_square H_1el_local;
  H_1el_local = T_local + V_local;
  // cout << "done" << endl;
  return H_1el_local;
}


// void state1d::set_local_basis(const vector<tensor_square>& locbas){
//   assert(locbas.size() == locbas1d_HF_orthonorm.n_cells);
//   locbas1d_HF_orthonorm.basis = locbas;
// }


void state1d::calc_and_set_local_basis(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  // for the moment:
  assert(occupied_mos.size()==1 && occupied_mos[0]==0);
  const tensor_1 mo_0 = mo_coeffs_ro().get_slice(0, XX);

  cout << locbas1d_HF_orthonorm.basis[15] << endl;
  for(int c=0; c<sys.n_cells; c++){
    const size_t index_begin = c*(sys.n_funcs_per_cell-sys.overlapping);
    const size_t index_end = c*(sys.n_funcs_per_cell-sys.overlapping) + sys.n_funcs_per_cell;
    // cout << g_begin << ", " << g_end << endl;
    // cout << "cell " << c << endl;
    const tensor_1 cell_basis({sys.n_funcs_per_cell}, mo_0.begin()+index_begin, mo_0.begin()+index_end);
    // cout << cell_basis << endl;
    const bool debug_print = true;
    // locbas1d_HF_orthonorm.basis[c] = orthonormalise_cell(sys, cell_basis, c, debug_print);
    locbas1d_HF_orthonorm.basis[c] = orthonormalise_cell_wo_boundary(sys, cell_basis, c, debug_print);
  }
  cout << locbas1d_HF_orthonorm.basis[15] << endl;
  // cout << new_basis << endl;

  auto t1 = Clock::now();
  times.append_item("calc and set local basis", t1 - t0);
  times.ind_down();
}


void state1d::calc_and_set_local_mo_coeffs(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  for(int c1=0; c1<sys.n_cells; c1++){
    tensor_2 mo_coeffs_cell({sys.n_gridpoints, sys.n_funcs_per_cell});
    for(int f1=0; f1<sys.n_gridpoints; f1++){
      for(int f2=0; f2<sys.n_funcs_per_cell; f2++){
        const int ind1 = f1;
        const int ind2 = c1*(sys.n_funcs_per_cell - sys.overlapping) + f2;
        mo_coeffs_cell(f1,f2) = mo_coeffs(ind1, ind2);
      }
    }
    tensor_2 local_mo_coeffs_cell({sys.n_gridpoints, sys.n_funcs_per_cell});

    char NTchar1 = 'N', NTchar2 = 'N';
    int N1 = locbas1d_HF_orthonorm.basis_inv[c1].get_dims()[0];
    int N2 = locbas1d_HF_orthonorm.basis_inv[c1].get_dims()[1];
    int N3 = mo_coeffs_cell.get_dims()[0];
    int N4 = mo_coeffs_cell.get_dims()[1];
    double alpha = 1.0, beta = 0.0;

    dgemm_(&NTchar1, &NTchar2, &N2, &N3, &N1,
           &alpha, locbas1d_HF_orthonorm.basis_inv[c1].begin(), &N2,
           mo_coeffs_cell.begin(), &N4,
           &beta, local_mo_coeffs_cell.begin(), &N2);

    // local_mo_coeffs_cell = matmul(locbas1d_HF_orthonorm.basis_inv[c1], mo_coeffs_cell, 'N', 'N');
    for(int f1=0; f1<sys.n_gridpoints; f1++){
      for(int f2=0; f2<sys.n_funcs_per_cell; f2++){
        const int ind1 = f1;
        const int ind2 = c1*sys.n_funcs_per_cell + f2;
        local_mo_coeffs(ind1,ind2) = local_mo_coeffs_cell(f1,f2);
      }
    }
  }
  has_local_mo_coeffs = true;

  auto t1 = Clock::now();
  times.append_item("calc and set local mo coeffs", t1 - t0);
  times.ind_down();
}


// assuming doubly occupied orbitals
// assuming a normalised local basis
void state1d::calc_and_set_local_density_from_local_mo_coeffs(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  for(int a=0; a<occupied_mos.size(); a++){
    const size_t occ = occupied_mos[a];
    for(int i=0; i<sys.n_cells*sys.n_funcs_per_cell; i++){
      const double norm_i = 1.0;
      for(int j=0; j<sys.n_cells*sys.n_funcs_per_cell; j++){
        const double norm_j = 1.0;
        local_density(i,j) = 2 * local_mo_coeffs(occ,i) * norm_i * local_mo_coeffs(occ,j) * norm_j;
      }
    }
  }

  double trace = 0;
  for(int i=0; i< local_density.get_dims()[0]; i++) trace += local_density(i,i);
  cout << "trace, local: " << trace << endl;
  const double eps = 1.e-6;
  assert(abs(trace - 2*occupied_mos.size()) < eps); // trace should be the number of electrons

  has_local_density = true;

  auto t1 = Clock::now();
  times.append_item("update local density from local mo coeffs", t1 - t0);
  times.ind_down();
}


void state1d::calc_and_set_local_density_from_density(const system1d& sys, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  local_density = tensor_square({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell});

  for(int c1=0; c1<sys.n_cells; c1++){
    for(int c2=0; c2<sys.n_cells; c2++){
      // make copy of one cell of the density
      tensor_square dens_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
      for(int i=0; i<sys.n_funcs_per_cell; i++){
        const size_t ind1 = c1 * (sys.n_funcs_per_cell - sys.overlapping) + i;
        for(int j=0; j<sys.n_funcs_per_cell; j++){
          const size_t ind2 = c2 * (sys.n_funcs_per_cell - sys.overlapping) + j;
          dens_cell(i,j) = density(ind1, ind2);
        }
      }

      const tensor_square& trafo_inv1 = locbas1d_HF_orthonorm.basis_inv[c1];
      const tensor_square& trafo_inv2 = locbas1d_HF_orthonorm.basis_inv[c2];
      tensor_square local_dens_cell = matmul(trafo_inv2, matmul(dens_cell, trafo_inv1, 'N', 'T'), 'N', 'N');

// if(c1==30 && c2==30){
// cout << "dc1 = " << dens_cell << endl;
// cout << "ti1 = " << trafo_inv1 << endl;
// // cout << "t1 = " << trafo_1 << endl;
// cout << "ti2 = " << trafo_inv2 << endl;
// cout << "ldc = " << local_dens_cell << endl;
// }

      for(int i=0; i<sys.n_funcs_per_cell; i++){
        const size_t ind1 = c1 * sys.n_funcs_per_cell + i;
        for(int j=0; j<sys.n_funcs_per_cell; j++){
          const size_t ind2 = c2 * sys.n_funcs_per_cell + j;
          local_density(ind1, ind2) = local_dens_cell(i,j);
        }
      }
    }
  }
  // cout << "local dens: " << local_density << endl;

  double trace = 0;
  for(int i=0; i< local_density.get_dims()[0]; i++) trace += local_density(i,i);
  cout << "loc trace: " << trace << endl;
  const double eps = 1.e-6;
  assert(abs(trace - 2*occupied_mos.size()) < eps); // trace should be the number of electrons

  has_local_density = true;

  auto t1 = Clock::now();
  times.append_item("update local density from density", t1 - t0);
  times.ind_down();
}


string state1d::mo_to_mathematica(const system1d& sys, const size_t mo) const {
  stringstream ss;
  ss << "Plot[Piecewise[{";
  for(int c=0; c<sys.n_cells; c++){
    ss << "{";
    for(int f=0; f<sys.n_funcs_per_cell; f++){
      const size_t gi = c*(sys.n_funcs_per_cell - sys.overlapping) + f;
      const double coeff = mo_coeffs(mo,gi);
      const double shift = sys.grid_begin+sys.cell_width*(c+0.5);
      ss << sys.legendre_polynomials[f].shift(shift) * coeff << endl;
      if(f<sys.n_funcs_per_cell-1) ss << "+";
    }
    ss << "," << sys.grid_begin+sys.cell_width*c << "<x<" << sys.grid_begin+sys.cell_width*(c+1);
    ss << "}";
    if(c<sys.n_cells-1) ss << ",";
  }
  ss << "}],{x," << sys.grid_begin << "," << sys.grid_begin+sys.cell_width*sys.n_cells << "}]";
  return ss.str();
}


