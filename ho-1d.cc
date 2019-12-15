
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <utility>

#include "aux_math.hh"
#include "auxiliary.hh"
#include "grid_1d.hh"
#include "integrate.hh"
#include "polynomial.hh"
#include "mo_1d.hh"
#include "state_1d.hh"
#include "system_1d.hh"
#include "timing.hh"


using namespace std;


int main(int ac, char **av) {

  timing times;
  auto t0 = Clock::now();

  if (ac != 6) {
    cout << "n_cells, cell_width, system_width, n_funcs_per_cell, interaction_scaling" << endl;
    cout << "try ./ho-1d 256 0.125 32 7 1.0" << endl;
    cout << "try ./ho-1d 64 0.25 16 5 1.0" << endl;
    cout << "try ./ho-1d 32 0.5 16 5 1.0" << endl;
    abort();
  }

  system1d sys;
  sys.n_cells = stol(av[1], 0, 0);
  sys.cell_width = stod(av[2], 0);
  const double system_width = stod(av[3], 0);
  sys.n_funcs_per_cell = stol(av[4], 0, 0);
  const double interaction_scaling  = stod(av[5], 0);

  // sys.n_cells = 400;
  // sys.cell_width = 0.1;
  sys.grid_begin = -system_width/2.0;
  // sys.n_funcs_per_cell = 10;
  sys.n_gp_per_cell = sys.n_funcs_per_cell;
  sys.n_gp_fine_per_cell = sys.n_funcs_per_cell + 1;
  const vector<double> bounds({-sys.cell_width/2.0, sys.cell_width/2.0});
#if 0 // not useful here, because one can neither diagonalise nor run Helmholtz
  sys.grid = grid1d(bounds, vector<size_t>({sys.n_funcs_per_cell}), GAUSSIAN_LIN);
  sys.overlapping = false;
#else
  // sys.grid = grid1d(bounds, vector<size_t>({sys.n_funcs_per_cell}), EQUIDISTANT);
  sys.grid = grid1d(bounds, vector<size_t>({sys.n_gp_per_cell}), GAUSS_LOBATTO);
  sys.grid_fine = grid1d(bounds, vector<size_t>({sys.n_gp_fine_per_cell}), GAUSS_LOBATTO);
  sys.overlapping = true;
#endif
  sys.n_gridpoints = sys.n_cells * (sys.grid.size() - sys.overlapping) + sys.overlapping;
  sys.n_gridpoints_fine = sys.n_cells * (sys.grid_fine.size() - sys.overlapping) + sys.overlapping;
  // sys.tgrid = grid1d(vector<double>({0.0, 2.0, 500.0}), vector<size_t>({20,16}), TGRID);
  sys.pot_coeffs = tensor_1({sys.n_gridpoints});
  for(int c=0; c<sys.n_cells; c++){
    for(int p=0; p<sys.n_funcs_per_cell; p++){
      // lets use V = x^2
      const double x = sys.grid_begin + (c+0.5)*sys.cell_width + sys.grid.p[p];
      const int gindex = c*(sys.n_funcs_per_cell - sys.overlapping) + p;
      sys.pot_coeffs[gindex] = 0.5*x*x - 100.0;
      // sys.pot_coeffs[gindex] = 225.0/2401.0*x*x*x*x - 150.0/49.0*x*x - 75.0;
      // sys.pot_coeffs[gindex] = 0.0173*x*x*x*x - 0.832*x*x - 90.0;
      // sys.pot_coeffs[gindex] = 0.0150*x*x*x*x - 0.547*x*x - 95.0;
      // sys.pot_coeffs[gindex] = 0.0138*x*x*x*x - 0.406*x*x - 97.0;
      // sys.pot_coeffs[gindex] = 0.0130*x*x*x*x - 0.323*x*x - 98.0;
      // sys.pot_coeffs[gindex] = 0.0121*x*x*x*x - 0.220*x*x - 99.0;
      // sys.pot_coeffs[gindex] = 0.0121*x*x*x*x - 0.220*x*x + 0.2*x - 99.0;
      // sys.pot_coeffs[gindex] = 0.01*x*x*x*x - 100.0;
    }
  }
  sys.calc_and_set_polynomials();
  sys.calc_and_set_derivatives();
  // sys.calc_and_set_normalised_polynomials();
  // sys.calc_and_set_normalised_derivatives();

  sys.calc_and_set_plain_basis(times);
  // sys.calc_and_set_normalised_basis();

  cout << "system: " << sys << endl;
  auto t1 = Clock::now();
  times.ind_up();
  times.append_item("setup", t1 - t0);
  times.ind_down();

  vector<size_t> occ(vector<size_t>({0})); // first alpha and first beta are occupied orbital
  state1d S(sys, occ);

// cout << "S = " << S.locbas1d.overlap << endl;

// tensor_square_sym S_invsqrt({sys.n_gridpoints, sys.n_gridpoints}), S_sqrt({sys.n_gridpoints, sys.n_gridpoints});
// tie(S_invsqrt, S_sqrt) = S.locbas1d.get_overlap_invsqrt_sqrt();

  cout << "getting T ..." << endl;
  tensor_square_sym T = get_kin_lip(sys, times); // the kinetic energy contribution to H
  // cout << "T= " << T << endl;

  cout << "getting V ..." << endl;
  tensor_square_sym V = get_pot_lip(sys, times); // the potential energy contribution to H
  // cout << "V= " << V << endl;

  cout << "getting H ..." << endl;
  tensor_square_sym H1_prime(T + V);
  // cout << "H'= " << H1_prime << endl;

  tensor_square_sym S_invsqrt(sys.globbas1d_plain.S_invsqrt);
  // cout << S_invsqrt << endl;
  tensor_square_sym H1 = matmul(S_invsqrt, matmul(H1_prime, S_invsqrt, 'N', 'N'), 'N', 'N');
  // cout << H1 << endl;

  cout << "diagonalising 1e H ..." << endl;
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = H1.eigensystem_sym(times);
  // get evecs of H1_prime
  evecs = matmul(S_invsqrt, evecs, 'N', 'N');
  // evecs = matmul(evecs, S_invsqrt, 'N', 'N');
  cout << "eigenvalues 1e: " << evals << endl;
  // cout << "eigenvectors: " << evecs << endl;


  double e0=0;
  for(int i=0; i<occ.size(); i++){
    e0 += evals[occ[i]];
  }

  S.set_mo_coeffs(evecs);
  S.set_mo_energies(evals);
  S.calc_and_set_4c2e_lip(sys, S, times);

  tensor_square_sym J_mat({sys.n_gridpoints, sys.n_gridpoints}); // coumlomb contribution to H
  tensor_square_sym K_mat({sys.n_gridpoints, sys.n_gridpoints}); // exchange contribution to H

  const double thr=1.e-8;
  const bool omit_coulomb = false;
  S.do_scf(sys, T, V, J_mat, K_mat, e0, thr, interaction_scaling, omit_coulomb, times);

  // cout << "dens= " << S.density << endl;
  // cout << "J=" << J_mat << endl;
  // cout << "K=" << K_mat << endl;

  for(int i=0; i<12; i++){
    S.normalise_mo_lip(sys, i);
  }

  for(int i=0; i<10; i++){ // for the first ten eigenvectors
    cout << setprecision(10) <<  S.mo_energies_ro()[i] << " -- " << S.mo_energies_ro()[i+1] - S.mo_energies_ro()[i] << endl;
  }

  cout << "store local basis ..." << endl;
  S.calc_and_set_local_basis(sys, times);
  // cout << "locbas= " << S.locbas1d_HF_orthonorm.basis << endl;
  cout << "done" << endl;

  cout << "store local basis inv ..." << endl;
  S.locbas1d_HF_orthonorm.calc_and_set_basis_inv();
  // cout << "locbasinv= " << S.locbas1d_HF_orthonorm.basis_inv << endl;
  cout << "done" << endl;

  cout << "store local mo_coeffs ..." << endl;
  S.calc_and_set_local_mo_coeffs(sys, times);
  // cout << "local_mo_coeffs= " << S.local_mo_coeffs << endl;
  cout << "done" << endl;

  cout << "store local density ..." << endl;
#if 0
  S.calc_and_set_local_density_from_density(sys, times);
  cout << "locdensfromd= " << S.local_density << endl;
#else
  S.calc_and_set_local_density_from_local_mo_coeffs(sys, times);
  // cout << "locdensfrommo= " << S.local_density << endl;
#endif
  cout << "done" << endl;

  // cout << "getting local T ..." << endl;
  // const tensor_square_sym T_local = get_kin_local(sys, S, times);
  // cout << "lT=" << T_local << endl;
  // cout << "done" << endl;

  // cout << "getting local V ..." << endl;
  // const tensor_square_sym V_local = get_pot_local(sys, S, times);
  // cout << "lV=" << V_local << endl;
  // cout << "done" << endl;

  // cout << "getting local J ..." << endl;
  // const tensor_square_sym coulomb_local = get_coulomb_local_precalc(sys, S, times);
  // cout << "lJ=" << coulomb_local << endl;
  // cout << "done" << endl;

  // cout << "getting local K ..." << endl;
  // const tensor_square_sym exchange_local = get_exchange_local(sys, S, times);
  // cout << "lK=" << exchange_local << endl;
  // cout << "done" << endl;

  // cout << "getting local Hprime ..." << endl;
  // tensor_square_sym H_prime_local;
  // H_prime_local = T_local + V_local + coulomb_local * 0.5;
  // cout << "lhp=" << H_prime_local << endl;
  // cout << "done" << endl;


  const double sci_energy = S.get_xci_energy(sys, "s", times);
  const double dci_energy = S.get_xci_energy(sys, "d", times);
  const double sdci_energy = S.get_xci_energy(sys, "sd", times);

//  const double ci_energy_apprx3 = S.get_ci_energy_apprx3(sys, "s", times);

  cout << "SCI (aka HF): " << sci_energy << endl;
  cout << "DCI: " << dci_energy << endl;
  cout << "SDCI: " << sdci_energy << endl;
//  cout << "appr-SDCI: " << ci_energy_apprx3 << endl;
//  cout << "appr-SCI: " << ci_energy_apprx2 << endl;

  auto t2 = Clock::now();
  times.append_item("total", t2 - t0);

  return 0;
}

