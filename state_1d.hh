#ifndef STATE_1D_HH
#define STATE_1D_HH

#include <vector>
#include <utility>

#include "timing.hh"
#include "local_basis_1d.hh"
#include "tensor.hh"
#include "system_1d.hh"

using namespace std;


class state1d{
public:
  // number of functions per MO in the global basis
  size_t n_funcs;
  // a local orthonormalised basis, based on the plain global basis of non-orthgonal + non-normalised basis functions
  localbasis1d locbas1d_HF_orthonorm;
  // number of mo * n_funcs, ie, n_funcs*n_funcs.  Each mo is stored contigously
  // until further notice the basis is non-normalised LIP
  tensor_square mo_coeffs;
  bool has_local_mo_coeffs;
  // number of mo * n_funcs_per_cell * n_cells, ie, n_funcs_per_cell * n_cells * n_funcs_per_cell * n_cells.  Each mo is stored contigously
  tensor_2 local_mo_coeffs;
  // symmetric, sum over occupied orbitals a: 2 * mo_coeff_a * mo_coeff_a
  // until further notice the basis is non-normalised LIP
  bool has_density;
  tensor_square_sym density;
  bool has_local_density;
  tensor_square_sym local_density;
  tensor_1 mo_energies;
  vector<size_t> occupied_mos; // list of occupied MO indices

  tensor_6 lip_4c2e;
  bool has_lip_4c2e;

// public:

  state1d(const system1d& sys, const vector<size_t>& occ):
      n_funcs(sys.n_gridpoints),
      locbas1d_HF_orthonorm(sys),
      mo_coeffs({n_funcs, n_funcs}),
      has_local_mo_coeffs(false),
      local_mo_coeffs({sys.n_gridpoints, sys.n_cells*sys.n_funcs_per_cell}),
      has_density(false),
      density({n_funcs, n_funcs}),
      has_local_density(false),
      local_density({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell}),
      occupied_mos(occ),
      lip_4c2e({sys.n_cells, sys.n_cells, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell}),
      has_lip_4c2e(false)
  { };
  state1d(const system1d& sys, const vector<size_t>& occ, const tensor_1& evals, const tensor_square& evecs):
      n_funcs(sys.n_gridpoints),
      locbas1d_HF_orthonorm(sys),
      mo_coeffs(evecs),
      has_local_mo_coeffs(false),
      local_mo_coeffs({sys.n_gridpoints, sys.n_cells*sys.n_funcs_per_cell}),
      has_density(false),
      density({n_funcs, n_funcs}),
      has_local_density(false),
      local_density({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell}),
      mo_energies(evals),
      occupied_mos(occ),
      lip_4c2e({sys.n_cells, sys.n_cells, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell, sys.n_funcs_per_cell}),
      has_lip_4c2e(false)
  {
    update_density(sys);
  };

  size_t n_funcs_ro() const {return n_funcs;}
  const localbasis1d& locbas1d_ro() const {return locbas1d_HF_orthonorm;}
  const tensor_square& mo_coeffs_ro() const {return mo_coeffs;}
  const tensor_square_sym& density_ro() const {return density;}
  const tensor_1& mo_energies_ro() const {return mo_energies;}

  const tensor_2& local_mo_coeffs_ro() const {return local_mo_coeffs;}
  const tensor_square_sym& local_density_ro() const {return local_density;}

  // guess occupations ... lets use sin(a*x) from 0 to pi, where a enumerates the mos
  void guess_occupied_mos(const system1d& sys);
  // update density from mo_coeffs and occupied_mo's
  void update_density(const system1d& sys);

  void set_mo_coeffs(const tensor_square& coeffs);
  void set_mo_energies(const tensor_1& en);

  void calc_and_set_local_basis(const system1d& sys, timing& times);
  void set_local_basis(const vector<tensor_square>& locbas);

  void calc_and_set_local_mo_coeffs(const system1d& sys, timing& times);

  void calc_and_set_local_density_from_local_mo_coeffs(const system1d& sys, timing& times);
  void calc_and_set_local_density_from_density(const system1d& sys, timing& times);

  // normalise the nth mo
  void normalise_mo_lip(const system1d& sys, const size_t n);
  void normalise_occupied_lip(const system1d& sys);
  void normalise_all_lip(const system1d& sys);

  void calc_and_set_4c2e_lip(const system1d& sys, const state1d& S, timing& times);

  void do_scf(const system1d& sys, const tensor_square_sym& T, const tensor_square_sym& V, tensor_square_sym& J_mat, tensor_square_sym& K_mat, double e0, const double thr, const double interaction_scaling, const bool omit_coulomb, timing& times);

  tensor_square get_H_1el_local(const system1d& sys, const operator_t op, timing& times) const;

  double get_xci_energy(const system1d& sys, const string& str, timing& times) const;
  // where bath and cell interact through a "1e-operator", and the final Hamiltonian is based on the same "1e-op"
  double get_ci_energy_apprx1(const system1d& sys, const string& str, timing& times) const;
  // where bath and cell interact through a "1e-operator", but the final H is based on an ordinary 2e-op
  double get_ci_energy_apprx2(const system1d& sys, const string& str, timing& times) const;
  // the version from 2019 01 22
  double get_ci_energy_apprx3(const system1d& sys, const string& str, timing& times) const;

  string mo_to_mathematica(const system1d& sys, const size_t mo) const;

};
  

tensor_square_sym get_kin_lip(const system1d& sys, timing& times);
tensor_square_sym get_pot_lip(const system1d& sys, timing& times);
tensor_square get_kin_local(const system1d& sys, const state1d& S, const operator_t op, timing& times);
tensor_square_sym get_pot_local(const system1d& sys, const state1d& S, timing& times);

double get_4c2eintegral_mo_local_ortho_grid(const system1d& sys, const state1d& state, size_t i, size_t j, size_t k, size_t l, timing& times);
double get_4c2eintegral_mo_precalc(const system1d& sys, const state1d& S, size_t mo_i, size_t mo_j, size_t mo_k, size_t mo_l, timing& times);
double get_4c2eintegral_1cexc_precalc(const system1d& sys, const state1d& S, const size_t cell1, const size_t cell2, const int mo_i, const int mo_j, const int mo_k, const int mo_l,
                                      const int mode, const bool normalised, int debug, timing& times);
double get_1cell_2electron_interaction(const system1d& sys, const state1d& S, const size_t cell1, const size_t cell2, const size_t mo_i, const size_t mo_j, const size_t mo_k, const size_t mo_l,
                                       const bool env, const bool normalised, timing& times);

tensor_square get_diag_single_cell_xci_mtx(const system1d& sys, const state1d& S, const tensor_square& H_1el_cell, const size_t cell, const string& str,
                                               const bool env, const bool normalised, timing& times);

tensor_square_sym get_coulomb_lip(const system1d& sys, const tensor_square_sym& dens, timing& times);
tensor_square_sym get_coulomb_lip_precalc(const system1d& sys, const state1d& state, const double scaling, timing& times);
tensor_square_sym get_coulomb_local(const system1d& sys, const state1d& state, timing& times);
tensor_square_sym get_coulomb_local_precalc(const system1d& sys, const state1d& state, timing& times);
tensor_square_sym get_exchange_lip(const system1d& sys, const tensor_square_sym& dens, timing& times);
tensor_square_sym get_exchange_local(const system1d& sys, const state1d& state, timing& times);

vector<vector<double>> tabulate_polynomial_coeffs(const vector<polynomial>& polys, const grid1d& grid);


#endif
  
