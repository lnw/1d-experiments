#include <vector>
#include <tuple>
#include <algorithm>

#include "timing.hh"
#include "tensor.hh"
#include "local_basis_1d.hh"
#include "system_1d.hh"
#include "state_1d.hh"

using namespace std;


vector<pair<size_t,size_t>> get_dets(const system1d& sys, const string& str){
  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell - n_occ;

  vector<pair<size_t, size_t>> dets;
  if(str=="s"){
    for(size_t s=0; s<n_occ+n_virt; s++){
      dets.push_back(make_pair(0,s));
    }
  }
  else if(str=="d"){
    for(size_t r=0; r<n_occ+n_virt; r++){
      for(size_t s=r; s<n_occ+n_virt; s++){
        if(r==0 && s!=0) continue;
        dets.push_back(make_pair(r,s));
      }
    }
  }
  else if(str=="sd"){
    for(size_t r=0; r<n_occ+n_virt; r++){
      for(size_t s=r; s<n_occ+n_virt; s++){
        dets.push_back(make_pair(r,s));
      }
    }
  }
  else
    assert(false);
  // cout << "dets in diag " << dets << endl;
  return dets;
}



// list of doubly excited determinants, assuming that the functions are in different cells
// ie, |ab> and |ba>
vector<pair<size_t,size_t>> get_dets_mixed(const system1d& sys){
  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell - n_occ;

  vector<pair<size_t, size_t>> dets;
  for(size_t r=1; r<n_occ+n_virt; r++){
    for(size_t s=1; s<n_occ+n_virt; s++){
      dets.push_back(make_pair(r,s));
    }
  }
  //cout << "dets in diag " << dets << endl;

  return dets;
}


// express the coulomb interaction with the other cells as a 1-electron op.
tensor_square_sym get_C1e_local(const system1d& sys, const state1d &S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  tensor_square_sym C1e_global({sys.n_cells*sys.n_funcs_per_cell, sys.n_cells*sys.n_funcs_per_cell});
  const tensor_6& lip_ints = S.lip_4c2e;
  const vector<tensor_square>& basis = S.locbas1d_HF_orthonorm.basis;
  // const vector<tensor_square>& basis_inv = S.locbas1d_HF_orthonorm.basis_inv;

  for(int cell=0; cell<sys.n_cells; cell++){

    // fill V with a x b entries (a b|lambda sigma)
    tensor_square_sym V_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});

    // first get coeffs of the lowest MO (= first local func), modify one cell later
    vector<tensor_1> D_lambda(sys.n_cells); // expansion of local funcs in LIP
    const int mo_0 = 0; // the occupied local orbital
    for(int c=0; c<sys.n_cells; c++){
      // cout << c << endl;
      tensor_1 C_lambda({sys.n_funcs_per_cell}); // expansion of the first mo (or a local func) in local funcs
      C_lambda[0] = S.local_mo_coeffs(mo_0, c*sys.n_funcs_per_cell+0); // and the rest are zero
      D_lambda[c] = basis[c] * C_lambda;
    }
    vector<tensor_1> D_sigma(D_lambda);

    // mu, nu: the single cell
    // lambda, sigma: all other cells
    // for(int cell_munu=0; cell_munu<sys.n_cells; cell_munu++){ // cell of mu and nu
    const int cell_munu = cell;
    for(int cell_lasi=0; cell_lasi<sys.n_cells; cell_lasi++){ // cell of lambda and sigma
      if(cell_lasi == cell_munu) continue;
      // { int cell_lasi = cell2;
      for(int mu=0; mu<sys.n_funcs_per_cell; mu++){ // mu, an LIP
        // double tmp_nu=0;
        for(int nu=0; nu<sys.n_funcs_per_cell; nu++){ // nu, an LIP // do only one triangle
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
          V_cell(mu,nu) += tmp_lambda;
          // tmp_nu += tmp_lambda; // * D_nu[cell_munu][nu];
        }
        // int_4c2e += tmp_nu; // * D_mu[cell_munu][mu];
      }
    }

    tensor_square_sym V_cell_local = matmul(basis[cell], matmul(V_cell, basis[cell], 'N', 'N'), 'T', 'N');

    for(int i=0; i<sys.n_funcs_per_cell; i++){
      const int gi = cell*sys.n_funcs_per_cell + i;
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        const int gj = cell*sys.n_funcs_per_cell + j;
        C1e_global(gi,gj) = V_cell_local(i,j);
      }
    }
  }

  auto t1 = Clock::now();
  times.append_item("C pot 1el (local)", t1 - t0);
  times.ind_down();
  return C1e_global;
}



// for a single cell, return the dci matrix where excitations happen in the local basis of that cell
// r,s,t,u are local functions
// (mu nu|lambda sigma) are LIP
tensor_square get_diag_single_cell_xci_mtx(const system1d& sys, const state1d& S, const tensor_square& H_1el_cell, const size_t cell, const string& str,
                                           const bool env, const bool normalised, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell-n_occ;
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  size_t n_dets;
  if(str=="s")
    n_dets = n_occ + n_singly_exc_dets;
  else if(str=="d")
    n_dets = n_occ + n_doubly_exc_dets;
  else if(str=="sd")
    n_dets = n_occ + n_singly_exc_dets + n_doubly_exc_dets;
  else
    assert(false);
  // cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store DCI matrix" << endl;

  cout << "getting 1 electron integrals for all local functions ..." << endl;
  tensor_square one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  double coeff;
  if(!normalised)
    coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0]; // coeff of the first local function
  else
    coeff = 1.0;

  for(int i=0; i<n_occ+n_virt; i++){
    for(int j=0; j<n_occ+n_virt; j++){
      one_el_integrals(i,j) = H_1el_cell(i,j) * coeff * coeff;
    }
  }
  // cout << "1e ints " << one_el_integrals << endl;
  cout << "done" << endl;

  vector<pair<size_t, size_t>> dets(get_dets(sys, str));
  // cout << "dets in diag " << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square sc_xci_matrix({n_dets, n_dets});
  // cout << "sc " << str << "ci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
  for(size_t det1=0; det1<n_dets; det1++){
    for(size_t det2=det1; det2<n_dets; det2++){
      const size_t r=dets[det1].first,
                   s=dets[det1].second,
                   t=dets[det2].first,
                   u=dets[det2].second;
      // if( (r+s+t+u)%2 ) continue; // (tu|rs) is zero for odd parity
      // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;

      if(det1==det2 && r==s){ // a
        const double entry = 2*one_el_integrals(r,r) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, r, r, r, env, normalised, times); // aka rtsu
        sc_xci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if(det1==det2 && r!=s){ // b
        const double entry = one_el_integrals(r,r) +
                             one_el_integrals(s,s) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, s, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, s, r, s, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        // cout << "b" << endl;
      }
      else if( r==s && t==u ){ // d
        const double entry = get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, t, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "d" << endl;
      }
      else if( (r==s && r==t)){ // f1
        const double el1_1 = sqrt(2) * one_el_integrals(s,u);
        const double el1_2 = sqrt(2) * one_el_integrals(u,s);
        const double el2 = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, r, r, u, env, normalised, times);
        // const double entry = term1 + term2;
        // if(cell==15 || cell==16) cout << r << s << t << u << ":" << entry << " = " << term1 << " + " << term2 << endl;
        sc_xci_matrix(det1, det2) = el1_1 + el2;
        sc_xci_matrix(det2, det1) = el1_2 + el2;
        // cout << "f1" << endl;
      }
      // else if( (t==u && t==r)){ // f2
      //   const double entry = sqrt(2) * one_el_integrals(r,s) +
      //                        sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, r, times);
      //   sc_xci_matrix(det1, det2) = entry;
      //   cout << "f2" << endl;
      // }
      else if( (s==u && s==t)){ // f3
        const double term1 = sqrt(2) * one_el_integrals(r,t);
        const double term2 = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, t, t, env, normalised, times);
        const double entry = term1 + term2;
        // if(cell==15 || cell==16) cout << r << s << t << u << ":" << entry << " = " << term1 << " + " << term2 << endl;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "f3" << endl;
      }
      // else if( (s==u && s==r)){ // f4
      //   const double entry = sqrt(2) * one_el_integrals(r,t) +
      //                        sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, r, times);
      //   sc_xci_matrix(det1, det2) = entry;
      //   cout << "f4" << endl;
      // }
      else if( r==t && s!=u ){ // c1
        const double entry = one_el_integrals(s,u) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, u, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, s, u, r, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "c1" << endl;
      }
      else if( (s==t && r!=u) ){ // c2
        const double entry = one_el_integrals(r,u) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, s, s, u, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, s, s, u, r, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "c2" << endl;
      }
      else if( (r!=t && s==u) ){ // c3
        const double entry = one_el_integrals(r,t) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, s, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, t, s, s, r, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "c3" << endl;
      }
      // else if( (s!=t && r==u) ){ // c4
      //   const double entry = one_el_integrals(s,t) +
      //                        get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, r, times) +
      //                        get_1cell_2electron_interaction(sys, S, cell, cell, t, s, r, r, times);
      //   sc_xci_matrix(det1, det2) = 11;//entry;
      //   cout << "c4" << endl;
      // }
      else if( r==s ){ // g
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, u, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "g1" << endl;
      }
      else if( t==u ){ // g
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, t, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "g2" << endl;
      }
      else if( r!=s && r!=t && r!=u && t!=u ){ // e
        const double entry = get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, u, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, t, s, u, r, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "e" << endl;
      }
      else
        assert(false);

    } // det2
  } // det1

  auto t1 = Clock::now();
  times.append_item("sc xci matrix", t1 - t0);
  times.ind_down();

  return sc_xci_matrix;
}



// for a single cell, return the xci matrix where excitations happen from the HF MO-0 to the local basis of that cell
// r,s,t,u are local functions
// (mu nu|lambda sigma) are LIP
// excitations are from the entire HF MO-0, to the virtual space in the chosen cell
tensor_square_sym get_diag_single_cell_xci_mtx_v2(const system1d& sys, const state1d& S, const tensor_square& H_1el_local, const size_t cell,
                                                  const string& str, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const bool env = false;

  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell - n_occ;
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  size_t n_dets;
  if(str=="s")
    n_dets = n_occ + n_singly_exc_dets;
  else if(str=="d")
    n_dets = n_occ + n_doubly_exc_dets;
  else if(str=="sd")
    n_dets = n_occ + n_singly_exc_dets + n_doubly_exc_dets;
  else
    assert(false);
  // cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store DCI matrix" << endl;

  // write down one-el-integrals for each cell, keeping in mind, that local-function-0 is part of the HF MO-0 and therefore not normalised, while all virtual functions are normalised.
  cout << "getting 1 electron integrals for all local functions ..." << endl;
  vector<tensor_square> one_el_integrals;
  for(size_t cell1=0; cell1<sys.n_cells; cell1++){
    tensor_square tmp({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    double coeff_i, coeff_j;
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      if(i==0) coeff_i = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0]; // coeff of the first local function
      else coeff_i = 1.0;
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        if(j==0) coeff_j = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0]; // coeff of the first local function
        else coeff_j = 1.0;
        tmp(i,j) = H_1el_local(cell1*sys.n_funcs_per_cell+i, cell1*sys.n_funcs_per_cell+j) * coeff_i * coeff_j;
      }
    }
    one_el_integrals.push_back(tmp);
  }
  cout << "1e ints " << one_el_integrals[cell] << endl;
  cout << "done" << endl;

  vector<pair<size_t, size_t>> dets(get_dets(sys, str));
  // cout << "dets in diag " << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym sc_xci_matrix({n_dets, n_dets});
  // cout << "sc " << str << "ci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
  for(size_t det1=0; det1<n_dets; det1++){
    for(size_t det2=det1; det2<n_dets; det2++){
      const size_t r=dets[det1].first,
                   s=dets[det1].second,
                   t=dets[det2].first,
                   u=dets[det2].second;
      // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;

      if(det1==det2 && r==s){ // r, s, t, u can be 0
        if(r==0){
          const bool normalised = false;
          double el1=0, el2=0;
          for(size_t cell1=0; cell1<sys.n_cells; cell1++){
            el1 += 2*one_el_integrals[cell1](r,r);
            for(size_t cell2=0; cell2<sys.n_cells; cell2++){
              el2 += get_1cell_2electron_interaction(sys, S, cell1, cell2, r, r, r, r, env, normalised, times); // aka rtsu
            }
          }
          // cout << "2 x 1el, 2el, sum: " << el1 << ", " << el2 << ", " << el1 + el2 << endl;
          sc_xci_matrix(det1, det2) += el1 + el2;
        }
        else{
          const bool normalised = true;
          const double entry = 2*one_el_integrals[cell](r,r) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, r, r, r, env, normalised, times); // aka rtsu
          sc_xci_matrix(det1, det2) = entry;
        }
        cout << "a" << endl;
      }
      else if(det1==det2 && r!=s){ // r, t can be 0
        if(r==0){
          const bool normalised = true;
          double el1_glob=0, el1_loc=0, el2_j=0, el2_k=0;
          const double coeff_sq = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0] * S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
          for(size_t cell1=0; cell1<sys.n_cells; cell1++){
            el1_glob += one_el_integrals[cell1](r,r);
            el2_j += get_1cell_2electron_interaction(sys, S, cell1, cell, r, r, s, s, env, normalised, times) * coeff_sq;
          }
          el1_loc += one_el_integrals[cell](s,s);
          el2_k += get_1cell_2electron_interaction(sys, S, cell, cell, r, s, r, s, env, normalised, times) * coeff_sq;
          // cout << "1g, 1l, 2j, 2k: " <<  el1_glob <<", "<< el1_loc <<", "<< el2_j <<", "<< el2_k << endl;
          sc_xci_matrix(det1, det2) = el1_glob + el1_loc + el2_j + el2_k;
        }
        else{
          const bool normalised = true;
          const double entry = one_el_integrals[cell](r,r) +
                               one_el_integrals[cell](s,s) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, s, env, normalised, times) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, s, r, s, env, normalised, times);
          sc_xci_matrix(det1, det2) = entry;
        }
        cout << "b" << endl;
      }
      else if(r==s && t==u){ // r, s can be 0
        double coeff_sq;
        if(r==0) coeff_sq = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0] * S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
        else coeff_sq = 1.0;
        const bool normalised = true;
        const double entry = get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, t, env, normalised, times) * coeff_sq;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "d" << endl;
      }
      else if(r==s && r==t){ // r, s, t can be 0
        if(r==0){
          const bool normalised = false;
          double el1 = sqrt(2) * one_el_integrals[cell](s,u), el2=0;
          for(size_t cell1=0; cell1<sys.n_cells; cell1++){
            el2 += sqrt(2) * get_1cell_2electron_interaction(sys, S, cell1, cell, r, r, r, u, env, normalised, times) / S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
          }
          // cout << "1el, 2el, sum: " << el1 << ", " << el2 << ", " << el1+el2 << endl;
          sc_xci_matrix(det1, det2) = el1+el2;
          sc_xci_matrix(det2, det1) = el1+el2;
        }
        else {
          const bool normalised = true;
          const double term1 = sqrt(2) * one_el_integrals[cell](s,u);
          const double term2 = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, r, r, u, env, normalised, times);
          const double entry = term1 + term2;
          sc_xci_matrix(det1, det2) = entry;
          sc_xci_matrix(det2, det1) = entry;
        }
        cout << "f1" << endl;
      }
//      // else if( (t==u && t==r)){ // f2
//      //   const double entry = sqrt(2) * one_el_integrals(r,s) +
//      //                        sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, r, times);
//      //   sc_xci_matrix(det1, det2) = entry;
//      //   cout << "f2" << endl;
//      // }
      else if(s==u && s==t){ // r can be 0
        const double term1 = sqrt(2) * one_el_integrals[cell](r,t);
        double coeff;
        if(r==0) coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
        else coeff = 1.0;
        const bool normalised = true;
        const double term2 = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, t, t, env, normalised, times) * coeff;
        const double entry = term1 + term2;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
        // cout << "f3" << endl;
      }
//      // else if( (s==u && s==r)){ // f4
//      //   const double entry = sqrt(2) * one_el_integrals(r,t) +
//      //                        sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, r, times);
//      //   sc_xci_matrix(det1, det2) = entry;
//      //   cout << "f4" << endl;
//      // }
      else if(r==t && s!=u){ // r, t can be 0
        if(r==0){
          const bool normalised = true;
          const double coeff_sq = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0] * S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
          double entry = one_el_integrals[cell](s,u);
          for(size_t cell1=0; cell1<sys.n_cells; cell1++){
            entry += get_1cell_2electron_interaction(sys, S, cell1, cell, r, r, s, u, env, normalised, times) * coeff_sq;
          }
          entry += get_1cell_2electron_interaction(sys, S, cell, cell, r, s, u, r, env, normalised, times) * coeff_sq;
          sc_xci_matrix(det1, det2) = entry;
          sc_xci_matrix(det2, det1) = entry;
        }
        else{
          const bool normalised = true;
          const double entry = one_el_integrals[cell](s,u) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, r, s, u, env, normalised, times) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, s, u, r, env, normalised, times);
          sc_xci_matrix(det1, det2) = entry;
          sc_xci_matrix(det2, det1) = entry;
        }
        cout << "c1" << endl;
      }
      else if( s==t && r!=u ){ // r can be 0
        // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;
        if(r==0){
          const bool normalised = true;
          const double coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
          double entry = one_el_integrals[cell](r,u) +
                         get_1cell_2electron_interaction(sys, S, cell, cell, r, s, s, u, env, normalised, times) * coeff +
                         get_1cell_2electron_interaction(sys, S, cell, cell, s, s, u, r, env, normalised, times) * coeff;
          sc_xci_matrix(det1, det2) = entry;
          sc_xci_matrix(det2, det1) = entry;
        }
        else{
          const bool normalised = true;
          const double entry = one_el_integrals[cell](r,u) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, r, s, s, u, env, normalised, times) +
                               get_1cell_2electron_interaction(sys, S, cell, cell, s, s, u, r, env, normalised, times);
          sc_xci_matrix(det1, det2) = entry;
          sc_xci_matrix(det2, det1) = entry;
        }
        // cout << "c2" << endl;
      }
      else if( r!=t && s==u ){ // r can be 0
        const bool normalised = true;
        double coeff;
        if(r==0) coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
        else coeff = 1.0;
        const double entry = one_el_integrals[cell](r,t) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, s, env, normalised, times) * coeff +
                             get_1cell_2electron_interaction(sys, S, cell, cell, t, s, s, r, env, normalised, times) * coeff;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
          // cout << "c3" << endl;
      }
//      // else if( (s!=t && r==u) ){ // c4
//      //   const double entry = one_el_integrals(s,t) +
//      //                        get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, r, times) +
//      //                        get_1cell_2electron_interaction(sys, S, cell, cell, t, s, r, r, times);
//      //   sc_xci_matrix(det1, det2) = 11;//entry;
//      //   cout << "c4" << endl;
//      // }
      else if( r==s ){ // r and s can be 0
        const bool normalised = true;
        double coeff_sq;
        if(r==0) coeff_sq = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0] * S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0] ;
        else coeff_sq = 1.0;
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, u, env, normalised, times) * coeff_sq;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
//        // cout << "g1" << endl;
      }
      else if( t==u ){ // r can be 0
        const bool normalised = true;
        double coeff;
        if(r==0) coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
        else coeff = 1.0;
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, t, env, normalised, times) * coeff;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
//        // cout << "g2" << endl;
      }
      else if( r!=s && r!=t && r!=u && t!=u ){ // r can be 0
        const bool normalised = true;
        double coeff;
        if(r==0) coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
        else coeff = 1.0;
        const double entry = get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, u, env, normalised, times) * coeff +
                             get_1cell_2electron_interaction(sys, S, cell, cell, t, s, u, r, env, normalised, times) * coeff;
        sc_xci_matrix(det1, det2) = entry;
        sc_xci_matrix(det2, det1) = entry;
//        // cout << "e" << endl;
      }
     else
       assert(false);

    } // det2
  } // det1

  auto t1 = Clock::now();
  times.append_item("sc xci matrix", t1 - t0);
  times.ind_down();

  return sc_xci_matrix;
}



// for a combination of two cell1 != cell2, calculate the XCI matrix.
// the CI matrix contains only 2-electron contributions because the two cells don't overlap.
// r,s,t,u are local functions
// unlike the symmetric case, no further normalisation is required
// env: include some kind of MO_0 environment
tensor_square get_offdiag_single_cell_xci_mtx(const system1d& sys, const state1d& S, const size_t cell1, const size_t cell2, const string& str,
                                              const bool env, const bool normalised, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell-n_occ;
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  size_t n_dets;
  if(str=="s")
    n_dets = n_occ + n_singly_exc_dets;
  else if(str=="d")
    n_dets = n_occ + n_doubly_exc_dets;
  else if(str=="sd")
    n_dets = n_occ + n_singly_exc_dets + n_doubly_exc_dets;
  else
    assert(false);
  // cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store DCI matrix" << endl;

  vector<pair<size_t, size_t>> dets;
  if(str=="s"){
    for(size_t s=0; s<n_occ+n_virt; s++){
      dets.push_back(make_pair(0,s));
    }
  }
  else if(str=="d"){
    for(size_t r=0; r<n_occ+n_virt; r++){
      for(size_t s=r; s<n_occ+n_virt; s++){
        if(r==0 && s!=0) continue;
        dets.push_back(make_pair(r,s));
      }
    }
  }
  else if(str=="sd"){
    for(size_t r=0; r<n_occ+n_virt; r++){
      for(size_t s=r; s<n_occ+n_virt; s++){
        dets.push_back(make_pair(r,s));
      }
    }
  }
  else
    assert(false);
  // cout << "dets in off-diag ... " << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym sc_xci_matrix({n_dets, n_dets});

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
  for(size_t det1=0; det1<n_dets; det1++){
    for(size_t det2=0; det2<n_dets; det2++){ // not triangular because off-diagonal CI mtx is not symmetric
      const size_t r=dets[det1].first,
                   s=dets[det1].second,
                   t=dets[det2].first,
                   u=dets[det2].second;
      // if( (r+s+t+u)%2 ) continue; // (tu|rs) is zero for odd parity
      // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;

      if(r==s && t==u){ // a
        const double entry = get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, r, t, env, normalised, times); // aka rtsu
        sc_xci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if( (r!=s && t==u) || (r==s && t!=u) ){ // b
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, s, u, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if(r!=s && t!=u){ // c
        const double entry = get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, s, u, env, normalised, times) +
                             get_1cell_2electron_interaction(sys, S, cell1, cell2, t, s, u, r, env, normalised, times);
        sc_xci_matrix(det1, det2) = entry;
        // cout << "b" << endl;
      }
      else
        assert(false);

    } // det2
  } // det1

  // cout << "offdiag scimtx " << cell1 << "," << cell2 << " : " << sc_xci_matrix << endl;

//   const double coeff1 = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0]; // n electrons in cell1
//   const double coeff2 = S.local_mo_coeffs[cell2*sys.n_funcs_per_cell + 0]; // n electrons in cell2
// if(cell1==0 && cell2==1) cout << "cb " << sc_xci_matrix << endl;
//   sc_xci_matrix = sc_xci_matrix * (1.0/(coeff1*coeff2)); // renormalise each ci matrix (in each cell) to 1.0
// if(cell1==0 && cell2==1) cout << "ca " << sc_xci_matrix << endl;

  auto t1 = Clock::now();
  times.append_item("sc offdiag xci matrix", t1 - t0);
  times.ind_down();

  return sc_xci_matrix;
}



// for a single cell, return the dci matrix where excitations happen in the local basis of that cell
// r,s,t,u are local functions
// (mu nu|lambda sigma) are LIP
tensor_1 get_offdiag_single_cell_xci_mtx_v2(const system1d& sys, const state1d& S, const tensor_square& H_1el_cell, const size_t cell, const string& str,
                                            const bool env, const bool normalised, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = sys.n_funcs_per_cell-n_occ;
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  size_t n_dets;
  if(str=="s")
    n_dets = n_occ + n_singly_exc_dets;
  else if(str=="d")
    n_dets = n_occ + n_doubly_exc_dets;
  else if(str=="sd")
    n_dets = n_occ + n_singly_exc_dets + n_doubly_exc_dets;
  else
    assert(false);
  // cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store DCI matrix" << endl;

  cout << "getting 1 electron integrals for all local functions ..." << endl;
  tensor_square_sym one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  double coeff;
  if(!normalised)
    coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0]; // coeff of the first local function
  else
    coeff = 1.0;

  for(int i=0; i<n_occ+n_virt; i++){
    for(int j=0; j<n_occ+n_virt; j++){
      one_el_integrals(i,j) = H_1el_cell(i,j) * coeff * coeff;
    }
  }
  // cout << "1e ints " << one_el_integrals << endl;
  cout << "done" << endl;

  vector<pair<size_t, size_t>> dets(get_dets(sys, str));
  // cout << "dets in offdiag " << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_1 sc_xci_matrix({n_dets});
  // cout << "sc " << str << "ci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
  for(size_t det1=0; det1<n_dets; det1++){
    const size_t r=dets[det1].first,
                 s=dets[det1].second,
                 t=0,
                 u=0;
      // if( (r+s+t+u)%2 ) continue; // (tu|rs) is zero for odd parity
      // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;

      if(r==0 && r==s){ // |00>
        const double entry = 2*one_el_integrals(r,r) +
                             get_1cell_2electron_interaction(sys, S, cell, cell, r, r, r, r, env, normalised, times); // aka rtsu
        sc_xci_matrix(det1) = entry;
        // cout << "a" << endl;
      }
      else if(r==0 && r!=s){ // |0a>
        const double entry = sqrt(2) * one_el_integrals(r,s) +
                             sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, u, env, normalised, times);
        sc_xci_matrix(det1) = entry;
        // cout << "b" << endl;
      }
      else if(r==s){ // |aa>
        const double entry = get_1cell_2electron_interaction(sys, S, cell, cell, r, t, r, t, env, normalised, times);
        sc_xci_matrix(det1) = entry;
        // cout << "d" << endl;
      }
      else if(r!=s){ // |ab>
        const double entry = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell, cell, r, t, s, u, env, normalised, times);
        sc_xci_matrix(det1) = entry;
        // cout << "f1" << endl;
      }
      else
        assert(false);

  } // det1

  // cout << "offdiag matrix: " << sc_xci_matrix << endl;

  auto t1 = Clock::now();
  times.append_item("sc xci matrix", t1 - t0);
  times.ind_down();

  return sc_xci_matrix;
}



// where CI matrices are one per cell,
// and C matrices contain the first n CI vectors for each cell.  obtained by diagonalising CI matrices
// two2one: include a 1e-coulomb op
// env: include some 2e-coulomb environment based on mo_0
tuple<vector<tensor_2>, vector<tensor_square>> get_C_CI_matrices(const system1d& sys, const state1d& S, size_t n_ci_funcs_per_cell, const string& str,
                                                                 const bool two2one, const bool env, timing& times) {
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);
  assert(S.has_lip_4c2e);
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  vector<tensor_2> C_matrices(sys.n_cells);
  vector<tensor_square> diag_CI_matrices(sys.n_cells);

  cout << "getting 1 electron Hamiltonian ..." << endl;
  const tensor_square_sym T_local = get_kin_local(sys, S, NABLA2_APPR, times);
  const tensor_square_sym V_local = get_pot_local(sys, S, times);
  tensor_square_sym H_1el_local;
  if(two2one){
    const tensor_square_sym C1e_local = get_C1e_local(sys, S, times);
    H_1el_local = T_local + V_local + C1e_local;
    cout << "t = " << T_local << endl;
    cout << "v = " << V_local << endl;
    cout << "c1e = " << C1e_local << endl;
  }
  else{
    H_1el_local = T_local + V_local;
  }
  cout << "done" << endl;

  // for each cell, do fci in the field of all other cells, using the local
  // orbitals in each cells.  So, as in the global MO case, the first local
  // function is occupied, the rest is unoccupied.  Orbitals are orthogonal and
  // normalised but really there is a fractional (but constant) number of
  // electrons in a cell.
  // iterate?

  for(size_t cell=0; cell<sys.n_cells; cell++){

    // slice local copy of cell from 1e-H
    tensor_square H_1el_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        H_1el_cell(i,j) = H_1el_local(cell*sys.n_funcs_per_cell+i, cell*sys.n_funcs_per_cell+j);
      }
    }

    const bool normalised = true;
    tensor_square_sym sc_xci_mtx(get_diag_single_cell_xci_mtx(sys, S, H_1el_cell, cell, str, env, normalised, times));
    // cout << str << "ci matrix: " << sc_xci_mtx << endl;

    const size_t n_dets = sc_xci_mtx.get_dims()[0];

    vector<double> diagonals;
    for(size_t i=0; i<sc_xci_mtx.get_dims()[0]; i++){
      diagonals.push_back(sc_xci_mtx(i,i));
    }
    sort(diagonals.begin(), diagonals.end());

    // diagonalise xci-hamiltonian
    tensor_1 evals;
    tensor_square evecs;
    tie(evals, evecs) = sc_xci_mtx.eigensystem_sym(times);

    cout << "evals(" << cell << "): " << evals << endl;
    cout << "HF levels (sorted): " << diagonals << endl;
    cout << "correlation in this cell: " << evals[0] - diagonals[0] << endl;
    // cout << evecs << endl;
    // cout << str << "ci coeffs of first state: " << evecs.get_slice(0, XX) << endl << endl;

    tensor_2 C_matrix({n_ci_funcs_per_cell, n_dets});
    for(size_t n=0; n<n_ci_funcs_per_cell; n++){ // 1, 2, ...
      const tensor_1 slice(evecs.get_slice(n, XX));
      for(size_t i=0; i<n_dets; i++){ // 16, ...
        C_matrix(n, i) = slice[i];
      }
    }
    // const double coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
    // sc_xci_mtx = sc_xci_mtx * (1.0/(coeff*coeff)); // renormalise each ci matrix (in each cell) to 1.0

    // cout << "cmtx" << C_matrix << endl;
    // C_matrix = tensor_2({n_ci_funcs_per_cell, n_dets}); C_matrix(0,0) = 1.0; // DEBUG, FIXME
    // cout << "cmtx" << C_matrix << endl;
    C_matrices[cell] = C_matrix;
    diag_CI_matrices[cell] = sc_xci_mtx;
  }

  auto t1 = Clock::now();
  times.append_item("get C/CI matrices", t1 - t0);
  times.ind_down();

  return make_tuple(C_matrices, diag_CI_matrices);
}



// where CI matrices are one per cell,
// the considered functions are the HF ground state MO + all local functions in one cell (without the one that contributes to the first MO)
tuple<vector<tensor_1>, vector<tensor_2>, vector<tensor_square_sym>> get_C_CI_matrices_v2(const system1d& sys, const state1d& S, size_t n_ci_funcs_per_cell,
                                                                                      const tensor_square& H_1el_local,
                                                                                      const string& str, timing& times) {
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);
  assert(S.has_lip_4c2e);
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  vector<tensor_1> energies(sys.n_cells);
  vector<tensor_2> C_matrices(sys.n_cells);
  vector<tensor_square_sym> diag_CI_matrices(sys.n_cells);

  for(size_t cell=0; cell<sys.n_cells; cell++){

    tensor_square_sym sc_xci_mtx(get_diag_single_cell_xci_mtx_v2(sys, S, H_1el_local, cell, str, times));
    cout << str << "cimatrix" << cell << " = " << sc_xci_mtx << endl;

    const size_t n_dets = sc_xci_mtx.get_dims()[0];

    vector<double> diagonals;
    for(size_t i=0; i<sc_xci_mtx.get_dims()[0]; i++){
      diagonals.push_back(sc_xci_mtx(i,i));
    }
    sort(diagonals.begin(), diagonals.end());

    // diagonalise xci-hamiltonian
    tensor_1 evals;
    tensor_square evecs;
    tie(evals, evecs) = sc_xci_mtx.eigensystem_sym(times);

    cout << "evals(" << cell << "): " << evals << endl;
    cout << "HF levels (sorted): " << diagonals << endl;
    cout << "correlation in this cell: " << evals[0] - diagonals[0] << endl;
    // cout << evecs << endl;
    cout << str << "ci coeffs of first state: " << evecs.get_slice(0, XX) << endl;
    // cout << str << "ci coeffs of first state (right):" << evecs_r.get_slice(0, XX) << endl << endl;

    tensor_2 C_matrix({n_ci_funcs_per_cell, n_dets});
    tensor_1 energy({n_ci_funcs_per_cell});
    for(size_t n=0; n<n_ci_funcs_per_cell; n++){ // 1, 2, ...
      energy[n] = evals[n];
      const tensor_1 slice(evecs.get_slice(n, XX));
      for(size_t i=0; i<n_dets; i++){ // 16, ...
        C_matrix(n, i) = slice[i];
      }
    }
    energies[cell] = energy;
    C_matrices[cell] = C_matrix;
    diag_CI_matrices[cell] = sc_xci_mtx;

    cout << endl;
  }

  auto t1 = Clock::now();
  times.append_item("get C/CI matrices v2", t1 - t0);
  times.ind_down();

  return make_tuple(energies, C_matrices, diag_CI_matrices);
}



// where CI matrices are one per cell,
// the considered functions are the HF ground state MO + all local functions in one cell (without the one that contributes to the first MO)
double get_exact_twocell_e(const system1d& sys, const state1d& S,
                           const tensor_square_sym& H_1el_local,
                           const vector<pair<size_t, size_t>>& dets, const vector<pair<size_t, size_t>>& dets_mixed,
                           const tensor_square& CI_mat_1, const tensor_square& CI_mat_2, const size_t cell1, const size_t cell2,
                           const string& str, timing& times) {
  times.ind_up();
  auto t0 = Clock::now();

  const bool env = false;

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);
  assert(S.has_lip_4c2e);
  assert(cell1 != cell2);

  vector<tensor_1> energies(sys.n_cells);
  vector<tensor_2> C_matrices(sys.n_cells);
  vector<tensor_square> diag_CI_matrices(sys.n_cells);

  // write down one-el-integrals for each cell, keeping in mind, that local-function-0 is part of the HF MO-0 and therefore not normalised, while all virtual functions are normalised.
  cout << "getting 1 electron integrals for all local functions ..." << endl;
  vector<tensor_square> one_el_integrals;
  for(size_t c1=0; c1<sys.n_cells; c1++){
    tensor_square tmp({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    double coeff_i, coeff_j;
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      if(i==0) coeff_i = S.local_mo_coeffs[c1*sys.n_funcs_per_cell + 0]; // coeff of the first local function
      else coeff_i = 1.0;
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        if(j==0) coeff_j = S.local_mo_coeffs[c1*sys.n_funcs_per_cell + 0]; // coeff of the first local function
        else coeff_j = 1.0;
        tmp(i,j) = H_1el_local(c1*sys.n_funcs_per_cell+i, c1*sys.n_funcs_per_cell+j) * coeff_i * coeff_j;
      }
    }
    one_el_integrals.push_back(tmp);
  }
  // cout << "1e ints " << one_el_integrals << endl;
  cout << "done" << endl;

  const size_t n_dets_in = dets.size();
  const size_t n_dets_mixed = dets_mixed.size();
  const size_t n_dets_out = 2*n_dets_in - 1 + n_dets_mixed;

  tensor_square_sym sc_xci_mtx({n_dets_out, n_dets_out});

  // copy all reusable elements from the one-cell matrices
  for(size_t i=0; i<n_dets_in; i++){
    for(size_t j=i; j<n_dets_in; j++){
      sc_xci_mtx(i,j) = CI_mat_1(i,j);
      sc_xci_mtx(j,i) = CI_mat_1(i,j);
    }
  }
  for(size_t i=1; i<n_dets_in; i++){
    sc_xci_mtx(0, n_dets_in + i - 1) = CI_mat_2(0,i);
    sc_xci_mtx(n_dets_in + i - 1, 0) = CI_mat_2(0,i);
  }
  for(size_t i=1; i<n_dets_in; i++){
    for(size_t j=i; j<n_dets_in; j++){
      sc_xci_mtx(n_dets_in + i - 1, n_dets_in + j - 1) = CI_mat_2(i,j);
      sc_xci_mtx(n_dets_in + j - 1, n_dets_in + i - 1) = CI_mat_2(i,j);
    }
  }

  // a in cell1, b in cell2
  // calculate remaining elements <0a||0b>
  for(size_t det1=1; det1<n_dets_in; det1++){
    for(size_t det2=1; det2<n_dets_in; det2++){
      const size_t r=dets[det1].first,
                   s=dets[det1].second,
                   t=dets[det2].first,
                   u=dets[det2].second;
      if(r != 0 || t != 0) continue; // only |0a> contributes

      const bool normalised = true;
      double coeff1 = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0];
      double coeff2 = S.local_mo_coeffs[cell2*sys.n_funcs_per_cell + 0];
      // det2 in cell2
      const double entry1 = get_1cell_2electron_interaction(sys, S, cell1, cell2, r, s, t, u, env, normalised, times) * coeff1 * coeff2;
      sc_xci_mtx(det1, n_dets_in + det2 - 1) = entry1;
      sc_xci_mtx(n_dets_in + det2 - 1, det1) = entry1;
      // cout << "a" << endl;
    }
  }

  // calculate remaining elements <00||ab>
  for(size_t detm=0; detm<n_dets_mixed; detm++){
    const size_t r=0,
                 s=0,
                 t=dets_mixed[detm].first,
                 u=dets_mixed[detm].second;

    const bool normalised = true;
    double coeff1 = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0];
    double coeff2 = S.local_mo_coeffs[cell2*sys.n_funcs_per_cell + 0];
    // det2 in cell2
    const double entry1 = sqrt(2) * get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, s, u, env, normalised, times) * coeff1 * coeff2;
    sc_xci_mtx(0, 2*n_dets_in - 1 + detm) = entry1;
    sc_xci_mtx(2*n_dets_in - 1 + detm, 0) = entry1;
    // cout << "a" << endl;
  }

  // calculate remaining elements <0a||ab> and <0b||ab>
  for(size_t det1=1; det1<n_dets_in; det1++){ // |0a>
    const size_t r=dets[det1].first,
                 s=dets[det1].second;
    if(r != 0) continue; // only |0a> contributes
    for(size_t detm=0; detm<n_dets_mixed; detm++){ // |ab>
      const size_t t=dets_mixed[detm].first,
                   u=dets_mixed[detm].second;

      const bool normalised = true;
      double coeff1 = S.local_mo_coeffs[cell1*sys.n_funcs_per_cell + 0];
      double coeff2 = S.local_mo_coeffs[cell2*sys.n_funcs_per_cell + 0];
      // <0a||ab>
      double entry1 = get_1cell_2electron_interaction(sys, S, cell2, cell1, r, u, s, t, env, normalised, times) * coeff2;
      if(s==t) entry1 += one_el_integrals[cell2](r,u);
      sc_xci_mtx(det1, 2*n_dets_in - 1 + detm) = entry1;
      sc_xci_mtx(2*n_dets_in - 1 + detm, det1) = entry1;
      // <0b||ab>
      double entry2 = get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, s, u, env, normalised, times) * coeff1;
      if(s==u) entry2 += one_el_integrals[cell1](r,t);
      sc_xci_mtx(n_dets_in - 1 + det1, 2*n_dets_in - 1 + detm) = entry2;
      sc_xci_mtx(2*n_dets_in - 1 + detm, n_dets_in - 1 + det1) = entry2;
      // cout << "a" << endl;
    }
  }

  // calculate remaining elements <ab||ab>
  for(size_t detm1=0; detm1<n_dets_mixed; detm1++){
    for(size_t detm2=0; detm2<n_dets_mixed; detm2++){
      const size_t r=dets_mixed[detm1].first,
                   s=dets_mixed[detm1].second,
                   t=dets_mixed[detm2].first,
                   u=dets_mixed[detm2].second;

      const bool normalised = true;
      // det2 in cell2
      double entry = get_1cell_2electron_interaction(sys, S, cell1, cell2, r, t, s, u, env, normalised, times);
      if(s==u) entry += one_el_integrals[cell1](r,t);
      if(r==t) entry += one_el_integrals[cell2](s,u);
      sc_xci_mtx(2*n_dets_in - 1 + detm1, 2*n_dets_in - 1 + detm2) = entry;
      sc_xci_mtx(2*n_dets_in - 1 + detm2, 2*n_dets_in - 1 + detm1) = entry;
      // cout << "a" << endl;
    }
  }
  cout << cell1 << cell2 << ", " << str << "ci matrix: " << sc_xci_mtx << endl;

//  vector<double> diagonals;
//  for(size_t i=0; i<sc_xci_mtx.get_dims()[0]; i++){
//    diagonals.push_back(sc_xci_mtx(i,i));
//  }
//  sort(diagonals.begin(), diagonals.end());

  // diagonalise xci-hamiltonian
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, std::ignore) = sc_xci_mtx.eigensystem_sym(times);

  // cout << "evals(" << cell1 << cell2 << "): " << evals << endl;
  // cout << "HF levels (sorted): " << diagonals << endl;
  // cout << "correlation in this cell: " << evals[0] - diagonals[0] << endl;
  // // cout << evecs << endl;
  // cout << str << "ci coeffs of first state: " << evecs.get_slice(0, XX) << endl << endl;

//    tensor_2 C_matrix({n_ci_funcs_per_cell, n_dets});
//    tensor_1 energy({n_ci_funcs_per_cell});
//    for(size_t n=0; n<n_ci_funcs_per_cell; n++){ // 1, 2, ...
//      energy[n] = evals[n];
//      const tensor_1 slice(evecs.get_slice(n, XX));
//      for(size_t i=0; i<n_dets; i++){ // 16, ...
//        C_matrix(n, i) = slice[i];
//      }
//    }
//    energies[cell] = energy;
//    C_matrices[cell] = C_matrix;
//    diag_CI_matrices[cell] = sc_xci_mtx;


  auto t1 = Clock::now();
  times.append_item("get C/CI matrices v2", t1 - t0);
  times.ind_down();

  return evals[0];
}



// return a vector with 'cell' entries of tensor_1 with 'n_ci_funcs_per_cell' entries
// the tensor entries are <\Psi^00|H|\Psi^CIx>, where |\Psi^00> is the first HF MO, and |\Psi^CIx> is one of the CI states in a specific cell.
// while the CI state is normalised, the slice of the MO is not

vector<tensor_1> get_offdiag_CI_matrices(const system1d& sys, const state1d& S, size_t n_ci_funcs_per_cell, const string& str, timing& times) {
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);
  assert(S.has_lip_4c2e);
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  vector<tensor_1> offdiag_CI_matrices(sys.n_cells);

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));

  // for each cell, do fci in the field of all other cells, using the local
  // orbitals in each cells.  So, as in the global MO case, the first local
  // function is occupied, the rest is unoccupied.  Orbitals are orthogonal and
  // normalised but really there is a fractional (but constant) number of
  // electrons in a cell.
  // iterate?

  for(size_t cell=0; cell<sys.n_cells; cell++){

    // slice local copy of cell from 1e-H
    tensor_square H_1el_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    for(int i=0; i<sys.n_funcs_per_cell; i++){
      for(int j=0; j<sys.n_funcs_per_cell; j++){
        H_1el_cell(i,j) = H_1el_local(cell*sys.n_funcs_per_cell+i, cell*sys.n_funcs_per_cell+j);
      }
    }

    const bool normalised = true;
    const bool env = false;
    tensor_1 sc_xci_mtx(get_offdiag_single_cell_xci_mtx_v2(sys, S, H_1el_cell, cell, str, env, normalised, times));
    // cout << str << "ci matrix: " << sc_xci_mtx << endl;

    const double coeff = S.local_mo_coeffs[cell*sys.n_funcs_per_cell + 0];
    offdiag_CI_matrices[cell] = sc_xci_mtx * coeff; // |\Psi^00> was treated as normalised, but it is not
  }

  auto t1 = Clock::now();
  times.append_item("get offdiagCI matrices", t1 - t0);
  times.ind_down();

  return offdiag_CI_matrices;
}


#if 0

// two2one: emulate the 2e operator by a 1e operator
// use for getting the local CI functions, as well as the Hamiltonian for the SI matrix
double state1d::get_ci_energy_apprx1(const system1d& sys, const string& str, timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  assert(has_local_mo_coeffs);
  assert(has_local_density);
  assert(has_lip_4c2e);
  const size_t n_ci_funcs_per_cell = 1;
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  // get C matrices which transform "excited local determinants" to "local CI functions", ie the are the CI coeffs of the first n CI functions
  vector<tensor_2> C_matrices;
  vector<tensor_square> diag_CI_matrices; // which are either on the diagonal or not
  cout << "get C(I) matrices ..." << endl;
  {
    const bool two2one = true; // express 2el-op as 1el-op
    const bool env = false;
    tie(C_matrices, diag_CI_matrices) = get_C_CI_matrices(sys, *this, n_ci_funcs_per_cell, str, two2one, env, times);
    // WRONG, the CI matrices are almost certainly wrong; too low values on diagonal
  }
  cout << "get C(I) matrices ... done" << endl;

  tensor_square_sym SI_matrix({n_ci_funcs_per_cell*sys.n_cells, n_ci_funcs_per_cell*sys.n_cells});
  for(size_t cell1=0; cell1<sys.n_cells; cell1++){
    const tensor_2& C_mtx1( C_matrices[cell1] );
    // cout << "C mtx " << cell1 << " : " << C_mtx1 << endl;
    for(size_t cell2=cell1; cell2<sys.n_cells; cell2++){
      const tensor_2& C_mtx2( C_matrices[cell2] );
      // here we need to calc and write n_ci_funcs_per_cell*n_ci_funcs_per_cell values to the SI_matrix

      tensor_square CI_mtx;
      if(cell1==cell2)
        CI_mtx = diag_CI_matrices[cell1];
      else{
        const bool env = false;
        const bool normalised = true;
        CI_mtx = get_offdiag_single_cell_xci_mtx(sys, *this, cell1, cell2, str, env, normalised, times);
      }

      // cout << CI_mtx.get_dims() << endl;
      // cout << C_mtx2.get_dims() << endl;
      // cout << C_mtx1.get_dims() << endl;

      // if(cell1==31 && cell2==32){
      //   cout << CI_mtx << endl;
      //   cout << C_mtx1 << endl;
      // }

      // FIXME check order of cell1 and cell2, fix later, seems not to matter
      tensor_square block( matmul(C_mtx2, matmul(CI_mtx, C_mtx1, 'N', 'N'), 'T', 'N') );

      // cout << cell1 << ", " << cell2 << ": " << block;
      // if(cell1==cell2) cout << " should be less than zero" << endl;
      // else cout << " should be greater than zero" << endl;

      if (cell1==cell2) cout << "cell, cell, block: " << cell1 << ", " << cell2 << ", " << block << endl;

      for(size_t i=0; i<n_ci_funcs_per_cell; i++){
        for(size_t j=0; j<n_ci_funcs_per_cell; j++){
          SI_matrix(cell1*n_ci_funcs_per_cell+i, cell2*n_ci_funcs_per_cell+j) = block(i, j);
          SI_matrix(cell2*n_ci_funcs_per_cell+j, cell1*n_ci_funcs_per_cell+i) = block(i, j); // fill the whole mtx for now
        }
      }

    }
  }
  cout << "SI= " << SI_matrix << endl;

  // diag SI_matrix etc
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = SI_matrix.eigensystem_sym(times);

  cout << "state interaction evals: " << evals << endl;
  cout << "evecs[0] = " << evecs.get_slice(0, XX) << endl;

  tensor_1 lmc_l = local_mo_coeffs.get_slice(0, XX);
  tensor_1 lmc_s({sys.n_cells});
  for(int i=0; i< sys.n_cells; i++)
    lmc_s[i] = lmc_l[i*sys.n_funcs_per_cell];
  cout << "lmcs = " << lmc_s << endl;

  cout << "si . evecs[0] = " << SI_matrix * evecs.get_slice(0, XX) << endl;
  cout << "si . lmcs = " << SI_matrix * lmc_s << endl;

  auto t1 = Clock::now();
  times.append_item("approx1 xci energy", t1 - t0);
  times.ind_down();

  return evals[0];
}
#endif



double hf_energy(const system1d& sys, const state1d& S, timing& times){
  assert(S.has_local_mo_coeffs);
  assert(S.has_local_density);
  assert(S.has_lip_4c2e);

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));
  const tensor_1 mos_l = S.local_mo_coeffs_ro().get_slice(0, XX);
  const double one_el_integral = mos_l * (H_1el_local * mos_l);
  // cout << "1ei " << one_el_integrals << endl;

  const double hf_energy = 2*one_el_integral +
                           get_4c2eintegral_mo_precalc(sys, S, 0, 0, 0, 0, times); // aka rtsu
  // cout << sci_matrix << endl;

  return hf_energy;
}



// two2one: emulate the 2e operator by a 1e operator, but don't use that Hamiltonian afterwards, only the locally obtained CI functions
double state1d::get_ci_energy_apprx2(const system1d& sys, const string& str, timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  assert(has_local_mo_coeffs);
  assert(has_local_density);
  assert(has_lip_4c2e);
  const size_t n_ci_funcs_per_cell = 1;
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  // get C matrices which transform "excited local determinants" to "local CI functions", ie the are the CI coeffs of the first n CI functions
  vector<tensor_2> C_matrices;
  vector<tensor_square> diag_CI_matrices; // which are square and on the diagonal
  vector<tensor_1> offdiag_CI_matrices; // which are tensor_1 and in the first row/column
  cout << "get C(I) matrices ..." << endl;
  {
    bool two2one = true; // express 2el-op as 1el-op
    const bool env = false;
    // for now, the CI expansion takes into account the surrounding field while the final SI entries do not
    tie(C_matrices, std::ignore) = get_C_CI_matrices(sys, *this, n_ci_funcs_per_cell, str, two2one, env, times); // FIXME: check, unite to a single call
    two2one = false;
    tie(std::ignore, diag_CI_matrices) = get_C_CI_matrices(sys, *this, n_ci_funcs_per_cell, str, two2one, env, times);
    offdiag_CI_matrices = get_offdiag_CI_matrices(sys, *this, n_ci_funcs_per_cell, str, times);
  }
  cout << "get C(I) matrices ... done" << endl;

  // one square for each cell, and one entry-row for the HF ground state
  tensor_square_sym SI_matrix({1 + sys.n_cells*n_ci_funcs_per_cell, 1 + sys.n_cells*n_ci_funcs_per_cell});

  // the first entry
  SI_matrix(0, 0) = hf_energy(sys, *this, times);

  // the first row
  for(size_t cell1=0; cell1<sys.n_cells; cell1++){
    const tensor_2& C_mtx1( C_matrices[cell1] );
    // cout << "C mtx " << cell1 << " : " << C_mtx1 << endl;

    tensor_1 CI_mtx = offdiag_CI_matrices[cell1];

    tensor_1 block({C_mtx1.get_size()});
    for(size_t i=0; i<n_ci_funcs_per_cell; i++)
      block[i] = CI_mtx[i] * C_mtx1[i];

    for(size_t i=0; i<n_ci_funcs_per_cell; i++){
      SI_matrix(0, 1 + cell1*n_ci_funcs_per_cell+i) = block(i);
      SI_matrix(1 + cell1*n_ci_funcs_per_cell+i, 0) = block(i); // fill the whole mtx for now
    }

  }

  // only the CI-functions local to cells
  for(size_t cell1=0; cell1<sys.n_cells; cell1++){
    const tensor_2& C_mtx1( C_matrices[cell1] );
    // cout << "C mtx " << cell1 << " : " << C_mtx1 << endl;

    tensor_square CI_mtx;
    CI_mtx = diag_CI_matrices[cell1];

    tensor_square block( matmul(C_mtx1, matmul(CI_mtx, C_mtx1, 'N', 'N'), 'T', 'N') );

    for(size_t i=0; i<n_ci_funcs_per_cell; i++){
      for(size_t j=i; j<n_ci_funcs_per_cell; j++){// since we are on the diagonal
        SI_matrix(1 + cell1*n_ci_funcs_per_cell+i, 1 + cell1*n_ci_funcs_per_cell+j) = block(i, j);
        SI_matrix(1 + cell1*n_ci_funcs_per_cell+j, 1 + cell1*n_ci_funcs_per_cell+i) = block(i, j); // fill the whole mtx for now
      }
    }
  }

  cout << "SI= " << SI_matrix << endl;

  // diag SI_matrix etc
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = SI_matrix.eigensystem_sym(times);

//  cout << "state interaction evals: " << evals << endl;
//  cout << "evecs[0] = " << evecs.get_slice(0, XX) << endl;
//
//  tensor_1 lmc_l = local_mo_coeffs.get_slice(0, XX);
//  tensor_1 lmc_s({sys.n_cells});
//  for(int i=0; i< sys.n_cells; i++)
//    lmc_s[i] = lmc_l[i*sys.n_funcs_per_cell];
//  cout << "lmcs = " << lmc_s << endl;
//
//  cout << "si . evecs[0] = " << SI_matrix * evecs.get_slice(0, XX) << endl;
//  cout << "si . lmcs = " << SI_matrix * lmc_s << endl;

  auto t1 = Clock::now();
  times.append_item("approx1 xci energy", t1 - t0);
  times.ind_down();

  return evals[0];
}



// excite from entire HF MO0 to local orbitals in one cell, without 1e Coulomb ...
double state1d::get_ci_energy_apprx3(const system1d& sys, const string& str, timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  assert(has_local_mo_coeffs);
  assert(has_local_density);
  assert(has_lip_4c2e);
  const size_t n_ci_funcs_per_cell = 1;
  assert(n_ci_funcs_per_cell <= sys.n_funcs_per_cell);

  const tensor_square H_1el_local(this->get_H_1el_local(sys, NABLA2_APPR, times));
  const vector<pair<size_t, size_t>> dets(get_dets(sys, str));
  const vector<pair<size_t, size_t>> dets_mixed(get_dets_mixed(sys));

  // get C matrices which transform "excited local determinants" to "local CI functions", ie the are the CI coeffs of the first n CI functions
  vector<tensor_1> energies;
  vector<tensor_2> C_matrices;
  vector<tensor_square_sym> diag_CI_matrices; // which are square and on the diagonal
  cout << "get C(I) matrices ..." << endl;
  {
    tie(energies, std::ignore, diag_CI_matrices) = get_C_CI_matrices_v2(sys, *this, n_ci_funcs_per_cell, H_1el_local, str, times);
  }
  cout << "get C(I) matrices ... done" << endl;

  // one square for each cell, and one entry-row for the HF ground state
  tensor_square_sym SI_matrix({ sys.n_cells, sys.n_cells });

  // the first row
  for(size_t cell1=0; cell1<sys.n_cells-1; cell1++){
    SI_matrix(cell1, cell1) = energies[cell1][0];
   //  for(size_t cell2=cell1+1; cell2<sys.n_cells; cell2++){

   //    const double two_cell_e = get_exact_twocell_e(sys, *this, H_1el_local, dets, dets_mixed, diag_CI_matrices[cell1], diag_CI_matrices[cell2], cell1, cell2, str, times);
   //    double off_diag_entry = 0;
   //    if( (energies[cell1][0] - two_cell_e)>1.e-12 && (energies[cell2][0] - two_cell_e)>1.e-12 )
   //      off_diag_entry = sqrt((energies[cell1][0] - two_cell_e)*(energies[cell2][0] - two_cell_e));
   //    // cout << energies[cell1][0] << ", " << energies[cell2][0] << ", " << two_cell_e << ", " << off_diag_entry << endl;
   //    SI_matrix(cell1, cell2) = off_diag_entry;
   //    SI_matrix(cell2, cell1) = off_diag_entry;
   //  }
  }
  cout << "SI= " << SI_matrix << endl;

  // diag SI_matrix etc
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = SI_matrix.eigensystem_sym(times);

  cout << "HF (for camparison) = " << hf_energy(sys, *this, times) << endl;

  auto t1 = Clock::now();
  times.append_item("approx1 xci energy", t1 - t0);
  times.ind_down();

  {
    const size_t cell=16;
    const double coeff16 = local_mo_coeffs(0, cell*sys.n_funcs_per_cell+0);
    cout << "coeff16 " << coeff16 << endl;
    const polynomial mo16 = (sys.legendre_polynomials[0] * locbas1d_HF_orthonorm.basis[cell](0,0) +
                            sys.legendre_polynomials[1] * locbas1d_HF_orthonorm.basis[cell](0,1) +
                            sys.legendre_polynomials[2] * locbas1d_HF_orthonorm.basis[cell](0,2) +
                            sys.legendre_polynomials[3] * locbas1d_HF_orthonorm.basis[cell](0,3) +
                            sys.legendre_polynomials[4] * locbas1d_HF_orthonorm.basis[cell](0,4) ) * coeff16;
    const polynomial dmo16 = mo16.derive();
    cout << "mo16 " << mo16 << endl;
    cout << "dmo16 " << dmo16 << endl;
    cout << "dmo16(-0.25) " << dmo16(-0.25) << endl;
    cout << "dmo16(0.25) " << dmo16(0.25) << endl;
  }
  {
    const size_t cell=17;
    const double coeff17 = local_mo_coeffs(0, cell*sys.n_funcs_per_cell+0);
    cout << "coeff17 " << coeff17 << endl;
    const polynomial mo17 = (sys.legendre_polynomials[0] * locbas1d_HF_orthonorm.basis[cell](0,0) +
                            sys.legendre_polynomials[1] * locbas1d_HF_orthonorm.basis[cell](0,1) +
                            sys.legendre_polynomials[2] * locbas1d_HF_orthonorm.basis[cell](0,2) +
                            sys.legendre_polynomials[3] * locbas1d_HF_orthonorm.basis[cell](0,3) +
                            sys.legendre_polynomials[4] * locbas1d_HF_orthonorm.basis[cell](0,4) ) * coeff17;
    const polynomial dmo17 = mo17.derive();
    cout << "mo17 " << mo17 << endl;
    cout << "dmo17 " << dmo17 << endl;
    cout << "dmo17(-0.25) " << dmo17(-0.25) << endl;
    cout << "dmo17(0.25) " << dmo17(0.25) << endl;
  }
  {
    const size_t cell=18;
    const double coeff18 = local_mo_coeffs(0, cell*sys.n_funcs_per_cell+0);
    cout << "coeff18 " << coeff18 << endl;
    const polynomial mo18 = (sys.legendre_polynomials[0] * locbas1d_HF_orthonorm.basis[cell](0,0) +
                            sys.legendre_polynomials[1] * locbas1d_HF_orthonorm.basis[cell](0,1) +
                            sys.legendre_polynomials[2] * locbas1d_HF_orthonorm.basis[cell](0,2) +
                            sys.legendre_polynomials[3] * locbas1d_HF_orthonorm.basis[cell](0,3) +
                            sys.legendre_polynomials[4] * locbas1d_HF_orthonorm.basis[cell](0,4) ) * coeff18;
    const polynomial dmo18 = mo18.derive();
    cout << "mo18 " << mo18 << endl;
    cout << "dmo18 " << dmo18 << endl;
    cout << "dmo18(-0.25) " << dmo18(-0.25) << endl;
    cout << "dmo18(0.25) " << dmo18(0.25) << endl;
  }
  {
    const size_t cell=19;
    const double coeff19 = local_mo_coeffs(0, cell*sys.n_funcs_per_cell+0);
    cout << "coeff19 " << coeff19 << endl;
    const polynomial mo19 = (sys.legendre_polynomials[0] * locbas1d_HF_orthonorm.basis[cell](0,0) +
                            sys.legendre_polynomials[1] * locbas1d_HF_orthonorm.basis[cell](0,1) +
                            sys.legendre_polynomials[2] * locbas1d_HF_orthonorm.basis[cell](0,2) +
                            sys.legendre_polynomials[3] * locbas1d_HF_orthonorm.basis[cell](0,3) +
                            sys.legendre_polynomials[4] * locbas1d_HF_orthonorm.basis[cell](0,4) ) * coeff19;
    const polynomial dmo19 = mo19.derive();
    cout << "mo19 " << mo19 << endl;
    cout << "dmo19 " << dmo19 << endl;
    cout << "dmo19(-0.25) " << dmo19(-0.25) << endl;
    cout << "dmo19(0.25) " << dmo19(0.25) << endl;
  }


  return evals[0];
}

