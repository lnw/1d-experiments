#include <vector>
#include <tuple>
#include <algorithm>

#include "timing.hh"
#include "tensor.hh"
#include "local_basis_1d.hh"
#include "system_1d.hh"
#include "state_1d.hh"

using namespace std;


tensor_square_sym get_sci_matrix(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  assert(S.has_local_mo_coeffs);

  const size_t n_occ = 1;
  const size_t n_virt = 10; //n*(m-1);
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_dets = n_occ + n_singly_exc_dets;
  cout << "number of distinct singly excited determinents: " << n_singly_exc_dets << endl;
  cout << "total number of distinct determinents: " << n_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store SCI matrix" << endl;

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));

  cout << "getting 1-electron integrals ..." << endl;
  tensor_square one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  vector<tensor_1> mos_l({n_occ+n_virt});
  for(int i=0; i<n_occ+n_virt; i++){
    mos_l[i] = S.local_mo_coeffs_ro().get_slice(i, XX);
  }
  for(int i=0; i<n_occ+n_virt; i++){
    for(int j=0; j<n_occ+n_virt; j++){
      one_el_integrals(i,j) = mos_l[i] * (H_1el_local * mos_l[j]);
    }
  }
  // cout << "1ei " << one_el_integrals << endl;
  cout << "done" << endl;

  // order of determinants
  vector<pair<size_t, size_t>> dets;
  for(size_t r=0; r<n_occ+n_virt; r++){
    dets.push_back(make_pair(0,r));
  }
  cout << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym sci_matrix({n_dets, n_dets});
  cout << "sci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
#pragma omp parallel for default(none) shared(sci_matrix, one_el_integrals, dets, times)
  for(size_t det1=0; det1<n_dets; det1++){
    for(size_t det2=det1; det2<n_dets; det2++){
      const size_t r=dets[det1].first,
                   s=dets[det1].second,
                   t=dets[det2].first,
                   u=dets[det2].second;
      // if( (r+s+t+u)%2 ) continue; // (tu|rs) is zero for odd parity
      // cout << "rstu: " << r << ", " << s << ", " << t << ", " << u << endl;
      if(det1==det2 && r==s){
        const double entry = 2*one_el_integrals(r,r) +
                             get_4c2eintegral_mo_precalc(sys, S, r, r, r, r, times); // aka rtsu
        // cout << entry << endl;
        sci_matrix(det1, det2) = entry;
      }
      else if(det1==det2 && r!=s){
        const double entry = one_el_integrals(r,r) +
                             one_el_integrals(s,s) +
                             get_4c2eintegral_mo_precalc(sys, S, r, r, s, s, times) +
                             get_4c2eintegral_mo_precalc(sys, S, r, s, r, s, times);
        // cout << entry << endl;
        sci_matrix(det1, det2) = entry;
      }
      else if( (r==s && t!=u) || (r!=s && t==u) ){ // tt rs  || tu rr
        const double entry = sqrt(2) * one_el_integrals(s,u) +
                             sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, s, u, times);
        // cout << entry << endl;
        sci_matrix(det1, det2) = entry;
      }
      else /*if( r!=s && t!=u )*/ { // tt rs  || tu rr
        const double entry = one_el_integrals(s,u) +
                             get_4c2eintegral_mo_precalc(sys, S, r, t, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, t, s, u, r, times);
        // cout << entry << endl;
        sci_matrix(det1, det2) = entry;
      }
    } // det2
  } // det1

  // cout << sci_matrix << endl;
  auto t1 = Clock::now();
  times.append_item("sci matrix", t1 - t0);
  times.ind_down();

  return sci_matrix;
}


tensor_square_sym get_dci_matrix(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = 10; //n*(m-1);
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  const size_t n_dets = n_occ + n_doubly_exc_dets;
  cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  // cout << "will need " << n_dets*n_dets *8.0/1000000.0 << " MB to store DCI matrix" << endl;

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));

  cout << "getting 1 electron integrals ..." << endl;
  tensor_square_sym one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  vector<tensor_1> mos({n_occ+n_virt});
  for(int i=0; i<n_occ+n_virt; i++){
    mos[i] = S.local_mo_coeffs_ro().get_slice(i, XX);
  }
  for(int i=0; i<n_occ+n_virt; i++){
    for(int j=0; j<n_occ+n_virt; j++){
      one_el_integrals(i,j) = mos[i] * (H_1el_local * mos[j]);
    }
  }
  // cout << "1ei " << one_el_integrals << endl;
  cout << "done" << endl;

  // order of determinants
  vector<pair<size_t, size_t>> dets;
  for(size_t r=0; r<n_occ+n_virt; r++){
    for(size_t s=r; s<n_occ+n_virt; s++){
      if(r==0 && s!=0) continue;
      dets.push_back(make_pair(r,s));
    }
  }
  // const size_t n_dets=dets.size();
  cout << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym dci_matrix({n_doubly_exc_dets+1, n_doubly_exc_dets+1});
  cout << "dci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
#pragma omp parallel for default(none) shared(dci_matrix, one_el_integrals, dets, times)
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
                             get_4c2eintegral_mo_precalc(sys, S, r, r, r, r, times); // aka rtsu
        dci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if(det1==det2 && r!=s){ // b
        const double entry = one_el_integrals(r,r) +
                             one_el_integrals(s,s) +
                             get_4c2eintegral_mo_precalc(sys, S, r, r, s, s, times) +
                             get_4c2eintegral_mo_precalc(sys, S, r, s, r, s, times);
        dci_matrix(det1, det2) = entry;
        // cout << "b" << endl;
      }
      else if( r==s && t==u ){ // d
        const double entry = get_4c2eintegral_mo_precalc(sys, S, r, t, r, t, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "d" << endl;
      }
      else if( (r==s && r==t)){ // f1
        const double entry = sqrt(2) * one_el_integrals(s,u) +
                             sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, r, r, u, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "f1" << endl;
      }
      // else if( (t==u && t==r)){ // f2
      //   const double entry = sqrt(2) * one_el_integrals(r,s) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, r, s, r, times);
      //   dci_matrix(det1, det2) = 5;//entry;
      //   cout << "f2" << endl;
      // }
      else if( (s==u && s==t)){ // f3
        const double entry = sqrt(2) * one_el_integrals(r,t) +
                             sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, t, t, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "f3" << endl;
      }
      // else if( (s==u && s==r)){ // f4
      //   const double entry = sqrt(2) * one_el_integrals(r,t) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, r, r, times);
      //   dci_matrix(det1, det2) = 7;//entry;
      //   cout << "f4" << endl;
      // }
      else if( r==t && s!=u ){ // c1
        const double entry = one_el_integrals(s,u) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, r, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, r, s, u, r, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "c1" << endl;
      }
      else if( s==t && r!=u ){ // c2
        const double entry = one_el_integrals(r,u) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, s, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, s, s, u, r, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "c2" << endl;
      }
      else if( r!=t && s==u ){ // c3
        const double entry = one_el_integrals(r,t) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, t, s, s, times) +
                             get_4c2eintegral_mo_precalc(sys, S, t, s, s, r, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "c3" << endl;
      }
      // else if( (s!=t && r==u) ){ // c4
      //   const double entry = one_el_integrals(s,t) + 
      //                        get_4c2eintegral_mo_precalc(sys, S, r, t, s, r, times) +
      //                        get_4c2eintegral_mo_precalc(sys, S, t, s, r, r, times);
      //   dci_matrix(det1, det2) = 11;//entry;
      //   cout << "c4" << endl;
      // }
      else if( r==s ){ // g
        const double entry = sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, r, u, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "g1" << endl;
      }
      else if( t==u ){ // g
        const double entry = sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, s, t, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "g2" << endl;
      }
      else if( r!=s && r!=t && r!=u && t!=u ){ // e
        const double entry = get_4c2eintegral_mo_precalc(sys, S, r, t, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, t, s, u, r, times);
        dci_matrix(det1, det2) = entry;
        dci_matrix(det2, det1) = entry;
        // cout << "e" << endl;
      }
      else
        assert(false);
    } // det2
  } // det1

  auto t1 = Clock::now();
  times.append_item("dci matrix", t1 - t0);
  times.ind_down();

  return dci_matrix;
}

#if 1

tensor_square_sym get_sdci_matrix(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = 10; //n*(m-1); // 20 is no better
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  const size_t n_dets = 1 + n_singly_exc_dets + n_doubly_exc_dets;
  cout << "number of distinct singly excited determinents: " << n_singly_exc_dets << endl;
  cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  cout << "total number of distinct determinents: " << n_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store SDCI matrix" << endl;

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));

  cout << "getting 1-electron integrals ..." << endl;
  tensor_square one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  vector<tensor_1> mos({n_occ+n_virt});
  for(int i=0; i<n_occ+n_virt; i++){
    mos[i] = S.local_mo_coeffs_ro().get_slice(i, XX);
  }
  for(int i=0; i<n_occ+n_virt; i++){
    for(int j=0; j<n_occ+n_virt; j++){
      one_el_integrals(i,j) = mos[i] * (H_1el_local * mos[j]);
    }
  }
  // cout << "1ei " << one_el_integrals << endl;
  cout << "done" << endl;

  // order of determinants
  vector<pair<size_t, size_t>> dets;
  for(size_t r=0; r<n_occ+n_virt; r++){
    for(size_t s=r; s<n_occ+n_virt; s++){
      dets.push_back(make_pair(r,s));
    }
  }
  cout << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym sdci_matrix({n_dets, n_dets});
  cout << "sdci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
#pragma omp parallel for default(none) shared(sdci_matrix, one_el_integrals, dets, times)
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
                             get_4c2eintegral_mo_precalc(sys, S, r, r, r, r, times); // aka rtsu
        sdci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if(det1==det2 && r!=s){ // b
        const double entry = one_el_integrals(r,r) +
                             one_el_integrals(s,s) +
                             get_4c2eintegral_mo_precalc(sys, S, r, r, s, s, times) +
                             get_4c2eintegral_mo_precalc(sys, S, r, s, r, s, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "b" << endl;
      }
      else if( r==s && t==u ){ // d
        const double entry = get_4c2eintegral_mo_precalc(sys, S, r, t, r, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "d" << endl;
      }
      else if( (r==s && r==t)){ // f1
        const double entry = sqrt(2) * one_el_integrals(s,u) +
                             sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, r, r, u, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "f1" << endl;
      }
      // else if( (t==u && t==r)){ // f2
      //   const double entry = sqrt(2) * one_el_integrals(r,s) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, r, s, r, times);
      //   sdci_matrix(det1, det2) = 5;//entry;
      //   cout << "f2" << endl;
      // }
      else if( (s==u && s==t)){ // f3
        const double entry = sqrt(2) * one_el_integrals(r,t) +
                             sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, t, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "f3" << endl;
      }
      // else if( (s==u && s==r)){ // f4
      //   const double entry = sqrt(2) * one_el_integrals(r,t) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, r, r, times);
      //   sdci_matrix(det1, det2) = 7;//entry;
      //   cout << "f4" << endl;
      // }
      else if( (r==t && s!=u) ){ // c1
        const double entry = one_el_integrals(s,u) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, r, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, r, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c1" << endl;
      }
      else if( (s==t && r!=u) ){ // c2
        const double entry = one_el_integrals(r,u) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, s, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, s, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c2" << endl;
      }
      else if( (r!=t && s==u) ){ // c3
        const double entry = one_el_integrals(r,t) + 
                             get_4c2eintegral_mo_precalc(sys, S, r, t, s, s, times) +
                             get_4c2eintegral_mo_precalc(sys, S, t, s, s, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c3" << endl;
      }
      // else if( (s!=t && r==u) ){ // c4
      //   const double entry = one_el_integrals(s,t) + 
      //                        get_4c2eintegral_mo_precalc(sys, S, r, t, s, r, times) +
      //                        get_4c2eintegral_mo_precalc(sys, S, t, s, r, r, times);
      //   sdci_matrix(det1, det2) = 11;//entry;
      //   cout << "c4" << endl;
      // }
      else if( r==s ){ // g
        const double entry = sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, r, u, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "g1" << endl;
      }
      else if( t==u ){ // g
        const double entry = sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, s, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "g2" << endl;
      }
      else if( r!=s && r!=t && r!=u && t!=u ){ // e
        const double entry = get_4c2eintegral_mo_precalc(sys, S, r, t, s, u, times) +
                             get_4c2eintegral_mo_precalc(sys, S, t, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "e" << endl;
      }
    } // det2
  } // det1

  // cout << sdci_matrix << endl;
  auto t1 = Clock::now();
  times.append_item("sdci matrix", t1 - t0);
  times.ind_down();

  return sdci_matrix;
}

#else

tensor_square_sym get_sdci_matrix(const system1d& sys, const state1d& S, timing& times){
  times.ind_up();
  auto t0 = Clock::now();

  const size_t n_occ = 1;
  const size_t n_virt = 2;
  const size_t n_singly_exc_dets = n_virt;
  const size_t n_doubly_exc_dets = (n_virt*(n_virt+1))/2;
  const size_t n_dets = 1 + n_singly_exc_dets + n_doubly_exc_dets;
  cout << "number of distinct singly excited determinents: " << n_singly_exc_dets << endl;
  cout << "number of distinct doubly excited determinents: " << n_doubly_exc_dets << endl;
  cout << "total number of distinct determinents: " << n_dets << endl;
  // cout << "will need " << n_dets*n_dets*8.0/1000000.0 << " MB to store SDCI matrix" << endl;

  const tensor_square H_1el_local(S.get_H_1el_local(sys, NABLA2_APPR, times));

  cout << "getting 1-electron integrals ..." << endl;
  tensor_square one_el_integrals({n_occ+n_virt, n_occ+n_virt});
  one_el_integrals(0,0) = -3.75206 + 1.83943;
  one_el_integrals(0,1) = -1.80427 + 1.6335;
  one_el_integrals(0,2) =  1.56579 + -1.51649;
  one_el_integrals(1,0) = -1.80427 + 1.6335;
  one_el_integrals(1,1) = -3.29297 + 4.96518;
  one_el_integrals(1,2) =  2.01592 + -1.87675;
  one_el_integrals(2,0) =  1.56579 + -1.51649;
  one_el_integrals(2,1) =  2.01592 + -1.87675;
  one_el_integrals(2,2) = -5.31424 + 14.4726;
  cout << "done" << endl;

  tensor_4 two_el_integrals({n_occ+n_virt, n_occ+n_virt, n_occ+n_virt, n_occ+n_virt});
  two_el_integrals(0,0,0,0) = 1.22506    ;
  two_el_integrals(0,0,0,1) = 0.318769   ;
  two_el_integrals(0,0,0,2) = -0.0862558 ;
  two_el_integrals(0,0,1,0) = 0.318769   ;
  two_el_integrals(0,0,1,1) = 1.06428    ;
  two_el_integrals(0,0,1,2) = -0.259811  ;
  two_el_integrals(0,0,2,0) = -0.0862558 ;
  two_el_integrals(0,0,2,1) = -0.259811  ;
  two_el_integrals(0,0,2,2) = 1.28354    ;
  two_el_integrals(0,1,0,0) = 0.318769   ;
  two_el_integrals(0,1,0,1) = 0.22808    ;
  two_el_integrals(0,1,0,2) = -0.113595  ;
  two_el_integrals(0,1,1,0) = 0.22808    ;
  two_el_integrals(0,1,1,1) = 0.241054   ;
  two_el_integrals(0,1,1,2) = -0.1814    ;
  two_el_integrals(0,1,2,0) = -0.113595  ;
  two_el_integrals(0,1,2,1) = -0.1814    ;
  two_el_integrals(0,1,2,2) = 0.387033   ;
  two_el_integrals(0,2,0,0) = -0.0862558 ;
  two_el_integrals(0,2,0,1) = -0.113595  ;
  two_el_integrals(0,2,0,2) = 0.144492   ;
  two_el_integrals(0,2,1,0) = -0.113595  ;
  two_el_integrals(0,2,1,1) = -0.0921472 ;
  two_el_integrals(0,2,1,2) = 0.127557   ;
  two_el_integrals(0,2,2,0) = 0.144492   ;
  two_el_integrals(0,2,2,1) = 0.127557   ;
  two_el_integrals(0,2,2,2) = -0.177838  ;
  two_el_integrals(1,0,0,0) = 0.318769   ;
  two_el_integrals(1,0,0,1) = 0.22808    ;
  two_el_integrals(1,0,0,2) = -0.113595  ;
  two_el_integrals(1,0,1,0) = 0.22808    ;
  two_el_integrals(1,0,1,1) = 0.241054   ;
  two_el_integrals(1,0,1,2) = -0.1814    ;
  two_el_integrals(1,0,2,0) = -0.113595  ;
  two_el_integrals(1,0,2,1) = -0.1814    ;
  two_el_integrals(1,0,2,2) = 0.387033   ;
  two_el_integrals(1,1,0,0) = 1.06428    ;
  two_el_integrals(1,1,0,1) = 0.241054   ;
  two_el_integrals(1,1,0,2) = -0.0921472 ;
  two_el_integrals(1,1,1,0) = 0.241054   ;
  two_el_integrals(1,1,1,1) = 0.956731   ;
  two_el_integrals(1,1,1,2) = -0.208711  ;
  two_el_integrals(1,1,2,0) = -0.0921472 ;
  two_el_integrals(1,1,2,1) = -0.208711  ;
  two_el_integrals(1,1,2,2) = 1.11977    ;
  two_el_integrals(1,2,0,0) = -0.259811  ;
  two_el_integrals(1,2,0,1) = -0.1814    ;
  two_el_integrals(1,2,0,2) = 0.127557   ;
  two_el_integrals(1,2,1,0) = -0.1814    ;
  two_el_integrals(1,2,1,1) = -0.208711  ;
  two_el_integrals(1,2,1,2) = 0.177322   ;
  two_el_integrals(1,2,2,0) = 0.127557   ;
  two_el_integrals(1,2,2,1) = 0.177322   ;
  two_el_integrals(1,2,2,2) = -0.352371  ;
  two_el_integrals(2,0,0,0) = -0.0862558 ;
  two_el_integrals(2,0,0,1) = -0.113595  ;
  two_el_integrals(2,0,0,2) = 0.144492   ;
  two_el_integrals(2,0,1,0) = -0.113595  ;
  two_el_integrals(2,0,1,1) = -0.0921472 ;
  two_el_integrals(2,0,1,2) = 0.127557   ;
  two_el_integrals(2,0,2,0) = 0.144492   ;
  two_el_integrals(2,0,2,1) = 0.127557   ;
  two_el_integrals(2,0,2,2) = -0.177838  ;
  two_el_integrals(2,1,0,0) = -0.259811  ;
  two_el_integrals(2,1,0,1) = -0.1814    ;
  two_el_integrals(2,1,0,2) = 0.127557   ;
  two_el_integrals(2,1,1,0) = -0.1814    ;
  two_el_integrals(2,1,1,1) = -0.208711  ;
  two_el_integrals(2,1,1,2) = 0.177322   ;
  two_el_integrals(2,1,2,0) = 0.127557   ;
  two_el_integrals(2,1,2,1) = 0.177322   ;
  two_el_integrals(2,1,2,2) = -0.352371  ;
  two_el_integrals(2,2,0,0) = 1.28354    ;
  two_el_integrals(2,2,0,1) = 0.387033   ;
  two_el_integrals(2,2,0,2) = -0.177838  ;
  two_el_integrals(2,2,1,0) = 0.387033   ;
  two_el_integrals(2,2,1,1) = 1.11977    ;
  two_el_integrals(2,2,1,2) = -0.352371  ;
  two_el_integrals(2,2,2,0) = -0.177838  ;
  two_el_integrals(2,2,2,1) = -0.352371  ;
  two_el_integrals(2,2,2,2) = 1.43178    ;


  // order of determinants
  vector<pair<size_t, size_t>> dets;
  for(size_t r=0; r<n_occ+n_virt; r++){
    for(size_t s=r; s<n_occ+n_virt; s++){
      dets.push_back(make_pair(r,s));
    }
  }
  cout << dets << endl;

  // shape of the matrix:
  // psi_0, psi_00^11, ..., psi_00^1n, psi_00^22, psi_00^23, ..., psi_00^2n, ..., ..., psi_00^nn-1, psi_00^nn
  tensor_square_sym sdci_matrix({n_dets, n_dets});
  cout << "sdci matrix allocated" << endl;

  // get matrix elements <\psi_aa^tu|H|\psi_aa^rs>
#pragma omp parallel for default(none) shared(sdci_matrix, one_el_integrals, dets, times)
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
                             two_el_integrals(r,r,r,r); 
                             //get_4c2eintegral_mo_precalc(sys, S, r, r, r, r, times); // aka rtsu
        sdci_matrix(det1, det2) = entry;
        // cout << "a" << endl;
      }
      else if(det1==det2 && r!=s){ // b
        const double entry = one_el_integrals(r,r) +
                             one_el_integrals(s,s) +
                             two_el_integrals(r,r,s,s) +
                             two_el_integrals(r,s,r,s); 
                             // get_4c2eintegral_mo_precalc(sys, S, r, r, s, s, times) +
                             // get_4c2eintegral_mo_precalc(sys, S, r, s, r, s, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "b" << endl;
      }
      else if( r==s && t==u ){ // d
        const double entry = two_el_integrals(r,t,r,t);//get_4c2eintegral_mo_precalc(sys, S, r, t, r, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "d" << endl;
      }
      else if( (r==s && r==t)){ // f1
        const double entry = sqrt(2) * one_el_integrals(s,u) +
                             sqrt(2) * two_el_integrals(r,r,r,u); // get_4c2eintegral_mo_precalc(sys, S, r, r, r, u, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "f1" << endl;
      }
      // else if( (t==u && t==r)){ // f2
      //   const double entry = sqrt(2) * one_el_integrals(r,s) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, r, s, r, times);
      //   sdci_matrix(det1, det2) = 5;//entry;
      //   cout << "f2" << endl;
      // }
      else if( (s==u && s==t)){ // f3
        const double entry = sqrt(2) * one_el_integrals(r,t) +
                             sqrt(2) * two_el_integrals(r,t,t,t); // get_4c2eintegral_mo_precalc(sys, S, r, t, t, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "f3" << endl;
      }
      // else if( (s==u && s==r)){ // f4
      //   const double entry = sqrt(2) * one_el_integrals(r,t) +
      //                        sqrt(2) * get_4c2eintegral_mo_precalc(sys, S, r, t, r, r, times);
      //   sdci_matrix(det1, det2) = 7;//entry;
      //   cout << "f4" << endl;
      // }
      else if( (r==t && s!=u) ){ // c1
        const double entry = one_el_integrals(s,u) + 
                             two_el_integrals(r,r,s,u) + // get_4c2eintegral_mo_precalc(sys, S, r, r, s, u, times) +
                             two_el_integrals(r,s,u,r); // get_4c2eintegral_mo_precalc(sys, S, r, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c1" << endl;
      }
      else if( (s==t && r!=u) ){ // c2
        const double entry = one_el_integrals(r,u) + 
                             two_el_integrals(r,s,s,u) + //get_4c2eintegral_mo_precalc(sys, S, r, s, s, u, times) +
                             two_el_integrals(s,s,u,r); //get_4c2eintegral_mo_precalc(sys, S, s, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c2" << endl;
      }
      else if( (r!=t && s==u) ){ // c3
        const double entry = one_el_integrals(r,t) + 
                             two_el_integrals(r,t,s,s)+ //get_4c2eintegral_mo_precalc(sys, S, r, t, s, s, times) +
                             two_el_integrals(t,s,s,r); //get_4c2eintegral_mo_precalc(sys, S, t, s, s, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "c3" << endl;
      }
      // else if( (s!=t && r==u) ){ // c4
      //   const double entry = one_el_integrals(s,t) + 
      //                        get_4c2eintegral_mo_precalc(sys, S, r, t, s, r, times) +
      //                        get_4c2eintegral_mo_precalc(sys, S, t, s, r, r, times);
      //   sdci_matrix(det1, det2) = 11;//entry;
      //   cout << "c4" << endl;
      // }
      else if( r==s ){ // g
        const double entry = sqrt(2) * two_el_integrals(r,t,r,u); // get_4c2eintegral_mo_precalc(sys, S, r, t, r, u, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "g1" << endl;
      }
      else if( t==u ){ // g
        const double entry = sqrt(2) * two_el_integrals(r,t,s,t); // get_4c2eintegral_mo_precalc(sys, S, r, t, s, t, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "g2" << endl;
      }
      else if( r!=s && r!=t && r!=u && t!=u ){ // e
        const double entry = two_el_integrals(r,t,s,u) + //get_4c2eintegral_mo_precalc(sys, S, r, t, s, u, times) +
                             two_el_integrals(t,s,u,r); //get_4c2eintegral_mo_precalc(sys, S, t, s, u, r, times);
        sdci_matrix(det1, det2) = entry;
        // cout << "e" << endl;
      }
    } // det2
  } // det1

  // cout << sdci_matrix << endl;
  auto t1 = Clock::now();
  times.append_item("sdci matrix", t1 - t0);
  times.ind_down();

  return sdci_matrix;
}

#endif

double state1d::get_xci_energy(const system1d& sys, const string& str, timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  assert(has_lip_4c2e);

  tensor_square_sym xci_matrix;
  if(str=="s")
    xci_matrix = get_sci_matrix(sys, *this, times);
  else if(str=="d")
    xci_matrix = get_dci_matrix(sys, *this, times);
  else if(str=="sd")
    xci_matrix = get_sdci_matrix(sys, *this, times);
  else
    assert(false);
  cout << str << "mat = " << xci_matrix << endl;

  // diagonalise xci-hamiltonian
  tensor_1 evals;
  tensor_square evecs;
  tie(evals, evecs) = xci_matrix.eigensystem_sym(times);

  vector<double> diagonals;
  for(int i=0; i<xci_matrix.get_dims()[0]; i++){
    diagonals.push_back(xci_matrix(i,i));
  }
  sort(diagonals.begin(), diagonals.end());

  cout << "HF levels (sorted): " << diagonals << endl;
  cout << "correlated levels: " << evals << endl;
  // cout << evecs << endl;
  cout << "xci coeffs of first state: " << evecs.get_slice(0, XX) << endl;

  const double xci_energy = evals[0];
  cout << "HF energy: " << diagonals[0] << endl;
  cout << "xci energy: " << xci_energy << endl;
  cout << "correlation energy: " << xci_energy - diagonals[0] << endl;

  auto t1 = Clock::now();
  times.append_item("xci energy", t1 - t0);
  times.ind_down();

  return xci_energy;
}


