
#include <vector>
#include <chrono>

#include <omp.h>

#include "aux_math.hh"
#include "integrate.hh"
#include "polynomial.hh"
#include "system_1d.hh"
#include "system_3d.hh"
#include "tensor.hh"

typedef std::chrono::high_resolution_clock Clock;


/////////////////////////////////////////////////////////
// one function integrals
/////////////////////////////////////////////////////////

// double integrate_funcs(const system1d& sys, const size_t c, const size_t f1, operator_t op){
//   double integral = 0.0;
//   const polynomial& p1 = sys.legendre_polynomials[f1];
//   if(op==NONE){
//     integral += sys.grid.w[f];
//   }
//  return integral;
// }

double integrate_cell(const system1d& sys, const tensor_1& coeffs1, operator_t op){
  double integral = 0.0;
  if(op==NONE){
    for(size_t f=0; f<sys.n_funcs_per_cell; f++){
      integral += coeffs1[f] * sys.grid.w[f];
    }
  }
  else{assert(false);}
  return integral;
}


double integrate_system(const system1d& sys, const tensor_1& coeffs1, operator_t op){
  double integral = 0.0;
  if(op==NONE){
    for(int c=0; c<sys.n_cells; c++){
      const size_t first_coeff = c*(sys.n_funcs_per_cell-sys.overlapping);
      const size_t last_coeff = (c+1)*(sys.n_funcs_per_cell-sys.overlapping) - (!sys.overlapping);
      const tensor_1 coeff_frame({last_coeff-first_coeff+1}, coeffs1.begin()+first_coeff, coeffs1.begin()+last_coeff+1);
      integral += integrate_cell(sys, coeff_frame, op);
    }
  }
  else{assert(false);}
  return integral;
}


///////////////////////////////////////////////////////////////
// two function integrals
///////////////////////////////////////////////////////////////

// integrate the product of two LIPs in one cell
double integrate_funcs(const system1d& sys, const size_t c, const size_t f1, const size_t f2, operator_t op){
  // assert(sys.legendre_polynomials_normed.size()>0);
  double integral = 0.0;
  const polynomial& p1 = sys.legendre_polynomials[f1];
  const polynomial& p2 = sys.legendre_polynomials[f2];
  if(op==NONE){
#if 1 // analytical / integrate on grid
    const polynomial p1p2 = p1*p2;
    integral = p1p2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
#else
#endif
  }
  else if(op==NABLA2_APPR){
// Justified only if the function is zero at the border of the outer two cells.
// Then the integrals are wrong in every cell, but the error cancels out.  It
// does not have to be zero in every cell.
    // const polynomial dp1 = p1.derive();
    const polynomial dp1 = sys.d_legendre_polynomials[f1];
    const polynomial dp2 = sys.d_legendre_polynomials[f2];
#if 1 // analytical / integrate on grid
    const polynomial dp1dp2 = dp1*dp2;
    integral = -dp1dp2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
#else
    for(int p=0; p<sys.n_funcs_per_cell; p++){
      const double x = sys.grid.p[p];
      const double w = sys.grid.w[p];
      integral += -dp1(x) * dp2(x) * w;
    }
#endif
  }
  else if(op==NABLA2_REAL){
    const polynomial dp2 = sys.d_legendre_polynomials[f2];
    const polynomial ddp2 = dp2.derive();
#if 1 // analytical / integrate on grid
    const polynomial p1ddp2 = p1*ddp2;
    integral = p1ddp2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
#else
    for(int p=0; p<sys.n_funcs_per_cell; p++){
      const double x = sys.grid.p[p];
      const double w = sys.grid.w[p];
      integral += p1(x) * ddp2(x) * w;
    }
#endif
  }
  else if(op==X){
#if 1 // analytical / integrate on grid
    assert(false); // currently wrong
    const polynomial p1p2 = p1*p2;
    integral = p1p2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
#else
    for(int p=0; p<sys.n_funcs_per_cell; p++){ // sum over grid points
      const double x = sys.grid.p[p];
      const double w = sys.grid.w[p];
      const double gx = sys.grid_begin + (c+0.5)*sys.cell_width + x;
      // assuming w and p1p2 are centred at 0, and an offset potential
      // integral += w * p1p2(x) * gx;
      integral += w * p1(x) * p2(x) * gx;
    }
#endif
  }
  else if(op==POTENTIAL_E){
#if 1 // analytical / integrate on grid
    const polynomial p1p2 = p1*p2;
    for(size_t f=0; f<sys.n_funcs_per_cell; f++){ // we dont have an analytic potential, hence sum over polynomials
      const size_t gi = c*(sys.n_funcs_per_cell-sys.overlapping) + f;
// cout << gi << endl;
      const polynomial& p3 = sys.legendre_polynomials[f];
      const polynomial p1p2p3 = p1p2*p3;
// cout << p1p2p3 << endl;
      integral += sys.pot_coeffs[gi] * p1p2p3.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
    }
#else
#endif
  }
  else{assert(false);}
  return integral;
}


// 
// // integrate the product of two LIPs in one cell
// double integrate_funcs_custom_pot(const system1d& sys, const size_t c, const size_t f1, const size_t f2, const tensor_1 pot, operator_t op){
//   // assert(sys.legendre_polynomials_normed.size()>0);
//   double integral = 0.0;
//   const polynomial& p1 = sys.legendre_polynomials[f1];
//   const polynomial& p2 = sys.legendre_polynomials[f2];
// 
//   if(op==POTENTIAL_E){
// #if 0 // analytical / integrate on grid
// #else
//     for(size_t f=0; f<sys.n_gp_fine_per_cell; f++){ // integrate f1*f2*pot in the fine grid
//       integral += sys.grid_fine.w[f] * p1(sys.grid_fine.p[f]) * p2(sys.grid_fine.p[f]) * pot[f];
//     }
// #endif
//   }
//   else{assert(false);}
//   return integral;
// }
// 


// integrate the product of two functions, expanded in LIPs
double integrate_cell(const system1d& sys, const size_t c, const tensor_1& coeffs1, const tensor_1 &coeffs2, operator_t op){
  double integral = 0.0;
  if(op==NONE || op==NABLA2_REAL || op==NABLA2_APPR || op==X || op==POTENTIAL_E){
    // c not used in NABLA2
    for(size_t f1=0; f1<sys.n_funcs_per_cell; f1++){
      for(size_t f2=0; f2<sys.n_funcs_per_cell; f2++){
        integral += coeffs1[f1] * coeffs2[f2] * integrate_funcs(sys, c, f1, f2, op);
      }
    }
  }
  else{assert(false);}
  return integral;
}


double integrate_system(const system1d& sys, const tensor_1& coeffs1, const tensor_1& coeffs2, operator_t op){
  double integral = 0.0;
  if(op==NONE || op==X || op==POTENTIAL_E){
    for(int c=0; c<sys.n_cells; c++){
      const size_t first_coeff = c*(sys.n_funcs_per_cell-sys.overlapping);
      const size_t last_coeff = (c+1)*(sys.n_funcs_per_cell-sys.overlapping) - (!sys.overlapping);
// cout << c << ": " << first_coeff << ", " << last_coeff << endl;
      const tensor_1 coeff1_frame({last_coeff-first_coeff+1}, coeffs1.get_data_ptr(first_coeff), coeffs1.get_data_ptr(last_coeff+1));
      const tensor_1 coeff2_frame({last_coeff-first_coeff+1}, coeffs2.get_data_ptr(first_coeff), coeffs2.get_data_ptr(last_coeff+1));
// cout << coeff1_frame <<", "<< coeff2_frame << endl;
      integral += integrate_cell(sys, c, coeff1_frame, coeff2_frame, op);
// cout << integral << endl;
    }
  }
  else if(op==NABLA2_REAL || op==NABLA2_APPR){
    tensor_square single_cell({sys.n_funcs_per_cell, sys.n_funcs_per_cell});
    for(int f1=0; f1<sys.n_funcs_per_cell; f1++){
      for(int f2=0; f2<sys.n_funcs_per_cell; f2++){ // single_cell is not symmetric for NABLA2_REAL
        const int cell = 0; // does not matter
        const double cell_integral = integrate_funcs(sys, cell, f1, f2, op);
        single_cell(f1,f2) = cell_integral;
      }
    }
    for(int c=0; c<sys.n_cells; c++){
      const size_t first_coeff = c*(sys.n_funcs_per_cell-sys.overlapping);
      const size_t last_coeff = (c+1)*(sys.n_funcs_per_cell-sys.overlapping) - (!sys.overlapping);
      const tensor_1 coeff1_frame({last_coeff-first_coeff+1}, coeffs1.get_data_ptr(first_coeff), coeffs1.get_data_ptr(last_coeff+1));
      const tensor_1 coeff2_frame({last_coeff-first_coeff+1}, coeffs2.get_data_ptr(first_coeff), coeffs2.get_data_ptr(last_coeff+1));
      for(int i=0; i<sys.n_funcs_per_cell; i++){
        double tmp=0;
        for(int j=0; j<sys.n_funcs_per_cell; j++){
          tmp += coeff1_frame[j] * single_cell(i,j);
        }
        integral += tmp * coeff2_frame[i];
      }
    }
  }
  else{assert(false);}
  return integral;
}


////////////////////////////////////////////////////////
// 3D, one function integrals
////////////////////////////////////////////////////////


double integrate_cell(const system3d& sys, const size_t cz, const size_t cy, const size_t cx, const tensor_3& coeffs1, operator_t op){
  double integral = 0.0;
  if(op==NONE){
    int f=0;
    for(size_t fz=0; fz<sys.n_funcs_per_cell1d; fz++){
    for(size_t fy=0; fy<sys.n_funcs_per_cell1d; fy++){
    for(size_t fx=0; fx<sys.n_funcs_per_cell1d; fx++){
      integral += coeffs1[f++] * sys.grid.w[fz] * sys.grid.w[fy] * sys.grid.w[fx];
    } } }
  }
  else if(op==POTENTIAL_E){
    int f=0;
    for(size_t fz=0; fz<sys.n_funcs_per_cell1d; fz++){
    for(size_t fy=0; fy<sys.n_funcs_per_cell1d; fy++){
    for(size_t fx=0; fx<sys.n_funcs_per_cell1d; fx++){
      integral += coeffs1[f++] * sys.grid.w[fz] * sys.grid.w[fy] * sys.grid.w[fx] * sys.pot_coeffs(cz,cy,cx,fz,fy,fx);
    } } }
  }
  else{assert(false);}
  return integral;
}


double integrate_system(const system3d& sys, const tensor_6& coeffs, operator_t op){
#if defined(_OPENMP)
  const int np = omp_get_num_procs();
  omp_set_num_threads(np);
#endif
  double integral = 0.0;
  if(op==NONE || op==POTENTIAL_E){
#pragma omp parallel for default(none) shared(sys,coeffs,op) reduction(+:integral)
    for(int cz=0; cz<sys.n_cells1d; cz++){
    for(int cy=0; cy<sys.n_cells1d; cy++){
    for(int cx=0; cx<sys.n_cells1d; cx++){
      const int c = index_3(cz, sys.n_cells1d, cy, sys.n_cells1d, cx);
      const int first_coeff = c * sys.n_funcs_per_cell3d;
      const int last_coeff = (c+1) * sys.n_funcs_per_cell3d - 1;
      const tensor_3 coeff_frame({sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d}, coeffs.begin()+first_coeff, coeffs.begin()+last_coeff+1);
      integral += integrate_cell(sys, cz, cy, cx, coeff_frame, op);
    } } }
  }
  else{assert(false);}
  return integral;
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


// integrate a pair of functions in one cell
double integrate_funcs(const system3d& sys, const size_t cz, const size_t cy, const size_t cx,
                       const size_t fz1, const size_t fy1, const size_t fx1,
                       const size_t fz2, const size_t fy2, const size_t fx2, operator_t op){
  double integral = 0.0;
  // const polynomial& pz1 = sys.legendre_polynomials[fz1];
  // const polynomial& py1 = sys.legendre_polynomials[fy1];
  // const polynomial& px1 = sys.legendre_polynomials[fx1];
  // const polynomial& pz2 = sys.legendre_polynomials[fz2];
  // const polynomial& py2 = sys.legendre_polynomials[fy2];
  // const polynomial& px2 = sys.legendre_polynomials[fx2];
  if(op==NABLA2_APPR){
// Justified only if the function is zero at the border of the outer two cells.
// Then the integrals are wrong in every cell, but the error cancels out.  It
// does not have to be zero in every cell.
    const polynomial& dpz1 = sys.d_legendre_polynomials[fz1];
    const polynomial& dpy1 = sys.d_legendre_polynomials[fy1];
    const polynomial& dpx1 = sys.d_legendre_polynomials[fx1];
    const polynomial& dpz2 = sys.d_legendre_polynomials[fz2];
    const polynomial& dpy2 = sys.d_legendre_polynomials[fy2];
    const polynomial& dpx2 = sys.d_legendre_polynomials[fx2];
#if 0 // analytical / integrate on grid
    const polynomial dp1dp2 = dp1*dp2;
    integral = -dp1dp2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
#else
    for(int pz=0; pz<sys.n_funcs_per_cell1d; pz++){
      const double z = sys.grid.p[pz];
      const double wz = sys.grid.w[pz];
      const double z_contrib = dpz1(z) * dpz2(z) * wz;
      for(int py=0; py<sys.n_funcs_per_cell1d; py++){
        const double y = sys.grid.p[py];
        const double wy = sys.grid.w[py];
        const double y_contrib = dpy1(y) * dpy2(y) * wy;
        for(int px=0; px<sys.n_funcs_per_cell1d; px++){
          const double x = sys.grid.p[px];
          const double wx = sys.grid.w[px];
          const double x_contrib = dpx1(x) * dpx2(x) * wx;
          integral += - z_contrib * y_contrib * x_contrib;
    } } }
#endif
  }
  else if(op==NABLA2_REAL){
    assert(false);
//    const polynomial dpx2 = px2.derive();
//    const polynomial ddpx2 = dpx2.derive();
// #if 0 // analytical / integrate on grid
//     const polynomial p1ddp2 = p1*ddp2;
//     integral = p1ddp2.integrate_limits(-sys.cell_width/2.0, sys.cell_width/2.0);
// #else
//     for(int p=0; p<sys.n_funcs_per_cell; p++){
//       const double x = sys.grid.p[p];
//       const double w = sys.grid.w[p];
//       integral += p1(x) * ddp2(x) * w;
//     }
// #endif
  }
  else{assert(false);}
  return integral;
}


double integrate_cell(const system3d& sys, const size_t cz, const size_t cy, const size_t cx, const tensor_3& coeffs1, const tensor_3& coeffs2, operator_t op){
  // auto start = Clock::now();
  double integral = 0.0;
  if(op==NONE){
    for(size_t fz=0; fz<sys.n_funcs_per_cell1d; fz++){
    for(size_t fy=0; fy<sys.n_funcs_per_cell1d; fy++){
    for(size_t fx=0; fx<sys.n_funcs_per_cell1d; fx++){
      const int f = index_3(fz, sys.n_funcs_per_cell1d, fy, sys.n_funcs_per_cell1d, fx);
      integral += coeffs1[f] * coeffs2[f] * sys.grid.w[fz] * sys.grid.w[fy] * sys.grid.w[fx];
    } } }
  }
  else if(op==NABLA2_REAL || op==NABLA2_APPR || op==X){
    // c not used in NABLA2
    int f1 = 0;
    for(size_t f1z=0; f1z<sys.n_funcs_per_cell1d; f1z++){
    for(size_t f1y=0; f1y<sys.n_funcs_per_cell1d; f1y++){
    for(size_t f1x=0; f1x<sys.n_funcs_per_cell1d; f1x++){
      int f2 = 0;
      for(size_t f2z=0; f2z<sys.n_funcs_per_cell1d; f2z++){
      for(size_t f2y=0; f2y<sys.n_funcs_per_cell1d; f2y++){
      for(size_t f2x=0; f2x<sys.n_funcs_per_cell1d; f2x++){
        integral += coeffs1[f1] * coeffs2[f2] * integrate_funcs(sys, cz, cy, cx, f1z, f1y, f1x, f2z, f2y, f2x, op);
        f2++;
      } } }
      f1++;
    } } }
  }
  // assuming two instances of the same function which is non-zero in one grid point
  else if(op==POTENTIAL_E){
    int f1 = 0;
    int gi = index_6(cz, sys.n_cells1d, cy, sys.n_cells1d, cx, sys.n_funcs_per_cell1d, 0, sys.n_funcs_per_cell1d, 0, sys.n_funcs_per_cell1d, 0);
    for(size_t fz=0; fz<sys.n_funcs_per_cell1d; fz++){
    for(size_t fy=0; fy<sys.n_funcs_per_cell1d; fy++){
    for(size_t fx=0; fx<sys.n_funcs_per_cell1d; fx++){
      integral += coeffs1[f1] * coeffs2[f1] * sys.grid.w[fz] * sys.grid.w[fy] * sys.grid.w[fx] * sys.pot_coeffs[gi];
      f1++;
      gi++;
    } } }
  }
  else{assert(false);}
  // auto stop = Clock::now();
  // cout << "integrate cell: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
  return integral;
}


double integrate_system(const system3d& sys, const tensor_6& coeffs1, const tensor_6& coeffs2, operator_t op){
  // auto start = Clock::now();
#if defined(_OPENMP)
  const int np = omp_get_num_procs();
  omp_set_num_threads(np);
#endif
  double integral = 0.0;
  if(op==NONE || op==POTENTIAL_E){
#pragma omp parallel for default(none) shared(sys,coeffs1,coeffs2,op) reduction(+:integral)
    for(int cz=0; cz<sys.n_cells1d; cz++){
    for(int cy=0; cy<sys.n_cells1d; cy++){
    for(int cx=0; cx<sys.n_cells1d; cx++){
      const int c = index_3(cz, sys.n_cells1d, cy, sys.n_cells1d, cx);
      const int first_coeff = c * sys.n_funcs_per_cell3d;
      const int last_coeff = (c+1) * sys.n_funcs_per_cell3d - 1;
      const tensor_3 coeff1_frame({sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d}, coeffs1.begin()+first_coeff, coeffs1.begin()+last_coeff+1);
      const tensor_3 coeff2_frame({sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d}, coeffs2.begin()+first_coeff, coeffs2.begin()+last_coeff+1);
      integral += integrate_cell(sys, cz, cy, cx, coeff1_frame, coeff2_frame, op);
    } } }
  }
  else if(op==NABLA2_REAL || op==NABLA2_APPR){
    tensor_square single_cell({sys.n_funcs_per_cell3d, sys.n_funcs_per_cell3d});
#pragma omp parallel for default(none) shared(sys,op,single_cell)
    for(int fz1=0; fz1<sys.n_funcs_per_cell1d; fz1++){
    for(int fy1=0; fy1<sys.n_funcs_per_cell1d; fy1++){
    for(int fx1=0; fx1<sys.n_funcs_per_cell1d; fx1++){
      const int f1 = index_3(fz1, sys.n_funcs_per_cell1d, fy1, sys.n_funcs_per_cell1d, fx1);
      int f2 = 0;
      for(int fz2=0; fz2<sys.n_funcs_per_cell1d; fz2++){ // single_cell is not symmetric for NABLA2_REAL
      for(int fy2=0; fy2<sys.n_funcs_per_cell1d; fy2++){ // single_cell is not symmetric for NABLA2_REAL
      for(int fx2=0; fx2<sys.n_funcs_per_cell1d; fx2++){ // single_cell is not symmetric for NABLA2_REAL
        const int cell = 0; // does not matter
        const double cell_integral = integrate_funcs(sys, cell, cell, cell, fz1, fy1, fx1, fz2, fy2, fx2, op);
        // const int f2 = index_3(fz2, sys.n_funcs_per_cell1d, fy2, sys.n_funcs_per_cell1d, fx2);
        single_cell(f1,f2) = cell_integral;
        f2++;
      } } }
    } } }
    for(int c=0; c<sys.n_cells3d; c++){
      const int first_coeff = c * sys.n_funcs_per_cell3d;
      // const int last_coeff = (c+1) * sys.n_funcs_per_cell3d - 1;
      // const vector<double> coeff1_frame(coeffs1.begin()+first_coeff, coeffs1.begin()+last_coeff+1);
      // const vector<double> coeff2_frame(coeffs2.begin()+first_coeff, coeffs2.begin()+last_coeff+1);
      for(int i=0; i<sys.n_funcs_per_cell3d; i++){
        double tmp=0;
        for(int j=0; j<sys.n_funcs_per_cell3d; j++){
          // tmp += coeff1_frame[j] * single_cell(i,j);
          tmp += coeffs1[first_coeff+j] * single_cell(i,j);
        }
        //integral += tmp * coeff2_frame[i];
        integral += tmp * coeffs2[first_coeff+i];
      }
    }
  }
  else{assert(false);}
  // auto stop = Clock::now();
  // cout << "integrate system("<<op<<"): " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
  return integral;
}


