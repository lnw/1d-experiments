
#include <algorithm>
#include <cmath>
#include <chrono>
#include <vector>

#include <omp.h>

#include "aux_math.hh"
#include "lina.hh"
#include "integrate.hh"
#include "mo_1d.hh"
#include "mo_3d.hh"
#include "system_1d.hh"
#include "tensor.hh"

using namespace std;
typedef std::chrono::high_resolution_clock Clock;


#pragma omp declare reduction(tp: tensor_1: omp_out += omp_in) //initializer( omp_priv=vector<double>() )


void mo1d::guess(const system1d& sys, const localbasis1d& lo1d){
  auto start = Clock::now();
  orbital_e = -1.0;
  // lets put some uniform occupation in the middle third or the whole vector if it's short
  if(dim < 10){
    fill(orbital_coeffs.begin(), orbital_coeffs.end(), 1.0);
  }
  else{
    const int middle_begin = dim/5.0, middle_end = dim - middle_begin;
    fill(orbital_coeffs.begin()+middle_begin, orbital_coeffs.begin()+middle_end, 1.0);
  }
  normalise(sys, lo1d);
  auto stop = Clock::now();
  cout << "guess: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
}


#if 0

void mo1d::do_helmholtz_scf_step(const system1d& sys){
  auto start = Clock::now();
cout << "--- begin helmholtz ---" << endl;
  const double k = sqrt(-2.0*orbital_e);
  const double eps = 1.e-10;
  const double t_f = 1000.0; // that's the upper boundary of the log grid
  tensor_1 orbital_x_potential({sys.n_cells, sys.n_funcs_per_cell});
  for (int i=0; i<dim; i++) orbital_x_potential[i] = orbital_coeffs[i] * sys.pot_coeffs[i];
  tensor_1 coeffs_new({sys.n_cells, sys.n_funcs_per_cell});
  for(int c_r=0; c_r<sys.n_cells; c_r++){ // cells for r
    for(int p_r=0; p_r<sys.n_funcs_per_cell; p_r++){ // points as far as r is concerned
      const int index_r = c_r*(sys.n_funcs_per_cell-sys.overlapping) + p_r;
      double coeff_new = 0.0;
      tensor_1 t_integrals({sys.n_cells, sys.n_funcs_per_cell});
      for(int c_rp=0; c_rp<sys.n_cells; c_rp++){ // cells for r prime
        for(int p_rp=0; p_rp<sys.n_funcs_per_cell; p_rp++){ // points for r prime
          double t_integral = 0.0;
          const int index_rp = c_rp*(sys.n_funcs_per_cell-sys.overlapping) + p_rp;
// cout << "index_rp: " << index_rp << endl;
          const double delta_r = fabs((c_r - c_rp)*sys.cell_width + (sys.grid.p[p_r] - sys.grid.p[p_rp]));
          for(size_t t=0; t<sys.tgrid.size(); t++){
            t_integral += sys.tgrid.w[t] * exp(-k*k/(4.0*sys.tgrid.p[t]*sys.tgrid.p[t])) * exp(-sys.tgrid.p[t]*sys.tgrid.p[t]*delta_r*delta_r);
// cout <<  sys.tgrid.w[t] <<", "<< exp(-k*k/(4*sys.tgrid.p[t]*sys.tgrid.p[t])) << ", " << exp(-sys.tgrid.p[t]*sys.tgrid.p[t]*delta_r*delta_r) << endl;
          }
          if(delta_r < eps) t_integral += M_PI/(t_f*t_f);
          t_integrals[index_rp] = t_integral;
        }
      }
      coeff_new += integrate_system(sys, t_integrals, orbital_x_potential, NONE);
//cout << "index_r: " << index_r << ": " <<  coeff_new << ", " << sys.grid.w[p_r] << endl;
      coeffs_new[index_r] = -pow(M_PI, -1.5) * coeff_new;
// if(index_r == 5 || index_r == 6 || index_r == 7) {
 cout << "t integral: " << t_integrals << endl;
// cout << "orbital x pot: " << orbital_x_potential << endl;
//cout << "orbital: " << orbital_coeffs << endl;
//cout << "pot: " << sys.pot_coeffs << endl;
// }
// cout << "new coeff of the first function: " << coeff_new << endl;
// assert(false);
    }
  }
  orbital_coeffs = coeffs_new;
  cout << "new coeffs: " << coeffs_new << endl;
  auto stop = Clock::now();
  cout << "helmholtz: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
}
#endif


void mo1d::update_e(const system1d& sys, const localbasis1d& lo1d){
  auto start = Clock::now();
  cout << "--- begin update e ---" << endl;

  const double kinE = -0.5 * integrate_system(sys, orbital_coeffs, orbital_coeffs, NABLA2_APPR);
  cout << "kinE: " << kinE << endl;

  const double potE = integrate_system(sys, orbital_coeffs, orbital_coeffs, POTENTIAL_E);
  cout << "potE: " << potE << endl;
  orbital_e = potE + kinE;

  auto stop = Clock::now();
  cout << "update e: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
}


void mo1d::normalise(const system1d& sys, const localbasis1d& lo1d){
  auto start = Clock::now();
  cout << "--- begin normalise ---" << endl;
  const double norm = integrate_system(sys, orbital_coeffs, NONE);
  const double norm2 = integrate_system(sys, orbital_coeffs, orbital_coeffs, NONE);
  for(double &d: orbital_coeffs) d /= norm;
  orbital_e /= norm2;
  cout << "norm: " << norm << endl;
  cout << "norm2: " << norm2 << endl;
  auto stop = Clock::now();
  cout << "normalise: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
}



