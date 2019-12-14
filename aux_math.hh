#ifndef AUX_MATH
#define AUX_MATH

#include <vector>

#include "polynomial.hh"
#include "grid_1d.hh"
#include "system_1d.hh"

using namespace std;

int binomial(double n, double k);
double kahan_sum(const vector<double>& summands);
double sum_x_over_n_except_j(const size_t x, const vector<size_t> except, const vector<double> &vals);

double integral_gauss_polynomial_64bit(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);
double integral_gauss_polynomial_64bit_desc(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);
double integral_gauss_polynomial_64bit_asc(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);
double integral_gauss_polynomial_128bit(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);
double integral_gauss_polynomial_128bit_desc(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);
double integral_gauss_polynomial_128bit_asc(const double p, const double pos_r, const double c1, const double c2, const polynomial& poly);

inline size_t index_2(const size_t y, const size_t max_z, const size_t z){return y*max_z + z;}
inline size_t index_3(const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return x*max_y*max_z + index_2(y,max_z,z);}
inline size_t index_4(const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return w*max_x*max_y*max_z + index_3(x,max_y,y,max_z,z);}
inline size_t index_5(const size_t v, const size_t max_w, const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return v*max_w*max_x*max_y*max_z + index_4(w,max_x,x,max_y,y,max_z,z);}
inline size_t index_6(const size_t u, const size_t max_v, const size_t v, const size_t max_w, const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return u*max_v*max_w*max_x*max_y*max_z + index_5(v,max_w,w,max_x,x,max_y,y,max_z,z);}
inline size_t index_7(const size_t t, const size_t max_u, const size_t u, const size_t max_v, const size_t v, const size_t max_w, const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return t*max_u*max_v*max_w*max_x*max_y*max_z + index_6(u,max_v,v,max_w,w,max_x,x,max_y,y,max_z,z);}
inline size_t index_8(const size_t s, const size_t max_t, const size_t t, const size_t max_u, const size_t u, const size_t max_v, const size_t v, const size_t max_w, const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return s*max_t*max_u*max_v*max_w*max_x*max_y*max_z + index_7(t,max_u,u,max_v,v,max_w,w,max_x,x,max_y,y,max_z,z);}
inline size_t index_9(const size_t r, const size_t max_s, const size_t s, const size_t max_t, const size_t t, const size_t max_u, const size_t u, const size_t max_v, const size_t v, const size_t max_w, const size_t w, const size_t max_x, const size_t x, const size_t max_y, const size_t y, const size_t max_z, const size_t z) {return r*max_s*max_t*max_u*max_v*max_w*max_x*max_y*max_z + index_8(s,max_t,t,max_u,u,max_v,v,max_w,w,max_x,x,max_y,y,max_z,z);}

// in one cell with N functions, input one occupied orbital and return N-1 orthonormal functions
tensor_2 orthonormalise_cell(const system1d& sys, const tensor_1& occ, int c, bool debug_print);
// like above but omitting the first and the last polynomial from the virtual space
tensor_2 orthonormalise_cell_wo_boundary(const system1d& sys, const tensor_1& occ, int c, bool debug_print);

#endif 

