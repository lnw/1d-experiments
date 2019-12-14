#ifndef INTEGRATE_HH
#define INTEGRATE_HH

#include <vector>

// #include "local_basis_1d.hh"
#include "polynomial.hh"
#include "tensor.hh"
// #include "system_1d.hh"
// #include "system_3d.hh"

typedef enum {NONE, POTENTIAL_E, NABLA2_REAL, NABLA2_APPR, X} operator_t;

class localbasis1d;
class system1d;
class system3d;

//
// 1D
//
// this is for the construction of operators and does not include coefficients
double integrate_funcs(const system1d& sys, const polynomial& p1, operator_t op);
double integrate_cell(const system1d& sys, const tensor_1& coeffs1, operator_t op);
double integrate_system(const system1d& sys, const tensor_1& ev1, operator_t op);

// c is the global index on the cell, f1 and f2 are indices of functions in cell
// this is for the construction of operators and does not include coefficients
double integrate_funcs(const system1d& sys, const size_t c, const size_t f1, const size_t f2, operator_t op);
// double integrate_funcs_custom_pot(const system1d& sys, const size_t c, const size_t f1, const size_t f2,
//                                   const tensor_1 pot, operator_t op);
// coeffs1 and coeffs2 are coefficents of functions in the chosen cell
double integrate_cell(const system1d& sys, const size_t c, const tensor_1& coeffs1, const tensor_1& coeffs2, operator_t op);
double integrate_system(const system1d& sys, const tensor_1& coeffs1, const tensor_1& coeffs2, operator_t op);


//
// 3D
//
double integrate_system(const system3d& sys, const tensor_6& coeffs1, operator_t op);

double integrate_funcs(const system3d& sys, const size_t cz, const size_t cy, const size_t cx, const size_t fx1, const size_t fx2, const size_t fy1, const size_t fy2, const size_t fz1, const size_t fz2, operator_t op);
double integrate_system(const system3d& sys, const tensor_6& coeffs1, const tensor_6& coeffs2, operator_t op);

#endif
