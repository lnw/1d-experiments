
#include <cmath>
#include <fstream>  // for ofstream
#include <utility>  // for pair
#include <vector>

#include "auxiliary.hh"
#include "aux_math.hh"
#include "polynomial.hh"

using namespace std;


double polynomial::operator()(const double x) const {
#ifdef KAHAN
  vector<double> terms;
  double x_power = 1.0;
  for(int i=0; i<coeffs.size(); i++){
    terms.push_back(coeffs[i] * x_power);
    x_power *= x;
  }
  return kahan_sum(terms);
#else
  double result = 0.0;
  double x_power = 1.0;
  for(int i=0; i<coeffs.size(); i++){
    result += coeffs[i] * x_power;
    x_power *= x;
  }
  return result;
#endif
}

bool polynomial::operator==(const vector<int>& V){
  const double eps = 1.e-14;
  if(size() != V.size()) return false;
  for(int i=0; i<size(); i++){
    if( abs((*this)[i] - V[i])>eps ) return false;
  }
  return true;
}

bool polynomial::operator==(const polynomial& P){
  const double eps = 1.e-14;
  if(size() != P.size()) return false;
  for(int i=0; i<size(); i++){
    if( abs((*this)[i] - P[i])>eps ) return false;
  }
  return true;
}

// shift x by a to the right, x --> (x-a)
polynomial polynomial::shift(double a) const {
  polynomial p;
  p.resize(size());
  for(int i=0; i<size(); i++){
    for(int j=0; j<=i; j++){
      p[j] += binomial(i,j) * (*this)[i] * pow(-a, i-j); 
    }
  }
  return p;
}

polynomial polynomial::derive() const {
  polynomial p(*this);
  for(int i=1; i<p.size(); i++){
    p.coeffs[i-1] = i*p.coeffs[i];
  }
  if(size()>1) p.coeffs.resize(size()-1);
  else p.coeffs[0] = 0.0;
  return p;
}

polynomial polynomial::integrate() const {
  polynomial p(*this);
  p.coeffs.resize(size()+1);
  for(int i=p.size()-1; i>0; i--){
    p.coeffs[i] = p.coeffs[i-1]/i;
  }
  p.coeffs.front() = 0.0;
  return p;
}

double polynomial::integrate_limits(double x1, double x2) const {
  const polynomial p(this->integrate());
  return p(x2) - p(x1);
}

polynomial& polynomial::trim(){
  const double eps = 1.e-20;
  while(size() > 1 && fabs(coeffs.back()) < eps)
    coeffs.resize(size()-1);
  return *this;
}

