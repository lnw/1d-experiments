#ifndef POLYNOMIAL_HH
#define POLYNOMIAL_HH

#include <cmath>
#include <fstream>  // for ofstream
#include <utility>  // for pair
#include <vector>

#include "auxiliary.hh"

using namespace std;


class polynomial{
public:
  vector<double> coeffs; // index = exponent


  polynomial(const vector<double> v): coeffs(v){}
  polynomial(const initializer_list<int> il) {for(int i: il) coeffs.push_back(double(i)); }
  polynomial(): coeffs(){}

  polynomial operator*(const double& d)     const { return polynomial(*this) *= d; }
  polynomial operator*(const polynomial& p) const { return polynomial(*this) *= p; }
  polynomial operator/(const double& d)     const { return polynomial(*this) /= d; }
  polynomial operator+(const polynomial& p) const { return polynomial(*this) += p; }
  polynomial operator-(const polynomial& p) const { return polynomial(*this) -= p; }
  polynomial& operator*=(const double& d){ for(size_t i=0; i<size(); i++) coeffs[i]*=d; return *this; }
  polynomial& operator*=(const polynomial& p){ 
    polynomial pp;
    pp.coeffs.resize((size()-1)+(p.size()-1)+1, 0.0);
    for(size_t i=0; i<size(); i++){
      for(size_t j=0; j<p.size(); j++){
        pp[i+j] += coeffs[i]*p[j];
      }
    }
    *this = pp;
    return *this;
  }
  polynomial& operator/=(const double& d){ for(size_t i=0; i<size(); i++) coeffs[i]/=d; return *this; }
  polynomial& operator+=(const polynomial& p){
    coeffs.resize(max(size(),p.size()));
    for(size_t i=0; i<p.size(); i++)
      coeffs[i] += p[i];
    return *this;
  } 
  polynomial& operator-=(const polynomial& p){
    coeffs.resize(max(size(),p.size()));
    for(size_t i=0; i<p.size(); i++)
      coeffs[i] -= p[i];
    return *this;
  } 

  bool operator==(const vector<int>& V);
  bool operator==(const polynomial& P);
  double& operator[](unsigned int i){ return coeffs[i]; }
  double  operator[](unsigned int i) const { return coeffs[i]; }
  size_t size() const {return coeffs.size();}

  polynomial& trim();
  polynomial& resize(const size_t k){coeffs.resize(k, 0.0); return *this;}

  double operator()(const double x) const;

  polynomial derive() const;
  polynomial integrate() const;
  double integrate_limits(double, double) const;

  polynomial shift(double x) const; // shift by x to the right

  
  friend ostream& operator<<(ostream& S, const polynomial& p) {
    // S << p.coeffs;
    S << "{";
    for (int i=0; i<p.coeffs.size(); i++){
      S << setprecision(12) <<  p.coeffs[i] << "*x^" << i << (i==p.coeffs.size()-1?"":" + ");
    }
    S << "}";
    return S;
  }
};

#endif
