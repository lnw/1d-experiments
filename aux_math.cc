
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
// #include <iomanip>

#include "tensor.hh"
#include "aux_math.hh"
#include "polynomial.hh"
#include "integrate.hh"


using namespace std;

// sum numbers in a vector without loosing too much accuracy
double kahan_sum(const vector<double>& summands){
  double sum = 0.0;
  double corr = 0.0;                 // A running compensation for lost low-order bits.
  for(int i=0; i<summands.size(); i++){
    double y = summands[i] - corr;    // So far, so good: c is zero.
    double t = sum + y;         // Alas, sum is big, y small, so low-order digits of y are lost.
    corr = (t - sum) - y;       // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
    sum = t;                 // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
  }                   // Next time around, the lost low part will be added to y in a fresh attempt.
  return sum;
}

template<typename T>
bool abs_less(T a, T b) {return abs(a) < abs(b);}

template<typename T>
bool abs_greater(T a, T b) {return abs(a) > abs(b);}

template<typename T>
T sort_desc_and_accumulate(vector<T>& vals){
  sort(vals.begin(), vals.end(), abs_greater<T>);
  T sum = 0;
  for(const T &val: vals) sum += val;
  return sum;
}

template<typename T>
T sort_asc_and_accumulate(vector<T>& vals){
  sort(vals.begin(), vals.end(), abs_less<T>);
  T sum = 0;
  for(const T &val: vals) sum += val;
  return sum;
}


// used to get one coefficient of the polynomial that is the product of (x-x1)(x-x2)...(x-xn)
// where some x_i are missing
// there are vals.size() possible terms, except are excluded, x are used.
// the sum is over all permutations of used over available terms
double sum_x_over_n_except_j(const size_t x, const vector<size_t> except, const vector<double>& vals) {
  // sorted bit-array with zeros and ones (for included and excluded values)
  // cout << "x,j,vals " << x << ", " << j << ", "<< vals << endl;
  vector<short> subset(vals.size()-x,0);
  subset.resize(vals.size()-except.size(),1);

  double sum=0.0;
  do {
    // cout << subset << endl;
    double product=1.0;
    for(int i=0, k=0;i<vals.size(); i++, k++){
      while(find(except.begin(), except.end(), k)!=except.end()) k++; // k in except ?
      if(k==vals.size()) break;
      bool bit(subset[i]);
      if(bit){
        product *= -vals[k]; // because a 'zero' at x_i enters as (x - x_i)
        // cout << "factor : " << -vals[k] << endl;
      }
      // cout << "prod: " << product << endl;
    }
    sum += product;
    // cout << "product : " << product << endl;
  } while(std::next_permutation(subset.begin(), subset.end()));
  // cout << "sum:" << sum << endl;

  return sum;
}


int binomial(double n, double k){
  if(k>n) return 0;
  if(2*k > n) k=n-k;
  if(int(k+0.5)==0) return 1;
  return int(n/k*binomial(n-1, k-1) + 0.5);
}


// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_64bit(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(double)>=8);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const double t2 = t*t;
  const double t3 = t2*t;
  const double t4 = t2*t2;
  const double t6 = t3*t3;
  const double t8 = t4*t4;
  const double p2 = p*p;
  const double p3 = p2*p;
  const double p4 = p2*p2;
  const double p5 = p3*p2;
  const double p6 = p4*p2;
  const double gauss_c1 = exp(-t2*(p-c1)*(p-c1));
  const double gauss_c2 = exp(-t2*(p-c2)*(p-c2));
  const double sqrt_api = t*sqrt(M_PI);
  const double sqrt_pi  = sqrt(M_PI);
  const double erf_c1   = erf(t*(c1-p));
  const double erf_c2   = erf(t*(c2-p));

//  cout << t << "," << p << "," << c1 << "," << c2 << poly << endl;

  const double term_0 = poly[0] * sqrt(M_PI)/(2*t) * (erf_c2 - erf_c1);
  integral += term_0;
cout << "(64)  0 " << term_0 << endl;

  if(p_deg >= 1){
    const double term_1 = poly[1]/(2*t2) * (-gauss_c2 + gauss_c1 + sqrt_api*p*(erf_c2 - erf_c1));
    integral += term_1;
cout << "(64)  1 " << term_1 << endl;
  }

  if(p_deg >= 2){
    const double aux_1 = -1/(2*t2) * (gauss_c2*(p+c2)-gauss_c1*(p+c1));
    const double aux_2 = sqrt_pi/(4*t3) * (1+2*t2*p2) * (erf_c2-erf_c1);
    const double term_2 = poly[2] * (aux_1 + aux_2);
    integral += term_2;
cout << "(64)  2 " << term_2 << endl;
  }

  if(p_deg >= 3){
    const double aux_11 = -gauss_c2 + gauss_c1;
    const double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    const double aux_13 = (-gauss_c2*c2 + gauss_c1*c1)*t2*p;
    const double aux_14 = (-gauss_c2*c2*c2 + gauss_c1*c1*c1)*t2;
    const double aux_2 = 0.5 * sqrt_api * p * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const double term_3 = poly[3] / (2*t4) * (aux_11 + aux_12 + aux_13 + aux_14 + aux_2);
    integral += term_3;
 cout << "(64)  3 " << term_3 << endl;
  }

  if(p_deg >= 4){
    const double aux_11 = -gauss_c2*(5*p + 2*t2*p3 + 3*c2 + 2*t2*p2*c2 + 2*t2*p*c2*c2 + 2*t2*c2*c2*c2);
    const double aux_12 = -gauss_c1*(5*p + 2*t2*p3 + 3*c1 + 2*t2*p2*c1 + 2*t2*p*c1*c1 + 2*t2*c1*c1*c1);
    const double aux_2 = sqrt_pi/(2*t) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const double term_4 = poly[4]/(4*t4)* (aux_11 - aux_12 + aux_2);
    integral += term_4;
 cout << "(64)  4 " << term_4 << endl;
  }

  if(p_deg >= 5){
    const double aux_11 = 8*(-gauss_c2 + gauss_c1);
    const double aux_12 = 18*t2*p2*(-gauss_c2 + gauss_c1);
    const double aux_13 = 14*t2*p*(-gauss_c2*c2 + gauss_c1*c1);
    const double aux_14 = 8*t2*(-gauss_c2*c2*c2 + gauss_c1*c1*c1);
    const double aux_15 = 4*t4*p4*(-gauss_c2 + gauss_c1);
    const double aux_16 = 4*t4*p3*(-gauss_c2*c2 + gauss_c1*c1);
    const double aux_17 = 4*t4*p2*(-gauss_c2*c2*c2 + gauss_c1*c1*c1);
    const double aux_18 = 4*t4*p*(-gauss_c2*c2*c2*c2 + gauss_c1*c1*c1*c1);
    const double aux_19 = 4*t4*(-gauss_c2*c2*c2*c2*c2 + gauss_c1*c1*c1*c1*c1);
    const double aux_2 = sqrt_api * p * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const double term_5 = poly[5]/(8*t6) * (aux_11 + aux_12 + aux_13 + aux_14 + aux_15 + aux_16 + aux_17 + aux_18 + aux_19 + aux_2);
    // const double aux_11 = -2*gauss_c2*(4/t4 + 1/t2*(9*p2 + 7*p*c2 + 4*c2*c2) + 2*(p4 + p3*c2 + p2*c2*c2 + p*c2*c2*c2 + c2*c2*c2*c2));
    // const double aux_12 = -2*gauss_c1*(4/t4 + 1/t2*(9*p2 + 7*p*c1 + 4*c1*c1) + 2*(p4 + p3*c1 + p2*c1*c1 + p*c1*c1*c1 + c1*c1*c1*c1));
    // const double aux_2 = sqrt_api * p * (15/t4 + 20*p2/t2 + 4*p4) * (erf_c2 - erf_c1);
    // const double term_5 = poly[5]/(8*t2) * (aux_11 - aux_12 + aux_2);
    integral += term_5;
 cout << "(64)  5 " << term_5 << endl;
  }

  if(p_deg >= 6){
    const double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*c2 + 4*t2*p2*c2*(6 + t2*c2*c2) + 4*t2*p3*(7 + t2*c2*c2) + c2*(15 + 10*t2*c2*c2 + 4*t4*c2*c2*c2*c2) + p*(33 + 18*t2*c2*c2 + 4*t4*c2*c2*c2*c2));
    const double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*c1 + 4*t2*p2*c1*(6 + t2*c1*c1) + 4*t2*p3*(7 + t2*c1*c1) + c1*(15 + 10*t2*c1*c1 + 4*t4*c1*c1*c1*c1) + p*(33 + 18*t2*c1*c1 + 4*t4*c1*c1*c1*c1));
    const double aux_2 = sqrt_pi/(2*t) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);
    const double term_6 = poly[6]/(8*t6) * (aux_11 - aux_12 + aux_2);
    integral += term_6;
 cout << "(64)  6 " << term_6 << endl;
  }

  if(p_deg >= 7){
    const double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*p*c2 + 8*c2*c2) + 2*t4*(20*p4 + 18*p3*c2 + 15*p2*c2*c2 + 11*p*c2*c2*c2 + 6*c2*c2*c2*c2) + 4*t6*(p6 + p5*c2 + p4*c2*c2 + p3*c2*c2*c2 + p2*c2*c2*c2*c2 + p*c2*c2*c2*c2*c2 + c2*c2*c2*c2*c2*c2));
    const double aux_12 = -2*gauss_c1*(24 + 3*t2*(29*p2 + 19*p*c1 + 8*c1*c1) + 2*t4*(20*p4 + 18*p3*c1 + 15*p2*c1*c1 + 11*p*c1*c1*c1 + 6*c1*c1*c1*c1) + 4*t6*(p6 + p5*c1 + p4*c1*c1 + p3*c1*c1*c1 + p2*c1*c1*c1*c1 + p*c1*c1*c1*c1*c1 + c1*c1*c1*c1*c1*c1));
    const double aux_2 = sqrt_api * p * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);
    const double term_7 = poly[7]/(16*t8) * (aux_11 - aux_12 + aux_2);
    integral += term_7;
 cout << "(64)  7 " << term_7 << endl;
  }

  if(p_deg > 7){
    assert(false);
  }

  return integral;
}


// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_128bit(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(long double)>=16);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const long double lt = t;
  const long double t2 = lt*lt;
  const long double t3 = t2*lt;
  const long double t4 = t2*t2;
  const long double t6 = t3*t3;
  const long double t8 = t4*t4;
  const long double lp = p;
  const long double p2 = lp*lp;
  const long double p3 = p2*lp;
  const long double p4 = p2*p2;
  const long double p5 = p3*p2;
  const long double p6 = p4*p2;
  const long double lc1 = c1;
  const long double lc2 = c2;
  const long double gauss_c1 = exp(-t2*(lp-lc1)*(lp-lc1));
  const long double gauss_c2 = exp(-t2*(lp-lc2)*(lp-lc2));
  const long double sqrt_api = lt*sqrt(M_PIl);
  const long double sqrt_pi  = sqrt(M_PIl);
  const long double erf_c1   = erf(lt*(lc1-lp));
  const long double erf_c2   = erf(lt*(lc2-lp));

//  cout << t << "," << p << "," << lc1 << "," << lc2 << poly << endl;

  {
    const long double term_0 = poly[0] * sqrt(M_PIl)/(2*lt) * (erf_c2 - erf_c1);
    integral += term_0;
    //cout << "  0 " << term_0 << endl;
    if(p_deg < 1) return integral;
  }

  {
    const long double term_1 = poly[1]/(2*t2) * (-gauss_c2 + gauss_c1 + sqrt_api*lp*(erf_c2 - erf_c1));
    integral += term_1;
    //cout << "  1 " << term_1 << endl;
    if(p_deg < 2) return integral;
  }

  {
    const long double aux_1 = -1/(2*t2) * (gauss_c2*(lp+lc2) - gauss_c1*(lp+lc1));
    const long double aux_2 = sqrt_pi/(4*t3) * (1 + 2*t2*p2) * (erf_c2 - erf_c1);
    const long double term_2 = poly[2] * (aux_1 + aux_2);
    integral += term_2;
    //cout << "  2 " << term_2 << endl;
    if(p_deg < 3) return integral;
  }

  {
    const long double aux_11 = -gauss_c2 + gauss_c1;
    const long double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    const long double aux_13 = (-gauss_c2*lc2 + gauss_c1*lc1)*t2*lp;
    const long double aux_14 = (-gauss_c2*lc2*lc2  + gauss_c1*lc1*lc1)*t2;
    // const long double aux_12 = -gauss_c2*t2*(p2 + lp*lc2 + lc2*lc2);
    // const long double aux_13 =  gauss_c1*t2*(p2 + lp*lc1 + lc1*lc1);
    const long double aux_2 = 0.5 * sqrt_api * lp * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const long double term_3 = poly[3] / (2*t4) * (aux_11 + aux_12 + aux_13 + aux_14 + aux_2);
    integral += term_3;
    //cout << "  3 " << term_3 << endl;
    if(p_deg < 4) return integral;
  }

  {
    const long double aux_11 = -gauss_c2*(5*lp + 2*t2*p3 + 3*lc2 + 2*t2*p2*lc2 + 2*t2*lp*lc2*lc2 + 2*t2*lc2*lc2*lc2);
    const long double aux_12 = gauss_c1*(5*lp + 2*t2*p3 + 3*lc1 + 2*t2*p2*lc1 + 2*t2*lp*lc1*lc1 + 2*t2*lc1*lc1*lc1);
    const long double aux_2 = sqrt_pi/(2*lt) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const long double term_4 = poly[4]/(4*t4)* (aux_11 + aux_12 + aux_2);
    integral += term_4;
    //cout << "  4 " << term_4 << endl;
    if(p_deg < 5) return integral;
  }

  {
    const long double aux_11 = 8*(-gauss_c2 + gauss_c1);
    const long double aux_12 = 18*t2*p2*(-gauss_c2 + gauss_c1);
    const long double aux_13 = 14*t2*lp*(-gauss_c2*lc2 + gauss_c1*lc1);
    const long double aux_14 = 8*t2*(-gauss_c2*lc2*lc2 + gauss_c1*lc1*lc1);
    const long double aux_15 = 4*t4*p4*(-gauss_c2 + gauss_c1);
    const long double aux_16 = 4*t4*p3*(-gauss_c2*lc2 + gauss_c1*lc1);
    const long double aux_17 = 4*t4*p2*(-gauss_c2*lc2*lc2 + gauss_c1*lc1*lc1);
    const long double aux_18 = 4*t4*lp*(-gauss_c2*lc2*lc2*lc2 + gauss_c1*lc1*lc1*lc1);
    const long double aux_19 = 4*t4*(-gauss_c2*lc2*lc2*lc2*lc2 + gauss_c1*lc1*lc1*lc1*lc1);
    const long double aux_2 = sqrt_api * lp * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const long double term_5 = poly[5]/(8*t6) * (aux_11 + aux_12 + aux_13 + aux_14 + aux_15 + aux_16 + aux_17 + aux_18 + aux_19 + aux_2);
    integral += term_5;
  // cout << "  5 " << term_5 << endl;
    if(p_deg < 6) return integral;
  }

  {
    const long double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*lc2 + 4*t2*p2*lc2*(6 + t2*lc2*lc2) + 4*t2*p3*(7 + t2*lc2*lc2) + lc2*(15 + 10*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2) + lp*(33 + 18*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2));
    const long double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*lc1 + 4*t2*p2*lc1*(6 + t2*lc1*lc1) + 4*t2*p3*(7 + t2*lc1*lc1) + lc1*(15 + 10*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1) + lp*(33 + 18*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1));
    const long double aux_2 = sqrt_pi/(2*lt) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);
    const long double term_6 = poly[6]/(8*t6) * (aux_11 - aux_12 + aux_2);
    integral += term_6;
    // cout << "  6 " << term_6 << endl;
    if(p_deg < 7) return integral;
  }

  {
    const long double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*lp*lc2 + 8*lc2*lc2) + 2*t4*(20*p4 + 18*p3*lc2 + 15*p2*lc2*lc2 + 11*lp*lc2*lc2*lc2 + 6*lc2*lc2*lc2*lc2) + 4*t6*(p6 + p5*lc2 + p4*lc2*lc2 + p3*lc2*lc2*lc2 + p2*lc2*lc2*lc2*lc2 + lp*lc2*lc2*lc2*lc2*lc2 + lc2*lc2*lc2*lc2*lc2*lc2));
    const long double aux_12 = -2*gauss_c1*(24 + 3*t2*(29*p2 + 19*lp*lc1 + 8*lc1*lc1) + 2*t4*(20*p4 + 18*p3*lc1 + 15*p2*lc1*lc1 + 11*lp*lc1*lc1*lc1 + 6*lc1*lc1*lc1*lc1) + 4*t6*(p6 + p5*lc1 + p4*lc1*lc1 + p3*lc1*lc1*lc1 + p2*lc1*lc1*lc1*lc1 + lp*lc1*lc1*lc1*lc1*lc1 + lc1*lc1*lc1*lc1*lc1*lc1));
    const long double aux_2 = sqrt_api * lp * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);
    const long double term_7 = poly[7]/(16*t8) * (aux_11 - aux_12 + aux_2);
    integral += term_7;
    // cout << "  7 " << term_7 << endl;
    if(p_deg < 7) return integral;
  }

  assert(false);
}


// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_128bit_desc(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(long double)>=16);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const long double lt = t;
  const long double t2 = lt*lt;
  const long double t3 = t2*lt;
  const long double t4 = t2*t2;
  const long double t6 = t3*t3;
  const long double t8 = t4*t4;
  const long double lp = p;
  const long double p2 = lp*lp;
  const long double p3 = p2*lp;
  const long double p4 = p2*p2;
  const long double p5 = p3*p2;
  const long double p6 = p4*p2;
  const long double lc1 = c1;
  const long double lc2 = c2;
  const long double gauss_c1 = exp(-t2*(lp-lc1)*(lp-lc1));
  const long double gauss_c2 = exp(-t2*(lp-lc2)*(lp-lc2));
  const long double sqrt_api = lt*sqrt(M_PIl);
  const long double sqrt_pi  = sqrt(M_PIl);
  const long double erf_c1   = erf(lt*(lc1-lp));
  const long double erf_c2   = erf(lt*(lc2-lp));

//  cout << t << "," << p << "," << lc1 << "," << lc2 << poly << endl;

  {
    const long double term_0 = poly[0] * sqrt(M_PIl)/(2*lt) * (erf_c2 - erf_c1);
    integral += term_0;
    //cout << "  0 " << term_0 << endl;
    if(p_deg < 1) return integral;
  }

  {
    vector<long double> terms(4);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    terms[2] = sqrt_api*lp*erf_c2;
    terms[3] = -sqrt_api*lp*erf_c1;
    const long double sum = sort_desc_and_accumulate(terms);
// cout << setprecision(15) << terms << endl;
    integral += poly[1]/(2*t2) * sum;
    //cout << "  1 " << term_1 << endl;
    if(p_deg < 2) return integral;
  }

  {
    vector<long double> terms(8);
    terms[0] = -1/(2*t2) * gauss_c2*lp;
    terms[1] = -1/(2*t2) * gauss_c2*lc2;
    terms[2] = 1/(2*t2) * gauss_c1*lp;
    terms[3] = 1/(2*t2) * gauss_c1*lc1;
    terms[4] = sqrt_pi/(4*t3) * erf_c2;
    terms[5] = - sqrt_pi/(4*t3) * erf_c1;
    terms[6] = sqrt_pi/(4*t3) * 2*t2*p2 * erf_c2;
    terms[7] = - sqrt_pi/(4*t3) *  2*t2*p2 * erf_c1;
    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[2] * sum;
    //cout << "  2 " << term_2 << endl;
    if(p_deg < 3) return integral;
  }

  {
    vector<long double> terms(12);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    // const long double aux_11 = -gauss_c2 + gauss_c1;
    terms[2] = -gauss_c2*t2*p2;
    terms[3] = gauss_c1*t2*p2;
    // const long double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    terms[4] = -gauss_c2*lc2 *t2*lp;
    terms[5] = gauss_c1*lc1*t2*lp;
    terms[6] = -gauss_c2*lc2*lc2*t2;
    terms[7] = gauss_c1*lc1*lc1*t2;
    // const long double aux_13 = (-gauss_c2*lc2 + gauss_c1*lc1)*t2*lp;
    // const long double aux_14 = (-gauss_c2*lc2*lc2  + gauss_c1*lc1*lc1)*t2;
    // const long double aux_12 = -gauss_c2*t2*(p2 + lp*lc2 + lc2*lc2);
    // const long double aux_13 =  gauss_c1*t2*(p2 + lp*lc1 + lc1*lc1);
    terms[8] = 0.5 * sqrt_api * lp * 3 * erf_c2;
    terms[9] = -0.5 * sqrt_api * lp * 3 * erf_c1;
    terms[10] = 0.5 * sqrt_api * lp * 2*t2*p2 * erf_c2;
    terms[11] = -0.5 * sqrt_api * lp * 2*t2*p2 * erf_c1;
    // const long double aux_2 = 0.5 * sqrt_api * lp * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[3] / (2*t4) * sum;
    //cout << "  3 " << term_3 << endl;
    if(p_deg < 4) return integral;
  }

  {
    vector<long double> terms(18);
    terms[0] = -gauss_c2* 5*lp;
    terms[1] = -gauss_c2* 2*t2*p3;
    terms[2] = -gauss_c2* 3*lc2;
    terms[3] = -gauss_c2* 2*t2*p2*lc2;
    terms[4] = -gauss_c2* 2*t2*lp*lc2*lc2;
    terms[5] = -gauss_c2* 2*t2*lc2*lc2*lc2;
    // const long double aux_11 = -gauss_c2*(5*lp + 2*t2*p3 + 3*lc2 + 2*t2*p2*lc2 + 2*t2*lp*lc2*lc2 + 2*t2*lc2*lc2*lc2);
    terms[6] =  gauss_c1* 5*lp;
    terms[7] =  gauss_c1* 2*t2*p3;
    terms[8] =  gauss_c1* 3*lc1;
    terms[9] =  gauss_c1* 2*t2*p2*lc1;
    terms[10] = gauss_c1* 2*t2*lp*lc1*lc1;
    terms[11] = gauss_c1* 2*t2*lc1*lc1*lc1;
    // const long double aux_12 = gauss_c1*(5*lp + 2*t2*p3 + 3*lc1 + 2*t2*p2*lc1 + 2*t2*lp*lc1*lc1 + 2*t2*lc1*lc1*lc1);
    terms[12] = -sqrt_pi/(2*lt) * 3 * erf_c1;
    terms[13] = sqrt_pi/(2*lt) * 3 * erf_c2;
    terms[14] = -sqrt_pi/(lt) * 6*t2*p2 * erf_c1;
    terms[15] = sqrt_pi/(lt) * 6*t2*p2 * erf_c2;
    terms[16] = -sqrt_pi/(lt) * 2*t4*p4 * erf_c1;
    terms[17] = sqrt_pi/(lt) * 2*t4*p4 * erf_c2;
    // const long double aux_2 = sqrt_pi/(2*lt) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[4]/(4*t4) * sum;
    //cout << "  4 " << term_4 << endl;
    if(p_deg < 5) return integral;
  }

  {
    vector<long double> terms(24);
    terms[0] = -14*t2*lp*gauss_c2*lc2;
    terms[1] = 14*t2*lp*gauss_c1*lc1;
    terms[2] = -18*t2*p2*gauss_c2;
    terms[3] = 18*t2*p2*gauss_c1;
    terms[4] = -4*t4*gauss_c2*lc2*lc2*lc2*lc2;
    terms[5] = 4*t4*gauss_c1*lc1*lc1*lc1*lc1;
    terms[6] = -4*t4*lp*gauss_c2*lc2*lc2*lc2;
    terms[7] = 4*t4*lp*gauss_c1*lc1*lc1*lc1;
    terms[8] = -4*t4*p2*gauss_c2*lc2*lc2;
    terms[9] = 4*t4*p2*gauss_c1*lc1*lc1;
    terms[10] = -4*t4*p3*gauss_c2*lc2;
    terms[11] = 4*t4*p3*gauss_c1*lc1;
    terms[12] = -4*t4*p4*gauss_c2;
    terms[13] = 4*t4*p4*gauss_c1;
    terms[14] = -8*gauss_c2;
    terms[15] = 8*gauss_c1;
    terms[16] = -8*t2*gauss_c2*lc2*lc2;
    terms[17] = 8*t2*gauss_c1*lc1*lc1;
    terms[18] = sqrt_api * lp * 15.0 * erf_c2;
    terms[19] = -sqrt_api * lp * 15.0 * erf_c1;
    terms[20] = sqrt_api * lp * 20*p2*t2 * erf_c2;
    terms[21] = -sqrt_api * lp * 20*p2*t2 * erf_c1;
    terms[22] = sqrt_api * lp * 4*p4*t4 * erf_c2;
    terms[23] = -sqrt_api * lp * 4*p4*t4 * erf_c1;
    // const long double aux_2 = sqrt_api * lp * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[5]/(8*t6) * sum;
  // cout << "  5 " << term_5 << endl;
    if(p_deg < 6) return integral;
  }

  {
    vector<long double> terms(32);
    terms[0] = -gauss_c2* 4*t4*p5;
    terms[1] = -gauss_c2* 4*t4*p4*lc2;
    terms[2] = -gauss_c2* 24*t2*p2*lc2;
    terms[3] = -gauss_c2* 4*t2*p2*lc2* t2*lc2*lc2;
    terms[4] = -gauss_c2* 28*t2*p3;
    terms[5] = -gauss_c2* 4*t2*p3* t2*lc2*lc2;
    terms[6] = -gauss_c2* lc2*15;
    terms[7] = -gauss_c2* lc2* 10*t2*lc2*lc2;
    terms[8] = -gauss_c2* lc2* 4*t4*lc2*lc2*lc2*lc2;
    terms[9] = -gauss_c2* lp*33;
    terms[10] = -gauss_c2* lp*18*t2*lc2*lc2;
    terms[11] = -gauss_c2* lp* 4*t4*lc2*lc2*lc2*lc2;
    // const long double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*lc2 + 4*t2*p2*lc2*(6 + t2*lc2*lc2) + 4*t2*p3*(7 + t2*lc2*lc2) + lc2*(15 + 10*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2) + lp*(33 + 18*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2));

    terms[12] = 4* gauss_c1* t4*p5;
    terms[13] = 4* gauss_c1* t4*p4*lc1;
    terms[14] = 24*gauss_c1* t2*p2*lc1;
    terms[15] = 4* gauss_c1* t2*p2*lc1* t2*lc1*lc1;
    terms[16] = 28*gauss_c1* t2*p3;
    terms[17] = 4* gauss_c1* t2*p3* t2*lc1*lc1;
    terms[18] = 15*gauss_c1* lc1;
    terms[19] = 10*gauss_c1* lc1* t2*lc1*lc1;
    terms[20] = 4* gauss_c1* lc1* t4*lc1*lc1*lc1*lc1;
    terms[21] = 33*gauss_c1* lp;
    terms[22] = 18*gauss_c1* lp*t2*lc1*lc1;
    terms[23] = 4* gauss_c1* lp* t4*lc1*lc1*lc1*lc1;
    // const long double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*lc1 + 4*t2*p2*lc1*(6 + t2*lc1*lc1) + 4*t2*p3*(7 + t2*lc1*lc1) + lc1*(15 + 10*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1) + lp*(33 + 18*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1));

    terms[24] =  sqrt_pi/(2*lt) * 15 * (erf_c2);
    terms[25] = -sqrt_pi/(2*lt) * 15 * (erf_c1);
    terms[26] =  sqrt_pi/(lt) * (45*t2*p2) * (erf_c2);
    terms[27] = -sqrt_pi/(lt) * (45*t2*p2) * (erf_c1);
    terms[28] =  sqrt_pi/(lt) * (30*t4*p4) * (erf_c2);
    terms[29] = -sqrt_pi/(lt) * (30*t4*p4) * (erf_c1);
    terms[30] =  sqrt_pi/(lt) * (4*t6*p6) * (erf_c2);
    terms[31] = -sqrt_pi/(lt) * (4*t6*p6) * (erf_c1);
    // const long double aux_2 = sqrt_pi/(2*lt) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[6]/(8*t6) * sum;
    // cout << "  6 " << term_6 << endl;
    if(p_deg < 7) return integral;
  }

  {
    vector<long double> terms(40);
    terms[0] = -2*gauss_c2* 24;
    terms[1] = -2*gauss_c2* 3*t2*29*p2;
    terms[2] = -2*gauss_c2* 3*t2* 19*lp*lc2;
    terms[3] = -2*gauss_c2* 3*t2* 8*lc2*lc2;
    terms[4] = -2*gauss_c2* 2*t4*20*p4;
    terms[5] = -2*gauss_c2* 2*t4* 18*p3*lc2;
    terms[6] = -2*gauss_c2* 2*t4* 15*p2*lc2*lc2;
    terms[7] = -2*gauss_c2* 2*t4* 11*lp*lc2*lc2*lc2;
    terms[8] = -2*gauss_c2* 2*t4* 6*lc2*lc2*lc2*lc2;
    terms[9] = -2*gauss_c2* 4*t6*p6;
    terms[10] = -2*gauss_c2* 4*t6* p5*lc2;
    terms[11] = -2*gauss_c2* 4*t6* p4*lc2*lc2;
    terms[12] = -2*gauss_c2* 4*t6* p3*lc2*lc2*lc2;
    terms[13] = -2*gauss_c2* 4*t6* p2*lc2*lc2*lc2*lc2;
    terms[14] = -2*gauss_c2* 4*t6* lp*lc2*lc2*lc2*lc2*lc2;
    terms[15] = -2*gauss_c2* 4*t6* lc2*lc2*lc2*lc2*lc2*lc2;
    // const long double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*lp*lc2 + 8*lc2*lc2) + 2*t4*(20*p4 + 18*p3*lc2 + 15*p2*lc2*lc2 + 11*lp*lc2*lc2*lc2 + 6*lc2*lc2*lc2*lc2) + 4*t6*(p6 + p5*lc2 + p4*lc2*lc2 + p3*lc2*lc2*lc2 + p2*lc2*lc2*lc2*lc2 + lp*lc2*lc2*lc2*lc2*lc2 + lc2*lc2*lc2*lc2*lc2*lc2));

    terms[16] = 2*gauss_c1* 24;
    terms[17] = 2*gauss_c1* 3*t2*29*p2;
    terms[18] = 2*gauss_c1* 3*t2* 19*lp*lc1;
    terms[19] = 2*gauss_c1* 3*t2* 8*lc1*lc1;
    terms[20] = 2*gauss_c1* 2*t4*20*p4;
    terms[21] = 2*gauss_c1* 2*t4* 18*p3*lc1;
    terms[22] = 2*gauss_c1* 2*t4* 15*p2*lc1*lc1;
    terms[23] = 2*gauss_c1* 2*t4* 11*lp*lc1*lc1*lc1;
    terms[24] = 2*gauss_c1* 2*t4* 6*lc1*lc1*lc1*lc1;
    terms[25] = 2*gauss_c1* 4*t6*p6;
    terms[26] = 2*gauss_c1* 4*t6* p5*lc1;
    terms[27] = 2*gauss_c1* 4*t6* p4*lc1*lc1;
    terms[28] = 2*gauss_c1* 4*t6* p3*lc1*lc1*lc1;
    terms[29] = 2*gauss_c1* 4*t6* p2*lc1*lc1*lc1*lc1;
    terms[30] = 2*gauss_c1* 4*t6* lp*lc1*lc1*lc1*lc1*lc1;
    terms[31] = 2*gauss_c1* 4*t6* lc1*lc1*lc1*lc1*lc1*lc1;
    // const long double aux_12 = 2*gauss_c1*(24 + 3*t2*(29*p2 + 19*lp*lc1 + 8*lc1*lc1) + 2*t4*(20*p4 + 18*p3*lc1 + 15*p2*lc1*lc1 + 11*lp*lc1*lc1*lc1 + 6*lc1*lc1*lc1*lc1) + 4*t6*(p6 + p5*lc1 + p4*lc1*lc1 + p3*lc1*lc1*lc1 + p2*lc1*lc1*lc1*lc1 + lp*lc1*lc1*lc1*lc1*lc1 + lc1*lc1*lc1*lc1*lc1*lc1));

    terms[32] = -sqrt_api * lp * (105) * (erf_c1);
    terms[33] =  sqrt_api * lp * (105) * (erf_c2);
    terms[34] = -sqrt_api * lp * (210*t2*p2) * (erf_c1);
    terms[35] =  sqrt_api * lp * (210*t2*p2) * (erf_c2);
    terms[36] = -sqrt_api * lp * (84*t4*p4) * (erf_c1);
    terms[37] =  sqrt_api * lp * (84*t4*p4) * (erf_c2);
    terms[38] = -sqrt_api * lp * (8*t6*p6) * (erf_c1);
    terms[39] =  sqrt_api * lp * (8*t6*p6) * (erf_c2);
    // const long double aux_2 = sqrt_api * lp * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const long double sum = sort_desc_and_accumulate(terms);
    integral += poly[7]/(16*t8) * sum;
    // cout << "  7 " << term_7 << endl;
    if(p_deg < 7) return integral;
  }

  assert(false);
}


// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_128bit_asc(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(long double)>=16);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const long double lt = t;
  const long double t2 = lt*lt;
  const long double t3 = t2*lt;
  const long double t4 = t2*t2;
  const long double t6 = t3*t3;
  const long double t8 = t4*t4;
  const long double lp = p;
  const long double p2 = lp*lp;
  const long double p3 = p2*lp;
  const long double p4 = p2*p2;
  const long double p5 = p3*p2;
  const long double p6 = p4*p2;
  const long double lc1 = c1;
  const long double lc2 = c2;
  const long double gauss_c1 = exp(-t2*(lp-lc1)*(lp-lc1));
  const long double gauss_c2 = exp(-t2*(lp-lc2)*(lp-lc2));
  const long double sqrt_api = lt*sqrt(M_PIl);
  const long double sqrt_pi  = sqrt(M_PIl);
  const long double erf_c1   = erf(lt*(lc1-lp));
  const long double erf_c2   = erf(lt*(lc2-lp));

//  cout << t << "," << p << "," << lc1 << "," << lc2 << poly << endl;

  {
    const long double term_0 = poly[0] * sqrt(M_PIl)/(2*lt) * (erf_c2 - erf_c1);
    integral += term_0;
    //cout << "  0 " << term_0 << endl;
    if(p_deg < 1) return integral;
  }

  {
    vector<long double> terms(4);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    terms[2] = sqrt_api*lp*erf_c2;
    terms[3] = -sqrt_api*lp*erf_c1;
    const long double sum = sort_asc_and_accumulate(terms);
// cout << setprecision(15) << terms << endl;
    integral += poly[1]/(2*t2) * sum;
    //cout << "  1 " << term_1 << endl;
    if(p_deg < 2) return integral;
  }

  {
    vector<long double> terms(8);
    terms[0] = -1/(2*t2) * gauss_c2*lp;
    terms[1] = -1/(2*t2) * gauss_c2*lc2;
    terms[2] = 1/(2*t2) * gauss_c1*lp;
    terms[3] = 1/(2*t2) * gauss_c1*lc1;
    terms[4] = sqrt_pi/(4*t3) * erf_c2;
    terms[5] = - sqrt_pi/(4*t3) * erf_c1;
    terms[6] = sqrt_pi/(4*t3) * 2*t2*p2 * erf_c2;
    terms[7] = - sqrt_pi/(4*t3) *  2*t2*p2 * erf_c1;
    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[2] * sum;
    //cout << "  2 " << term_2 << endl;
    if(p_deg < 3) return integral;
  }

  {
    vector<long double> terms(12);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    // const long double aux_11 = -gauss_c2 + gauss_c1;
    terms[2] = -gauss_c2*t2*p2;
    terms[3] = gauss_c1*t2*p2;
    // const long double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    terms[4] = -gauss_c2*lc2 *t2*lp;
    terms[5] = gauss_c1*lc1*t2*lp;
    terms[6] = -gauss_c2*lc2*lc2*t2;
    terms[7] = gauss_c1*lc1*lc1*t2;
    // const long double aux_13 = (-gauss_c2*lc2 + gauss_c1*lc1)*t2*lp;
    // const long double aux_14 = (-gauss_c2*lc2*lc2  + gauss_c1*lc1*lc1)*t2;
    // const long double aux_12 = -gauss_c2*t2*(p2 + lp*lc2 + lc2*lc2);
    // const long double aux_13 =  gauss_c1*t2*(p2 + lp*lc1 + lc1*lc1);
    terms[8] = 0.5 * sqrt_api * lp * 3 * erf_c2;
    terms[9] = -0.5 * sqrt_api * lp * 3 * erf_c1;
    terms[10] = 0.5 * sqrt_api * lp * 2*t2*p2 * erf_c2;
    terms[11] = -0.5 * sqrt_api * lp * 2*t2*p2 * erf_c1;
    // const long double aux_2 = 0.5 * sqrt_api * lp * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[3] / (2*t4) * sum;
    //cout << "  3 " << term_3 << endl;
    if(p_deg < 4) return integral;
  }

  {
    vector<long double> terms(18);
    terms[0] = -gauss_c2* 5*lp;
    terms[1] = -gauss_c2* 2*t2*p3;
    terms[2] = -gauss_c2* 3*lc2;
    terms[3] = -gauss_c2* 2*t2*p2*lc2;
    terms[4] = -gauss_c2* 2*t2*lp*lc2*lc2;
    terms[5] = -gauss_c2* 2*t2*lc2*lc2*lc2;
    // const long double aux_11 = -gauss_c2*(5*lp + 2*t2*p3 + 3*lc2 + 2*t2*p2*lc2 + 2*t2*lp*lc2*lc2 + 2*t2*lc2*lc2*lc2);
    terms[6] =  gauss_c1* 5*lp;
    terms[7] =  gauss_c1* 2*t2*p3;
    terms[8] =  gauss_c1* 3*lc1;
    terms[9] =  gauss_c1* 2*t2*p2*lc1;
    terms[10] = gauss_c1* 2*t2*lp*lc1*lc1;
    terms[11] = gauss_c1* 2*t2*lc1*lc1*lc1;
    // const long double aux_12 = gauss_c1*(5*lp + 2*t2*p3 + 3*lc1 + 2*t2*p2*lc1 + 2*t2*lp*lc1*lc1 + 2*t2*lc1*lc1*lc1);
    terms[12] = -sqrt_pi/(2*lt) * 3 * erf_c1;
    terms[13] = sqrt_pi/(2*lt) * 3 * erf_c2;
    terms[14] = -sqrt_pi/(lt) * 6*t2*p2 * erf_c1;
    terms[15] = sqrt_pi/(lt) * 6*t2*p2 * erf_c2;
    terms[16] = -sqrt_pi/(lt) * 2*t4*p4 * erf_c1;
    terms[17] = sqrt_pi/(lt) * 2*t4*p4 * erf_c2;
    // const long double aux_2 = sqrt_pi/(2*lt) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[4]/(4*t4) * sum;
    //cout << "  4 " << term_4 << endl;
    if(p_deg < 5) return integral;
  }

  {
    vector<long double> terms(24);
    terms[0] = -14*t2*lp*gauss_c2*lc2;
    terms[1] = 14*t2*lp*gauss_c1*lc1;
    terms[2] = -18*t2*p2*gauss_c2;
    terms[3] = 18*t2*p2*gauss_c1;
    terms[4] = -4*t4*gauss_c2*lc2*lc2*lc2*lc2;
    terms[5] = 4*t4*gauss_c1*lc1*lc1*lc1*lc1;
    terms[6] = -4*t4*lp*gauss_c2*lc2*lc2*lc2;
    terms[7] = 4*t4*lp*gauss_c1*lc1*lc1*lc1;
    terms[8] = -4*t4*p2*gauss_c2*lc2*lc2;
    terms[9] = 4*t4*p2*gauss_c1*lc1*lc1;
    terms[10] = -4*t4*p3*gauss_c2*lc2;
    terms[11] = 4*t4*p3*gauss_c1*lc1;
    terms[12] = -4*t4*p4*gauss_c2;
    terms[13] = 4*t4*p4*gauss_c1;
    terms[14] = -8*gauss_c2;
    terms[15] = 8*gauss_c1;
    terms[16] = -8*t2*gauss_c2*lc2*lc2;
    terms[17] = 8*t2*gauss_c1*lc1*lc1;
    terms[18] = sqrt_api * lp * 15.0 * erf_c2;
    terms[19] = -sqrt_api * lp * 15.0 * erf_c1;
    terms[20] = sqrt_api * lp * 20*p2*t2 * erf_c2;
    terms[21] = -sqrt_api * lp * 20*p2*t2 * erf_c1;
    terms[22] = sqrt_api * lp * 4*p4*t4 * erf_c2;
    terms[23] = -sqrt_api * lp * 4*p4*t4 * erf_c1;
    // const long double aux_2 = sqrt_api * lp * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[5]/(8*t6) * sum;
  // cout << "  5 " << term_5 << endl;
    if(p_deg < 6) return integral;
  }

  {
    vector<long double> terms(32);
    terms[0] = -gauss_c2* 4*t4*p5;
    terms[1] = -gauss_c2* 4*t4*p4*lc2;
    terms[2] = -gauss_c2* 24*t2*p2*lc2;
    terms[3] = -gauss_c2* 4*t2*p2*lc2* t2*lc2*lc2;
    terms[4] = -gauss_c2* 28*t2*p3;
    terms[5] = -gauss_c2* 4*t2*p3* t2*lc2*lc2;
    terms[6] = -gauss_c2* lc2*15;
    terms[7] = -gauss_c2* lc2* 10*t2*lc2*lc2;
    terms[8] = -gauss_c2* lc2* 4*t4*lc2*lc2*lc2*lc2;
    terms[9] = -gauss_c2* lp*33;
    terms[10] = -gauss_c2* lp*18*t2*lc2*lc2;
    terms[11] = -gauss_c2* lp* 4*t4*lc2*lc2*lc2*lc2;
    // const long double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*lc2 + 4*t2*p2*lc2*(6 + t2*lc2*lc2) + 4*t2*p3*(7 + t2*lc2*lc2) + lc2*(15 + 10*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2) + lp*(33 + 18*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2));

    terms[12] = 4* gauss_c1* t4*p5;
    terms[13] = 4* gauss_c1* t4*p4*lc1;
    terms[14] = 24*gauss_c1* t2*p2*lc1;
    terms[15] = 4* gauss_c1* t2*p2*lc1* t2*lc1*lc1;
    terms[16] = 28*gauss_c1* t2*p3;
    terms[17] = 4* gauss_c1* t2*p3* t2*lc1*lc1;
    terms[18] = 15*gauss_c1* lc1;
    terms[19] = 10*gauss_c1* lc1* t2*lc1*lc1;
    terms[20] = 4* gauss_c1* lc1* t4*lc1*lc1*lc1*lc1;
    terms[21] = 33*gauss_c1* lp;
    terms[22] = 18*gauss_c1* lp*t2*lc1*lc1;
    terms[23] = 4* gauss_c1* lp* t4*lc1*lc1*lc1*lc1;
    // const long double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*lc1 + 4*t2*p2*lc1*(6 + t2*lc1*lc1) + 4*t2*p3*(7 + t2*lc1*lc1) + lc1*(15 + 10*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1) + lp*(33 + 18*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1));

    terms[24] =  sqrt_pi/(2*lt) * 15 * (erf_c2);
    terms[25] = -sqrt_pi/(2*lt) * 15 * (erf_c1);
    terms[26] =  sqrt_pi/(lt) * (45*t2*p2) * (erf_c2);
    terms[27] = -sqrt_pi/(lt) * (45*t2*p2) * (erf_c1);
    terms[28] =  sqrt_pi/(lt) * (30*t4*p4) * (erf_c2);
    terms[29] = -sqrt_pi/(lt) * (30*t4*p4) * (erf_c1);
    terms[30] =  sqrt_pi/(lt) * (4*t6*p6) * (erf_c2);
    terms[31] = -sqrt_pi/(lt) * (4*t6*p6) * (erf_c1);
    // const long double aux_2 = sqrt_pi/(2*lt) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[6]/(8*t6) * sum;
    // cout << "  6 " << term_6 << endl;
    if(p_deg < 7) return integral;
  }

  {
    vector<long double> terms(40);
    terms[0] = -2*gauss_c2* 24;
    terms[1] = -2*gauss_c2* 3*t2*29*p2;
    terms[2] = -2*gauss_c2* 3*t2* 19*lp*lc2;
    terms[3] = -2*gauss_c2* 3*t2* 8*lc2*lc2;
    terms[4] = -2*gauss_c2* 2*t4*20*p4;
    terms[5] = -2*gauss_c2* 2*t4* 18*p3*lc2;
    terms[6] = -2*gauss_c2* 2*t4* 15*p2*lc2*lc2;
    terms[7] = -2*gauss_c2* 2*t4* 11*lp*lc2*lc2*lc2;
    terms[8] = -2*gauss_c2* 2*t4* 6*lc2*lc2*lc2*lc2;
    terms[9] = -2*gauss_c2* 4*t6*p6;
    terms[10] = -2*gauss_c2* 4*t6* p5*lc2;
    terms[11] = -2*gauss_c2* 4*t6* p4*lc2*lc2;
    terms[12] = -2*gauss_c2* 4*t6* p3*lc2*lc2*lc2;
    terms[13] = -2*gauss_c2* 4*t6* p2*lc2*lc2*lc2*lc2;
    terms[14] = -2*gauss_c2* 4*t6* lp*lc2*lc2*lc2*lc2*lc2;
    terms[15] = -2*gauss_c2* 4*t6* lc2*lc2*lc2*lc2*lc2*lc2;
    // const long double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*lp*lc2 + 8*lc2*lc2) + 2*t4*(20*p4 + 18*p3*lc2 + 15*p2*lc2*lc2 + 11*lp*lc2*lc2*lc2 + 6*lc2*lc2*lc2*lc2) + 4*t6*(p6 + p5*lc2 + p4*lc2*lc2 + p3*lc2*lc2*lc2 + p2*lc2*lc2*lc2*lc2 + lp*lc2*lc2*lc2*lc2*lc2 + lc2*lc2*lc2*lc2*lc2*lc2));

    terms[16] = 2*gauss_c1* 24;
    terms[17] = 2*gauss_c1* 3*t2*29*p2;
    terms[18] = 2*gauss_c1* 3*t2* 19*lp*lc1;
    terms[19] = 2*gauss_c1* 3*t2* 8*lc1*lc1;
    terms[20] = 2*gauss_c1* 2*t4*20*p4;
    terms[21] = 2*gauss_c1* 2*t4* 18*p3*lc1;
    terms[22] = 2*gauss_c1* 2*t4* 15*p2*lc1*lc1;
    terms[23] = 2*gauss_c1* 2*t4* 11*lp*lc1*lc1*lc1;
    terms[24] = 2*gauss_c1* 2*t4* 6*lc1*lc1*lc1*lc1;
    terms[25] = 2*gauss_c1* 4*t6*p6;
    terms[26] = 2*gauss_c1* 4*t6* p5*lc1;
    terms[27] = 2*gauss_c1* 4*t6* p4*lc1*lc1;
    terms[28] = 2*gauss_c1* 4*t6* p3*lc1*lc1*lc1;
    terms[29] = 2*gauss_c1* 4*t6* p2*lc1*lc1*lc1*lc1;
    terms[30] = 2*gauss_c1* 4*t6* lp*lc1*lc1*lc1*lc1*lc1;
    terms[31] = 2*gauss_c1* 4*t6* lc1*lc1*lc1*lc1*lc1*lc1;
    // const long double aux_12 = 2*gauss_c1*(24 + 3*t2*(29*p2 + 19*lp*lc1 + 8*lc1*lc1) + 2*t4*(20*p4 + 18*p3*lc1 + 15*p2*lc1*lc1 + 11*lp*lc1*lc1*lc1 + 6*lc1*lc1*lc1*lc1) + 4*t6*(p6 + p5*lc1 + p4*lc1*lc1 + p3*lc1*lc1*lc1 + p2*lc1*lc1*lc1*lc1 + lp*lc1*lc1*lc1*lc1*lc1 + lc1*lc1*lc1*lc1*lc1*lc1));

    terms[32] = -sqrt_api * lp * (105) * (erf_c1);
    terms[33] =  sqrt_api * lp * (105) * (erf_c2);
    terms[34] = -sqrt_api * lp * (210*t2*p2) * (erf_c1);
    terms[35] =  sqrt_api * lp * (210*t2*p2) * (erf_c2);
    terms[36] = -sqrt_api * lp * (84*t4*p4) * (erf_c1);
    terms[37] =  sqrt_api * lp * (84*t4*p4) * (erf_c2);
    terms[38] = -sqrt_api * lp * (8*t6*p6) * (erf_c1);
    terms[39] =  sqrt_api * lp * (8*t6*p6) * (erf_c2);
    // const long double aux_2 = sqrt_api * lp * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const long double sum = sort_asc_and_accumulate(terms);
    integral += poly[7]/(16*t8) * sum;
    // cout << "  7 " << term_7 << endl;
    if(p_deg < 7) return integral;
  }

  assert(false);
}



// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_64bit_desc(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(double)>=8);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const double lt = t;
  const double t2 = lt*lt;
  const double t3 = t2*lt;
  const double t4 = t2*t2;
  const double t6 = t3*t3;
  const double t8 = t4*t4;
  const double lp = p;
  const double p2 = lp*lp;
  const double p3 = p2*lp;
  const double p4 = p2*p2;
  const double p5 = p3*p2;
  const double p6 = p4*p2;
  const double lc1 = c1;
  const double lc2 = c2;
  const double gauss_c1 = exp(-t2*(lp-lc1)*(lp-lc1));
  const double gauss_c2 = exp(-t2*(lp-lc2)*(lp-lc2));
  const double sqrt_api = lt*sqrt(M_PIl);
  const double sqrt_pi  = sqrt(M_PIl);
  const double erf_c1   = erf(lt*(lc1-lp));
  const double erf_c2   = erf(lt*(lc2-lp));

//  cout << t << "," << p << "," << lc1 << "," << lc2 << poly << endl;

  {
    const double term_0 = poly[0] * sqrt(M_PIl)/(2*lt) * (erf_c2 - erf_c1);
    integral += term_0;
    //cout << "  0 " << term_0 << endl;
    if(p_deg < 1) return integral;
  }

  {
    vector<double> terms(4);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    terms[2] = sqrt_api*lp*erf_c2;
    terms[3] = -sqrt_api*lp*erf_c1;
    const double sum = sort_desc_and_accumulate(terms);
// cout << setprecision(15) << terms << endl;
    integral += poly[1]/(2*t2) * sum;
    //cout << "  1 " << term_1 << endl;
    if(p_deg < 2) return integral;
  }

  {
    vector<double> terms(8);
    terms[0] = -1/(2*t2) * gauss_c2*lp;
    terms[1] = -1/(2*t2) * gauss_c2*lc2;
    terms[2] = 1/(2*t2) * gauss_c1*lp;
    terms[3] = 1/(2*t2) * gauss_c1*lc1;
    terms[4] = sqrt_pi/(4*t3) * erf_c2;
    terms[5] = - sqrt_pi/(4*t3) * erf_c1;
    terms[6] = sqrt_pi/(4*t3) * 2*t2*p2 * erf_c2;
    terms[7] = - sqrt_pi/(4*t3) *  2*t2*p2 * erf_c1;
    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[2] * sum;
    //cout << "  2 " << term_2 << endl;
    if(p_deg < 3) return integral;
  }

  {
    vector<double> terms(12);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    // const double aux_11 = -gauss_c2 + gauss_c1;
    terms[2] = -gauss_c2*t2*p2;
    terms[3] = gauss_c1*t2*p2;
    // const double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    terms[4] = -gauss_c2*lc2 *t2*lp;
    terms[5] = gauss_c1*lc1*t2*lp;
    terms[6] = -gauss_c2*lc2*lc2*t2;
    terms[7] = gauss_c1*lc1*lc1*t2;
    // const double aux_13 = (-gauss_c2*lc2 + gauss_c1*lc1)*t2*lp;
    // const double aux_14 = (-gauss_c2*lc2*lc2  + gauss_c1*lc1*lc1)*t2;
    // const double aux_12 = -gauss_c2*t2*(p2 + lp*lc2 + lc2*lc2);
    // const double aux_13 =  gauss_c1*t2*(p2 + lp*lc1 + lc1*lc1);
    terms[8] = 0.5 * sqrt_api * lp * 3 * erf_c2;
    terms[9] = -0.5 * sqrt_api * lp * 3 * erf_c1;
    terms[10] = 0.5 * sqrt_api * lp * 2*t2*p2 * erf_c2;
    terms[11] = -0.5 * sqrt_api * lp * 2*t2*p2 * erf_c1;
    // const double aux_2 = 0.5 * sqrt_api * lp * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[3] / (2*t4) * sum;
    //cout << "  3 " << term_3 << endl;
    if(p_deg < 4) return integral;
  }

  {
    vector<double> terms(18);
    terms[0] = -gauss_c2* 5*lp;
    terms[1] = -gauss_c2* 2*t2*p3;
    terms[2] = -gauss_c2* 3*lc2;
    terms[3] = -gauss_c2* 2*t2*p2*lc2;
    terms[4] = -gauss_c2* 2*t2*lp*lc2*lc2;
    terms[5] = -gauss_c2* 2*t2*lc2*lc2*lc2;
    // const double aux_11 = -gauss_c2*(5*lp + 2*t2*p3 + 3*lc2 + 2*t2*p2*lc2 + 2*t2*lp*lc2*lc2 + 2*t2*lc2*lc2*lc2);
    terms[6] =  gauss_c1* 5*lp;
    terms[7] =  gauss_c1* 2*t2*p3;
    terms[8] =  gauss_c1* 3*lc1;
    terms[9] =  gauss_c1* 2*t2*p2*lc1;
    terms[10] = gauss_c1* 2*t2*lp*lc1*lc1;
    terms[11] = gauss_c1* 2*t2*lc1*lc1*lc1;
    // const double aux_12 = gauss_c1*(5*lp + 2*t2*p3 + 3*lc1 + 2*t2*p2*lc1 + 2*t2*lp*lc1*lc1 + 2*t2*lc1*lc1*lc1);
    terms[12] = -sqrt_pi/(2*lt) * 3 * erf_c1;
    terms[13] = sqrt_pi/(2*lt) * 3 * erf_c2;
    terms[14] = -sqrt_pi/(lt) * 6*t2*p2 * erf_c1;
    terms[15] = sqrt_pi/(lt) * 6*t2*p2 * erf_c2;
    terms[16] = -sqrt_pi/(lt) * 2*t4*p4 * erf_c1;
    terms[17] = sqrt_pi/(lt) * 2*t4*p4 * erf_c2;
    // const double aux_2 = sqrt_pi/(2*lt) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[4]/(4*t4) * sum;
    //cout << "  4 " << term_4 << endl;
    if(p_deg < 5) return integral;
  }

  {
    vector<double> terms(24);
    terms[0] = -14*t2*lp*gauss_c2*lc2;
    terms[1] = 14*t2*lp*gauss_c1*lc1;
    terms[2] = -18*t2*p2*gauss_c2;
    terms[3] = 18*t2*p2*gauss_c1;
    terms[4] = -4*t4*gauss_c2*lc2*lc2*lc2*lc2;
    terms[5] = 4*t4*gauss_c1*lc1*lc1*lc1*lc1;
    terms[6] = -4*t4*lp*gauss_c2*lc2*lc2*lc2;
    terms[7] = 4*t4*lp*gauss_c1*lc1*lc1*lc1;
    terms[8] = -4*t4*p2*gauss_c2*lc2*lc2;
    terms[9] = 4*t4*p2*gauss_c1*lc1*lc1;
    terms[10] = -4*t4*p3*gauss_c2*lc2;
    terms[11] = 4*t4*p3*gauss_c1*lc1;
    terms[12] = -4*t4*p4*gauss_c2;
    terms[13] = 4*t4*p4*gauss_c1;
    terms[14] = -8*gauss_c2;
    terms[15] = 8*gauss_c1;
    terms[16] = -8*t2*gauss_c2*lc2*lc2;
    terms[17] = 8*t2*gauss_c1*lc1*lc1;
    terms[18] = sqrt_api * lp * 15.0 * erf_c2;
    terms[19] = -sqrt_api * lp * 15.0 * erf_c1;
    terms[20] = sqrt_api * lp * 20*p2*t2 * erf_c2;
    terms[21] = -sqrt_api * lp * 20*p2*t2 * erf_c1;
    terms[22] = sqrt_api * lp * 4*p4*t4 * erf_c2;
    terms[23] = -sqrt_api * lp * 4*p4*t4 * erf_c1;
    // const double aux_2 = sqrt_api * lp * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[5]/(8*t6) * sum;
  // cout << "  5 " << term_5 << endl;
    if(p_deg < 6) return integral;
  }

  {
    vector<double> terms(32);
    terms[0] = -gauss_c2* 4*t4*p5;
    terms[1] = -gauss_c2* 4*t4*p4*lc2;
    terms[2] = -gauss_c2* 24*t2*p2*lc2;
    terms[3] = -gauss_c2* 4*t2*p2*lc2* t2*lc2*lc2;
    terms[4] = -gauss_c2* 28*t2*p3;
    terms[5] = -gauss_c2* 4*t2*p3* t2*lc2*lc2;
    terms[6] = -gauss_c2* lc2*15;
    terms[7] = -gauss_c2* lc2* 10*t2*lc2*lc2;
    terms[8] = -gauss_c2* lc2* 4*t4*lc2*lc2*lc2*lc2;
    terms[9] = -gauss_c2* lp*33;
    terms[10] = -gauss_c2* lp*18*t2*lc2*lc2;
    terms[11] = -gauss_c2* lp* 4*t4*lc2*lc2*lc2*lc2;
    // const double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*lc2 + 4*t2*p2*lc2*(6 + t2*lc2*lc2) + 4*t2*p3*(7 + t2*lc2*lc2) + lc2*(15 + 10*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2) + lp*(33 + 18*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2));

    terms[12] = 4* gauss_c1* t4*p5;
    terms[13] = 4* gauss_c1* t4*p4*lc1;
    terms[14] = 24*gauss_c1* t2*p2*lc1;
    terms[15] = 4* gauss_c1* t2*p2*lc1* t2*lc1*lc1;
    terms[16] = 28*gauss_c1* t2*p3;
    terms[17] = 4* gauss_c1* t2*p3* t2*lc1*lc1;
    terms[18] = 15*gauss_c1* lc1;
    terms[19] = 10*gauss_c1* lc1* t2*lc1*lc1;
    terms[20] = 4* gauss_c1* lc1* t4*lc1*lc1*lc1*lc1;
    terms[21] = 33*gauss_c1* lp;
    terms[22] = 18*gauss_c1* lp*t2*lc1*lc1;
    terms[23] = 4* gauss_c1* lp* t4*lc1*lc1*lc1*lc1;
    // const double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*lc1 + 4*t2*p2*lc1*(6 + t2*lc1*lc1) + 4*t2*p3*(7 + t2*lc1*lc1) + lc1*(15 + 10*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1) + lp*(33 + 18*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1));

    terms[24] =  sqrt_pi/(2*lt) * 15 * (erf_c2);
    terms[25] = -sqrt_pi/(2*lt) * 15 * (erf_c1);
    terms[26] =  sqrt_pi/(lt) * (45*t2*p2) * (erf_c2);
    terms[27] = -sqrt_pi/(lt) * (45*t2*p2) * (erf_c1);
    terms[28] =  sqrt_pi/(lt) * (30*t4*p4) * (erf_c2);
    terms[29] = -sqrt_pi/(lt) * (30*t4*p4) * (erf_c1);
    terms[30] =  sqrt_pi/(lt) * (4*t6*p6) * (erf_c2);
    terms[31] = -sqrt_pi/(lt) * (4*t6*p6) * (erf_c1);
    // const double aux_2 = sqrt_pi/(2*lt) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[6]/(8*t6) * sum;
    // cout << "  6 " << term_6 << endl;
    if(p_deg < 7) return integral;
  }

  {
    vector<double> terms(40);
    terms[0] = -2*gauss_c2* 24;
    terms[1] = -2*gauss_c2* 3*t2*29*p2;
    terms[2] = -2*gauss_c2* 3*t2* 19*lp*lc2;
    terms[3] = -2*gauss_c2* 3*t2* 8*lc2*lc2;
    terms[4] = -2*gauss_c2* 2*t4*20*p4;
    terms[5] = -2*gauss_c2* 2*t4* 18*p3*lc2;
    terms[6] = -2*gauss_c2* 2*t4* 15*p2*lc2*lc2;
    terms[7] = -2*gauss_c2* 2*t4* 11*lp*lc2*lc2*lc2;
    terms[8] = -2*gauss_c2* 2*t4* 6*lc2*lc2*lc2*lc2;
    terms[9] = -2*gauss_c2* 4*t6*p6;
    terms[10] = -2*gauss_c2* 4*t6* p5*lc2;
    terms[11] = -2*gauss_c2* 4*t6* p4*lc2*lc2;
    terms[12] = -2*gauss_c2* 4*t6* p3*lc2*lc2*lc2;
    terms[13] = -2*gauss_c2* 4*t6* p2*lc2*lc2*lc2*lc2;
    terms[14] = -2*gauss_c2* 4*t6* lp*lc2*lc2*lc2*lc2*lc2;
    terms[15] = -2*gauss_c2* 4*t6* lc2*lc2*lc2*lc2*lc2*lc2;
    // const double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*lp*lc2 + 8*lc2*lc2) + 2*t4*(20*p4 + 18*p3*lc2 + 15*p2*lc2*lc2 + 11*lp*lc2*lc2*lc2 + 6*lc2*lc2*lc2*lc2) + 4*t6*(p6 + p5*lc2 + p4*lc2*lc2 + p3*lc2*lc2*lc2 + p2*lc2*lc2*lc2*lc2 + lp*lc2*lc2*lc2*lc2*lc2 + lc2*lc2*lc2*lc2*lc2*lc2));

    terms[16] = 2*gauss_c1* 24;
    terms[17] = 2*gauss_c1* 3*t2*29*p2;
    terms[18] = 2*gauss_c1* 3*t2* 19*lp*lc1;
    terms[19] = 2*gauss_c1* 3*t2* 8*lc1*lc1;
    terms[20] = 2*gauss_c1* 2*t4*20*p4;
    terms[21] = 2*gauss_c1* 2*t4* 18*p3*lc1;
    terms[22] = 2*gauss_c1* 2*t4* 15*p2*lc1*lc1;
    terms[23] = 2*gauss_c1* 2*t4* 11*lp*lc1*lc1*lc1;
    terms[24] = 2*gauss_c1* 2*t4* 6*lc1*lc1*lc1*lc1;
    terms[25] = 2*gauss_c1* 4*t6*p6;
    terms[26] = 2*gauss_c1* 4*t6* p5*lc1;
    terms[27] = 2*gauss_c1* 4*t6* p4*lc1*lc1;
    terms[28] = 2*gauss_c1* 4*t6* p3*lc1*lc1*lc1;
    terms[29] = 2*gauss_c1* 4*t6* p2*lc1*lc1*lc1*lc1;
    terms[30] = 2*gauss_c1* 4*t6* lp*lc1*lc1*lc1*lc1*lc1;
    terms[31] = 2*gauss_c1* 4*t6* lc1*lc1*lc1*lc1*lc1*lc1;
    // const double aux_12 = 2*gauss_c1*(24 + 3*t2*(29*p2 + 19*lp*lc1 + 8*lc1*lc1) + 2*t4*(20*p4 + 18*p3*lc1 + 15*p2*lc1*lc1 + 11*lp*lc1*lc1*lc1 + 6*lc1*lc1*lc1*lc1) + 4*t6*(p6 + p5*lc1 + p4*lc1*lc1 + p3*lc1*lc1*lc1 + p2*lc1*lc1*lc1*lc1 + lp*lc1*lc1*lc1*lc1*lc1 + lc1*lc1*lc1*lc1*lc1*lc1));

    terms[32] = -sqrt_api * lp * (105) * (erf_c1);
    terms[33] =  sqrt_api * lp * (105) * (erf_c2);
    terms[34] = -sqrt_api * lp * (210*t2*p2) * (erf_c1);
    terms[35] =  sqrt_api * lp * (210*t2*p2) * (erf_c2);
    terms[36] = -sqrt_api * lp * (84*t4*p4) * (erf_c1);
    terms[37] =  sqrt_api * lp * (84*t4*p4) * (erf_c2);
    terms[38] = -sqrt_api * lp * (8*t6*p6) * (erf_c1);
    terms[39] =  sqrt_api * lp * (8*t6*p6) * (erf_c2);
    // const double aux_2 = sqrt_api * lp * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const double sum = sort_desc_and_accumulate(terms);
    integral += poly[7]/(16*t8) * sum;
    // cout << "  7 " << term_7 << endl;
    if(p_deg < 7) return integral;
  }

  assert(false);
}


// int_c1^c2 e^(-t^2(p-x)^2) * poly dx
double integral_gauss_polynomial_64bit_asc(const double t, const double p, const double c1, const double c2, const polynomial& poly){
  assert(sizeof(double)>=8);
  double integral = 0;

  const size_t p_deg = poly.size() - 1;
  const double lt = t;
  const double t2 = lt*lt;
  const double t3 = t2*lt;
  const double t4 = t2*t2;
  const double t6 = t3*t3;
  const double t8 = t4*t4;
  const double lp = p;
  const double p2 = lp*lp;
  const double p3 = p2*lp;
  const double p4 = p2*p2;
  const double p5 = p3*p2;
  const double p6 = p4*p2;
  const double lc1 = c1;
  const double lc2 = c2;
  const double gauss_c1 = exp(-t2*(lp-lc1)*(lp-lc1));
  const double gauss_c2 = exp(-t2*(lp-lc2)*(lp-lc2));
  const double sqrt_api = lt*sqrt(M_PIl);
  const double sqrt_pi  = sqrt(M_PIl);
  const double erf_c1   = erf(lt*(lc1-lp));
  const double erf_c2   = erf(lt*(lc2-lp));

//  cout << t << "," << p << "," << lc1 << "," << lc2 << poly << endl;

  {
    const double term_0 = poly[0] * sqrt(M_PIl)/(2*lt) * (erf_c2 - erf_c1);
    integral += term_0;
    //cout << "  0 " << term_0 << endl;
    if(p_deg < 1) return integral;
  }

  {
    vector<double> terms(4);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    terms[2] = sqrt_api*lp*erf_c2;
    terms[3] = -sqrt_api*lp*erf_c1;
    const double sum = sort_asc_and_accumulate(terms);
// cout << setprecision(15) << terms << endl;
    integral += poly[1]/(2*t2) * sum;
    //cout << "  1 " << term_1 << endl;
    if(p_deg < 2) return integral;
  }

  {
    vector<double> terms(8);
    terms[0] = -1/(2*t2) * gauss_c2*lp;
    terms[1] = -1/(2*t2) * gauss_c2*lc2;
    terms[2] = 1/(2*t2) * gauss_c1*lp;
    terms[3] = 1/(2*t2) * gauss_c1*lc1;
    terms[4] = sqrt_pi/(4*t3) * erf_c2;
    terms[5] = - sqrt_pi/(4*t3) * erf_c1;
    terms[6] = sqrt_pi/(4*t3) * 2*t2*p2 * erf_c2;
    terms[7] = - sqrt_pi/(4*t3) *  2*t2*p2 * erf_c1;
    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[2] * sum;
    //cout << "  2 " << term_2 << endl;
    if(p_deg < 3) return integral;
  }

  {
    vector<double> terms(12);
    terms[0] = -gauss_c2;
    terms[1] = gauss_c1;
    // const double aux_11 = -gauss_c2 + gauss_c1;
    terms[2] = -gauss_c2*t2*p2;
    terms[3] = gauss_c1*t2*p2;
    // const double aux_12 = (-gauss_c2 + gauss_c1)*t2*p2;
    terms[4] = -gauss_c2*lc2 *t2*lp;
    terms[5] = gauss_c1*lc1*t2*lp;
    terms[6] = -gauss_c2*lc2*lc2*t2;
    terms[7] = gauss_c1*lc1*lc1*t2;
    // const double aux_13 = (-gauss_c2*lc2 + gauss_c1*lc1)*t2*lp;
    // const double aux_14 = (-gauss_c2*lc2*lc2  + gauss_c1*lc1*lc1)*t2;
    // const double aux_12 = -gauss_c2*t2*(p2 + lp*lc2 + lc2*lc2);
    // const double aux_13 =  gauss_c1*t2*(p2 + lp*lc1 + lc1*lc1);
    terms[8] = 0.5 * sqrt_api * lp * 3 * erf_c2;
    terms[9] = -0.5 * sqrt_api * lp * 3 * erf_c1;
    terms[10] = 0.5 * sqrt_api * lp * 2*t2*p2 * erf_c2;
    terms[11] = -0.5 * sqrt_api * lp * 2*t2*p2 * erf_c1;
    // const double aux_2 = 0.5 * sqrt_api * lp * (3 + 2*t2*p2) * (erf_c2 - erf_c1);
    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[3] / (2*t4) * sum;
    //cout << "  3 " << term_3 << endl;
    if(p_deg < 4) return integral;
  }

  {
    vector<double> terms(18);
    terms[0] = -gauss_c2* 5*lp;
    terms[1] = -gauss_c2* 2*t2*p3;
    terms[2] = -gauss_c2* 3*lc2;
    terms[3] = -gauss_c2* 2*t2*p2*lc2;
    terms[4] = -gauss_c2* 2*t2*lp*lc2*lc2;
    terms[5] = -gauss_c2* 2*t2*lc2*lc2*lc2;
    // const double aux_11 = -gauss_c2*(5*lp + 2*t2*p3 + 3*lc2 + 2*t2*p2*lc2 + 2*t2*lp*lc2*lc2 + 2*t2*lc2*lc2*lc2);
    terms[6] =  gauss_c1* 5*lp;
    terms[7] =  gauss_c1* 2*t2*p3;
    terms[8] =  gauss_c1* 3*lc1;
    terms[9] =  gauss_c1* 2*t2*p2*lc1;
    terms[10] = gauss_c1* 2*t2*lp*lc1*lc1;
    terms[11] = gauss_c1* 2*t2*lc1*lc1*lc1;
    // const double aux_12 = gauss_c1*(5*lp + 2*t2*p3 + 3*lc1 + 2*t2*p2*lc1 + 2*t2*lp*lc1*lc1 + 2*t2*lc1*lc1*lc1);
    terms[12] = -sqrt_pi/(2*lt) * 3 * erf_c1;
    terms[13] = sqrt_pi/(2*lt) * 3 * erf_c2;
    terms[14] = -sqrt_pi/(lt) * 6*t2*p2 * erf_c1;
    terms[15] = sqrt_pi/(lt) * 6*t2*p2 * erf_c2;
    terms[16] = -sqrt_pi/(lt) * 2*t4*p4 * erf_c1;
    terms[17] = sqrt_pi/(lt) * 2*t4*p4 * erf_c2;
    // const double aux_2 = sqrt_pi/(2*lt) * (3 + 12*t2*p2 + 4*t4*p4) * (erf_c2 - erf_c1);
    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[4]/(4*t4) * sum;
    //cout << "  4 " << term_4 << endl;
    if(p_deg < 5) return integral;
  }

  {
    vector<double> terms(24);
    terms[0] = -14*t2*lp*gauss_c2*lc2;
    terms[1] = 14*t2*lp*gauss_c1*lc1;
    terms[2] = -18*t2*p2*gauss_c2;
    terms[3] = 18*t2*p2*gauss_c1;
    terms[4] = -4*t4*gauss_c2*lc2*lc2*lc2*lc2;
    terms[5] = 4*t4*gauss_c1*lc1*lc1*lc1*lc1;
    terms[6] = -4*t4*lp*gauss_c2*lc2*lc2*lc2;
    terms[7] = 4*t4*lp*gauss_c1*lc1*lc1*lc1;
    terms[8] = -4*t4*p2*gauss_c2*lc2*lc2;
    terms[9] = 4*t4*p2*gauss_c1*lc1*lc1;
    terms[10] = -4*t4*p3*gauss_c2*lc2;
    terms[11] = 4*t4*p3*gauss_c1*lc1;
    terms[12] = -4*t4*p4*gauss_c2;
    terms[13] = 4*t4*p4*gauss_c1;
    terms[14] = -8*gauss_c2;
    terms[15] = 8*gauss_c1;
    terms[16] = -8*t2*gauss_c2*lc2*lc2;
    terms[17] = 8*t2*gauss_c1*lc1*lc1;
    terms[18] = sqrt_api * lp * 15.0 * erf_c2;
    terms[19] = -sqrt_api * lp * 15.0 * erf_c1;
    terms[20] = sqrt_api * lp * 20*p2*t2 * erf_c2;
    terms[21] = -sqrt_api * lp * 20*p2*t2 * erf_c1;
    terms[22] = sqrt_api * lp * 4*p4*t4 * erf_c2;
    terms[23] = -sqrt_api * lp * 4*p4*t4 * erf_c1;
    // const double aux_2 = sqrt_api * lp * (15.0 + 20*p2*t2 + 4*p4*t4) * (erf_c2 - erf_c1);
    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[5]/(8*t6) * sum;
  // cout << "  5 " << term_5 << endl;
    if(p_deg < 6) return integral;
  }

  {
    vector<double> terms(32);
    terms[0] = -gauss_c2* 4*t4*p5;
    terms[1] = -gauss_c2* 4*t4*p4*lc2;
    terms[2] = -gauss_c2* 24*t2*p2*lc2;
    terms[3] = -gauss_c2* 4*t2*p2*lc2* t2*lc2*lc2;
    terms[4] = -gauss_c2* 28*t2*p3;
    terms[5] = -gauss_c2* 4*t2*p3* t2*lc2*lc2;
    terms[6] = -gauss_c2* lc2*15;
    terms[7] = -gauss_c2* lc2* 10*t2*lc2*lc2;
    terms[8] = -gauss_c2* lc2* 4*t4*lc2*lc2*lc2*lc2;
    terms[9] = -gauss_c2* lp*33;
    terms[10] = -gauss_c2* lp*18*t2*lc2*lc2;
    terms[11] = -gauss_c2* lp* 4*t4*lc2*lc2*lc2*lc2;
    // const double aux_11 = -gauss_c2*(4*t4*p5 + 4*t4*p4*lc2 + 4*t2*p2*lc2*(6 + t2*lc2*lc2) + 4*t2*p3*(7 + t2*lc2*lc2) + lc2*(15 + 10*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2) + lp*(33 + 18*t2*lc2*lc2 + 4*t4*lc2*lc2*lc2*lc2));

    terms[12] = 4* gauss_c1* t4*p5;
    terms[13] = 4* gauss_c1* t4*p4*lc1;
    terms[14] = 24*gauss_c1* t2*p2*lc1;
    terms[15] = 4* gauss_c1* t2*p2*lc1* t2*lc1*lc1;
    terms[16] = 28*gauss_c1* t2*p3;
    terms[17] = 4* gauss_c1* t2*p3* t2*lc1*lc1;
    terms[18] = 15*gauss_c1* lc1;
    terms[19] = 10*gauss_c1* lc1* t2*lc1*lc1;
    terms[20] = 4* gauss_c1* lc1* t4*lc1*lc1*lc1*lc1;
    terms[21] = 33*gauss_c1* lp;
    terms[22] = 18*gauss_c1* lp*t2*lc1*lc1;
    terms[23] = 4* gauss_c1* lp* t4*lc1*lc1*lc1*lc1;
    // const double aux_12 = -gauss_c1*(4*t4*p5 + 4*t4*p4*lc1 + 4*t2*p2*lc1*(6 + t2*lc1*lc1) + 4*t2*p3*(7 + t2*lc1*lc1) + lc1*(15 + 10*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1) + lp*(33 + 18*t2*lc1*lc1 + 4*t4*lc1*lc1*lc1*lc1));

    terms[24] =  sqrt_pi/(2*lt) * 15 * (erf_c2);
    terms[25] = -sqrt_pi/(2*lt) * 15 * (erf_c1);
    terms[26] =  sqrt_pi/(lt) * (45*t2*p2) * (erf_c2);
    terms[27] = -sqrt_pi/(lt) * (45*t2*p2) * (erf_c1);
    terms[28] =  sqrt_pi/(lt) * (30*t4*p4) * (erf_c2);
    terms[29] = -sqrt_pi/(lt) * (30*t4*p4) * (erf_c1);
    terms[30] =  sqrt_pi/(lt) * (4*t6*p6) * (erf_c2);
    terms[31] = -sqrt_pi/(lt) * (4*t6*p6) * (erf_c1);
    // const double aux_2 = sqrt_pi/(2*lt) * (15 + 90*t2*p2 + 60*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[6]/(8*t6) * sum;
    // cout << "  6 " << term_6 << endl;
    if(p_deg < 7) return integral;
  }

  {
    vector<double> terms(40);
    terms[0] = -2*gauss_c2* 24;
    terms[1] = -2*gauss_c2* 3*t2*29*p2;
    terms[2] = -2*gauss_c2* 3*t2* 19*lp*lc2;
    terms[3] = -2*gauss_c2* 3*t2* 8*lc2*lc2;
    terms[4] = -2*gauss_c2* 2*t4*20*p4;
    terms[5] = -2*gauss_c2* 2*t4* 18*p3*lc2;
    terms[6] = -2*gauss_c2* 2*t4* 15*p2*lc2*lc2;
    terms[7] = -2*gauss_c2* 2*t4* 11*lp*lc2*lc2*lc2;
    terms[8] = -2*gauss_c2* 2*t4* 6*lc2*lc2*lc2*lc2;
    terms[9] = -2*gauss_c2* 4*t6*p6;
    terms[10] = -2*gauss_c2* 4*t6* p5*lc2;
    terms[11] = -2*gauss_c2* 4*t6* p4*lc2*lc2;
    terms[12] = -2*gauss_c2* 4*t6* p3*lc2*lc2*lc2;
    terms[13] = -2*gauss_c2* 4*t6* p2*lc2*lc2*lc2*lc2;
    terms[14] = -2*gauss_c2* 4*t6* lp*lc2*lc2*lc2*lc2*lc2;
    terms[15] = -2*gauss_c2* 4*t6* lc2*lc2*lc2*lc2*lc2*lc2;
    // const double aux_11 = -2*gauss_c2*(24 + 3*t2*(29*p2 + 19*lp*lc2 + 8*lc2*lc2) + 2*t4*(20*p4 + 18*p3*lc2 + 15*p2*lc2*lc2 + 11*lp*lc2*lc2*lc2 + 6*lc2*lc2*lc2*lc2) + 4*t6*(p6 + p5*lc2 + p4*lc2*lc2 + p3*lc2*lc2*lc2 + p2*lc2*lc2*lc2*lc2 + lp*lc2*lc2*lc2*lc2*lc2 + lc2*lc2*lc2*lc2*lc2*lc2));

    terms[16] = 2*gauss_c1* 24;
    terms[17] = 2*gauss_c1* 3*t2*29*p2;
    terms[18] = 2*gauss_c1* 3*t2* 19*lp*lc1;
    terms[19] = 2*gauss_c1* 3*t2* 8*lc1*lc1;
    terms[20] = 2*gauss_c1* 2*t4*20*p4;
    terms[21] = 2*gauss_c1* 2*t4* 18*p3*lc1;
    terms[22] = 2*gauss_c1* 2*t4* 15*p2*lc1*lc1;
    terms[23] = 2*gauss_c1* 2*t4* 11*lp*lc1*lc1*lc1;
    terms[24] = 2*gauss_c1* 2*t4* 6*lc1*lc1*lc1*lc1;
    terms[25] = 2*gauss_c1* 4*t6*p6;
    terms[26] = 2*gauss_c1* 4*t6* p5*lc1;
    terms[27] = 2*gauss_c1* 4*t6* p4*lc1*lc1;
    terms[28] = 2*gauss_c1* 4*t6* p3*lc1*lc1*lc1;
    terms[29] = 2*gauss_c1* 4*t6* p2*lc1*lc1*lc1*lc1;
    terms[30] = 2*gauss_c1* 4*t6* lp*lc1*lc1*lc1*lc1*lc1;
    terms[31] = 2*gauss_c1* 4*t6* lc1*lc1*lc1*lc1*lc1*lc1;
    // const double aux_12 = 2*gauss_c1*(24 + 3*t2*(29*p2 + 19*lp*lc1 + 8*lc1*lc1) + 2*t4*(20*p4 + 18*p3*lc1 + 15*p2*lc1*lc1 + 11*lp*lc1*lc1*lc1 + 6*lc1*lc1*lc1*lc1) + 4*t6*(p6 + p5*lc1 + p4*lc1*lc1 + p3*lc1*lc1*lc1 + p2*lc1*lc1*lc1*lc1 + lp*lc1*lc1*lc1*lc1*lc1 + lc1*lc1*lc1*lc1*lc1*lc1));

    terms[32] = -sqrt_api * lp * (105) * (erf_c1);
    terms[33] =  sqrt_api * lp * (105) * (erf_c2);
    terms[34] = -sqrt_api * lp * (210*t2*p2) * (erf_c1);
    terms[35] =  sqrt_api * lp * (210*t2*p2) * (erf_c2);
    terms[36] = -sqrt_api * lp * (84*t4*p4) * (erf_c1);
    terms[37] =  sqrt_api * lp * (84*t4*p4) * (erf_c2);
    terms[38] = -sqrt_api * lp * (8*t6*p6) * (erf_c1);
    terms[39] =  sqrt_api * lp * (8*t6*p6) * (erf_c2);
    // const double aux_2 = sqrt_api * lp * (105 + 210*t2*p2 + 84*t4*p4 + 8*t6*p6) * (erf_c2 - erf_c1);

    const double sum = sort_asc_and_accumulate(terms);
    integral += poly[7]/(16*t8) * sum;
    // cout << "  7 " << term_7 << endl;
    if(p_deg < 7) return integral;
  }

  assert(false);
}


// in one cell with N functions, input one occupied orbital and return N
// orthonormal functions where the first is the input orbital
tensor_2 orthonormalise_cell(const system1d& sys, const tensor_1& occ, int c, bool debug_print){
  //const double eps = 1.e-30;
  const size_t N = sys.grid.size();
  tensor_2 orthn_orbitals({sys.n_funcs_per_cell, sys.n_funcs_per_cell});

  orthn_orbitals.set_slice(0, XX, occ);
  const bool adjust_order = false;
  if(c<16 || !adjust_order){
    for(int i=0; i<N-1; i++){
      tensor_1 orb({N});
      orb[i] = 1;
      orthn_orbitals.set_slice(i+1, XX, orb); // start with orthn_orbitals which are equal to the LIP
    }
    // cout << orthn_orbitals << endl;
  }
  else{
    for(int i=0; i<N-1; i++){
      tensor_1 orb({N});
      orb[N-1-i] = 1;
      orthn_orbitals.set_slice(i+1, XX, orb); // start with orthn_orbitals which are equal to the LIP
    }
  }
  // cout << orthn_orbitals << endl;


  // orthogonalise
  // iterate over local orbitals
  for(int i=0; i<N; i++){
    tensor_1 orb_work = orthn_orbitals.get_slice(i, XX);
    // and project out all other ones with larger indices
    for(int j=i-1; j>=0; j--){
      //const size_t c=0; // irrelevant
      const tensor_1 frame_i = orthn_orbitals.get_slice(i, XX);
      const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
      // cout << frame_i << endl;
      // cout << frame_j << endl;
      const double integral_ij = integrate_cell(sys, c, frame_i, frame_j, NONE);
      const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);

      //if(integral_jj > eps) {
        orb_work -= frame_j * (integral_ij/integral_jj);
      //}
    }
    orthn_orbitals.set_slice(i, XX, orb_work);
  }
  // cout << orthn_orbitals << endl;

  // normalise
  for(int j=0; j<N; j++){
    //const size_t c=0; // irrelevant
    tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
    const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);
    const double integral_jj_sqrt = sqrt(integral_jj);
    frame_j /= integral_jj_sqrt;
    orthn_orbitals.set_slice(j, XX, frame_j);
    // for(double &d: orthn_orbitals[j]) d/= integral_jj_sqrt;
  }

//  if(c==16){
//    cout << orthn_orbitals << endl;
//    tensor_1 frame_2 = orthn_orbitals.get_slice(1, XX);
//    tensor_1 frame_3 = orthn_orbitals.get_slice(2, XX);
//    tensor_1 frame_2new( (frame_2 + frame_3)*(1/sqrt(2)) );
//    tensor_1 frame_3new( (frame_2 - frame_3)*(1/sqrt(2)) );
//    orthn_orbitals.set_slice(1, XX, frame_2new);
//    orthn_orbitals.set_slice(2, XX, frame_3new);
//    cout << orthn_orbitals << endl;
//  }

  if(debug_print){
    cout << "normalised?" << endl;
    for(int j=0; j<N; j++){
      //const size_t c=0; // irrelevant
      const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
      const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);
      cout << integral_jj << endl;
    }
  
    cout << "orthogonal?" << endl;
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        if(i==j)continue;
        //const size_t c=0; // irrelevant
        const tensor_1 frame_i = orthn_orbitals.get_slice(i, XX);
        const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
        const double integral_ij = integrate_cell(sys, c, frame_i, frame_j, NONE);
        cout << integral_ij << endl;
      }
    }
  }

  return orthn_orbitals;
}


// in one cell with N functions, input one occupied orbital and return N-1
// orthonormal functions where the first is the input orbital, but from the
// virtual space the first and last polynomial are omitted
tensor_2 orthonormalise_cell_wo_boundary(const system1d& sys, const tensor_1& occ, int c, bool debug_print){
  //const double eps = 1.e-30;
  const size_t N = sys.grid.size();
  tensor_2 orthn_orbitals({sys.n_funcs_per_cell, sys.n_funcs_per_cell}); // 1 occ + n-3 virt

  orthn_orbitals.set_slice(0, XX, occ);
  tensor_1 dummy1({N}), dummy2({N});
  dummy1[0] = 1; dummy2[N-1] = 1;
  orthn_orbitals.set_slice(1, XX, dummy1);
  orthn_orbitals.set_slice(2, XX, dummy2);
  // orthn_orbitals(0,0)=0;
  // orthn_orbitals(0,N-1)=0;
  for(int i=0; i<N-3; i++){ // n-3 virt
    tensor_1 orb({N});
    orb[i+1] = 1;
    orthn_orbitals.set_slice(i+3, XX, orb); // start with orthn_orbitals which are equal to the LIP
  }
  cout << "create: " << orthn_orbitals << endl;


  // orthogonalise
  // iterate over local orbitals
  for(int i=3; i<N; i++){ // 1 occ + n-3 virt
    tensor_1 orb_work = orthn_orbitals.get_slice(i, XX);
    // and project out all other ones with larger indices
    for(int j=i-1; j>=0; j--){
      //const size_t c=0; // irrelevant
      const tensor_1 frame_i = orthn_orbitals.get_slice(i, XX);
      const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
      // cout << frame_i << endl;
      // cout << frame_j << endl;
      const double integral_ij = integrate_cell(sys, c, frame_i, frame_j, NONE);
      const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);

      //if(integral_jj > eps) {
        orb_work -= frame_j * (integral_ij/integral_jj);
      //}
    }
    orthn_orbitals.set_slice(i, XX, orb_work);
  }
  // orthn_orbitals.set_slice(0, XX, occ);
  cout << "ortho: " << orthn_orbitals << endl;

  // normalise
  for(int j=0; j<N; j++){
    //const size_t c=0; // irrelevant
    tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
    const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);
    const double integral_jj_sqrt = sqrt(integral_jj);
    frame_j /= integral_jj_sqrt;
    orthn_orbitals.set_slice(j, XX, frame_j);
    // for(double &d: orthn_orbitals[j]) d/= integral_jj_sqrt;
  }
  cout << "normalise " << orthn_orbitals << endl;

  if(debug_print){
    cout << "normalised?" << endl;
    for(int j=0; j<N; j++){
      //const size_t c=0; // irrelevant
      const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
      const double integral_jj = integrate_cell(sys, c, frame_j, frame_j, NONE);
      cout << integral_jj << endl;
    }
  
    cout << "orthogonal?" << endl;
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        if(i==j)continue;
        //const size_t c=0; // irrelevant
        const tensor_1 frame_i = orthn_orbitals.get_slice(i, XX);
        const tensor_1 frame_j = orthn_orbitals.get_slice(j, XX);
        const double integral_ij = integrate_cell(sys, c, frame_i, frame_j, NONE);
        cout << integral_ij << endl;
      }
    }
  }

  return orthn_orbitals;
}


