#ifndef TENSOR_HH
#define TENSOR_HH

#include <cmath>
#include <cstring>
#include <cassert>
#include <fstream>  // for ofstream
#include <sstream>
#include <utility>  // for pair
#include <vector>

#include "auxiliary.hh"
#include "timing.hh"

using namespace std;


typedef enum {XX, YY, ZZ} dir_t;

class system3d;
// class system1d;


class tensor{
  vector<size_t> dims;
  size_t size; // product of dims
  double* t;

public:
  explicit tensor(const initializer_list<size_t> il): dims(il), size(product(il)), t(new double[size]) { memset(t, 0, size*sizeof(double)); }
  tensor(): dims(1), size(0), t(nullptr) {}
  tensor(const tensor& T): dims(T.get_dims()), size(T.get_size()), t(new double[size]) { memcpy(t, T.t, size*sizeof(double)); }
  tensor(initializer_list<size_t> il, double* ptr1, double* ptr2): dims(il), size(product(il)), t(new double[size]) {assert(size==ptr2-ptr1); memcpy(t, ptr1, (ptr2-ptr1)*sizeof(double));}

  ~tensor(){delete [] t;}

  tensor operator=(const tensor& T);

  // void resize(initializer_list<size_t> il) {dims=il; size(product(il)); delete [] t; t(new double[size]);  memset(t, 0, size*sizeof(double));}

  double& operator[](size_t i) {return t[i];}
  double  operator[](size_t i) const {return t[i];}

  // yes, I know, it's not const
  double* begin() const {return t;}
  double* end() const  {return t + size;}
  double* get_data_ptr(size_t d) const {return t + d;} // a pointer into the array

  vector<size_t> get_dims() const {return dims;}
  size_t get_dims_size() const {return dims.size();}
  size_t get_size() const {return size;}

  void set_dims(const initializer_list<size_t> il){dims = vector<size_t>(il); size = product(il);}
  void set_dims(const vector<size_t> v){dims = v; size = product(v);}
  void set_data(size_t sz, double* source) {delete [] t; t = new double[sz]; memcpy(t, source, sz*sizeof(double)); }

  bool operator==(const tensor & T) const;

  tensor  operator-() const;
  tensor  operator+(const tensor& T) const { return tensor(*this) += T; }
  tensor& operator+=(const tensor& T);
  tensor  operator-(const tensor& T) const { return tensor(*this) -= T; }
  tensor& operator-=(const tensor& T);
  tensor  operator*(const double d) const { return tensor(*this) *= d; }
  tensor& operator*=(const double d);
  tensor  operator/(const double d) const { return tensor(*this) /= d; }
  tensor& operator/=(const double d);
//  tensor  operator*(const tensor& T) const { return tensor(*this) *= T; }
//  tensor& operator*=(const tensor& T);

////  friend vector<double> operator*(const vector<double>&V, const matrix& M);
////  friend vector<double>& operator*=(const vector<double>&V, const matrix& M);
////  friend double dot(const vector<double>& V1, const vector<double>& V2);
////  friend matrix outer(const vector<double>& V1, const vector<double>& V2);


  friend ostream& operator<<(ostream& S, const tensor& T){
    const size_t I(T.dims.size());
    for(int i=0; i<T.get_size(); i++){
      // opening braces
      int p=1;
      for(int k=I-1; k>=0; k--){
        p *= T.dims[k];
        if(i%p == 0) S << "{";
        else break;
      }
      // the actual content
      S << setprecision(10) << T.t[i];
      // commas
      if( (i+1)%T.dims.back()!=0 ) S << ", ";
      // closing braces
      p=1;
      for(int k=I-1; k>=0; k--){
        p *= T.dims[k];
        if((i+1)%p == 0) S << "}";
        else break;
      }
      // commas
      if( (i+1)%T.dims.back()==0 && i+1!=T.get_size()) S << ", ";
    }
    return S;
  }

};



class tensor_1: public tensor{
  public:
  explicit tensor_1(const initializer_list<size_t> il): tensor(il) { assert(get_dims_size()==1); }
  tensor_1(): tensor() { }
  tensor_1(const tensor& T): tensor(T) { assert(get_dims_size()==1); }
  tensor_1(const tensor_1& T): tensor(T) { assert(get_dims_size()==1); }
  tensor_1(initializer_list<size_t> il, double* ptr1, double* ptr2): tensor(il) { assert(get_size()==ptr2-ptr1); memcpy(begin(), ptr1, (ptr2-ptr1)*sizeof(double));}
  // explicit tensor_1(const vector<double> v): tensor({v.size()}) { memcpy( begin(), v.begin(), (v.size())*sizeof(double));}

  double operator*(const tensor_1& T) const;
  // because shadowing happens despite the different signature
  using tensor::operator*;
  using tensor::operator*=;
  // using tensor::operator/;
  // using tensor::operator/=;

  double& operator()(size_t i);
  double  operator()(size_t i) const;
};

class tensor_2: public tensor{
  public:
  explicit tensor_2(const initializer_list<size_t> il): tensor(il) { assert(get_dims_size()==2); }
  tensor_2(const initializer_list<size_t> il, double* ptr1, double* ptr2): tensor_2(il) { assert(get_size()==ptr2-ptr1); memcpy(begin(), ptr1, (ptr2-ptr1)*sizeof(double));}
  tensor_2(): tensor() {}
  tensor_2(const tensor& T): tensor(T) { assert(get_dims_size()==2); }
  tensor_2(const tensor_2& T): tensor(T) { assert(get_dims_size()==2); }

  double& operator()(size_t i, size_t j);
  double  operator()(size_t i, size_t j) const;
  tensor_1 get_slice(const size_t pos, const dir_t var) const;
  void set_slice(const size_t pos, const dir_t var, const tensor_1& sl);

  tensor_2 transpose() const;
};

class tensor_square: public tensor_2{
  public:
  explicit tensor_square(const initializer_list<size_t> il): tensor_2(il) { assert(get_dims()[0]==get_dims()[1]); }
  tensor_square(const initializer_list<size_t> il, double* ptr1, double* ptr2): tensor_square(il) {
    assert(get_size()==ptr2-ptr1);
    //assert(get_dims()[0]==get_dims()[1]);
    memcpy(begin(), ptr1, (ptr2-ptr1)*sizeof(double));}
  tensor_square(): tensor_2() { }
  tensor_square(const tensor& T): tensor_2(T) { assert(get_dims()[0]==get_dims()[1]); }
  tensor_square(const tensor_2& T): tensor_2(T) { assert(get_dims()[0]==get_dims()[1]); }
  tensor_square(const tensor_square& T): tensor_2(T) { }

  tensor_1 operator*(const tensor_1& T) const;

  bool is_symmetric() const;
  bool is_point_symmetric() const;

  tensor_square inverse() const;


  tuple<tensor_1, tensor_square, tensor_square> eigensystem_gen(timing& times) const;
};

class tensor_square_sym: public tensor_square{
  public:
  explicit tensor_square_sym(const initializer_list<size_t> il): tensor_square(il) { }
  tensor_square_sym(): tensor_square() { }
  tensor_square_sym(const tensor_square& T): tensor_square(T) { assert(is_symmetric()); }
  tensor_square_sym(const tensor_square_sym& T): tensor_square(T) { }

  tensor_square_sym operator=(const tensor& T);
  tensor_square_sym operator=(const tensor_square_sym& T);

  tensor_1 operator*(const tensor_1& T) const;
  tensor_square_sym operator*(const double d) const;

  tuple<tensor_1, tensor_square> eigensystem_sym(timing& times) const;
};

class tensor_3: public tensor{
  public:
  explicit tensor_3(const initializer_list<size_t> il): tensor(il) { assert(get_dims_size()==3); }
  tensor_3(initializer_list<size_t> il, double* ptr1, double* ptr2): tensor(il) { assert(get_size()==ptr2-ptr1); memcpy(begin(), ptr1, (ptr2-ptr1)*sizeof(double));}

  double& operator()(size_t i, size_t j, size_t k);
  double  operator()(size_t i, size_t j, size_t k) const;
};

class tensor_cube: public tensor_3{
  public:
  double& operator()(size_t i, size_t j, size_t k);
  double  operator()(size_t i, size_t j, size_t k) const;
  tensor_2 get_slice(const size_t c, const size_t p, dir_t direction) const;
};

class tensor_4: public tensor{
  public:
  explicit tensor_4(const initializer_list<size_t> il): tensor(il) { assert(get_dims_size()==4); }

  double& operator()(size_t i, size_t j, size_t k, size_t l);
  double  operator()(size_t i, size_t j, size_t k, size_t l) const;
  string to_ListPlot3D(const system3d& sys) const;
  string to_listplot(const system3d& sys) const;

};

class tensor_6: public tensor{
  public:
  explicit tensor_6(const initializer_list<size_t> il): tensor(il) { assert(get_dims_size()==6); }
  tensor_6(): tensor() { }

  double& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);
  double  operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;
  tensor_2 get_slice_flat(const size_t c, const size_t p, dir_t direction) const;

};

// tensor_square operator*(const tensor_square& T) const;
// tensor_square operator*(const tensor_square& T) const;
tensor_2 matmul(const tensor_2& T1, const tensor_2& T2, const char nt1, const char nt2);

tensor_2 matmul(const tensor_square& T1, const tensor_2& T2, const char nt1, const char nt2);
tensor_2 matmul(const tensor_2& T1, const tensor_square& T2, const char nt1, const char nt2);

tensor_square matmul(const tensor_square& T1, const tensor_square& T2, const char nt1, const char nt2);
tensor_square matmul(const tensor_square_sym& T1, const tensor_square& T2, const char nt1, const char nt2);
tensor_square matmul(const tensor_square& T1, const tensor_square_sym& T2, const char nt1, const char nt2);
tensor_square_sym matmul(const tensor_square_sym& T1, const tensor_square_sym& T2, const char nt1, const char nt2);

#endif


