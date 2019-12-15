
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>  // for ofstream
#include <sstream>
#include <utility>  // for pair
#include <vector>

#include "lina.hh"
#include "system_3d.hh"
#include "tensor.hh"

using namespace std;


double& tensor_1::operator()(size_t i)       { return (*this)[i];}
double  tensor_1::operator()(size_t i) const { return (*this)[i];}
double& tensor_2::operator()(size_t i, size_t j)       { return (*this)[i*get_dims()[1]+j];}
double  tensor_2::operator()(size_t i, size_t j) const { return (*this)[i*get_dims()[1]+j];}
double& tensor_3::operator()(size_t i, size_t j, size_t k)       { return (*this)[(i*get_dims()[1]+j)*get_dims()[2]+k];}
double  tensor_3::operator()(size_t i, size_t j, size_t k) const { return (*this)[(i*get_dims()[1]+j)*get_dims()[2]+k];}
double& tensor_4::operator()(size_t i, size_t j, size_t k, size_t l)       { return (*this)[((i*get_dims()[1]+j)*get_dims()[2]+k)*get_dims()[3]+l];}
double  tensor_4::operator()(size_t i, size_t j, size_t k, size_t l) const { return (*this)[((i*get_dims()[1]+j)*get_dims()[2]+k)*get_dims()[3]+l];}
// double& tensor_5::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)       { return (*this)[(((i*dims[1]+j)*dims[2]+k)*dims[3]+l)*dims[4]+m];}
// double  tensor_5::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const { return (*this)[(((i*dims[1]+j)*dims[2]+k)*dims[3]+l)*dims[4]+m];}
double& tensor_6::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)       { return (*this)[((((i*get_dims()[1]+j)*get_dims()[2]+k)*get_dims()[3]+l)*get_dims()[4]+m)*get_dims()[5]+n];}
double  tensor_6::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const { return (*this)[((((i*get_dims()[1]+j)*get_dims()[2]+k)*get_dims()[3]+l)*get_dims()[4]+m)*get_dims()[5]+n];}


tensor tensor::operator=(const tensor& T){
  if(dims != T.dims){
    delete [] t;
    dims = T.dims;
    size = T.size;
    t = new double[size];
  }
  memcpy(t, T.t, size*sizeof(double));
  return *this;
}


tensor tensor::operator-() const {
  tensor t_new;
  t_new.dims = dims;
  t_new.size = size;
  delete [] t_new.t;
  t_new.t = new double[t_new.size];
  for(size_t i=0; i<size; i++) t_new[i] = -t[i];
  return t_new;
}

tensor& tensor::operator+=(const tensor& T){
  if(size==0){
    *this = T;
    return *this;
  }
  assert(dims==T.dims);
  for(int i=0; i<size; i++) t[i] += T[i];
  return *this;
}

tensor& tensor::operator-=(const tensor& T){
  assert(dims==T.dims);
  for(int i=0; i<size; i++) t[i] -= T[i];
  return *this;
}

tensor& tensor::operator*=(const double d){
  for(int i=0; i<size; i++) t[i] *= d;
  return *this;
}

tensor& tensor::operator/=(const double d){
  for(int i=0; i<size; i++) t[i] /= d;
  return *this;
}

bool tensor::operator==(const tensor& T) const {
  const vector<size_t> Tdims = T.get_dims();
  if(dims.size() != Tdims.size()) return false;
  for(int i=0; i<dims.size(); i++){
    if(dims[i] != Tdims[i]) return false;
  }
  for(int i=0; i<size; i++){
    if(t[i] != T[i]) return false;
  }
  return true;
}


double tensor_1::operator*(const tensor_1& t1) const { 
  double p = 0;
  const double dim = get_dims()[0];
  for(int i=0; i<dim; i++)
    p += (*this)[i] * t1[i];
  return p;
}

// template<int I>
// template<int J>
// tensor<J> tensor<I>::reshape(initializer_list<size_t> il) const {
//   assert(product(il) == size);
//   tensor<J> newT;
//   newT.set_dims(il);
//   newT.set_data(size, t); // size, source
//   return newT;
// }

// tensor<J> tensor<I>::get_slice(const system3d& sys, const size_t cc, const size_t pp, dir_t normal) const {
//   tensor<J> sl({sys.n_cells1d, sys.n_cells1d, sys.n_funcs_per_cell1d, sys.n_funcs_per_cell1d});
//   int index=0;
//   for(int c1=0; c1<sys.n_cells1d; c1++){
//   for(int c2=0; c2<sys.n_cells1d; c2++){
//     for(int p1=0; p1<sys.n_funcs_per_cell1d; p1++){
//     for(int p2=0; p2<sys.n_funcs_per_cell1d; p2++){
//       switch(normal){
//       case XX:
//         sl[index++] = (*this)(c1, c2, cc, p1, p2, pp);
//         break;
//       case YY:
//         sl[index++] = (*this)(c1, cc, c2, p1, pp, p2);
//         break;
//       case ZZ:
//         sl[index++] = (*this)(cc, c1, c2, pp, p1, p2);
//         break;
//       }
//     } }
//   } }
//   return sl;
// }
//
//
// tensor<J> tensor<I>::get_slice_flat(const system3d& sys, const size_t cc, const size_t pp, dir_t normal) const {
//   cerr << "slice at " << cc << ", " << pp << endl;
//   tensor<J> sl({sys.n_cells1d * sys.n_funcs_per_cell1d, sys.n_cells1d * sys.n_funcs_per_cell1d});
//   int index=0;
//   for(int c1=0; c1<sys.n_cells1d; c1++){
//   for(int p1=0; p1<sys.n_funcs_per_cell1d; p1++){
//     for(int c2=0; c2<sys.n_cells1d; c2++){
//     for(int p2=0; p2<sys.n_funcs_per_cell1d; p2++){
//       switch(normal){
//       case XX:
//         sl[index++] = (*this)(c1, c2, cc, p1, p2, pp);
//         break;
//       case YY:
//         sl[index++] = (*this)(c1, cc, c2, p1, pp, p2);
//         break;
//       case ZZ:
//         sl[index++] = (*this)(cc, c1, c2, pp, p1, p2);
//         break;
//       }
//     } }
//   } }
//   return sl;
// }


// for a tensor of rank I, pos is a list of I values where I-J of them define the slice and the other J are irrelevant
// variable is a list of I values where I-J are zero and J are one, specifying which dimensions are included in the slice
tensor_1 tensor_2::get_slice(const size_t pos, const dir_t normal) const {
  switch(normal){
  case XX:{
    tensor_1 slice({get_dims()[1]});
    for(size_t y=0; y<get_dims()[1]; y++){
      slice[y] = (*this)(pos, y);
    }
    return slice;
    }
  case YY:{
    tensor_1 slice({get_dims()[0]});
    for(size_t x=0; x<get_dims()[0]; x++){
      slice[x] = (*this)(x, pos);
    }
    return slice;
    }
  default:
    assert(false);
  }
}


//  tensor<1> operator[](size_t i) const {
//   size_t b = dim*i;
//   size_t e = dim*(i+1);
//   return tensor<1>({e-b}, m+b, m+e);
//  }
//  void set_row(size_t i, const tensor<1>& v){ // assign v to ith eigenvector
//    assert(v.get_size()==dim);
//    size_t b = dim*i;
//    memcpy(m+b, v.begin(), dim*sizeof(double)); // dest, source, n-bytes
//  }



// for a tensor of rank I, pos is a list of I values where I-J of them define the slice and the other J are irrelevant
// variable is a list of I values where I-J are zero and J are one, specifying which dimensions are included in the slice
void tensor_2::set_slice(const size_t pos, const dir_t normal, const tensor_1& slice) {
  switch(normal){
  case XX:{
    for(size_t y=0; y<get_dims()[1]; y++){
      (*this)(pos, y) = slice[y];
    }
    break;
    }
  case YY:{
    for(size_t x=0; x<get_dims()[0]; x++){
      (*this)(x, pos) = slice[x];
    }
    break;
    }
  default:
    assert(false);
  }
}


tensor_2 tensor_2::transpose() const {
  const size_t M = get_dims()[0];
  const size_t N = get_dims()[1];
  tensor_2 T({N,M});
  for(int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      T(j,i) = (*this)(i,j);
    }
  }
  return T;
}


string tensor_4::to_ListPlot3D(const system3d& sys) const {
  ostringstream s;
  s << "ListPlot3D[";
  s << to_listplot(sys);
  s << "];" << endl;
  return s.str();
}


// {{x1,y1,z1},{x2,y2,z2},{x3,y3,z3},...}
// assume the system was centred at 0,0
string tensor_4::to_listplot(const system3d& sys) const {
  ostringstream s;
  const double x_origin = -0.5 * sys.n_cells1d * sys.cell_width;
  const double y_origin = -0.5 * sys.n_cells1d * sys.cell_width;
//cout << x_origin << ", " << y_origin << endl;
      s << "{";
  for(int xc=0; xc<sys.n_cells1d; xc++){
  for(int yc=0; yc<sys.n_cells1d; yc++){
    for(int xp=0; xp<sys.n_funcs_per_cell1d; xp++){
    for(int yp=0; yp<sys.n_funcs_per_cell1d; yp++){
      // const size_t xi = xc*sys.n_funcs_per_cell1d + xp;
      // const size_t yi = yc*sys.n_funcs_per_cell1d + yp;
      const double x_coord = x_origin + (xc+0.5)*sys.cell_width + sys.grid.p[xp];
      const double y_coord = y_origin + (yc+0.5)*sys.cell_width + sys.grid.p[yp];
      s << "{" << x_coord << ", " << y_coord << "," << (*this)(xc, xp, yc, yp) << "}";
      if(yp<sys.n_funcs_per_cell1d-1) s << ",";
    }
    if(xp<sys.n_funcs_per_cell1d-1) s << ",";
    }
  if(yc<sys.n_cells1d-1) s << ",";
  }
  if(xc<sys.n_cells1d-1) s << ",";
  }
  s << "}";
  return s.str();
}


tensor_2 matmul(const tensor_2& T1, const tensor_2& T2, const char nt1, const char nt2){
  assert(nt1 == 'N' || nt1 == 'T');
  assert(nt2 == 'N' || nt2 == 'T');
  char NTchar1 = nt1, NTchar2 = nt2;
  int N1, N2, N3, N4;
  int N1_init = T1.get_dims()[1],
      N3_init = T2.get_dims()[1];
  if(nt1=='N'){
    N1 = T1.get_dims()[1]; N2 = T1.get_dims()[0];
  }
  else { // if(nt1=='T'){
    N1 = T1.get_dims()[0]; N2 = T1.get_dims()[1];
  }
  if(nt2=='N'){
    N3 = T2.get_dims()[1]; N4 = T2.get_dims()[0];
  }
  else { // if(nt2=='T'){
    N3 = T2.get_dims()[0]; N4 = T2.get_dims()[1];
  }

  assert( N2 == N3 );
  size_t N1_l = N1, N2_l = N4;
  double alpha = 1.0, beta = 0.0;
  tensor_2 new_t2({N2_l, N1_l});

   dgemm_(&NTchar1, &NTchar2, &N1, &N4, &N2,
          &alpha, T1.begin(), &N1_init,
          T2.begin(), &N3_init,
          &beta, new_t2.begin(), &N4);

  return new_t2;
}


