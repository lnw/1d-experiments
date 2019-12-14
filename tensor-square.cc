
#include <tuple>
#include <cassert>
#include <fstream>  // for ofstream
#include <utility>  // for pair
#include <vector>

#include "timing.hh"
#include "auxiliary.hh"
#include "tensor.hh"
#include "lina.hh"

using namespace std;

tensor_square_sym tensor_square_sym::operator=(const tensor_square_sym& T){
  if(get_dims() != T.get_dims()){
    set_dims(T.get_dims());
  }
  set_data(T.get_size(), T.begin());
  return *this;
}

tensor_square_sym tensor_square_sym::operator=(const tensor& T){
  if(get_dims() != T.get_dims()){
    set_dims(T.get_dims());
  }
  set_data(T.get_size(), T.begin());
  return *this;
}

// matrix matrix::operator=(const matrix& M){
//   if(dim != M.size()){
//     delete [] m;
//     dim = M.size();
//     m = new double[dim*dim];
//   }
//   memcpy(m, M.m, dim*dim*sizeof(double));
//   return *this;
// }
//
// matrix& matrix::operator+=(const matrix& M){
//   assert(dim == M.size());
//   for(int i=0; i<dim; i++){
//     for(int j=0; j<dim; j++){
//       (*this)(i,j) += M(i,j);
//     }
//   }
//   return *this;
// }

tensor_square_sym tensor_square_sym::operator*(const double d) const {
  tensor_square_sym t_new(*this);
  const size_t size_l = get_size();
  for(int i=0; i<size_l; i++) t_new[i] = (*this)[i] * d;
  return t_new;
}

// matrix matrix::operator-() const {
//   matrix mat(dim);
//   for(int i=0; i<dim; i++){
//     for(int j=0; j<dim; j++){
//       mat(i,j) = -(*this)(i,j);
//     }
//   }
//   return mat;
// }



// two permutations of tensor_2/tensor_square
tensor_2 matmul(const tensor_2& T1, const tensor_square& T2, const char nt1, const char nt2){
  char NTchar1 = nt1, NTchar2 = nt2;
  int N1 = T1.get_dims()[0],
      N2 = T1.get_dims()[1],
      N3 = T2.get_dims()[0];
  double alpha = 1.0, beta = 0.0;

  size_t N1_l, N2_l = T2.get_dims()[0]; // square anyway
  if(nt1=='N'){
    assert( N1==N3 );
    N1_l = T1.get_dims()[1];
  }
  else if(nt1=='T'){
    assert( N2==N3 );
    N1_l = T1.get_dims()[0];
  }
  else assert(false);
  tensor_2 new_t2({N2_l, N1_l});

  if(nt1=='N'){
    dgemm_(&NTchar1, &NTchar2, &N2, &N3, &N3,
           &alpha, T1.begin(), &N2,
           T2.begin(), &N3,
           &beta, new_t2.begin(), &N2);
  }
  else if(nt1=='T'){
    dgemm_(&NTchar1, &NTchar2, &N1, &N3, &N3,
           &alpha, T1.begin(), &N2,
           T2.begin(), &N3,
           &beta, new_t2.begin(), &N1);
  }
  return new_t2;
}
tensor_2 matmul(const tensor_square& T1, const tensor_2& T2, const char nt1, const char nt2){
  char NTchar1 = nt1, NTchar2 = nt2;
  int N1 = T1.get_dims()[0],
      N3 = T2.get_dims()[0],
      N4 = T2.get_dims()[1];
  double alpha = 1.0, beta = 0.0;

  size_t N1_l = T1.get_dims()[0], N2_l;
  if(nt2=='N'){
    assert( N1==N4 );
    N2_l = T2.get_dims()[0];
  }
  else if(nt2=='T'){
    assert( N1==N3 );
    N2_l = T2.get_dims()[1];
  }
  else assert(false);
  tensor_2 new_t2({N2_l, N1_l});

  if(nt2=='N'){
    dgemm_(&NTchar1, &NTchar2, &N1, &N3, &N4,
           &alpha, T1.begin(), &N1,
           T2.begin(), &N4,
           &beta, new_t2.begin(), &N4);
  }
  else if(nt2=='T'){
    dgemm_(&NTchar1, &NTchar2, &N1, &N4, &N3,
           &alpha, T1.begin(), &N1,
           T2.begin(), &N4,
           &beta, new_t2.begin(), &N3);
  }
  return new_t2;
}



// four permutations of square_sym/square

tensor_square matmul(const tensor_square& T1, const tensor_square& T2, const char nt1, const char nt2){
  assert(T1.get_dims()==T2.get_dims());
  size_t N_s = T1.get_dims()[0];
  tensor_square new_tsq({N_s,N_s});

  char NTchar1 = nt1, NTchar2 = nt2;
  int N = T1.get_dims()[0];
  double alpha = 1.0, beta = 0.0;

  dgemm_(&NTchar1, &NTchar2, &N, &N, &N,
         &alpha, T1.begin(), &N,
         T2.begin(), &N,
         &beta, new_tsq.begin(), &N);

  return new_tsq;
}
tensor_square matmul(const tensor_square_sym& T1, const tensor_square& T2, const char nt1, const char nt2){
  assert(T1.get_dims()==T2.get_dims());
  size_t N_s = T1.get_dims()[0];
  tensor_square new_tsq({N_s,N_s});

  char NTchar1 = nt1, NTchar2 = nt2;
  int N = T1.get_dims()[0];
  double alpha = 1.0, beta = 0.0;

  dgemm_(&NTchar1, &NTchar2, &N, &N, &N,
         &alpha, T1.begin(), &N,
         T2.begin(), &N,
         &beta, new_tsq.begin(), &N);

  return new_tsq;
}
tensor_square matmul(const tensor_square& T1, const tensor_square_sym& T2, const char nt1, const char nt2){
  assert(T1.get_dims()==T2.get_dims());
  size_t N_s = T1.get_dims()[0];
  tensor_square new_tsq({N_s,N_s});

  char NTchar1 = nt1, NTchar2 = nt2;
  int N = T1.get_dims()[0];
  double alpha = 1.0, beta = 0.0;

  dgemm_(&NTchar1, &NTchar2, &N, &N, &N,
         &alpha, T1.begin(), &N,
         T2.begin(), &N,
         &beta, new_tsq.begin(), &N);

  return new_tsq;
}
tensor_square_sym matmul(const tensor_square_sym& T1, const tensor_square_sym& T2, const char nt1, const char nt2){
  assert(T1.get_dims()==T2.get_dims());
  size_t N_s = T1.get_dims()[0];
  tensor_square_sym new_tsq({N_s,N_s});

  char NTchar1 = nt1, NTchar2 = nt2;
  int N = T1.get_dims()[0];
  double alpha = 1.0, beta = 0.0;

  dgemm_(&NTchar1, &NTchar2, &N, &N, &N,
         &alpha, T1.begin(), &N,
         T2.begin(), &N,
         &beta, new_tsq.begin(), &N);

  return new_tsq;
}

// tensor_square tensor_square_sym::operator*(const tensor_square& T) const {
//   assert(get_dims()==T.get_dims());
//   size_t N_s = get_dims()[0];
//   tensor_square new_tsq({N_s,N_s});
// 
//   char Nchar = 'N';
//   int N = get_dims()[0];
//   double alpha = 1.0, beta = 0.0;
// 
//   dgemm_(&Nchar, &Nchar, &N, &N, &N,
//          &alpha, begin(), &N,
//          T.begin(), &N,
//          &beta, new_tsq.begin(), &N);
// 
//   return new_tsq;
// }


tensor_1 tensor_square::operator*(const tensor_1& t1) const {
  assert(get_dims()[0]==t1.get_dims()[0]);
  size_t N_s = get_dims()[0];
  tensor_1 tp({N_s});

  int N = get_dims()[0];
  char Nchar = 'N';
  double zero_d = 0.0;
  int one_i = 1.0;
  double one_d = 1.0;
  dgemv_(&Nchar, &N, &N, &one_d, begin(), &N, t1.begin(), &one_i, &zero_d, tp.begin(), &one_i);
  return tp;
}


tensor_1 tensor_square_sym::operator*(const tensor_1& t1) const {
  assert(get_dims()[0]==t1.get_dims()[0]);
  size_t N_s = get_dims()[0];
  tensor_1 tp({N_s});

  int N = get_dims()[0];
  char Nchar = 'N';
  double zero_d = 0.0;
  int one_i = 1.0;
  double one_d = 1.0;
  dgemv_(&Nchar, &N, &N, &one_d, begin(), &N, t1.begin(), &one_i, &zero_d, tp.begin(), &one_i);
  return tp;
}


tuple<tensor_1, tensor_square, tensor_square> tensor_square::eigensystem_gen(timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  char Vchar = 'V';
  int N = get_dims()[0];
  const size_t dim = get_dims()[0];
  tensor_square input(*this);
  tensor_square evecs_l({dim,dim}); // left
  tensor_square evecs_r({dim,dim}); // right

  tensor_1 evals_r({dim}); // real
  tensor_1 evals_i({dim}); // imaginary
  int info;

  int lwork = -1;
  double *work = new double[1];
  dgeev_(&Vchar, &Vchar, &N, input.begin(), &N, evals_r.begin(), evals_i.begin(), evecs_l.begin(), &N, evecs_r.begin(), &N, work, &lwork, &info);
  lwork = work[0];
  delete [] work;
  work = new double[lwork];
  dgeev_(&Vchar, &Vchar, &N, input.begin(), &N, evals_r.begin(), evals_i.begin(), evecs_l.begin(), &N, evecs_r.begin(), &N, work, &lwork, &info);
  delete [] work;

  // cout << "r: " << evals_r << endl;
  // cout << "i: " << evals_i << endl;

  auto t1 = Clock::now();
  times.append_item("solve eigensystem (gen)", t1 - t0);
  times.ind_down();

  return make_tuple(evals_r, evecs_l, evecs_r);
}


tuple<tensor_1, tensor_square> tensor_square_sym::eigensystem_sym(timing& times) const {
  times.ind_up();
  auto t0 = Clock::now();

  // assert(is_symmetric());

  char Vchar = 'V';
  char Uchar = 'L';
  int N = get_dims()[0];
  const size_t dim = get_dims()[0];
  tensor_square evecs(*this);
  tensor_1 evals({dim});
  int info;

  int lwork = -1;
  double *work = new double[1];
  dsyev_(&Vchar, &Uchar, &N, evecs.begin(), &N, evals.begin(), work, &lwork, &info);
  lwork = work[0];
  delete [] work;
  work = new double[lwork];
  dsyev_(&Vchar, &Uchar, &N, evecs.begin(), &N, evals.begin(), work, &lwork, &info);
  delete [] work;

  auto t1 = Clock::now();
  times.append_item("solve eigensystem (sym)", t1 - t0);
  times.ind_down();

  return make_tuple(evals, evecs);
}


bool tensor_square::is_symmetric() const{
  const double eps = 1.e-10;
  for(int i=0; i<get_dims()[0]; i++){
    for(int j=i; j<get_dims()[0]; j++){
      if( abs( (*this)(i,j) - (*this)(j,i) ) > eps ){
        if( abs( (*this)(i,j) / (*this)(j,i) - 1.) > eps ){
          cout << i << "," << j << ": " << setprecision(15) << (*this)(i,j) << " != " <<  (*this)(j,i) << endl;
          return false;
        }
      }
    }
  }
  return true;
}

bool tensor_square::is_point_symmetric() const{
  const double eps = 1.e-10;
  const size_t s = get_dims()[0];
  for(int i=0; i<s/2; i++){
    for(int j=0; j<s; j++){
      if( abs( (*this)(i,j) - (*this)(s-i-1,s-j-1) ) > eps ){
        if( abs( (*this)(i,j) / (*this)(s-i-1,s-j-1) - 1.) > eps ){
          cout << i << "," << j << ": " << setprecision(15) << (*this)(i,j) << " != " <<  (*this)(s-i-1,s-j-1) << endl;
          return false;
        }
      }
    }
  }
  return true;
}


tensor_square tensor_square::inverse() const {
  int N=get_dims()[0];
  tensor_square input(*this);
  int *ipiv = new int[N]; // pivots 
  int info;
  dgetrf_(&N, &N, input.begin(), &N, ipiv, &info);
  assert(info==0);

  int lwork = N*N;
  double *work = new double[lwork];
  dgetri_(&N, input.begin(), &N, ipiv, work, &lwork, &info);
  assert(info==0);
  delete [] ipiv;
  delete [] work;

  return input;
}

