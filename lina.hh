#ifndef LINA_HH
#define LINA_HH

extern "C" {
  // Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix.
  // DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
  extern void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);

  // Computes the eigenvalues and left and right eigenvectors of a general matrix.
  // DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
  extern void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);

  // DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
  // y := alpha*A*x + beta*y
  // A and X are unchanged on exit
  extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

  // DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  extern void dgemm_(char*, char*, int*, int*, int*,
                     double*, double*, int*,
                     double*, int*,
                     double*, double*, int*);

  // DGETRF( M, N, A, LDA, IPIV, INFO )
  // LU decomoposition of a general matrix
  extern void dgetrf_(int*, int*, double*, int*, int*, int*);

  // DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
  // generate inverse of a matrix given its LU decomposition
  extern void dgetri_(int*, double*, int*, int*, double*, int*, int*);
}

#endif
