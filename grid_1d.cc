
#include <cassert>
#include <cmath>
#include <cfloat> // for DBL_EPSILON
#include <numeric>
#include <vector>

#include "auxiliary.hh"
#include "grid_1d.hh"
#include "polynomial.hh"

using namespace std;

// integrates exactly up to 2n-1
// borrowed from numerical recipes (and modified)
void gauss_legendre(const double x1, const double x2, vector<double>& points, vector<double>& weights, const int n) {
/*******************************************************************************
Given the lower and upper limits of integration x1 and x2, and given n, this
routine returns arrays x[1..n] and w[1..n] of length n, containing the abscissas
and weights of the Gauss-Legendre n-point quadrature formula.
*******************************************************************************/
  points.resize(n,0.0);
  weights.resize(n,0.0);
  const double eps = 1.0e-14; /* EPS is the relative precision. */
  double z1(0),z(0);
  const int m=(n+1)/2; /* The roots are symmetric, so we only find half of them. */
  const double xm=0.5*(x2+x1);
  const double xl=0.5*(x2-x1);
  for (int i=0;i<m;i++) { /* Loop over the desired roots. */
    z=cos(M_PI*(i+1-0.25)/(n+0.5));
    /* Starting with the above approximation to the ith root, we enter */
    /* the main loop of refinement by Newton's method.                 */
    double pp;
    do {
      double p1=1.0;
      double p2=0.0;
      for (int j=1;j<=n;j++) { /* Recurrence to get Legendre polynomial. */
        const double p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      /* p1 is now the desired Legendre polynomial. We next compute */
      /* pp, its derivative, by a standard relation involving also  */
      /* p2, the polynomial of one lower order.                     */
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp; /* Newton's method. */
    } while (fabs(z-z1) > eps);
    points[i]=xm-xl*z;      /* Scale the root to the desired interval, */
    points[n-i-1]=xm+xl*z;  /* and put in its symmetric counterpart.   */
    weights[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight             */
    weights[n-i-1]=weights[i];                 /* and its symmetric counterpart. */
  }
}

// and the resulting error is supposed to be h^(n+1) for even and h^(n+2) for odd numbers of points
void equidistant(const double x1, const double x2, vector<double>& points, vector<double>& weights, const int n) {
  points.resize(n,0.0);
  weights.resize(n,0.0);
  const double span = x2 - x1;
  const double h = span/(n-1.0);
  for (int i=0; i<n; i++){
    points[i] = x1 + i*h;
  }
  switch(n){
  case 2:
    weights = vector<double>({1,1});
    for(double &w: weights) w *= 1.0/2.0*h;
    break;
  case 3:
    weights = vector<double>({1,4,1});
    for(double &w: weights) w *= 2.0/6.0*h;
    break;
  case 4:
    weights = vector<double>({1,3,3,1});
    for(double &w: weights) w *= 3.0/8.0*h;
    break;
  case 5:
    weights = vector<double>({7,32,12,32,7});
    for(double &w: weights) w *= 4.0/90.0*h;
    break;
  case 6:
    weights = vector<double>({19,75,50,50,75,19});
    for(double &w: weights) w *= 5.0/288.0*h;
    break;
  case 7:
    weights = vector<double>({41,216,27,272,27,216,41});
    for(double &w: weights) w *= 6.0/840.0*h;
    break;
  case 8:
    weights = vector<double>({751,3577,1323,2989,2989,1323,3577,751});
    for(double &w: weights) w *= 7.0/17280.0*h;
    break;
  case 9:
    weights = vector<double>({989,5888,-928,10496,-4540,10496,-928,5888,989});
    for(double &w: weights) w *= 8.0/28350.0*h;
    break;
  case 10:
    weights = vector<double>({2857,15741,1080,19344,5778,5778,19344,1080,15741,2857});
    for(double &w: weights) w *= 9.0/89600.0*h;
    break;
  case 11:
    weights = vector<double>({16067,106300,-48525,272400,-260550,427368,-260550,272400,-48525,106300,16067});
    for(double &w: weights) w *= 10.0/598752.0*h;
    break;
  default:
    cerr << "choose a different number of grid points" << endl;
    assert(false);
  }
  
}

// exact for polynomials of degree 2n-3
// borrowed from erkale/Susi Lehtola (and modified)
void gauss_lobatto(const double x1, const double x2, vector<double>& points, vector<double>& weights, const int n) {
  assert(x2>x1);
  assert(n>=2);

  points.resize(n,0.0);
  weights.resize(n,0.0);
  const double span = x2 - x1;

  // cout << x1 << ", " << x2 << endl;

  // Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
  for ( int i = 0; i < n; i++ ) points[i] = cos( M_PI * ( double ) ( i ) / ( n - 1. ) );

  double points_old[n];
  double p[n*n];
  double error;
  const double tolerance = 100.0 * DBL_EPSILON;

  do {
    for ( int i = 0; i < n; i++ ) points_old[i] = points[i];
    for ( int i = 0; i < n; i++ ) p[i+0*n] = 1.0;
    for ( int i = 0; i < n; i++ ) p[i+1*n] = points[i];

    for ( int j = 2; j <= n-1; j++ ) {
      for ( int i = 0; i < n; i++) {
        // p[i+j*n] = ( ( double ) ( 2. * j - 1. ) * points[i] * p[i+(j-1)*n] + ( double ) ( - j + 1. ) * p[i+(j-2)*n] ) / ( double ) ( j );
        p[i+j*n] = ( ( 2.*j - 1. ) * points[i] * p[i+(j-1)*n] + ( -j + 1. ) * p[i+(j-2)*n] ) /  j ;
      }
    }

    for ( int i = 0; i < n; i++ ) {
      points[i] = points_old[i] - ( points[i] * p[i+(n-1)*n] - p[i+(n-2)*n] ) / ( ( double ) ( n ) * p[i+(n-1)*n] );
    }

    error = 0.0;
    for ( int i = 0; i < n; i++ ) {
      const double test = fabs(points[i] - points_old[i]);
      if(test>error) error=test;
    }
  } while ( error > tolerance );

  reverse(points.begin(), points.end());

  for ( int i = 0; i < n; i++ ) weights[i] = 2.0 / ( ( (n - 1.) * n ) * pow( p[i+(n-1)*n], 2 ) );

  // shift and scale to interval width
  for(double& pp: points) pp = (pp+1.0)*span/2.0 + x1;
  for(double& w: weights) w = w*span/2.0;

// cout << points << endl;
}



grid1d::grid1d(const vector<double> bounds, const vector<size_t> n, const grid_t g):
               p(vector<double>(accumulate(n.begin(), n.end(), 0))),
               w(vector<double>(accumulate(n.begin(), n.end(), 0))){
  switch(g){
  case EQUIDISTANT:
    assert(bounds.size()==2);
    assert(n.size()==1);
    equidistant(bounds[0], bounds[1], p, w, n[0]);
    break;
  case GAUSSIAN_LIN:
    assert(bounds.size()==2);
    assert(n.size()==1);
    gauss_legendre(bounds[0], bounds[1], p, w, n[0]);
    break;
  case GAUSSIAN_LOG:
    {
    assert(bounds.size()==2);
    assert(n.size()==1);
    gauss_legendre(-1.0, 1.0, p, w, n[0]);
    for(int i=0; i<p.size(); i++) p[i] = bounds[0] * pow(bounds[1]/bounds[0], 0.5*(p[i]+1.0));
    for(int i=0; i<w.size(); i++) w[i] = 0.5 * w[i] * p[i] * log(bounds[1]/bounds[0]);
    break;
    }
  case TGRID:
    {
    assert(bounds.size()==3);
    assert(n.size()==2);
    const double tstart(bounds[0]), tlin(bounds[1]), tlog(bounds[2]);
    grid1d glin(vector<double>({tstart, tlin}), vector<size_t>({n[0]}), GAUSSIAN_LIN);
    grid1d glog(vector<double>({tlin, tlog}), vector<size_t>({n[1]}), GAUSSIAN_LOG);
    *this = glin + glog;
    break;
    }
  case GAUSS_LOBATTO:
    assert(bounds.size()==2);
    assert(n.size()==1);
    gauss_lobatto(bounds[0], bounds[1], p, w, n[0]);
    break;
//  default:
//    cerr << "unknown grid type" << endl;
//    assert(false);
  }
}

