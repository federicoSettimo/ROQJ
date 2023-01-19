#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <iostream>
#include <complex>

using namespace std;

// A 2d vector for pure states
class state {
private:
  complex<double> elem[2];
public:
  complex<double> at (int i) const {return elem[i];}
  void set_elem (int i, complex<double> val) {elem[i] = val;}

  state(complex<double> v0 = 1, complex<double> v1 = 0, bool normalized = false) {elem[0] = v0; elem[1] = v1; if(normalized) this->normalize();}

  double stateNorm () const {return sqrt(norm(elem[0]) + norm(elem[1]));}

  void normalize () {state u = (*this)/stateNorm(); elem[0] = u.at(0); elem[1] = u.at(1);}

  bool isNormalized () {return this->stateNorm() == 1.;}

  void print () const {cout << "(" << real(elem[0]) << " + " << imag(elem[0]) << "i, " << real(elem[1]) << " + " << imag(elem[1]) << "i)\n";}

  // Overriding operators + and * (for scalars, matrices and states)
  state operator + (const state &u) const {return state(elem[0]+u.at(0), elem[1]+u.at(1));}
  void operator += (const state &u) {elem[0]+=u.at(0), elem[1]+=u.at(1);}

  state operator - (const state &u) const {return state(elem[0]-u.at(0), elem[1]-u.at(1));}
  void operator -= (const state &u) {elem[0]-=u.at(0), elem[1]-=u.at(1);}

  // Multiplication by a scalar: not wrong, just ugly. Used as v*alpha
  state operator * (complex<double> a) const {return state(elem[0] * a, elem[1] * a);}
  state operator * (double a) const {return state(elem[0] * a, elem[1] * a);}
  void operator *= (complex<double> a) {state u = (*this)*a; elem[0] = u.at(0); elem[1] = u.at(1);}
  void operator *= (double a) {state u = (*this)*a; elem[0] = u.at(0); elem[1] = u.at(1);}

  state operator / (complex<double> a) const {return (*this)*(1./a);}
  state operator / (double a) const {return (*this)*(1./a);}
  void operator /= (complex<double> a) {state u = (*this)/a; elem[0] = u.at(0); elem[1] = u.at(1);}
  void operator /= (double a) {state u = (*this)/a; elem[0] = u.at(0); elem[1] = u.at(1);}
};

void print (const state &psi) {psi.print();}

// Inner product <psi, phi>
complex<double> innerProduct (const state &psi, const state &phi) {return conj(psi.at(0))*phi.at(0) + conj(psi.at(1))*phi.at(1);}

// A Hermitian 2x2 matrix
class matrix {
private:
  // Matrix elements
  complex<double> elem[2][2];
public:
  // Matrix elements
  complex<double> at (int i, int j) const {return elem[i][j];}
  void set_elem (int i, int j, complex<double> val = 0) {elem[i][j] = val;}
  
  matrix(complex<double> a00 = 0, complex<double> a01 = 0, complex<double> a10 = 0, complex<double> a11 = 0) {elem[0][0] = a00; elem[0][1] = a01; elem[1][0] = a10; elem[1][1] = a11;}

  // Overriding operators + and * (for scalars, matrices and states)
  matrix operator + (matrix const &B) const {return matrix(elem[0][0]+B.at(0,0),elem[0][1]+B.at(0,1),elem[1][0]+B.at(1,0),elem[1][1]+B.at(1,1));}

  matrix operator - (matrix const &B) const {return matrix(elem[0][0]-B.at(0,0),elem[0][1]-B.at(0,1),elem[1][0]-B.at(1,0),elem[1][1]-B.at(1,1));}

  // Multiplication by a scalar: not wrong, just ugly. Used as A*alpha
  matrix operator * (complex<double> a) const {return matrix(elem[0][0]*a,elem[0][1]*a,elem[1][0]*a,elem[1][1]*a);}
  matrix operator * (double a) const {return matrix(elem[0][0]*a,elem[0][1]*a,elem[1][0]*a,elem[1][1]*a);}

  matrix operator * (matrix B) const;

  state operator * (state v) const {return state (elem[0][0] * v.at(0) + elem[0][1] * v.at(1), elem[1][0] * v.at(0) + elem[1][1] * v.at(1));}

  void print () const;

  // Hermitian conjugate matrix
  void dagger ();

  // Assuming that the matrix is Hermitian
  //double eigenvalue1 () const {return 0.5*( real(at(0,0) + at(1,1)) - sqrt( norm(at(0,0)-at(1,1)) + 4.*norm(at(0,1)) ) );}
  //double eigenvalue2 () const {return 0.5*( real(at(0,0) + at(1,1)) + sqrt( norm(at(0,0)-at(1,1)) + 4.*norm(at(0,1)) ) );}
  double eigenvalue1 () const {return 0.5*( real(at(0,0) + at(1,1)) - sqrt( norm(at(0,0)) + norm(at(1,1)) - 2.*real(at(0,0)*at(1,1)) + 4.*norm(at(0,1)) ) );}
  double eigenvalue2 () const {return 0.5*( real(at(0,0) + at(1,1)) + sqrt( norm(at(0,0)) + norm(at(1,1)) - 2.*real(at(0,0)*at(1,1)) + 4.*norm(at(0,1)) ) );}
  
  //state eigenvector1 () const {return state((at(0,0) - at(1,1) - sqrt( norm(at(0,0)-at(1,1)) + 4.*norm(at(0,1)) ))/(2.*at(1,0)), 1., true);}
  //state eigenvector2 () const {return state((at(0,0) - at(1,1) + sqrt( norm(at(0,0)-at(1,1)) + 4.*norm(at(0,1)) ))/(2.*at(1,0)), 1., true);}
  state eigenvector1 () const {return state((at(0,0) - at(1,1) - sqrt( norm(at(0,0)) + norm(at(1,1)) - 2.*real(at(0,0)*at(1,1)) + 4.*norm(at(0,1)) ))/(2.*at(1,0)), 1., true);}
  state eigenvector2 () const {return state((at(0,0) - at(1,1) + sqrt( norm(at(0,0)) + norm(at(1,1)) - 2.*real(at(0,0)*at(1,1)) + 4.*norm(at(0,1)) ))/(2.*at(1,0)), 1., true);}
};

// Definitions
complex<double> I(0,1); // Imag unity

// Useful operators
matrix sigma_x(0,1,1,0), sigma_y(0,-I,I,0), sigma_z(1,0,0,-1), sigma_p(0,1,0,0), sigma_m(0,0,1,0), Id(1,0,0,1);

// Ground, excited and |+>, |-> states 
state G(0,1), E(1,0), Plus(sqrt(.5), sqrt(.5)), Minus(sqrt(.5),-sqrt(.5));

matrix operator * (double alpha, const matrix &M) {return M*alpha;}
matrix operator * (complex<double> alpha, const matrix &M) {return M*alpha;}
state operator * (double alpha, const state &psi) {return psi*alpha;}
state operator * (complex<double> alpha, const state &psi) {return psi*alpha;}
state operator * (int alpha, const state &psi) {return ((double)alpha)*psi;}

// Commutator
matrix comm (const matrix &A, const matrix &B) {return A*B - B*A;}

// Anticommutator
matrix anticomm (const matrix &A, const matrix &B) {return A*B + B*A;}

void print (const matrix& M) {M.print();}

// Trace
complex<double> tr (const matrix &A) {return A.at(0,0) + A.at(1,1);}

// Given psi, returns |psi><psi|
matrix projector (const state &psi) {return matrix(norm(psi.at(0)), psi.at(0)*conj(psi.at(1)), conj(psi.at(1))*psi.at(0), norm(psi.at(1)));}

// Bloch vector to statistical operator
matrix BlochVector (double x, double y, double z) {return .5*(Id + x*sigma_x + y*sigma_y + z*sigma_z);}


// to move to matrix.cpp
void matrix::dagger () {
  elem[0][0] = conj(elem[0][0]);
  elem[1][1] = conj(elem[1][1]);
  complex<double> tmp = elem[0][1];
  elem[0][1] = conj(elem[1][0]);
  elem[1][0] = conj(tmp);
}

matrix matrix::operator* (matrix B) const {
  complex<double> c00, c01, c10, c11;
  c00 = at(0,0)*B.at(0,0) + at(0,1)*B.at(1,0);
  c01 = at(0,0)*B.at(0,1) + at(0,1)*B.at(1,1);
  c10 = at(1,0)*B.at(0,0) + at(1,1)*B.at(1,0);
  c11 = at(1,0)*B.at(0,1) + at(1,1)*B.at(1,1);
  return matrix(c00,c01,c10,c11);
  /*matrix C;
  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++)
      C.set_elem(i,j, elem[i][0]*B.at(0,j) + elem[i][1]*B.at(1,j));
  }
  return C;*/
}

void matrix::print () const {
  cout << real(at(0,0)) << " + " << imag(at(0,0)) << "i\t" << real(at(0,1)) << " + " << imag(at(0,1)) << "i\n" << real(at(1,0)) << " + " << imag(at(1,0)) << "i\t" << real(at(1,1)) << " + " << imag(at(1,1)) << "i\n";
}

#endif