#include "functions.h"

bool isNormalized (const VectorXcd &psi) {return psi.norm() == 1;}

MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B) {return A*B-B*A;}

MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B) {return A*B+B*A;}

MatrixXcd projector (const VectorXcd &psi) {return psi*psi.adjoint();}

Matrix2cd BlochToMatrix (double x, double y, double z) {
  double r = sqrt(x*x + y*y + z*z);
  if (r > 1.) {x /= r; y /= r; z /= r;}
  return .5*(id + x*sigma_x + y*sigma_y + z*sigma_z);
}

Matrix2cd BlochToMatrix (const Vector3d &r) {return BlochToMatrix(r(0),r(1),r(2));}

Vector3d MatrixToBloch (const Matrix2cd &rho) {
  Vector3d r;
  r << real((rho*sigma_x).trace()), real((rho*sigma_y).trace()), real((rho*sigma_z).trace());
  return r;
}

Matrix2cd tr_1(const Matrix4cd &rho) {
	Matrix2cd A(2,2);
  A << rho(0,0) + rho(2,2), rho(0,1) + rho(2,3), rho(1,0) + rho(3,2), rho(1,1) + rho(3,3);
	return A;
}

Matrix2cd tr_2(const Matrix4cd &rho) {
	Matrix2cd A(2,2);
  A << rho(0,0) + rho(1,1), rho(0,2) + rho(1,3), rho(2,0) + rho(3,1), rho(2,2) + rho(3,3);
  return A;
}

Matrix4cd tens (const Matrix2cd &A, const Matrix2cd &B) {
  Matrix4cd C = MatrixXcd::Zero(4,4);
  C.topLeftCorner(2,2) = A(0,0)*B;
  C.topRightCorner(2,2) = A(0,1)*B;
  C.bottomLeftCorner(2,2) = A(1,0)*B;
  C.bottomRightCorner(2,2) = A(1,1)*B;
  return C;
}


Vector4cd tens_state (const Vector2cd &psi1, const Vector2cd &psi2) {
  Vector4cd psi = VectorXcd::Zero(4);
  psi(0) = psi1(0)*psi2(0);
  psi(1) = psi1(0)*psi2(1);
  psi(2) = psi1(1)*psi2(0);
  psi(3) = psi1(1)*psi2(1);
  return psi;
}

double entropy (const Matrix2cd &rho) {
  ComplexEigenSolver<Matrix2cd> eigs;
  eigs.compute(rho);
  double p = real(eigs.eigenvalues()(0));
  return -p*log2(p) - (1.-p)*log2(1.-p);
}