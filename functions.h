#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>

using namespace std;
using namespace Eigen;

// Checks whether the vector is normalized
bool isNormalized (const VectorXcd &psi);

// Commutator
MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B);

// Anticommutator
MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B);

// Projector |psi><psi|
MatrixXcd projector (const VectorXcd &psi);

// Density operator from its Bloch vector representation
Matrix2cd BlochToMatrix (double x, double y, double z);
Matrix2cd BlochToMatrix (const Vector3d &r);

// Bloch vector from state
Vector3d MatrixToBloch (const Matrix2cd &rho);

// Partial trace over the qubit 1 and 2
Matrix2cd tr_1 (const Matrix4cd &rho);
Matrix2cd tr_2 (const Matrix4cd &rho);

// Tensor product between two 2x2 matrices
Matrix4cd tens (const Matrix2cd &A, const Matrix2cd &B);

// Tensor producto between two 2-dim vectors
Vector4cd tens_state (const Vector2cd &psi1, const Vector2cd &psi2);

// Entropy for a qubit
double entropy (const Matrix2cd &rho);





// ------------------------- PAULI MATRICES -------------------------
static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::MatrixXcd sigma_x {{0,1},{1,0}};
static Eigen::MatrixXcd sigma_y {{0,-I},{I,0}};
static Eigen::MatrixXcd sigma_z {{1,0},{0,-1}};
static Eigen::MatrixXcd sigma_p {{0,1},{0,0}};
static Eigen::MatrixXcd sigma_m {{0,0},{1,0}};
static Eigen::MatrixXcd id {{1,0},{0,1}};





// ------------------------- SOME STATES -------------------------
static Eigen::VectorXcd ground_state {{0.,1.}};
static Eigen::VectorXcd excited_state {{1.,0.}};
static Eigen::VectorXcd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::VectorXcd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};
#endif 