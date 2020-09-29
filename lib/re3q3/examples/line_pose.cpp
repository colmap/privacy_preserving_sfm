/*
 * line_pose.cpp
 * Solves absolute pose problem from six 2D line to 3D point correspondences.
 * Author: Viktor Larsson
 */
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include "re3q3.h"
#include "line_pose.h"

/* Homogeneous linear constraints on rotation matrix
     Rcoeffs*R(:) = 0 
  converted into 3q3 problem. */
void rotation_to_e3q3(Eigen::Matrix<double, 3, 9> Rcoeffs, Eigen::Matrix<double, 3, 10> *coeffs) {
	for (int k = 0; k < 3; k++) {
		(*coeffs)(k, 0) = Rcoeffs(k, 0) - Rcoeffs(k, 4) - Rcoeffs(k, 8);
		(*coeffs)(k, 1) = 2 * Rcoeffs(k, 1) + 2 * Rcoeffs(k, 3);
		(*coeffs)(k, 2) = 2 * Rcoeffs(k, 2) + 2 * Rcoeffs(k, 6);
		(*coeffs)(k, 3) = Rcoeffs(k, 4) - Rcoeffs(k, 0) - Rcoeffs(k, 8);
		(*coeffs)(k, 4) = 2 * Rcoeffs(k, 5) + 2 * Rcoeffs(k, 7);
		(*coeffs)(k, 5) = Rcoeffs(k, 8) - Rcoeffs(k, 4) - Rcoeffs(k, 0);
		(*coeffs)(k, 6) = 2 * Rcoeffs(k, 5) - 2 * Rcoeffs(k, 7);
		(*coeffs)(k, 7) = 2 * Rcoeffs(k, 6) - 2 * Rcoeffs(k, 2);
		(*coeffs)(k, 8) = 2 * Rcoeffs(k, 1) - 2 * Rcoeffs(k, 3);
		(*coeffs)(k, 9) = Rcoeffs(k, 0) + Rcoeffs(k, 4) + Rcoeffs(k, 8);
	}
}


/* Inhomogeneous linear constraints on rotation matrix
	 Rcoeffs*[R(:);1] = 0
  converted into 3q3 problem. */
void rotation_to_e3q3(Eigen::Matrix<double, 3, 10> Rcoeffs, Eigen::Matrix<double, 3, 10> *coeffs) {
	for (int k = 0; k < 3; k++) {
		(*coeffs)(k, 0) = Rcoeffs(k, 0) - Rcoeffs(k, 4) - Rcoeffs(k, 8) + Rcoeffs(k, 9);
		(*coeffs)(k, 1) = 2 * Rcoeffs(k, 1) + 2 * Rcoeffs(k, 3);
		(*coeffs)(k, 2) = 2 * Rcoeffs(k, 2) + 2 * Rcoeffs(k, 6);
		(*coeffs)(k, 3) = Rcoeffs(k, 4) - Rcoeffs(k, 0) - Rcoeffs(k, 8) + Rcoeffs(k, 9);
		(*coeffs)(k, 4) = 2 * Rcoeffs(k, 5) + 2 * Rcoeffs(k, 7);
		(*coeffs)(k, 5) = Rcoeffs(k, 8) - Rcoeffs(k, 4) - Rcoeffs(k, 0) + Rcoeffs(k, 9);
		(*coeffs)(k, 6) = 2 * Rcoeffs(k, 5) - 2 * Rcoeffs(k, 7);
		(*coeffs)(k, 7) = 2 * Rcoeffs(k, 6) - 2 * Rcoeffs(k, 2);
		(*coeffs)(k, 8) = 2 * Rcoeffs(k, 1) - 2 * Rcoeffs(k, 3);
		(*coeffs)(k, 9) = Rcoeffs(k, 0) + Rcoeffs(k, 4) + Rcoeffs(k, 8) + Rcoeffs(k, 9);
	}
}

void cayley_param(Eigen::Matrix<double, 3, 1> c, Eigen::Matrix<double, 3, 3> *R) {
	*R << c(0)*c(0) - c(1)*c(1) - c(2)*c(2) + 1,
		2 * c(0)*c(1) - 2 * c(2),
		2 * c(1) + 2 * c(0)*c(2),
		2 * c(2) + 2 * c(0)*c(1),
		c(1)*c(1) - c(0)*c(0) - c(2)*c(2) + 1,
		2 * c(1)*c(2) - 2 * c(0),
		2 * c(0)*c(2) - 2 * c(1),
		2 * c(0) + 2 * c(1)*c(2),
		c(2)*c(2) - c(1)*c(1) - c(0)*c(0) + 1;
	*R /= 1 + c(0)*c(0) + c(1)*c(1) + c(2)*c(2);
}

int pose_from_lines(Eigen::Matrix<double, 3, 6> lines,
	Eigen::Matrix<double, 3, 6> points, std::vector<Eigen::Matrix<double, 3, 4>> *output,
	bool normalize_input)
{
	// If we should perform some normalization on the coordinate systems.
	// (Usually not required.)
	if (normalize_input) {
		lines.colwise().normalize();
		Eigen::Vector3d m = points.rowwise().mean();
		points.colwise() -= m;
		double s = points.rowwise().norm().mean();
		points /= s;
		int sols = pose_from_lines(lines, points, output, false);
		Eigen::Matrix4d H;
		H << 1, 0, 0, -m(0), 0, 1, 0, -m(1), 0, 0, 1, -m(2), 0, 0, 0, s;

		// revert change of variables
		for (int i = 0; i < sols; i++) {
			(*output)[i] *= H;
		}
		return sols;
	}

	// l'*(RX+t) = 0
	// l'*t + kron(X',l')*r = 0
	Eigen::Matrix<double, 3, 9> tt;
	Eigen::Matrix<double, 3, 9> Rcoeffs;
	for (int i = 0; i < 3; i++) {
		tt(i, 0) = points(0, i) * lines(0, i);
		tt(i, 1) = points(0, i) * lines(1, i);
		tt(i, 2) = points(0, i) * lines(2, i);
		tt(i, 3) = points(1, i) * lines(0, i);
		tt(i, 4) = points(1, i) * lines(1, i);
		tt(i, 5) = points(1, i) * lines(2, i);
		tt(i, 6) = points(2, i) * lines(0, i);
		tt(i, 7) = points(2, i) * lines(1, i);
		tt(i, 8) = points(2, i) * lines(2, i);

		Rcoeffs(i, 0) = points(0, i + 3) * lines(0, i + 3);
		Rcoeffs(i, 1) = points(0, i + 3) * lines(1, i + 3);
		Rcoeffs(i, 2) = points(0, i + 3) * lines(2, i + 3);
		Rcoeffs(i, 3) = points(1, i + 3) * lines(0, i + 3);
		Rcoeffs(i, 4) = points(1, i + 3) * lines(1, i + 3);
		Rcoeffs(i, 5) = points(1, i + 3) * lines(2, i + 3);
		Rcoeffs(i, 6) = points(2, i + 3) * lines(0, i + 3);
		Rcoeffs(i, 7) = points(2, i + 3) * lines(1, i + 3);
		Rcoeffs(i, 8) = points(2, i + 3) * lines(2, i + 3);
	}
	
	// Make sure we have non-singular block for solving for t!
	double det_tt = std::fabs(lines.block<3,3>(0,0).determinant());
	if (det_tt < 1e-10) {
		// if not we do a random combination which (hopefully) should solve the problem
		Eigen::Matrix<double, 3, 3> A;
		A.setRandom();
		tt += A * Rcoeffs.block<3, 9>(0, 0);
		lines.block<3, 3>(0, 0) += lines.block<3, 3>(0, 3) * A.transpose();
	}

	// Express t in R
	tt = lines.block<3, 3>(0, 0).transpose().partialPivLu().solve(tt);
	Rcoeffs -= (lines.block<3, 3>(0, 3)).transpose()*tt;

	// Convert to 3Q3 problem
	Eigen::Matrix<double, 3, 10> coeffs;
	rotation_to_e3q3(Rcoeffs, &coeffs);

	Eigen::Matrix<double, 3, 8> solutions;
	int n_sols = re3q3(coeffs, &solutions);

	Eigen::Matrix<double, 3, 3> R;
	Eigen::Matrix<double, 3, 1> t;
	Eigen::Matrix<double, 3, 4> pose;
	for (int i = 0; i < n_sols; i++) {
		cayley_param(solutions.block<3, 1>(0, i), &R);
		t = tt * Eigen::Map<Eigen::Matrix<double, 9, 1>>(R.data());
		pose << R, -t;
		output->push_back(pose);
	}
	return n_sols;
}