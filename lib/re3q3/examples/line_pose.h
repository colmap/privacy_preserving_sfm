/*
 * line_pose.h
 * Solves absolute pose problem from six 2D line to 3D point correspondences.
 * Author: Viktor Larsson
 */
#pragma once
#include <Eigen/Dense>
#include <vector>

/* Homogeneous linear constraints on rotation matrix
     Rcoeffs*R(:) = 0 
  converted into 3q3 problem. */
void rotation_to_e3q3(Eigen::Matrix<double, 3, 9> Rcoeffs,
	Eigen::Matrix<double, 3, 10> *coeffs);

/* Inhomogeneous linear constraints on rotation matrix
	 Rcoeffs*[R(:);1] = 0
  converted into 3q3 problem. */
void rotation_to_e3q3(Eigen::Matrix<double, 3, 10> Rcoeffs,
	Eigen::Matrix<double, 3, 10> *coeffs);

/* Cayley parameterization of rotation */
void cayley_param(Eigen::Matrix<double, 3, 1> c,
	Eigen::Matrix<double, 3, 3> *R);

/* Absolute pose from six 2D line to 3D point correspondences */
int pose_from_lines(Eigen::Matrix<double, 3, 6> lines,
	Eigen::Matrix<double, 3, 6> points, 
	std::vector<Eigen::Matrix<double, 3, 4>> *output,
	bool normalize_input = false);