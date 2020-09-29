// Copyright (c) 2020, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Authors: Johannes L. Schoenberger (jsch-at-demuc-dot-de)
//          Viktor Larsson (viktor.larsson@inf.ethz.ch)
//          Marcel Geppert (marcel.geppert@inf.ethz.ch)

#include "estimators/absolute_pose.h"
#include <re3q3/re3q3.h>
#include "base/polynomial.h"
#include "estimators/utils.h"
#include "util/logging.h"

namespace colmap {
namespace {

/* Homogeneous linear constraints on rotation matrix
	 Rcoeffs*R(:) = 0
  converted into 3q3 problem. */
void rotation_to_e3q3(
    const Eigen::Matrix<double, 3, 9> &Rcoeffs,
    Eigen::Matrix<double, 3, 10> *coeffs)
{
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

void cayley_param(Eigen::Matrix<double, 3, 1> c, Eigen::Matrix<double, 3, 3> *R) {
  *R << c(0) * c(0) - c(1) * c(1) - c(2) * c(2) + 1,
      2 * c(0) * c(1) - 2 * c(2),
      2 * c(1) + 2 * c(0) * c(2),
      2 * c(2) + 2 * c(0) * c(1),
      c(1) * c(1) - c(0) * c(0) - c(2) * c(2) + 1,
      2 * c(1) * c(2) - 2 * c(0),
      2 * c(0) * c(2) - 2 * c(1),
      2 * c(0) + 2 * c(1) * c(2),
      c(2) * c(2) - c(1) * c(1) - c(0) * c(0) + 1;
  *R /= 1 + c(0) * c(0) + c(1) * c(1) + c(2) * c(2);
}

}  // namespace

std::vector<P6LEstimator::M_t> P6LEstimator::Estimate(const std::vector<P6LEstimator::X_t>& lines2D,
                                                      const std::vector<P6LEstimator::Y_t>& points3D)
{


  Eigen::Matrix<double, 3, 6> lines;
  Eigen::Matrix<double, 3, 6> points;

  bool all_lines_aligned = true;
  for(int i = 0; i < 6; ++i){
    lines.col(i) = lines2D[i].Line();
    points.col(i) = points3D[i];

    all_lines_aligned &= lines2D[i].IsAligned();
  }

  if (all_lines_aligned) {
    return std::vector<P6LEstimator::M_t>();
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
	double det_tt = std::abs(lines.block<3, 3>(0, 0).determinant());
	Eigen::Matrix3d B = lines.block<3, 3>(0, 0);
	if (det_tt < 1e-10) {
		// if not we do a random combination which (hopefully) should solve the problem
		Eigen::Matrix<double, 3, 3> A;
		A.setRandom();
		tt += A * Rcoeffs.block<3, 9>(0, 0);
		B += lines.block<3, 3>(0, 3) * A.transpose();
	}

	// Express t in R
	tt = B.transpose().partialPivLu().solve(tt);
	Rcoeffs -= lines.block<3, 3>(0, 3).transpose() * tt;

	// Convert to 3Q3 problem
	Eigen::Matrix<double, 3, 10> coeffs;
	rotation_to_e3q3(Rcoeffs, &coeffs);

	Eigen::Matrix<double, 3, 8> solutions;
	int n_sols = re3q3(coeffs, &solutions);


  std::vector<P6LEstimator::M_t> output;

	Eigen::Matrix3d R;
	Eigen::Vector3d t;	
	for (int i = 0; i < n_sols; i++) {
		cayley_param(solutions.block<3, 1>(0, i), &R);
		t = -tt * Eigen::Map<Eigen::Matrix<double, 9, 1>>(R.data());		

    P6LEstimator::M_t pose;
    pose.block<3,3>(0,0) = R;
    pose.block<3,1>(0,3) = t;
		output.push_back(pose);
	}
	return output;
}


void P6LEstimator::Residuals(const std::vector<P6LEstimator::X_t>& lines2D,
                             const std::vector<P6LEstimator::Y_t>& points3D,
                             const P6LEstimator::M_t& proj_matrix, std::vector<double>* residuals)
{
  std::vector<Eigen::Vector3d> line_params;
  for (const auto& line : lines2D) {
    line_params.emplace_back(line.Line());
  }
  ComputeSquaredLineReprojectionError(line_params, points3D, proj_matrix, residuals);
}


}  // namespace colmap
