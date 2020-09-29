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

#include "estimators/utils.h"

#include "util/logging.h"

namespace colmap {

void ComputeSquaredLineReprojectionError(
    const std::vector<Eigen::Vector3d>& lines2D,
    const std::vector<Eigen::Vector3d>& points3D,
    const Eigen::Matrix3x4d& proj_matrix, std::vector<double>* residuals) {
  CHECK_EQ(lines2D.size(), points3D.size());

  residuals->resize(lines2D.size());

  // Note that this code might not be as nice as Eigen expressions,
  // but it is significantly faster in various tests.

  const double P_00 = proj_matrix(0, 0);
  const double P_01 = proj_matrix(0, 1);
  const double P_02 = proj_matrix(0, 2);
  const double P_03 = proj_matrix(0, 3);
  const double P_10 = proj_matrix(1, 0);
  const double P_11 = proj_matrix(1, 1);
  const double P_12 = proj_matrix(1, 2);
  const double P_13 = proj_matrix(1, 3);
  const double P_20 = proj_matrix(2, 0);
  const double P_21 = proj_matrix(2, 1);
  const double P_22 = proj_matrix(2, 2);
  const double P_23 = proj_matrix(2, 3);

  for (size_t i = 0; i < lines2D.size(); ++i) {
    const double X_0 = points3D[i](0);
    const double X_1 = points3D[i](1);
    const double X_2 = points3D[i](2);

    // Project 3D point from world to camera.
    const double px_2 = P_20 * X_0 + P_21 * X_1 + P_22 * X_2 + P_23;

    // Check if 3D point is in front of camera.
    if (px_2 > std::numeric_limits<double>::epsilon()) {
      const double px_0 = P_00 * X_0 + P_01 * X_1 + P_02 * X_2 + P_03;
      const double px_1 = P_10 * X_0 + P_11 * X_1 + P_12 * X_2 + P_13;

      const double l_0 = lines2D[i](0);
      const double l_1 = lines2D[i](1);
      const double l_2 = lines2D[i](2);
      
      const double inv_px_2 = 1.0 / px_2;
      const double res = px_0 * l_0 * inv_px_2 + px_1 * l_1 * inv_px_2 + l_2;      

      (*residuals)[i] = res * res;
    } else {
      (*residuals)[i] = std::numeric_limits<double>::max();
    }
  }
}

}  // namespace colmap
