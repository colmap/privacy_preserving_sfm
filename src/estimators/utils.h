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

#ifndef COLMAP_SRC_ESTIMATORS_UTILS_H_
#define COLMAP_SRC_ESTIMATORS_UTILS_H_

#include <vector>

#include <Eigen/Core>

#include "util/alignment.h"
#include "util/types.h"

namespace colmap {

// Calculate the squared reprojection error in normalized coordinates given
// a set of 2D-3D line-point correspondences and a projection matrix.
// Returns DBL_MAX if a 3D point is behind the given camera.
//
// @param lines2D       Normalized 2D image lines. Assumes first 2x1 vector is unit
// @param points3D      3D world points.
// @param proj_matrix   3x4 projection matrix.
// @param residuals     Output vector of residuals.
void ComputeSquaredLineReprojectionError(
    const std::vector<Eigen::Vector3d>& lines2D,
    const std::vector<Eigen::Vector3d>& points3D,
    const Eigen::Matrix3x4d& proj_matrix, std::vector<double>* residuals);

}  // namespace colmap

#endif  // COLMAP_SRC_ESTIMATORS_UTILS_H_
