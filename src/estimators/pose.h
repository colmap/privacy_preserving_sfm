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

#ifndef COLMAP_SRC_ESTIMATORS_POSE_H_
#define COLMAP_SRC_ESTIMATORS_POSE_H_

#include <vector>

#include <Eigen/Core>

#include <ceres/ceres.h>

#include "base/camera.h"
#include "base/camera_models.h"
#include "feature/types.h"
#include "optim/loransac.h"
#include "util/alignment.h"
#include "util/logging.h"
#include "util/threading.h"
#include "util/types.h"

namespace colmap {

struct AbsolutePoseEstimationOptions {
  // Whether to estimate the focal length.
  bool estimate_focal_length = false;

  // Number of discrete samples for focal length estimation.
  size_t num_focal_length_samples = 30;

  // Minimum focal length ratio for discrete focal length sampling
  // around focal length of given camera.
  double min_focal_length_ratio = 0.2;

  // Maximum focal length ratio for discrete focal length sampling
  // around focal length of given camera.
  double max_focal_length_ratio = 5;

  // Number of threads for parallel estimation of focal length.
  int num_threads = ThreadPool::kMaxNumThreads;

  // Options used for P3P RANSAC.
  RANSACOptions ransac_options;

  void Check() const {
    CHECK_GT(num_focal_length_samples, 0);
    CHECK_GT(min_focal_length_ratio, 0);
    CHECK_GT(max_focal_length_ratio, 0);
    CHECK_LT(min_focal_length_ratio, max_focal_length_ratio);
    ransac_options.Check();
  }
};

struct AbsolutePoseRefinementOptions {
  // Convergence criterion.
  double gradient_tolerance = 1.0;

  // Maximum number of solver iterations.
  int max_num_iterations = 100;

  // Scaling factor determines at which residual robustification takes place.
  double loss_function_scale = 1.0;

  // Whether to refine the focal length parameter group.
  bool refine_focal_length = false;

  // Whether to refine the extra parameter group.
  bool refine_extra_params = false;

  // Whether to print final summary.
  bool print_summary = true;

  void Check() const {
    CHECK_GE(gradient_tolerance, 0.0);
    CHECK_GE(max_num_iterations, 0);
    CHECK_GE(loss_function_scale, 0.0);
  }
};

bool EstimateAbsolutePoseFromLines(const RANSACOptions& options,
                                   const FeatureLines& lines2D,
                                   const std::vector<Eigen::Vector3d>& points3D,
                                   Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                                   size_t* num_inliers,
                                   std::vector<char>* inlier_mask);

bool RefineAbsolutePoseFromLines(const AbsolutePoseRefinementOptions& options,
                                 const std::vector<char>& inlier_mask,
                                 const std::vector<Eigen::Vector3d>& lines2D,
                                 const std::vector<Eigen::Vector3d>& points3D,
                                 Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                                 Camera* camera);

}  // namespace colmap

#endif  // COLMAP_SRC_ESTIMATORS_POSE_H_
