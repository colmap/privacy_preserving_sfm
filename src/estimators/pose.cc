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

#include "estimators/pose.h"

#include "base/camera_models.h"
#include "base/cost_functions.h"
#include "base/pose.h"
#include "estimators/absolute_pose.h"
#include "optim/bundle_adjustment.h"
#include "util/matrix.h"
#include "util/misc.h"
#include "util/threading.h"

namespace colmap {
namespace {

typedef RANSAC<P6LEstimator> AbsolutePoseFromLinesRANSAC;

}  // namespace

bool EstimateAbsolutePoseFromLines(const RANSACOptions& options,
                                   const FeatureLines& lines2D,
                                   const std::vector<Eigen::Vector3d>& points3D,
                                   Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                                   size_t* num_inliers,
                                   std::vector<char>* inlier_mask) {
  options.Check();

  AbsolutePoseFromLinesRANSAC ransac(options);
  const auto report = ransac.Estimate(lines2D, points3D);
  *num_inliers = report.support.num_inliers;
  *inlier_mask = report.inlier_mask;

  if (*num_inliers == 0) {
    return false;
  }

  // If we used mostly aligned features for pose estimation it is likely that
  // the pose is off. Report no pose found in that case.
  int num_aligned_inliers = 0;
  const int num_lines = static_cast<int>(lines2D.size());
  for (int i = 0; i < num_lines; ++i) {
    if (inlier_mask->at(i) && lines2D.at(i).IsAligned()) {
      num_aligned_inliers += 1;
    }
  }

  // If we have more than 90% aligned lines as inliers something is likely
  // wrong.
  if (num_aligned_inliers > *num_inliers * 0.9) {
    return false;
  }

  // Extract pose parameters.
  *qvec = RotationMatrixToQuaternion(report.model.leftCols<3>());
  *tvec = report.model.rightCols<1>();

  if (IsNaN(*qvec) || IsNaN(*tvec)) {
    return false;
  }

  return true;
}

bool RefineAbsolutePoseFromLines(const AbsolutePoseRefinementOptions& options,
                                 const std::vector<char>& inlier_mask,
                                 const std::vector<Eigen::Vector3d>& lines2D,
                                 const std::vector<Eigen::Vector3d>& points3D,
                                 Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                                 Camera* camera) {
  CHECK_EQ(inlier_mask.size(), lines2D.size());
  CHECK_EQ(lines2D.size(), points3D.size());
  options.Check();

  ceres::LossFunction* loss_function =
          new ceres::CauchyLoss(options.loss_function_scale);

  double* camera_params_data = camera->ParamsData();
  double* qvec_data = qvec->data();
  double* tvec_data = tvec->data();

  std::vector<Eigen::Vector3d> points3D_copy = points3D;

  ceres::Problem problem;

  for (size_t i = 0; i < lines2D.size(); ++i) {
      // Skip outlier observations
      if (!inlier_mask[i]) {
          continue;
      }

      ceres::CostFunction* cost_function = nullptr;

      switch (camera->ModelId()) {
#define CAMERA_MODEL_CASE(CameraModel)                                             \
        case CameraModel::kModelId:                                                \
          cost_function =                                                          \
              BundleAdjustmentLineCostFunction<CameraModel>::Create(lines2D[i]);   \
          break;

          CAMERA_MODEL_SWITCH_CASES
#undef CAMERA_MODEL_CASE

      }

      problem.AddResidualBlock(cost_function, loss_function, qvec_data, tvec_data,
                               points3D_copy[i].data(), camera_params_data);
      problem.SetParameterBlockConstant(points3D_copy[i].data());
  }

  if (problem.NumResiduals() > 0) {
      // Quaternion parameterization.
      *qvec = NormalizeQuaternion(*qvec);
      ceres::LocalParameterization* quaternion_parameterization =
              new ceres::QuaternionParameterization;
      problem.SetParameterization(qvec_data, quaternion_parameterization);

      // Camera parameterization.
      if (!options.refine_focal_length && !options.refine_extra_params) {
          problem.SetParameterBlockConstant(camera->ParamsData());
      } else {
          // Always set the principal point as fixed.
          std::vector<int> camera_params_const;
          const std::vector<size_t>& principal_point_idxs =
                  camera->PrincipalPointIdxs();
          camera_params_const.insert(camera_params_const.end(),
                                     principal_point_idxs.begin(),
                                     principal_point_idxs.end());

          if (!options.refine_focal_length) {
              const std::vector<size_t>& focal_length_idxs =
                      camera->FocalLengthIdxs();
              camera_params_const.insert(camera_params_const.end(),
                                         focal_length_idxs.begin(),
                                         focal_length_idxs.end());
          }

          if (!options.refine_extra_params) {
              const std::vector<size_t>& extra_params_idxs =
                      camera->ExtraParamsIdxs();
              camera_params_const.insert(camera_params_const.end(),
                                         extra_params_idxs.begin(),
                                         extra_params_idxs.end());
          }

          if (camera_params_const.size() == camera->NumParams()) {
              problem.SetParameterBlockConstant(camera->ParamsData());
          } else {
              ceres::SubsetParameterization* camera_params_parameterization =
                      new ceres::SubsetParameterization(
                              static_cast<int>(camera->NumParams()), camera_params_const);
              problem.SetParameterization(camera->ParamsData(),
                                          camera_params_parameterization);
          }
      }
  }

  ceres::Solver::Options solver_options;
  solver_options.gradient_tolerance = options.gradient_tolerance;
  solver_options.max_num_iterations = options.max_num_iterations;
  solver_options.linear_solver_type = ceres::DENSE_QR;

  // The overhead of creating threads is too large.
  solver_options.num_threads = 1;
  #if CERES_VERSION_MAJOR < 2
  solver_options.num_linear_solver_threads = 1;
  #endif  // CERES_VERSION_MAJOR

  ceres::Solver::Summary summary;
  ceres::Solve(solver_options, &problem, &summary);

  if (solver_options.minimizer_progress_to_stdout) {
      std::cout << std::endl;
  }

  if (options.print_summary) {
      PrintHeading2("Pose refinement report");
      PrintSolverSummary(summary);
  }

  return summary.IsSolutionUsable();
}

}  // namespace colmap
