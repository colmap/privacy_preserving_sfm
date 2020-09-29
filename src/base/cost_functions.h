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

#ifndef COLMAP_SRC_BASE_COST_FUNCTIONS_H_
#define COLMAP_SRC_BASE_COST_FUNCTIONS_H_

#include <Eigen/Core>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

namespace colmap {

// Bundle adjustment cost function with feature lines for variable
// camera pose and calibration and point parameters.
template <typename CameraModel>
class BundleAdjustmentLineCostFunction {
public:
    explicit BundleAdjustmentLineCostFunction(const Eigen::Vector3d& line_2D)
            : a_(line_2D(0)), b_(line_2D(1)), c_(line_2D(2)) {
      const double norm = sqrt(a_ * a_ + b_ * b_);
      CHECK_NEAR(norm, 1.0, 1e-6);
    }

    static ceres::CostFunction* Create(const Eigen::Vector3d& line_2D) {
        return (new ceres::AutoDiffCostFunction<
                BundleAdjustmentLineCostFunction<CameraModel>, 2, 4, 3, 3,
                    CameraModel::kNumParams>(
                new BundleAdjustmentLineCostFunction(line_2D)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
                    const T* const point3D, const T* const camera_params,
                    T* residuals) const {
      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += tvec[0];
      projection[1] += tvec[1];
      projection[2] += tvec[2];

      // Project to image plane.
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      // Compute the closest point on the line
      T alpha = T(a_) * projection[0] + T(b_) * projection[1] + T(c_);

      T line_point[2];
      line_point[0] = projection[0] - alpha * T(a_);
      line_point[1] = projection[1] - alpha * T(b_);

      // Transform both points to pixel coordinates
      T im_projection[2];
      T im_line_point[2];

      // Distort and transform to pixel space.
      CameraModel::WorldToImage(camera_params, projection[0], projection[1],
                                &im_projection[0], &im_projection[1]);

      CameraModel::WorldToImage(camera_params, line_point[0], line_point[1],
                                &im_line_point[0], &im_line_point[1]);

      // Compute the error
      residuals[0] = im_projection[0] - im_line_point[0];
      residuals[1] = im_projection[1] - im_line_point[1];

      return true;
    }

private:
  const double a_;
  const double b_;
  const double c_;
};

// Bundle adjustment cost function for variable
// camera calibration and point parameters, and fixed camera pose with line constraints.
template <typename CameraModel>
class BundleAdjustmentConstantPoseLineCostFunction {
public:
    BundleAdjustmentConstantPoseLineCostFunction(const Eigen::Vector4d& qvec,
                                             const Eigen::Vector3d& tvec,
                                             const Eigen::Vector3d& line2D)
            : qw_(qvec(0)),
              qx_(qvec(1)),
              qy_(qvec(2)),
              qz_(qvec(3)),
              tx_(tvec(0)),
              ty_(tvec(1)),
              tz_(tvec(2)),
              a_(line2D(0)),
              b_(line2D(1)),
              c_(line2D(2)) {
      const double norm = sqrt(a_ * a_ + b_ * b_);
      CHECK_NEAR(norm, 1.0, 1e-6);
    }

    static ceres::CostFunction* Create(const Eigen::Vector4d& qvec,
                                       const Eigen::Vector3d& tvec,
                                       const Eigen::Vector3d& line2D) {
        return (new ceres::AutoDiffCostFunction<
                BundleAdjustmentConstantPoseLineCostFunction<CameraModel>, 2,
                    3, CameraModel::kNumParams>(
                new BundleAdjustmentConstantPoseLineCostFunction(qvec, tvec, line2D)));
    }

    template <typename T>
    bool operator()(const T* const point3D, const T* const camera_params,
                    T* residuals) const {
      const T qvec[4] = {T(qw_), T(qx_), T(qy_), T(qz_)};

      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += T(tx_);
      projection[1] += T(ty_);
      projection[2] += T(tz_);

      // Project to image plane.
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      // Compute the closest point on the line
      T alpha = T(a_) * projection[0] + T(b_) * projection[1] + T(c_);

      T line_point[2];
      line_point[0] = projection[0] - alpha * T(a_);
      line_point[1] = projection[1] - alpha * T(b_);

      // Transform both points to pixel coordinates
      T im_projection[2];
      T im_line_point[2];

      // Distort and transform to pixel space.
      CameraModel::WorldToImage(camera_params, projection[0], projection[1],
                                &im_projection[0], &im_projection[1]);

      CameraModel::WorldToImage(camera_params, line_point[0], line_point[1],
                                &im_line_point[0], &im_line_point[1]);

      // Compute the error
      residuals[0] = im_projection[0] - im_line_point[0];
      residuals[1] = im_projection[1] - im_line_point[1];

      return true;
    }

private:
    const double qw_;
    const double qx_;
    const double qy_;
    const double qz_;
    const double tx_;
    const double ty_;
    const double tz_;
    const double a_;
    const double b_;
    const double c_;
};

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_COST_FUNCTIONS_H_
