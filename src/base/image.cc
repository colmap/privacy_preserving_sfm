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

#include "base/image.h"

#include "base/pose.h"
#include "base/projection.h"

namespace colmap {
namespace {

static const double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

Image::Image()
    : image_id_(kInvalidImageId),
      name_(""),
      camera_id_(kInvalidCameraId),
      registered_(false),
      num_points3D_(0),
      num_observations_(0),
      num_correspondences_(0),
      num_visible_points3D_(0),
      qvec_(1.0, 0.0, 0.0, 0.0),
      tvec_(0.0, 0.0, 0.0),
      qvec_prior_(kNaN, kNaN, kNaN, kNaN),
      tvec_prior_(kNaN, kNaN, kNaN),
      gravity_(kNaN, kNaN, kNaN) {}

void Image::SetLines(const FeatureLines& lines) {
  num_correspondences_have_point3D_.resize(lines.size(), 0);
  lines_ = lines;
}

void Image::SetPoint3DForLine(const point2D_t line_idx,
                                 const point3D_t point3D_id) {
  CHECK_NE(point3D_id, kInvalidPoint3DId);
  class FeatureLine& line = lines_.at(line_idx);
  if (!line.HasPoint3D()) {
    num_points3D_ += 1;
  }
  line.SetPoint3DId(point3D_id);
}

void Image::ResetPoint3DForLine(const point2D_t line_idx) {
  class FeatureLine& line = lines_.at(line_idx);
  if (line.HasPoint3D()) {
    line.SetPoint3DId(-1);
    num_points3D_ -= 1;
  }
}

bool Image::HasPoint3D(const point3D_t point3D_id) const {
  return std::find_if(lines_.begin(), lines_.end(),
                      [point3D_id](const class FeatureLine& line) {
                        return line.Point3DId() == point3D_id;
                      }) != lines_.end();
}

void Image::IncrementCorrespondenceHasPoint3D(const point2D_t line_idx) {
  num_correspondences_have_point3D_[line_idx] += 1;
  if (num_correspondences_have_point3D_[line_idx] == 1) {
    num_visible_points3D_ += 1;
  }

  assert(num_visible_points3D_ <= num_observations_);
}

void Image::DecrementCorrespondenceHasPoint3D(const point2D_t point2D_idx) {
  num_correspondences_have_point3D_[point2D_idx] -= 1;
  if (num_correspondences_have_point3D_[point2D_idx] == 0) {
    num_visible_points3D_ -= 1;
  }

  assert(num_visible_points3D_ <= num_observations_);
}

void Image::NormalizeQvec() { qvec_ = NormalizeQuaternion(qvec_); }

Eigen::Matrix3x4d Image::ProjectionMatrix() const {
  return ComposeProjectionMatrix(qvec_, tvec_);
}

Eigen::Matrix3x4d Image::InverseProjectionMatrix() const {
  return InvertProjectionMatrix(ComposeProjectionMatrix(qvec_, tvec_));
}

Eigen::Matrix3d Image::RotationMatrix() const {
  return QuaternionToRotationMatrix(qvec_);
}

Eigen::Vector3d Image::ProjectionCenter() const {
  return ProjectionCenterFromPose(qvec_, tvec_);
}

Eigen::Vector3d Image::ViewingDirection() const {
  return RotationMatrix().row(2);
}

}  // namespace colmap
