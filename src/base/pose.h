// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
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
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#ifndef COLMAP_SRC_BASE_POSE_H_
#define COLMAP_SRC_BASE_POSE_H_

#include <vector>

#include <Eigen/Core>

#include "util/alignment.h"
#include "util/types.h"

namespace colmap {

// Convert 3D rotation matrix to Quaternion representation.
//
// @param rot_mat        3x3 rotation matrix.
//
// @return               Unit Quaternion rotation coefficients (w, x, y, z).
Eigen::Vector4d RotationMatrixToQuaternion(const Eigen::Matrix3d& rot_mat);

// Convert Quaternion representation to 3D rotation matrix.
//
// @param qvec           Unit Quaternion rotation coefficients (w, x, y, z).
//
// @return               3x3 rotation matrix.
Eigen::Matrix3d QuaternionToRotationMatrix(const Eigen::Vector4d& qvec);

// Compose the Quaternion vector corresponding to a  identity transformation.
inline Eigen::Vector4d ComposeIdentityQuaternion();

// Normalize Quaternion vector.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
//
// @return              Unit Quaternion rotation coefficients (w, x, y, z).
Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec);

// Invert Quaternion vector to return Quaternion of inverse rotation.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
//
// @return              Inverse Quaternion rotation coefficients (w, x, y, z).
Eigen::Vector4d InvertQuaternion(const Eigen::Vector4d& qvec);

// Transform point by quaternion rotation.
//
// @param qvec          Quaternion rotation coefficients (w, x, y, z).
// @param point         Point to rotate.
//
// @return              Rotated point.
Eigen::Vector3d QuaternionRotatePoint(const Eigen::Vector4d& qvec,
                                      const Eigen::Vector3d& point);

// Extract camera projection center from projection matrix, i.e. the projection
// center in world coordinates `-R^T t`.
Eigen::Vector3d ProjectionCenterFromMatrix(
    const Eigen::Matrix3x4d& proj_matrix);

// Extract camera projection center from projection parameters.
//
// @param qvec           Unit Quaternion rotation coefficients (w, x, y, z).
// @param tvec           3x1 translation vector.
//
// @return               3x1 camera projection center.
Eigen::Vector3d ProjectionCenterFromPose(const Eigen::Vector4d& qvec,
                                         const Eigen::Vector3d& tvec);

// Invert transformation of the pose.
// @param qvec, tvec          Input camera pose.
// @param inv_qvec, inv_tvec  Inverse camera pose.
void InvertPose(const Eigen::Vector4d& qvec, const Eigen::Vector3d& tvec,
                Eigen::Vector4d* inv_qvec, Eigen::Vector3d* inv_tvec);

// Linearly interpolate camera pose.
//
// @param qvec1, tvec1      Camera pose at t0 = 0.
// @param qvec2, tvec2      Camera pose at t1 = 1.
// @param t                 Interpolation time.
// @param qveci, tveci      Camera pose at time t.
void InterpolatePose(const Eigen::Vector4d& qvec1, const Eigen::Vector3d& tvec1,
                     const Eigen::Vector4d& qvec2, const Eigen::Vector3d& tvec2,
                     const double t, Eigen::Vector4d* qveci,
                     Eigen::Vector3d* tveci);

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

Eigen::Vector4d ComposeIdentityQuaternion() {
  return Eigen::Vector4d(1, 0, 0, 0);
}

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_POSE_H_
