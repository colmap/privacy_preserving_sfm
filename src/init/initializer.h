// Copyright (c) 2020, ETH Zurich.
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
//     * Neither the name of ETH Zurich nor the names of its contributors may be
//       used to endorse or promote products derived from this software without
//       specific prior written permission.
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

#ifndef COLMAP_SRC_INIT_INITIALIZER_H
#define COLMAP_SRC_INIT_INITIALIZER_H

#include <Eigen/Dense>

#include <vector>

#include "util/types.h"
#include "feature/types.h"

namespace colmap {
namespace init {

typedef Eigen::Matrix3x4d Pose;

struct InitOptions {
  // Minimum mean triangulation angle (in rad)
  double min_tri_angle = 0.1;

  // Minimum number of inliers to accept a model
  double min_num_inliers = 6;

  // Maximum reprojection error in ransac (in normalized coordinates)
  double max_error = 0.005;
};



class PlanarOffsetEstimator {
    public:
    struct Reconstruction {
        std::vector<Pose> cams;
        std::vector<Eigen::Vector3d> X;
    };
    typedef std::vector<Reconstruction> ReconstructionVector;

    PlanarOffsetEstimator(const std::vector<Pose> &poses, const std::vector<std::vector<Eigen::Vector3d>> &lines, const std::vector<Eigen::Matrix3d> Rg, const double inlier_threshold) : poses_(poses), lines_(lines), Rg_(Rg), inlier_threshold_(inlier_threshold) {}

    int min_sample_size() const {
        return 3;
    }

    int non_minimal_sample_size() const {
        return 20;
    }
    
    int num_data() const {
        return lines_[0].size();
    }
    
    
    int MinimalSolver(const std::vector<int>& sample,
                        ReconstructionVector* models) const;

    int NonMinimalSolver(const std::vector<int>& sample,
                        Reconstruction* model) const;
    
    double EvaluateModelOnPoint(const Reconstruction& model, int i) const;

    void LeastSquares(const std::vector<int>& sample, Reconstruction* model) const;

    private:
        std::vector<Pose> poses_;
        std::vector<std::vector<Eigen::Vector3d>> lines_;
        std::vector<Eigen::Matrix3d> Rg_;
        const double inlier_threshold_;
        
};

bool initialize_reconstruction(
    const std::vector<FeatureLines> &lines,
    const std::vector<Eigen::Vector3d> &gravity,
    const InitOptions& options,
    std::vector<Pose> *output,
    double* inlier_ratio);


} //namespace init
} //namespace colmap

#endif // COLMAP_SRC_INIT_INITIALIZER_H