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

#ifndef COLMAP_SRC_INIT_SFM_2D_H
#define COLMAP_SRC_INIT_SFM_2D_H

#include <vector>
#include <Eigen/Dense>

namespace colmap {
namespace init {

    typedef Eigen::Matrix<double, 2, 3> Pose2d;
    typedef Eigen::Matrix<double, 8, 1> TrifocalTensor;
    typedef std::vector<TrifocalTensor> TrifocalTensorVector;
    typedef std::vector<Pose2d> Pose2dVector;

    class FourView2dEstimator {
        public:

        struct Reconstruction {
            std::vector<Pose2d> cams;
            std::vector<Eigen::Vector2d> X;
        };

        typedef std::vector<Reconstruction> ReconstructionVector;

        FourView2dEstimator(const std::vector<Eigen::Vector2d> &x1, const std::vector<Eigen::Vector2d> &x2,
                            const std::vector<Eigen::Vector2d> &x3, const std::vector<Eigen::Vector2d> &x4,
                            const double inlier_threshold) :
            x1_(x1), x2_(x2), x3_(x3), x4_(x4), inlier_threshold_(inlier_threshold)
        {
            for(int i = 0; i < x1_.size(); ++i) {
                x1_[i].normalize();
                x2_[i].normalize();
                x3_[i].normalize();
                x4_[i].normalize();
            }
        }

        int min_sample_size() const {
            return 5;
        }
    
        int non_minimal_sample_size() const {
            return 2 * min_sample_size();
        }
        
        int num_data() const {
            return x1_.size();
        }
        
        
        int MinimalSolver(const std::vector<int>& sample,
                            ReconstructionVector* models) const;
        
        int AbsPoseSolver(const std::vector<int>& sample, const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &X_, Pose2d* model) const;

        int NonMinimalSolver(const std::vector<int>& sample,
                            Reconstruction* model) const;
        
        double EvaluateModelOnPoint(const Reconstruction& model, int i) const;

        void LeastSquares(const std::vector<int>& sample, Reconstruction* model) const;

        private:
            std::vector<Eigen::Vector2d> x1_, x2_, x3_, x4_;
            const double inlier_threshold_;
            
    };

    class AbsolutePose2dEstimator {
        public:

        AbsolutePose2dEstimator(const std::vector<Eigen::Vector2d> &x, 
            const std::vector<Eigen::Vector2d> &X) : x_(x), X_(X)
        {
            for(Eigen::Vector2d &xi : x_) {
                xi.normalize();
            }
        }

        int min_sample_size() const {
            return 3;
        }
    
        int non_minimal_sample_size() const {
            return 2 * min_sample_size();
        }
        
        int num_data() const {
            return x_.size();
        }
        
        
        int MinimalSolver(const std::vector<int>& sample,
                            Pose2dVector* models) const {
            Pose2d cam;
            
            NonMinimalSolver(sample, &cam);
            models->clear();
            models->push_back(cam);
            return 1;
        }

        int NonMinimalSolver(const std::vector<int>& sample,
                            Pose2d* model) const;
        
        double EvaluateModelOnPoint(const Pose2d& model, int i) const;

        void LeastSquares(const std::vector<int>& sample, Pose2d* model) const {
            NonMinimalSolver(sample, model);
        }

        private:
            std::vector<Eigen::Vector2d> x_, X_;
            
    };

} // namespace init
} // namespace colmap


#endif // COLMAP_SRC_INIT_SFM_2D_H