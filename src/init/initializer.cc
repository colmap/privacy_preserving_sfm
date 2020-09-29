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

#include "init/initializer.h"
#include "util/misc.h"
#include <iostream>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "init/sfm2d.h"
#include <RansacLib/ransac.h>

namespace colmap {
namespace init {

void lift_camera(const Pose2d &pose2d, Pose &pose3d) {
    pose3d.setZero();
    pose3d(0,0) = pose2d(0,0);
    pose3d(0,2) = pose2d(0,1);
    pose3d(2,0) = pose2d(1,0);
    pose3d(2,2) = pose2d(1,1);
    pose3d(1,1) = 1.0;

    pose3d(0,3) = pose2d(0,2);
    pose3d(2,3) = pose2d(1,2);
}

bool initialize_reconstruction(
    const std::vector<FeatureLines> &lines,
    const std::vector<Eigen::Vector3d> &gravity,
    const InitOptions& options,
    std::vector<Pose> *output,
    double* inlier_ratio) {
    // Step 1. Align with gravity

    *inlier_ratio = 0;

    std::vector<std::vector<Eigen::Vector2d>> x(4); // gravity aligned lines are converted to points for the 2d sfm
    std::vector<std::vector<Eigen::Vector3d>> lines_r(4); // random lines

    std::vector<Eigen::Matrix3d> Rg(4);

    for(int i = 0; i < 4; ++i) {
        Rg[i] = Eigen::Quaterniond::FromTwoVectors(gravity[i], Eigen::Vector3d{0.0, 1.0, 0.0});

        const Eigen::Vector3d e3{0.0, 0.0, 1.0};
      std::cout << (e3.transpose() - Rg[i].row(1).dot(e3) * Rg[i].row(1)).normalized() << " || " << Rg[i].row(2) << "\n";

        const int num_lines = lines[i].size();
        for(int j = 0; j < num_lines; ++j) {
            Eigen::Vector3d l =  lines[i][j].Line();
           
            if(lines[i][j].IsAligned()) {
                // We only pre-rotate the aligned lines.
                l = Rg[i] * l;

                CHECK_NEAR(l(1), 0.0, 1e-6);
                

                Eigen::Vector2d xl = Eigen::Vector2d{l(2),-l(0)};
                if( xl(1) < 0 )
                    xl *= -1.0;

                xl.normalize();
                x[i].push_back( xl );                
            } else {
                lines_r[i].push_back( l );
            }
        }
    }
    std::cout << "sizes: " << x[0].size() << ", " << x[1].size() << ", " << x[2].size() << ", " << x[3].size() << "\n";
    CHECK_EQ(x[0].size(), x[1].size());
    CHECK_EQ(x[0].size(), x[2].size());
    CHECK_EQ(x[0].size(), x[3].size());
    CHECK_EQ(lines_r[0].size(), lines_r[1].size());
    CHECK_EQ(lines_r[0].size(), lines_r[2].size());
    CHECK_EQ(lines_r[0].size(), lines_r[3].size());


    std::cout << x[0].size() << " aligned correspondences, " << lines_r[0].size() << " unaligned.\n";


    const double normalized_reproj_error_threshold = options.max_error;

    ransac_lib::LORansacOptions ransac_options;
    ransac_options.final_least_squares_ = true;
    ransac_options.min_num_iterations_ = 1000;
    ransac_options.squared_inlier_threshold_ = normalized_reproj_error_threshold;

    // Estimate a four view reconstruction
    FourView2dEstimator solver(x[0], x[1], x[2], x[3], ransac_options.squared_inlier_threshold_);
    ransac_lib::LocallyOptimizedMSAC<FourView2dEstimator::Reconstruction, FourView2dEstimator::ReconstructionVector, FourView2dEstimator> fourview_ransac;
    ransac_lib::RansacStatistics ransac_stats;
    FourView2dEstimator::Reconstruction rec;
    int inliers = fourview_ransac.EstimateModel(ransac_options, solver, &rec, &ransac_stats);


    if(inliers < options.min_num_inliers)
        return false;

    std::vector<Eigen::Vector2d> X_2d = rec.X;

    Pose2d P1, P2, P3, P4;
    P1 = rec.cams[0];
    P2 = rec.cams[1];
    P3 = rec.cams[2];
    P4 = rec.cams[3];
    
    double total_error = 0;
    for(int i = 0; i < ransac_stats.inlier_indices.size(); ++i) {
        Eigen::Vector2d z1 = (rec.cams[0] * rec.X[ransac_stats.inlier_indices[i]].homogeneous()).normalized();
        Eigen::Vector2d z2 = (rec.cams[1] * rec.X[ransac_stats.inlier_indices[i]].homogeneous()).normalized();
        Eigen::Vector2d z3 = (rec.cams[2] * rec.X[ransac_stats.inlier_indices[i]].homogeneous()).normalized();
        Eigen::Vector2d z4 = (rec.cams[3] * rec.X[ransac_stats.inlier_indices[i]].homogeneous()).normalized();    

        double err1 = 1.0 - z1.dot(x[0][ransac_stats.inlier_indices[i]]);
        double err2 = 1.0 - z2.dot(x[1][ransac_stats.inlier_indices[i]]);
        double err3 = 1.0 - z3.dot(x[2][ransac_stats.inlier_indices[i]]);
        double err4 = 1.0 - z4.dot(x[3][ransac_stats.inlier_indices[i]]);

        total_error += err1 + err2 + err3 + err4;

    }

    double angle_sum = 0;
    // Compute the mean triangulation angle. This should give an indication on how stable the triangulation is
    for(int i = 0; i < ransac_stats.inlier_indices.size(); ++i) {
      int inlier_index = ransac_stats.inlier_indices[i];

      // Find the maximum angle between views
      // We check the triangulation angles between the first three cameras because these are used for triangulation
      int min_cam1, min_cam2;
      double min_angle = std::numeric_limits<double>::max();
      for (int cam1_idx = 0; cam1_idx < 3; ++cam1_idx) {
        const Eigen::Vector2d cam1_center = -rec.cams[cam1_idx].leftCols<2>().transpose() * rec.cams[cam1_idx].rightCols<1>();
        for (int cam2_idx = cam1_idx+1; cam2_idx < 3; ++cam2_idx) {
          const Eigen::Vector2d cam2_center = -rec.cams[cam2_idx].leftCols<2>().transpose() * rec.cams[cam2_idx].rightCols<1>();
          const Eigen::Vector2d v1 = cam1_center - rec.X[inlier_index];
          const Eigen::Vector2d v2 = cam2_center - rec.X[inlier_index];
          const double angle = acos(v1.normalized().dot(v2.normalized()));

          if (angle < min_angle) {
            min_angle = angle;
            min_cam1 = cam1_idx;
            min_cam2 = cam2_idx;
          }
        }
      }
      angle_sum += min_angle;
    }

    const double mean_tri_angle = (angle_sum / ransac_stats.inlier_indices.size()) / M_PI * 180.0;
    std::cout << "Mean minimum triangulation angle for 2D model: " << mean_tri_angle << " (" << ransac_stats.inlier_indices.size() << " inliers)\n";

    if (mean_tri_angle < options.min_tri_angle) {
      return false;
    }

    // Convert 2d cameras to 3d cameras (only t_y is missing now)
    std::vector<Pose> poses(4);
    lift_camera(P1, poses[0]);
    lift_camera(P2, poses[1]);
    lift_camera(P3, poses[2]);
    lift_camera(P4, poses[3]);
        
    // Step 4. Estimate out-of-plane translations
    ransac_lib::LORansacOptions planar_offset_options;
    planar_offset_options.final_least_squares_ = true;
    planar_offset_options.min_num_iterations_ = 1000;    
    planar_offset_options.squared_inlier_threshold_ = normalized_reproj_error_threshold; // This error is in normalized image plane!

    PlanarOffsetEstimator planar_offset_solver(poses, lines_r, Rg, planar_offset_options.squared_inlier_threshold_);

    ransac_lib::LocallyOptimizedMSAC<PlanarOffsetEstimator::Reconstruction, PlanarOffsetEstimator::ReconstructionVector, PlanarOffsetEstimator> planar_offset_ransac;
    PlanarOffsetEstimator::Reconstruction rec3d;
    
    inliers = planar_offset_ransac.EstimateModel(planar_offset_options, planar_offset_solver, &rec3d, &ransac_stats);
    if(inliers < options.min_tri_angle)
        return false;

    *output = rec3d.cams;

    *inlier_ratio = ransac_stats.inlier_ratio;

    return inliers >= options.min_num_inliers;
}



void four_view_triangulate(const std::vector<Pose> &cams, const std::vector<std::vector<Eigen::Vector3d>> &lines, std::vector<Eigen::Vector3d> *X)
{
    Eigen::Matrix<double, 4,3> A;
    Eigen::Matrix<double, 4,1> b;
    const int num_points = static_cast<int>(lines[0].size());
    for(int i = 0; i < num_points; ++i) {
        
        for(int j = 0; j < 4; ++j) {
            A.row(j) = lines[j][i].transpose() * cams[j].block<3,3>(0,0);
            b(j) = -lines[j][i].transpose() * cams[j].block<3,1>(0,3);
        }
        X->push_back( A.colPivHouseholderQr().solve(b) );

    }
}


int PlanarOffsetEstimator::MinimalSolver(const std::vector<int>& sample, ReconstructionVector* models) const {
    Eigen::Matrix<double, Eigen::Dynamic, 3> A(sample.size(), 3);
    Eigen::Matrix<double, Eigen::Dynamic, 1> b(sample.size(), 1);

    const int sample_size = static_cast<int>(sample.size());
    for(int i = 0; i < sample_size; ++i) {
        Eigen::Matrix3d A0;
        Eigen::Matrix<double, 3, 4> B0;
        B0.setZero();

        for(int j = 1; j < 4; ++j) {
            Eigen::Vector3d lg = Rg_[j] * lines_[j][sample[i]];

            A0.row(j-1) = lg.transpose() * poses_[j].block<3,3>(0,0);
            B0(j-1, j-1) = lg(1);
            B0(j-1, 3) = lg(0) * poses_[j](0,3) + lg(2) * poses_[j](2,3);
        }

        B0 = Rg_[0].transpose() * A0.partialPivLu().solve(B0);

        A.row(i) = lines_[0][sample[i]].transpose() * B0.block<3,3>(0,0);
        b(i) = -lines_[0][sample[i]].transpose() * B0.block<3,1>(0,3);
    }

    Eigen::Vector3d tt = A.colPivHouseholderQr().solve(b);

    Reconstruction rec;
    rec.cams.resize(4);

    for(int i = 0; i < 4; ++i) {
        rec.cams[i] = poses_[i];
        if(i > 0)
            rec.cams[i](1,3) = tt(i-1);

        rec.cams[i] = Rg_[i].transpose() * rec.cams[i];
    }

    four_view_triangulate(rec.cams, lines_, &rec.X);

    models->clear();
    models->push_back(rec);


    return 1;

}

int PlanarOffsetEstimator::NonMinimalSolver(const std::vector<int>& sample, Reconstruction* model) const
{
    ReconstructionVector models;
    MinimalSolver(sample, &models);

    double best_score = std::numeric_limits<double>::max();

    for(int i = 0; i < models.size(); ++i) {
        double score = 0;
        const int kNumData = num_data();
        for(int j = 0; j < kNumData; ++j) {
          const double tmp_score = EvaluateModelOnPoint(models[i], j);
          score += std::min(tmp_score, inlier_threshold_);
        }

        if(score < best_score) {
            best_score = score;
            *model = models[i];
        }
    }
    if(models.size() > 0) {
        LeastSquares(sample, model);
        return 1;
    } else {
        return 0;
    }
}
    
double PlanarOffsetEstimator::EvaluateModelOnPoint(const Reconstruction& model, int i) const
{
    
    Eigen::Vector3d z1 = (model.cams[0] * model.X[i].homogeneous());
    Eigen::Vector3d z2 = (model.cams[1] * model.X[i].homogeneous());
    Eigen::Vector3d z3 = (model.cams[2] * model.X[i].homogeneous());
    Eigen::Vector3d z4 = (model.cams[3] * model.X[i].homogeneous());

    if(z1(2) < 0 || z2(2) < 0 || z3(2) < 0 || z4(2) < 0) {
        return 100000.0;
    }
    z1 /= z1(2);
    z2 /= z2(2);
    z3 /= z3(2);
    z4 /= z4(2);

    double err1 = std::abs(lines_[0][i].dot(z1) / lines_[0][i].block<2,1>(0,0).norm());
    double err2 = std::abs(lines_[1][i].dot(z2) / lines_[1][i].block<2,1>(0,0).norm());
    double err3 = std::abs(lines_[2][i].dot(z3) / lines_[2][i].block<2,1>(0,0).norm());
    double err4 = std::abs(lines_[3][i].dot(z4) / lines_[3][i].block<2,1>(0,0).norm());

    return std::max(err1, std::max(err2, std::max(err3, err4)));   
}



class SimpleBundleAdjustmentCostFunction {
public:
    explicit SimpleBundleAdjustmentCostFunction(const Eigen::Vector3d& line)
            : a_(line(0)), b_(line(1)), c_(line(2)) {
      
    }

    static ceres::CostFunction* Create(const Eigen::Vector3d& line) {
        return (new ceres::AutoDiffCostFunction<SimpleBundleAdjustmentCostFunction, 1, 4, 3, 3>(
                new SimpleBundleAdjustmentCostFunction(line)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
                    const T* const point3D, T* residuals) const {
      // Rotate and translate.
      
      Eigen::Matrix<T, 3, 1> projection;
      Eigen::Matrix<T, 3, 1> translation;      
      Eigen::Matrix<T, 3, 1> point;


      point << point3D[0], point3D[1], point3D[2];
      translation << tvec[0], tvec[1], tvec[2];

      Eigen::Quaternion<T> q;
      q.coeffs()(0) = qvec[0];
      q.coeffs()(1) = qvec[1];
      q.coeffs()(2) = qvec[2];
      q.coeffs()(3) = qvec[3];
      Eigen::Matrix<T, 3, 3> R = q.toRotationMatrix();

      projection = R * point + translation;

      // Project to image plane.
      projection(0) /= projection(2);
      projection(1) /= projection(2);

      // Compute the closest point on the line
      T alpha = T(a_) * projection(0) + T(b_) * projection(1) + T(c_);

      residuals[0] = alpha / ceres::sqrt(a_*a_ + b_*b_);
      return true;
    }

private:
  const double a_;
  const double b_;  
  const double c_;  

};

void bundle_adjust_simple(std::vector<Pose> &cams, const std::vector<std::vector<Eigen::Vector3d>> &lines, std::vector<Eigen::Vector3d> &X) {    

    if(lines.size() == 0)
        return;
    if(lines[0].size() < 17) // only bundle when there are enough points to make it worthwhile
        return;
    

    // Build the problem.
    ceres::Problem problem;

    ceres::LossFunction* loss = nullptr;

    std::vector<Eigen::Quaterniond> cam_qvec;
    std::vector<Eigen::Vector3d> cam_tvec;
    for(int i = 0; i < cams.size(); ++i) {      
        cam_qvec.push_back(Eigen::Quaterniond(cams[i].block<3,3>(0,0)));
        cam_tvec.push_back(cams[i].block<3,1>(0,3));
    }

    for(int i = 0; i < cams.size(); ++i) {
        for (int j = 0; j < lines[i].size(); ++j) {
            ceres::CostFunction *cost_function = SimpleBundleAdjustmentCostFunction::Create(lines[i][j]);
            problem.AddResidualBlock(cost_function, loss, cam_qvec[i].coeffs().data(), cam_tvec[i].data(), X[j].data());
        }
    }

    // Setup parameterization
    for(int i = 0; i < cams.size(); ++i) {
        problem.SetParameterization(cam_qvec[i].coeffs().data(), new ceres::EigenQuaternionParameterization());
    }
    // Fix scale
    problem.SetParameterBlockConstant(cam_qvec[0].coeffs().data());
    problem.SetParameterBlockConstant(cam_tvec[0].data());

    problem.SetParameterization(cam_tvec[1].data(), new ceres::HomogeneousVectorParameterization(3));


    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = false;
    //options.logging_type = ceres::LoggingType::SILENT;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-10;
    options.parameter_tolerance = 1e-10;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
	
    if(summary.termination_type == ceres::FAILURE) {			
        std::cout << summary.FullReport() << "\n";
    }

    for(int i = 0; i < cams.size(); ++i) {
        cams[i].block<3,3>(0,0) = cam_qvec[i].toRotationMatrix();
        cams[i].block<3,1>(0,3) = cam_tvec[i];
    }
}




 void PlanarOffsetEstimator::LeastSquares(const std::vector<int>& sample, Reconstruction* model) const {
  return;
    std::vector<std::vector<Eigen::Vector3d>> lines(4);
    std::vector<Eigen::Vector3d> X;
    for(int i = 0; i < sample.size(); ++i) {
        lines[0].push_back(lines_[0][sample[i]]);
        lines[1].push_back(lines_[1][sample[i]]);
        lines[2].push_back(lines_[2][sample[i]]);
        lines[3].push_back(lines_[3][sample[i]]);
        X.push_back(model->X[sample[i]]);
    }

    bundle_adjust_simple(model->cams, lines, X);

    for(int i = 0; i < sample.size(); ++i) {
        model->X[sample[i]] = X[i];
    } 
}

} // namespace init
} // namespace colmap
