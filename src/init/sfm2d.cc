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

#include "sfm2d.h"
#include <ceres/ceres.h>

namespace colmap {
namespace init {

namespace {

class BundleAdjustment2DCostFunction {
 public:
  explicit BundleAdjustment2DCostFunction(const Eigen::Vector2d& x)
      : x_(x(0)), y_(x(1)) {
    const double norm = std::sqrt(x_ * x_ + y_ * y_);
    CHECK_NEAR(norm, 1.0, 1e-6);
  }

  static ceres::CostFunction* Create(const Eigen::Vector2d& x) {
    return (new ceres::AutoDiffCostFunction<BundleAdjustment2DCostFunction, 1, 2, 2, 2>(
        new BundleAdjustment2DCostFunction(x)));
  }

  template <typename T>
  bool operator()(const T* const qvec, const T* const tvec,
                  const T* const point2D, T* residuals) const {
    // Rotate and translate.
    T projection[2];

    projection[0] = (qvec[0] * point2D[0]) - (qvec[1] * point2D[1]) + tvec[0];
    projection[1] = (qvec[1] * point2D[0]) + (qvec[0] * point2D[1]) + tvec[1];

    // Project to image plane.
    residuals[0] = (projection[0] / projection[1]) - (T(x_) / T(y_));

    return true;
  }

 private:
  const double x_;
  const double y_;
};

void optimize_points2d(std::vector<Pose2d> &cams, const std::vector<std::vector<Eigen::Vector2d>> &x, std::vector<Eigen::Vector2d> &X) {

  if(x.empty())
    return;

  // Build the problem.
  ceres::Problem problem;

  ceres::LossFunction* loss = nullptr;

  std::vector<Eigen::Vector2d> cam_qvec;
  std::vector<Eigen::Vector2d> cam_tvec;
  for (const auto& cam : cams) {
    cam_qvec.emplace_back(Eigen::Vector2d(cam(0,0), cam(1,0)));
    cam_tvec.emplace_back(Eigen::Vector2d(cam(0,2), cam(1,2)));
  }

  for(int i = 0; i < cams.size(); ++i) {
    for (int j = 0; j < x[i].size(); ++j) {
      ceres::CostFunction *cost_function = BundleAdjustment2DCostFunction::Create(x[i][j]);
      problem.AddResidualBlock(cost_function, loss, cam_qvec[i].data(), cam_tvec[i].data(), X[j].data());
    }
  }

  // Setup parameterization
  for(int i = 0; i < cams.size(); ++i) {
    problem.SetParameterBlockConstant(cam_qvec[i].data());
    problem.SetParameterBlockConstant(cam_tvec[i].data());
  }

  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  options.function_tolerance = 1e-10;
  options.gradient_tolerance = 1e-10;
  options.parameter_tolerance = 1e-10;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  if(summary.termination_type == ceres::FAILURE) {
    std::cout << summary.FullReport() << "\n";
  }
}

void bundle_adjust2d(std::vector<Pose2d> &cams, const std::vector<std::vector<Eigen::Vector2d>> &x, std::vector<Eigen::Vector2d> &X) {

  if(x.size() == 0)
    return;
  if(x[0].size() < 10) // only bundle when there are enough points to make it worthwhile
    return;

  // Build the problem.
  ceres::Problem problem;

  ceres::LossFunction* loss = nullptr;

  std::vector<Eigen::Vector2d> cam_qvec;
  std::vector<Eigen::Vector2d> cam_tvec;
  for (const auto& cam : cams) {
    cam_qvec.emplace_back(Eigen::Vector2d(cam(0,0), cam(1,0)));
    cam_tvec.emplace_back(Eigen::Vector2d(cam(0,2), cam(1,2)));
  }

  for(int i = 0; i < cams.size(); ++i) {
    for (int j = 0; j < x[i].size(); ++j) {
      ceres::CostFunction *cost_function = BundleAdjustment2DCostFunction::Create(x[i][j]);
      problem.AddResidualBlock(cost_function, loss, cam_qvec[i].data(), cam_tvec[i].data(), X[j].data());
    }
  }

  // Setup parameterization
  for(int i = 0; i < cams.size(); ++i) {
    problem.SetParameterization(cam_qvec[i].data(), new ceres::HomogeneousVectorParameterization(2));
  }
  // Fix scale and coordinate system
  problem.SetParameterBlockConstant(cam_qvec[0].data());
  problem.SetParameterBlockConstant(cam_tvec[0].data());
  problem.SetParameterization(cam_tvec[1].data(), new ceres::HomogeneousVectorParameterization(2));


  ceres::Solver::Options options;
  options.minimizer_progress_to_stdout = false;
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
    cams[i](0,0) = cam_qvec[i](0);
    cams[i](0,1) =-cam_qvec[i](1);
    cams[i](1,0) = cam_qvec[i](1);
    cams[i](1,1) = cam_qvec[i](0);
    cams[i](0,2) = cam_tvec[i](0);
    cams[i](1,2) = cam_tvec[i](1);
  }
}

/* Estimates a 3x3 transform such that P2 and P3 become calibrated 2d cameras. Assumes P1 is [eye(2) [0;0]] */
void metric_upgrade(const Pose2d &P2, const Pose2d &P3, Eigen::Matrix3d &H) {

  Eigen::Matrix<double, 4, 2> A;
  Eigen::Matrix<double, 4, 1> b;

  A << P2(0,2), -P2(1,2), P2(1,2), P2(0,2), P3(0,2), -P3(1,2), P3(1,2), P3(0,2);
  b << P2(1,1)-P2(0,0), -P2(0,1)-P2(1,0), P3(1,1)-P3(0,0), -P3(0,1)-P3(1,0);

  Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);

  H.setIdentity();
  H(2,0) = x(0);
  H(2,1) = x(1);
}


void three_view_triangulate2d(const Pose2d &P1, const Pose2d &P2, const Pose2d &P3, const std::vector<Eigen::Vector2d> &x1,
                              const std::vector<Eigen::Vector2d> &x2, const std::vector<Eigen::Vector2d> &x3, std::vector<Eigen::Vector2d> *X)
{
  Eigen::Matrix<double, 3,2> A;
  Eigen::Matrix<double, 3,1> b;
  const int num_points = static_cast<int>(x1.size());
  for(int i = 0; i < num_points; ++i) {
    A.row(0) = x1[i](0) * P1.block<1,2>(1,0) - x1[i](1) * P1.block<1,2>(0,0);
    b(0) = x1[i](1) * P1(0,2) - x1[i](0) * P1(1,2);

    A.row(1) = x2[i](0) * P2.block<1,2>(1,0) - x2[i](1) * P2.block<1,2>(0,0);
    b(1) = x2[i](1) * P2(0,2) - x2[i](0) * P2(1,2);

    A.row(2) = x3[i](0) * P3.block<1,2>(1,0) - x3[i](1) * P3.block<1,2>(0,0);
    b(2) = x3[i](1) * P3(0,2) - x3[i](0) * P3(1,2);

    X->push_back( A.colPivHouseholderQr().solve(b) );

  }
}

void trifocal_tensor_coord_change(const TrifocalTensor &T, const Eigen::Matrix2d &A1, const Eigen::Matrix2d &A2, const Eigen::Matrix2d &A3, TrifocalTensor *Tout) {
  (*Tout)(0) = A3(0)*(A2(0)*(A1(0)*T(0) + A1(1)*T(1)) + A2(1)*(A1(0)*T(2) + A1(1)*T(3))) + A3(1)*(A2(0)*(A1(0)*T(4) + A1(1)*T(5)) + A2(1)*(A1(0)*T(6) + A1(1)*T(7)));
  (*Tout)(1) = A3(0)*(A2(0)*(A1(2)*T(0) + A1(3)*T(1)) + A2(1)*(A1(2)*T(2) + A1(3)*T(3))) + A3(1)*(A2(0)*(A1(2)*T(4) + A1(3)*T(5)) + A2(1)*(A1(2)*T(6) + A1(3)*T(7)));
  (*Tout)(2) = A3(0)*(A2(2)*(A1(0)*T(0) + A1(1)*T(1)) + A2(3)*(A1(0)*T(2) + A1(1)*T(3))) + A3(1)*(A2(2)*(A1(0)*T(4) + A1(1)*T(5)) + A2(3)*(A1(0)*T(6) + A1(1)*T(7)));
  (*Tout)(3) = A3(0)*(A2(2)*(A1(2)*T(0) + A1(3)*T(1)) + A2(3)*(A1(2)*T(2) + A1(3)*T(3))) + A3(1)*(A2(2)*(A1(2)*T(4) + A1(3)*T(5)) + A2(3)*(A1(2)*T(6) + A1(3)*T(7)));
  (*Tout)(4) = A3(2)*(A2(0)*(A1(0)*T(0) + A1(1)*T(1)) + A2(1)*(A1(0)*T(2) + A1(1)*T(3))) + A3(3)*(A2(0)*(A1(0)*T(4) + A1(1)*T(5)) + A2(1)*(A1(0)*T(6) + A1(1)*T(7)));
  (*Tout)(5) = A3(2)*(A2(0)*(A1(2)*T(0) + A1(3)*T(1)) + A2(1)*(A1(2)*T(2) + A1(3)*T(3))) + A3(3)*(A2(0)*(A1(2)*T(4) + A1(3)*T(5)) + A2(1)*(A1(2)*T(6) + A1(3)*T(7)));
  (*Tout)(6) = A3(2)*(A2(2)*(A1(0)*T(0) + A1(1)*T(1)) + A2(3)*(A1(0)*T(2) + A1(1)*T(3))) + A3(3)*(A2(2)*(A1(0)*T(4) + A1(1)*T(5)) + A2(3)*(A1(0)*T(6) + A1(1)*T(7)));
  (*Tout)(7) = A3(2)*(A2(2)*(A1(2)*T(0) + A1(3)*T(1)) + A2(3)*(A1(2)*T(2) + A1(3)*T(3))) + A3(3)*(A2(2)*(A1(2)*T(4) + A1(3)*T(5)) + A2(3)*(A1(2)*T(6) + A1(3)*T(7)));
}


int factorize_trifocal_tensor(const TrifocalTensor& T, Pose2d P1[2], Pose2d P2[2], Pose2d P3[2]) {

  // This factorization method degenerates sometimes (e.g. pure rotation around y axis),
  // so we do a random projective change of variables in the images
  TrifocalTensor AT;
  Eigen::Matrix2d A1; A1.setRandom();
  Eigen::Matrix2d A2; A2.setRandom();
  Eigen::Matrix2d A3; A3.setRandom();
  trifocal_tensor_coord_change(T,A1,A2,A3,&AT);

  double alpha = AT(2) * AT(7) - AT(3) * AT(6);
  double beta = AT(1) * AT(6) + AT(3) * AT(4) - AT(0) * AT(7) - AT(2) * AT(5);
  double gamma = AT(0) * AT(5) - AT(1) * AT(4);

  double aa1[2];
  double b2m4ac = beta * beta - 4.0 * alpha * gamma;

  if(b2m4ac < 0) {
    return 0; // no real factorization, (or is there? dun dun duuun. No really, we should figure this out...)
  }

  double sq = std::sqrt(b2m4ac);
  // Choose sign to avoid cancellations
  aa1[0] = (beta > 0) ? (2 * gamma) / (-beta - sq) : (2.0 * gamma) / (-beta + sq);
  aa1[1] = gamma / (alpha * aa1[0]);

  Eigen::Matrix<double, 7, 6> G;

  int n_sols = 0;
  for (int i = 0; i < 2; ++i) {
    double a1 = aa1[i];
    double s = std::sqrt(1 + a1 * a1);
    a1 /= s;
    double a2 = 1 / s;

    double rho = -(AT(1) * a2 - AT(3) * a1) / (AT(2) * a1 - AT(0) * a2);
    double b1 = rho * a1;
    double b2 = rho * a2;
    double c1 = -a2;
    double c2 = a1;

    G << 0,         AT(7) * c2, -AT(0) * c1,       0,          AT(0)* b1,             -AT(7) * a2,
        0,         0,         -AT(1) * c1,       AT(7)* c2,   AT(1)* b1,             -AT(7) * b2,
        0,        -AT(7) * c1, -AT(2) * c1,       0,          AT(2)* b1,              AT(7)* a1,
        0,         0,         -AT(3) * c1,      -AT(7) * c1,  AT(3)* b1,              AT(7)* b1,
        -AT(7) * c2, 0,         -AT(4) * c1,       0,          AT(7) * a2 + AT(4) * b1, 0,
        0,         0,         -AT(5)*c1-AT(7)*c2, 0,          AT(7) * b2 + AT(5) * b1, 0,
        AT(7)* c1,  0,         -AT(6) * c1,       0,         -AT(7) * a1 + AT(6) * b1, 0;

    Eigen::JacobiSVD<decltype(G)> svd(G, Eigen::ComputeFullV);
    Eigen::Matrix<double, 6, 1> def = svd.matrixV().rightCols(1);

    P1[n_sols] << 1, 0, 0, 0, 1, 0;
    P2[n_sols] << a1, b1, c1, a2, b2, c2;
    P3[n_sols] << def(0), def(2), def(4), def(1), def(3), def(5);

    n_sols++;
  }

  // Revert change of coordinates
  for(int i = 0; i < n_sols; ++i) {
    // P1[i] = A1 * P1[i];
    P2[i] = A2 * P2[i];
    P3[i] = A3 * P3[i];

    // Transform first camera back to [I2 0]
    P2[i].block<2,2>(0,0) *= A1.inverse();
    P3[i].block<2,2>(0,0) *= A1.inverse();
  }

  return n_sols;
}

} // namespace

double FourView2dEstimator::EvaluateModelOnPoint(const Reconstruction& model, int i) const {
    Eigen::Vector2d z1 = (model.cams[0] * model.X[i].homogeneous());
    Eigen::Vector2d z2 = (model.cams[1] * model.X[i].homogeneous());
    Eigen::Vector2d z3 = (model.cams[2] * model.X[i].homogeneous());
    Eigen::Vector2d z4 = (model.cams[3] * model.X[i].homogeneous());

    if(z1(1) < 0 || z2(1) < 0 || z3(1) < 0 || z4(1) < 0)
        return 1000000.0;

    double err1 = (x1_[i].hnormalized() - z1.hnormalized()).norm();
    double err2 = (x2_[i].hnormalized() - z2.hnormalized()).norm();
    double err3 = (x3_[i].hnormalized() - z3.hnormalized()).norm();
    double err4 = (x4_[i].hnormalized() - z4.hnormalized()).norm();

    double err = std::max(err1, std::max(err2, std::max(err3, err4)));

    return err;
}

int FourView2dEstimator::AbsPoseSolver(const std::vector<int>& sample, const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &X_, Pose2d* model) const
{
    
    Eigen::Matrix<double, Eigen::Dynamic, 2> A(sample.size(), 2);
    Eigen::Matrix<double, Eigen::Dynamic, 2> B(sample.size(), 2);

    const int sample_size = static_cast<int>(sample.size());
    for(int i = 0; i < sample_size; ++i) {
        double x1 = x_[sample[i]](0);
        double x2 = x_[sample[i]](1);
        double X1 = X_[sample[i]](0);
        double X2 = X_[sample[i]](1);        

        A.row(i) << X1*x2 - X2*x1, - X1*x1 - X2*x2;
        B.row(i) << x2, -x1;        
    }

    // Solve for t w.r.t. (a,b)
    Eigen::Matrix<double, 2, 2> C = -(B.transpose() * B).inverse() * (B.transpose() * A);
        
    Eigen::JacobiSVD<decltype(A)> svd(A + B*C, Eigen::ComputeFullV);

    Eigen::Vector2d ab = svd.matrixV().col(1);
    ab.normalize();
    Eigen::Vector2d t = C * ab;

    (*model)(0,0) = ab(0);
    (*model)(1,1) = ab(0);
    (*model)(0,1) = -ab(1);
    (*model)(1,0) = ab(1);
    (*model)(0,2) = t(0);
    (*model)(1,2) = t(1);

    // Choose sign by checking chirality of first point

    //if( x_[sample[0]].dot((*model) * X_[sample[0]].homogeneous()) < 0 )
    if( model->row(1) * X_[sample[0]].homogeneous() < 0 )
        (*model) *= -1.0;

     return 1;
}

int FourView2dEstimator::MinimalSolver(const std::vector<int>& sample, ReconstructionVector* models) const {     
    Eigen::Matrix<double, Eigen::Dynamic, 6> A(sample.size(), 6);
    const int sample_size = static_cast<int>(sample.size());
    for(int i = 0; i < sample_size; ++i) {
        double x1_1 = x1_[sample[i]](0);
        double x1_2 = x1_[sample[i]](1);
        double x2_1 = x2_[sample[i]](0);
        double x2_2 = x2_[sample[i]](1);
        double x3_1 = x3_[sample[i]](0);
        double x3_2 = x3_[sample[i]](1);
        A.row(i) << x1_1*x2_2*x3_1 - x1_2*x2_1*x3_1, x1_1*x2_1*x3_1 + x1_2*x2_2*x3_1, x1_1*x2_1*x3_2 - x1_2*x2_1*x3_1, x1_1*x2_1*x3_1 + x1_2*x2_1*x3_2, x1_1*x2_1*x3_1 + x1_1*x2_2*x3_2, x1_2*x2_1*x3_1 + x1_2*x2_2*x3_2;        
    }
    Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeFullV);
    Eigen::Matrix<double, 6, 1> t = svd.matrixV().col(5);   
    TrifocalTensor tensor;
    tensor(0) = t(1)+t(3)+t(4);
    tensor(1) = -t(2)-t(0)+t(5);
    tensor.block<6,1>(2,0) = t;
 
    Pose2d P1[2], P2[2], P3[2];
    int n_fact = factorize_trifocal_tensor(tensor,P1,P2,P3);

    if(n_fact == 0)
        return 0;
    

    models->clear();

    // there are 32 possible factorizations
    Eigen::Matrix3d H;
    for(int fact_ind = 0; fact_ind < n_fact; ++fact_ind) {
        metric_upgrade(P2[fact_ind], P3[fact_ind], H);

        Pose2d P1f = P1[fact_ind];
        Pose2d P2f = P2[fact_ind];
        Pose2d P3f = P3[fact_ind];

        P2f = P2f * H;
        P3f = P3f * H;
        P2f /= P2f.col(0).norm();
        P3f /= P3f.col(0).norm();

        double s = P2f.col(2).norm();
        P2f.col(2) /= s;
        P3f.col(2) /= s;


        for(int flip1 = 0; flip1 < 2; ++flip1) {
            for(int flip2 = 0; flip2 < 2; ++flip2) {
                for(int flip3 = 0; flip3 < 2; ++flip3) {
                    Reconstruction rec;
                    rec.cams.resize(4);
                    rec.cams[0] = P1f; rec.cams[1] = P2f; rec.cams[2] = P3f;

                    rec.cams[2].col(2) /= rec.cams[1].col(2).norm();
                    rec.cams[1].col(2) /= rec.cams[1].col(2).norm();

                    if(flip1 == 1) {
                        rec.cams[1].col(2) *= -1.0;
                        rec.cams[2].col(2) *= -1.0;
                    }

                    if(flip2 == 1) {
                        rec.cams[1] *= -1.0;
                    }

                    if(flip3 == 1) {
                        rec.cams[2] *= -1.0;
                    }

                    three_view_triangulate2d(rec.cams[0], rec.cams[1], rec.cams[2], x1_, x2_, x3_, &rec.X);

                    AbsPoseSolver(sample, x4_, rec.X, &(rec.cams[3]));

                    models->push_back(rec);
                }
            }
        }
    }

    return models->size();
}

int FourView2dEstimator::NonMinimalSolver(const std::vector<int>& sample, Reconstruction* model) const
{
    ReconstructionVector models;
    MinimalSolver(sample, &models);

    double best_score = std::numeric_limits<double>::max();

    for(int i = 0; i < models.size(); ++i) {
        double score = 0;
        const int kNumData = num_data();
        for (int j = 0; j < kNumData; ++j) {
          const double tmp_score = EvaluateModelOnPoint(models[i], j);
          score += std::min(inlier_threshold_, tmp_score);
        }
        
        if(score < best_score) {
            best_score = score;
            *model = models[i];
        }
    }
    return models.size() > 0 ? 1 : 0;
}

void FourView2dEstimator::LeastSquares(const std::vector<int>& sample, Reconstruction* model) const {

    std::vector<std::vector<Eigen::Vector2d>> x(4);
    std::vector<Eigen::Vector2d> X;
    for(int i = 0; i < sample.size(); ++i) {
        x[0].push_back(x1_[sample[i]]);
        x[1].push_back(x2_[sample[i]]);
        x[2].push_back(x3_[sample[i]]);
        x[3].push_back(x4_[sample[i]]);
        X.push_back(model->X[sample[i]]);
    }

    bundle_adjust2d(model->cams, x, X);

    for (int i = 0; i < sample.size(); ++i) {
      model->X[sample[i]] = X.at(i);
    }

    std::vector<std::vector<Eigen::Vector2d>> all_x{x1_, x2_, x3_, x4_};
    optimize_points2d(model->cams, all_x, model->X);
}

int AbsolutePose2dEstimator::NonMinimalSolver(const std::vector<int>& sample, Pose2d* model) const
{

  Eigen::Matrix<double, Eigen::Dynamic, 4> A(sample.size(), 4);

  const int sample_size = static_cast<int>(sample.size());
  for(int i = 0; i < sample_size; ++i) {
    double x1 = x_[sample[i]](0);
    double x2 = x_[sample[i]](1);
    double X1 = X_[sample[i]](0);
    double X2 = X_[sample[i]](1);

    A.row(i) << X1*x2 - X2*x1, - X1*x1 - X2*x2, x2, -x1;
  }

  Eigen::JacobiSVD<decltype(A)> svd(A, Eigen::ComputeFullV);

  Eigen::Matrix<double, 4, 1> t = svd.matrixV().col(3);
  t /= t.block<2,1>(0,0).norm();

  (*model)(0,0) = t(0);
  (*model)(1,1) = t(0);
  (*model)(0,1) = -t(1);
  (*model)(1,0) = t(1);
  (*model)(0,2) = t(2);
  (*model)(1,2) = t(3);

  // Choose sign by checking cheirality of first point
  if( model->row(1) * X_[sample[0]].homogeneous() < 0 )
    (*model) *= -1.0;

  return 1;
}

double AbsolutePose2dEstimator::EvaluateModelOnPoint(const Pose2d& model, int i) const
{
  // Cosine error
  double err = 1.0 - x_[i].dot((model * X_[i].homogeneous()).normalized());
  return err;
}

} // namespace init
} // namespace colmap
