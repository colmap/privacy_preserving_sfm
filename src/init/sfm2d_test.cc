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

#define TEST_NAME "init/sfm2d_test"
#include "util/testing.h"
#include <RansacLib/ransac.h>
#include "init/sfm2d.h"
#include <iostream>
#include <random>
using namespace colmap;
using namespace colmap::init;

void set_random_pose2d(Pose2d *cam) {
  cam->setRandom();
  (*cam)(1,1) = (*cam)(0,0);
  (*cam)(0,1) = -(*cam)(1,0);
  *cam /= cam->col(0).norm();
}

void setup_plausible_scene(int Ncams, int Npoints, std::vector<Pose2d> &cams, std::vector<std::vector<Eigen::Vector2d>> &x, std::vector<Eigen::Vector2d> &X ) {
  bool done = false;
  x.resize(Ncams);

  while(!done) {
    // Setup cameras
    cams.clear();
    for(int i = 0; i < Ncams; ++i) {
      Pose2d cam;
      if(i == 0) {
        cam.setZero();
        cam.block<2,2>(0,0).setIdentity();
      } else {
        set_random_pose2d(&cam);
      }
      cams.push_back(cam);
      x[i].clear();
    }

    done = true;
    X.clear();
    for(int i = 0; i < Npoints; ++i) {
      Eigen::Vector2d Z;
      Z.setRandom();
      if(Z(1) < 0)
        Z(1) = -Z(1);
      X.push_back(Z);
      
      for(int j = 0; j < Ncams; ++j) {
        Eigen::Vector2d z = cams[j]*Z.homogeneous();
        if(z(1) < 0) {
          done = false;
          break;
        }
        x[j].push_back(z.normalized());          
      }

      if(!done)
        break;
    }
  }
}


void add_outliers(int Noutliers, std::vector<std::vector<Eigen::Vector2d>> &x) {
    int Ncams = x.size();
    int Npts = x[0].size();

    std::vector<int> ind;
    ind.resize(Npts);
    std::iota(ind.begin(), ind.end(), 0);
    std::shuffle(ind.begin(), ind.end(), std::mt19937{std::random_device{}()});

    for(int i = 0; i < Noutliers; ++i) {
      int cam_idx = rand() % Ncams;
      // std::cout << "Making point " << ind[i] << " outlier in camera " << cam_idx << " .\n";
      x[cam_idx][ind[i]].setRandom();
    }


}

BOOST_AUTO_TEST_CASE(AbsPoseSolver) {
  Pose2d P_gt, P;
  set_random_pose2d(&P_gt);
  
  std::vector<Eigen::Vector2d> x;
  std::vector<Eigen::Vector2d> X;
  
  for(int i = 0; i < 3; ++i) {
      Eigen::Vector3d Xh;
      Xh.setRandom();
      Xh(2) = 1.0;

      // Choose sign of P_gt such that first point is infront
      if(i == 0 && (P_gt.row(1)*Xh < 0))
        P_gt *= -1.0;

      x.push_back((P_gt*Xh).normalized());
      X.push_back(Xh.topRows(2));
  }

  std::vector<int> sample = {0,1,2};

  AbsolutePose2dEstimator estimator(x,X);

  estimator.NonMinimalSolver(sample, &P);

  BOOST_CHECK_SMALL( (P - P_gt).norm(), 1e-8 );
}

BOOST_AUTO_TEST_CASE(SceneGeneratorTest) {

  std::vector<std::vector<Eigen::Vector2d>> x;
  std::vector<Eigen::Vector2d> X;
  std::vector<Pose2d> cams;

  setup_plausible_scene(4, 10, cams, x, X);
  BOOST_CHECK(X.size() == 10);
  BOOST_CHECK(x.size() == 4);
  for(int i = 0; i < 4; ++i) {
    BOOST_CHECK(x[i].size() == 10);
  
    for(int j = 0; j < 10; j++) {
      Eigen::Vector2d z = cams[i] * X[j].homogeneous();

      BOOST_CHECK( z(1) > 0 );

      // Reprojection error
      BOOST_CHECK_SMALL( (z.normalized() - x[i][j]).norm(), 1e-8 );
    }
  }
}

BOOST_AUTO_TEST_CASE(RansacTestAbsolutePoseNoOutliers) {
  std::vector<std::vector<Eigen::Vector2d>> x;
  std::vector<Eigen::Vector2d> X;
  std::vector<Pose2d> cams;

  const int Npts = 10;
  const int Ncams = 4;

  // First camera is always identity...
  setup_plausible_scene(Ncams, Npts, cams, x, X);
  ransac_lib::LORansacOptions options;

  for(int i = 0; i < Ncams; ++i) {
    AbsolutePose2dEstimator solver(x[i], X);
    
    ransac_lib::LocallyOptimizedMSAC<Pose2d, Pose2dVector, AbsolutePose2dEstimator> lomsac;
    ransac_lib::RansacStatistics ransac_stats;

    Pose2d best_model;

    int inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);


    BOOST_CHECK( inliers == Npts );

    BOOST_CHECK_SMALL( (best_model - cams[i]).norm(), 1e-8 );
  }

}

BOOST_AUTO_TEST_CASE(RansacTestAbsolutePose) {
  std::vector<std::vector<Eigen::Vector2d>> x;
  std::vector<Eigen::Vector2d> X;
  std::vector<Pose2d> cams;

  const int Npts = 100;
  const int Ncams = 4;

  const int Noutliers = 20;

  // First camera is always identity...
  setup_plausible_scene(Ncams, Npts, cams, x, X);

  add_outliers(Noutliers, x);

  ransac_lib::LORansacOptions options;
  options.squared_inlier_threshold_ = 2e-5;
  

  for(int i = 0; i < Ncams; ++i) {
    AbsolutePose2dEstimator solver(x[i], X);
    
    ransac_lib::LocallyOptimizedMSAC<Pose2d, Pose2dVector, AbsolutePose2dEstimator> lomsac;
    ransac_lib::RansacStatistics ransac_stats;


    Pose2d best_model;

    int inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

    /*
    std::cout << "Found pose with " << inliers << " inliers. (" << ransac_stats.num_iterations << ")\n";
    std::cout << "P = \n" << best_model << "\n";
    std::cout << "Pgt = \n" << cams[i] << "\n";
    */
  

    BOOST_CHECK( inliers >= Npts-Noutliers );

    BOOST_CHECK_SMALL( (best_model - cams[i]).norm(), 1e-8 );
  }

}

BOOST_AUTO_TEST_CASE(RansacTestFourViewEstimator) {
  std::vector<std::vector<Eigen::Vector2d>> x;
  std::vector<Eigen::Vector2d> X;
  std::vector<Pose2d> cams;

  const int Npts = 100;
  const int Ncams = 4;
  const int Noutliers = 20;


  // First camera is always identity...
  setup_plausible_scene(Ncams, Npts, cams, x, X);

  add_outliers(Noutliers, x);

  ransac_lib::LORansacOptions options;
  options.squared_inlier_threshold_ = 2e-5; // approx 0.5 deg
  options.squared_inlier_threshold_ = 1e-7; // approx 0.5 deg
  
  FourView2dEstimator solver(x[0], x[1], x[2], x[3], options.squared_inlier_threshold_);
  

  ransac_lib::LocallyOptimizedMSAC<FourView2dEstimator::Reconstruction, FourView2dEstimator::ReconstructionVector, FourView2dEstimator> lomsac;
  ransac_lib::RansacStatistics ransac_stats;

  FourView2dEstimator::Reconstruction best_model;
  int inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);


  std::cout << "GT = \n" << cams[0] << "\n" << cams[1] << "\n" << cams[2] << "\n" << cams[3] << "\n";
  std::cout << "PP = \n" << best_model.cams[0] << "\n" << best_model.cams[1] << "\n" << best_model.cams[2] << "\n" << best_model.cams[3] << "\n";
  std::cout << "inliers = " << inliers << "\n";

  BOOST_CHECK( inliers >= Npts-Noutliers );
}