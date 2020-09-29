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

#define TEST_NAME "init/initializer_test"
#include "util/testing.h"
#include <RansacLib/ransac.h>
#include "init/initializer.h"
#include <iostream>
#include <random>
using namespace colmap;
using namespace colmap::init;

void set_random_camera(Pose *cam, bool upright = false) {

    cam->block<3,3>(0,0) = Eigen::Quaterniond::UnitRandom().toRotationMatrix();
    if(upright) {
        cam->block<3,3>(0,0) = Eigen::Quaterniond::FromTwoVectors(cam->col(1), Eigen::Vector3d{0.0, 1.0, 0.0}) * cam->block<3,3>(0,0);
    }
    cam->block<3,1>(0,3).setRandom();
}

void setup_plausible_scene(int Ncams, int Npoints, std::vector<Pose> &cams, std::vector<std::vector<Eigen::Vector2d>> &x, std::vector<Eigen::Vector3d> &X, bool upright = false ) {
  bool done = false;
  x.resize(Ncams);

  while(!done) {
    // Setup cameras
    cams.clear();
    for(int i = 0; i < Ncams; ++i) {
      Pose cam;
      if(i == 0) {
        cam.setZero();
        cam.block<3,3>(0,0).setIdentity();
      } else {
        set_random_camera(&cam, upright);
      }

      if(i == 1) {
          cam.col(3).normalize();
      }

      cams.push_back(cam);
      x[i].clear();
    }

    done = true;
    X.clear();
    for(int i = 0; i < Npoints; ++i) {
      Eigen::Vector3d Z;
      Z.setRandom();
      if(Z(2) < 0)
        Z(2) = -Z(2);
      X.push_back(Z);
      
      for(int j = 0; j < Ncams; ++j) {
        Eigen::Vector3d z = cams[j]*Z.homogeneous();
        if(z(2) < 0) {
          done = false;
          break;
        }
        x[j].push_back((z.hnormalized()));          
      }

      if(!done)
        break;
    }
  }
}



void setup_random_lines(int Ngravity_aligned, const std::vector<std::vector<Eigen::Vector2d>> &x, const std::vector<Pose> &cams, std::vector<FeatureLines> &lines, std::vector<Eigen::Vector3d> &gravity, double gravity_noise = 0.0) {
    int Ncams = cams.size();
    int Npts = x[0].size();

    std::vector<int> ind;
    ind.resize(Npts);
    std::iota(ind.begin(), ind.end(), 0);
    std::shuffle(ind.begin(), ind.end(), std::mt19937{std::random_device{}()});

    lines.resize(Ncams);
    for(int i = 0; i < Ncams; ++i) {
        gravity.push_back(cams[i].col(1));
        if(gravity_noise > 0) {
            Eigen::AngleAxisd aa;
            aa.axis().setRandom().normalize();
            aa.angle() = gravity_noise;
            gravity[i] = aa.toRotationMatrix() * gravity[i];

        }

        for(int j = 0; j < Npts; ++j) {
            FeatureLine line;
            line.SetAligned(ind[j] < Ngravity_aligned);

            if(line.IsAligned()) {
                line.SetLine(x[i][j].homogeneous().cross(gravity[i]).normalized());
            } else {
                Eigen::Vector3d n;
                n.setRandom(); // TODO: This is not uniformly distributed! (Biased towards corners)
                line.SetLine(x[i][j].homogeneous().cross(n).normalized());
            }
            
            lines[i].push_back(line);
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



BOOST_AUTO_TEST_CASE(SceneGeneratorTest) {

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    

    setup_plausible_scene(4, 10, cams, x, X);
    BOOST_REQUIRE(X.size() == 10);

    BOOST_REQUIRE(x.size() == 4);

    for(int i = 0; i < 4; ++i) {
        BOOST_REQUIRE(x[i].size() == 10);

        for(int j = 0; j < 10; j++) {
            Eigen::Vector3d z = cams[i] * X[j].homogeneous();

            BOOST_CHECK( z(2) > 0 );

            // Reprojection error
            BOOST_CHECK_SMALL( (z.hnormalized() - x[i][j]).norm(), 1e-6 );
        }
    }
}

BOOST_AUTO_TEST_CASE(PrivacyGeneratorTest) {

    const int Ncams = 4;
    const int Npts = 20;
    const int Ngrav = 10;

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X);

    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;

    setup_random_lines(Ngrav, x, cams, lines, gravity);
    BOOST_REQUIRE(lines.size() == Ncams);

    BOOST_REQUIRE(gravity.size() == Ncams);


    for(int i = 0; i < Ncams; ++i) {
        BOOST_REQUIRE(lines[i].size() == Npts);
        
        int n_grav = 0;

        for(int j = 0; j < Npts; j++) {


            Eigen::Vector3d z = cams[i] * X[j].homogeneous();

            BOOST_CHECK( z(2) > 0 );

            if(lines[i][j].IsAligned()) {
                BOOST_CHECK_SMALL( lines[i][j].Line().cast<double>().dot(gravity[i]) , 1e-6 );
                n_grav++;
            }

            // Reprojection error
            BOOST_CHECK_SMALL( lines[i][j].Line().cast<double>().dot(z) , 1e-6 );
        }

        BOOST_REQUIRE(n_grav == Ngrav);
    }
}




BOOST_AUTO_TEST_CASE(PlanarOffsetEstimatorNoOutliers) {

    const int Ncams = 4;
    const int Npts = 20;
    const int Ngrav = 0;

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X, true);

    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;

    setup_random_lines(Ngrav, x, cams, lines, gravity);

    std::vector<std::vector<Eigen::Vector3d>> lines_r(Ncams);
    Eigen::Vector3d t_gt;

    for(int i = 0; i < Ncams; ++i) {
        for(int j = 0; j < Npts; ++j) {
            if(!lines[i][j].IsAligned())
                lines_r[i].push_back(lines[i][j].Line().cast<double>());
        }
        if(i > 0) {
            t_gt(i-1) = cams[i](1,3);
            cams[i](1,3) = 0.0;
        }
    }
    std::vector<Eigen::Matrix3d> Rg;
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());


    ransac_lib::LORansacOptions options;
    PlanarOffsetEstimator solver(cams, lines_r, Rg, 0.005 * 0.005);

    ransac_lib::LocallyOptimizedMSAC<PlanarOffsetEstimator::Reconstruction, PlanarOffsetEstimator::ReconstructionVector, PlanarOffsetEstimator> lomsac;
    ransac_lib::RansacStatistics ransac_stats;

    PlanarOffsetEstimator::Reconstruction best_model;

    int inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

    std::cout << "PlanarOffsetEstimator: inliers " << inliers << "\n";

    BOOST_REQUIRE(inliers == Npts);

    
}


BOOST_AUTO_TEST_CASE(PlanarOffsetEstimatorWithOutliers) {

    const int Ncams = 4;
    const int Npts = 100;
    const int Ngrav = 0;
    const int Noutliers = 20;
    
    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X, true);
    add_outliers(Noutliers, x);

    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;

    setup_random_lines(Ngrav, x, cams, lines, gravity);

    std::vector<std::vector<Eigen::Vector3d>> lines_r(Ncams);
    Eigen::Vector3d t_gt;

    for(int i = 0; i < Ncams; ++i) {
        for(int j = 0; j < Npts; ++j) {
            if(!lines[i][j].IsAligned())
                lines_r[i].push_back(lines[i][j].Line().cast<double>());
        }
        if(i > 0) {
            t_gt(i-1) = cams[i](1,3);
            cams[i](1,3) = 0.0;
        }
    }

    std::vector<Eigen::Matrix3d> Rg;
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());
    Rg.push_back(Eigen::Matrix3d::Identity());

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 1e-6;
    PlanarOffsetEstimator solver(cams, lines_r, Rg, options.squared_inlier_threshold_);

    ransac_lib::LocallyOptimizedMSAC<PlanarOffsetEstimator::Reconstruction, PlanarOffsetEstimator::ReconstructionVector, PlanarOffsetEstimator> lomsac;
    ransac_lib::RansacStatistics ransac_stats;

    PlanarOffsetEstimator::Reconstruction best_model;

    int inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

    BOOST_REQUIRE(inliers >- Npts - Noutliers);
    //BOOST_REQUIRE_SMALL( (best_model - t_gt).norm(), 1e-6 );
}




BOOST_AUTO_TEST_CASE(InitializerNoOutliers) {

    const int Ncams = 4;
    const int Npts = 100;
    const int Ngrav = 50;

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X, true);

    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;

    setup_random_lines(Ngrav, x, cams, lines, gravity);

    std::vector<Pose> poses;
    double inlier_ratio;
    InitOptions init_options;
    initialize_reconstruction(lines, gravity, init_options, &poses, &inlier_ratio);

    

    BOOST_REQUIRE(poses.size() == Ncams);

    // normalize w.r.t. second pose
    double s = poses[1].col(3).norm();

    for(int i = 0; i < Ncams; ++i) {

        poses[i].col(3) /= s;

        BOOST_CHECK_SMALL( (poses[i] - cams[i]).norm(), 1e-6 );

    }


}



BOOST_AUTO_TEST_CASE(InitializerWithOutliers) {

    const int Ncams = 4;
    const int Npts = 100;
    const int Ngrav = 50;
    const int Noutliers = 10;

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X, true);

    add_outliers(Noutliers, x);


    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;

    setup_random_lines(Ngrav, x, cams, lines, gravity);

    std::vector<Pose> poses;
    double inlier_ratio;
    InitOptions init_options;
    initialize_reconstruction(lines, gravity, init_options, &poses, &inlier_ratio);

    

    BOOST_REQUIRE(poses.size() == Ncams);

    // normalize w.r.t. second pose
    double s = poses[1].col(3).norm();

    for(int i = 0; i < Ncams; ++i) {
        

        poses[i].col(3) /= s;

        std::cout << "Cam " << i << "\n";
        std::cout << "EST: \n" << cams[i] << "\n";
        std::cout << "GT : \n" << poses[i] << "\n\n";
        BOOST_CHECK_SMALL( (poses[i] - cams[i]).norm(), 1e-4 );

    }


}


BOOST_AUTO_TEST_CASE(InitializerWithOutliersAndNoisyGravity) {

    const int Ncams = 4;
    const int Npts = 100;
    const int Ngrav = 50;
    const int Noutliers = 10;

    std::vector<std::vector<Eigen::Vector2d>> x;
    std::vector<Eigen::Vector3d> X;
    std::vector<Pose> cams;
    
    setup_plausible_scene(Ncams, Npts, cams, x, X, true);

    add_outliers(Noutliers, x);


    std::vector<FeatureLines> lines;
    std::vector<Eigen::Vector3d> gravity;
    double gravity_noise = 1 * M_PI / 180.0;

    setup_random_lines(Ngrav, x, cams, lines, gravity, gravity_noise);

    std::vector<Pose> poses;
    double inlier_ratio;
    InitOptions init_options;
    initialize_reconstruction(lines, gravity, init_options, &poses, &inlier_ratio);

    

    BOOST_REQUIRE(poses.size() == Ncams);

    // normalize w.r.t. second pose
    double s = poses[1].col(3).norm();
    Eigen::Matrix3d R0 = poses[0].block<3,3>(0,0);

    for(int i = 0; i < Ncams; ++i) {
        poses[i].block<3,3>(0,0) = poses[i].block<3,3>(0,0) * R0.transpose();

        poses[i].col(3) /= s;

        std::cout << "Cam " << i << "\n";
        std::cout << "GT: \n" << cams[i] << "\n";
        std::cout << "EST : \n" << poses[i] << "\n\n";
        BOOST_CHECK_SMALL( (poses[i] - cams[i]).norm(), 0.05 );

    }


}

