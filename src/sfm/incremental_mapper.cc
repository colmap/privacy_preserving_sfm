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

#include "sfm/incremental_mapper.h"

#include <array>
#include <fstream>
#include <util/random.h>

#include "base/pose.h"
#include "base/projection.h"
#include "base/triangulation.h"
#include "estimators/pose.h"
#include "init/initializer.h"
#include "util/bitmap.h"

namespace colmap {
namespace {

void SortAndAppendNextImages(std::vector<std::pair<image_t, float>> image_ranks,
                             std::vector<image_t>* sorted_images_ids) {
  std::sort(image_ranks.begin(), image_ranks.end(),
            [](const std::pair<image_t, float>& image1,
               const std::pair<image_t, float>& image2) {
              return image1.second > image2.second;
            });

  sorted_images_ids->reserve(sorted_images_ids->size() + image_ranks.size());
  for (const auto& image : image_ranks) {
    sorted_images_ids->push_back(image.first);
  }

  image_ranks.clear();
}

float RankNextImageMaxVisiblePointsNum(const Image& image) {
  return static_cast<float>(image.NumVisiblePoints3D());
}

float RankNextImageMaxVisiblePointsRatio(const Image& image) {
  return static_cast<float>(image.NumVisiblePoints3D()) /
         static_cast<float>(image.NumObservations());
}

}  // namespace

bool IncrementalMapper::Options::Check() const {
  CHECK_OPTION_GT(init_min_num_inliers, 0);
  CHECK_OPTION_GT(init_max_error, 0.0);
  CHECK_OPTION_GE(init_min_tri_angle, 0.0);
  CHECK_OPTION_GT(abs_pose_max_error, 0.0);
  CHECK_OPTION_GT(abs_pose_min_num_inliers, 0);
  CHECK_OPTION_GE(abs_pose_min_inlier_ratio, 0.0);
  CHECK_OPTION_LE(abs_pose_min_inlier_ratio, 1.0);
  CHECK_OPTION_GE(local_ba_num_images, 2);
  CHECK_OPTION_GE(local_ba_min_tri_angle, 0.0);
  CHECK_OPTION_GE(min_focal_length_ratio, 0.0);
  CHECK_OPTION_GE(max_focal_length_ratio, min_focal_length_ratio);
  CHECK_OPTION_GE(max_extra_param, 0.0);
  CHECK_OPTION_GE(filter_max_reproj_error, 0.0);
  CHECK_OPTION_GE(filter_min_tri_angle, 0.0);
  CHECK_OPTION_GE(max_reg_trials, 1);
  return true;
}

IncrementalMapper::IncrementalMapper(const DatabaseCache* database_cache)
    : database_cache_(database_cache),
      reconstruction_(nullptr),
      triangulator_(nullptr),
      num_total_reg_images_(0),
      num_shared_reg_images_(0) {}

void IncrementalMapper::BeginReconstruction(Reconstruction* reconstruction) {
  CHECK(reconstruction_ == nullptr);
  reconstruction_ = reconstruction;
  reconstruction_->Load(*database_cache_);
  reconstruction_->SetUp(&database_cache_->CorrespondenceGraph());
  triangulator_.reset(new IncrementalTriangulator(
      &database_cache_->CorrespondenceGraph(), reconstruction));

  num_shared_reg_images_ = 0;
  num_reg_images_per_camera_.clear();
  for (const image_t image_id : reconstruction_->RegImageIds()) {
    RegisterImageEvent(image_id);
  }

  existing_image_ids_ =
      std::unordered_set<image_t>(reconstruction->RegImageIds().begin(),
                                  reconstruction->RegImageIds().end());

  filtered_images_.clear();
  num_reg_trials_.clear();
}

void IncrementalMapper::EndReconstruction(const bool discard) {
  CHECK_NOTNULL(reconstruction_);

  if (discard) {
    for (const image_t image_id : reconstruction_->RegImageIds()) {
      DeRegisterImageEvent(image_id);
    }
  }

  reconstruction_->TearDown();
  reconstruction_ = nullptr;
  triangulator_.reset();
}

std::vector<image_t> IncrementalMapper::FindNextImages(const Options& options) {
  CHECK_NOTNULL(reconstruction_);
  CHECK(options.Check());

  std::function<float(const Image&)> rank_image_func;
  switch (options.image_selection_method) {
    case Options::ImageSelectionMethod::MAX_VISIBLE_POINTS_NUM:
      rank_image_func = RankNextImageMaxVisiblePointsNum;
      break;
    case Options::ImageSelectionMethod::MAX_VISIBLE_POINTS_RATIO:
      rank_image_func = RankNextImageMaxVisiblePointsRatio;
      break;
  }

  std::vector<std::pair<image_t, float>> image_ranks;
  std::vector<std::pair<image_t, float>> other_image_ranks;

  // Append images that have not failed to register before.
  for (const auto& image : reconstruction_->Images()) {
    // Skip images that are already registered.
    if (image.second.IsRegistered()) {
      continue;
    }

    // Only consider images with a sufficient number of visible points.
    if (image.second.NumVisiblePoints3D() <
        static_cast<size_t>(options.abs_pose_min_num_inliers)) {
      continue;
    }

    // Only try registration for a certain maximum number of times.
    const size_t num_reg_trials = num_reg_trials_[image.first];
    if (num_reg_trials >= static_cast<size_t>(options.max_reg_trials)) {
      continue;
    }

    // If image has been filtered or failed to register, place it in the
    // second bucket and prefer images that have not been tried before.
    const float rank = rank_image_func(image.second);
    if (filtered_images_.count(image.first) == 0 && num_reg_trials == 0) {
      image_ranks.emplace_back(image.first, rank);
    } else {
      other_image_ranks.emplace_back(image.first, rank);
    }
  }

  std::vector<image_t> ranked_images_ids;
  SortAndAppendNextImages(image_ranks, &ranked_images_ids);
  SortAndAppendNextImages(other_image_ranks, &ranked_images_ids);

  return ranked_images_ids;
}

bool IncrementalMapper::RegisterInitialLineImages(const Options& options, const DatabaseCache& aligned_db_cache) {

  CHECK_NOTNULL(reconstruction_);
  CHECK_EQ(reconstruction_->NumRegImages(), 0);

  const CorrespondenceGraph& corr_graph = aligned_db_cache.CorrespondenceGraph();

  // Store all 4-image tracks that only consist of aligned or unaligned
  // features, respectively.
  // We store them sorted to make it later simpler to remove duplicates

  using ImageSet = std::array<image_t,4>;
  using FeatureIdxs = std::array<int,4>;

  // Comparison function for image sets.
  // This is required to store image sets in an std::map
  const auto image_set_before =
      [](const ImageSet& lhs, const ImageSet& rhs) {
        for (int i = 0; i < 4; ++i) {
          const image_t lid = lhs.at(i);
          const image_t rid = rhs.at(i);
          if (lid < rid) {
            return true;
          } else if (lid > rid) {
            return false;
          }
        }
        return false;
      };

  // Comparison function for feature indices.
  // This is required to store the feature indices for the selected images in a
  // std::set.
  const auto feature_idxs_before =
      [](const FeatureIdxs& lhs, const FeatureIdxs& rhs) {
        for (int i = 0; i < 4; ++i) {
          const int lidx = lhs.at(i);
          const int ridx = rhs.at(i);
          if (lidx < ridx) {
            return true;
          } else if (lidx > ridx) {
            return false;
          }
        }
        return false;
      };

  // Containers to collect all aligned or unaligned feature tracks,
  // respectively.
  // We store them for each image set separately, this way we can let sted::set
  // take care of removing duplicates.
  std::map<ImageSet,
      std::set<FeatureIdxs, decltype(feature_idxs_before)>,
      decltype(image_set_before)>
      all_aligned_tracks(image_set_before);
  std::map<ImageSet,
      std::set<FeatureIdxs, decltype(feature_idxs_before)>,
      decltype(image_set_before)>
      all_unaligned_tracks(image_set_before);

  // Comparison for image IDs inside correspondences.
  // We use this to get a consistent ordering of image sets when assembling
  // tracks from different reference images to avoid duplicates.
  const auto lower_image_id =
      [](const CorrespondenceGraph::Correspondence& lhs,
         const CorrespondenceGraph::Correspondence& rhs) {
        return lhs.image_id < rhs.image_id;
      };

  // For all correspondences to a reference feature, we generate all possible
  // 4-view tracks that include the reference.
  // We could also add the ones that are just connected by the reference but
  // do not include it, but the runtime will be horrible.
  const auto assemble_tracks =
      [lower_image_id](
          const CorrespondenceGraph::Correspondence& reference_corr,
          const std::vector<CorrespondenceGraph::Correspondence>& corrs) {

        std::vector<std::vector<CorrespondenceGraph::Correspondence>> tracks;
        const int num_corrs = corrs.size();
        for (int i = 0; i < num_corrs; ++i) {
          for (int j = i+1; j < num_corrs; ++j) {
            for (int k = j+1; k < num_corrs; ++k) {
              std::set<CorrespondenceGraph::Correspondence,
                  decltype(lower_image_id)>
                  track_candidate(lower_image_id);
              track_candidate.emplace(reference_corr);
              track_candidate.emplace(corrs.at(i));
              track_candidate.emplace(corrs.at(j));
              track_candidate.emplace(corrs.at(k));

              if (track_candidate.size() == 4) {
                tracks.emplace_back(track_candidate.begin(), track_candidate.end());
              }
            }
          }
        }
        return tracks;
      };

  // Again, only loop over images that actually have aligned features
  // We sort the IDs first, this way we might be able to do something smart
  // once we know that all images in a set were checked already.
  std::vector<image_t> image_ids;
  for (const auto& image : aligned_db_cache.Images()) {
    image_ids.emplace_back(image.first);
  }
  std::sort(image_ids.begin(), image_ids.end());

  // We randomly select some images to run the check
  const int num_images = static_cast<int>(image_ids.size());
  const int num_check_images = std::min(10, num_images);
  std::unordered_set<image_t> check_image_ids;
  SetPRNGSeed(std::chrono::system_clock::now().time_since_epoch().count());
  while(check_image_ids.size() < num_check_images) {
    const image_t check_image_id = image_ids.at(RandomInteger(0,num_images-1));
    check_image_ids.emplace(check_image_id);
  }

  for (const auto& image_id : check_image_ids) {

    const Image& image = aligned_db_cache.Image(image_id);
    std::cout << "Collecting tracks for image " << image_id << std::endl;

    // Search for correspondences for all features, aligned and unaligned.
    // Wen can split the tracks later
    const int num_features = image.Lines().size();
    for (int line_idx = 0; line_idx < num_features; ++line_idx) {

      const bool aligned_feature = image.Line(line_idx).IsAligned();
      const auto& corrs = corr_graph.FindCorrespondences(image_id, line_idx);

      // Only consider correspondences with the same alignment (either
      // gravity-aligned or random) as the reference.
      std::vector<CorrespondenceGraph::Correspondence> alignment_corrs;
      for (const auto& corr : corrs) {
        if (aligned_db_cache.Image(corr.image_id).Line(corr.line_idx).IsAligned() ==
            aligned_feature) {
          alignment_corrs.emplace_back(corr);
        }
      }

      const CorrespondenceGraph::Correspondence
          reference_corr(image_id, line_idx);

      if (alignment_corrs.size() >= 3) {

        const auto new_tracks = assemble_tracks(reference_corr, alignment_corrs);

        // Make sure we add the tracks to the correct container
        auto& track_container =
            aligned_feature ? all_aligned_tracks : all_unaligned_tracks;

        for (const auto& track : new_tracks) {
          ImageSet track_set;
          FeatureIdxs track_features;
          for (int i = 0; i < 4; ++i) {
            track_set[i] = track.at(i).image_id;
            track_features[i] = track.at(i).line_idx;
          }

          // We need to generate elements for new images sets explicitly since
          // they require the comparison function as argument.
          if (track_container.count(track_set) == 0) {
            track_container.emplace(
                std::make_pair(track_set, feature_idxs_before));
          }
          track_container.at(track_set).emplace(track_features);
        }
      }
    }
  }

  std::map<int,int> num_aligned_tracks;
  for (auto& image_set_tracks : all_aligned_tracks) {
    num_aligned_tracks[image_set_tracks.second.size()] += 1;
  }

  std::cout << "Number of image sets with N aligned tracks:\n";
  for (const auto& num : num_aligned_tracks) {
    std::cout << "N = " << num.first << ": " << num.second << std::endl;
  }

  std::map<int,int> num_unaligned_tracks;
  for (auto& image_set_tracks : all_unaligned_tracks) {
    num_unaligned_tracks[image_set_tracks.second.size()] += 1;
  }

  std::cout << "Number of image sets with N unaligned tracks:\n";
  for (const auto& num : num_unaligned_tracks) {
    std::cout << "N = " << num.first << ": " << num.second << std::endl;
  }

  // Collect the numbers of both track types to find the best image sets for initialization
  struct TrackNumbers {
    ImageSet image_set;
    size_t aligned_tracks;
    size_t unaligned_tracks;
  };
  std::vector<TrackNumbers> track_numbers;

  // We can enforce some minimum numbers of tracks here.
  // Make sure that we actually have both aligned and random tracks for the
  // image set.
  const size_t min_num_aligned_tracks = 20;
  const size_t min_num_random_tracks = 20;
  for (const auto& aligned_set_tracks : all_aligned_tracks) {
    if (aligned_set_tracks.second.size() >= min_num_aligned_tracks
        && all_unaligned_tracks.count(aligned_set_tracks.first) > 0
        && all_unaligned_tracks.at(aligned_set_tracks.first).size()
           >= min_num_random_tracks) {
      track_numbers.emplace_back(
          TrackNumbers{
              aligned_set_tracks.first,
              aligned_set_tracks.second.size(),
              all_unaligned_tracks.at(aligned_set_tracks.first).size()});
    }
  }

  if (track_numbers.empty()) {
    std::cerr << "Error - could not find sufficient tracks\n";
    return false;
  }

  // Sort the image sets according to their track numbers.
  // We put a higher weight on the aligned tracks since we usually have less.
  std::sort(track_numbers.begin(), track_numbers.end(),
            [](const TrackNumbers& lhs, const TrackNumbers& rhs) {
              const double unaligned_weight = 0.0;
              return lhs.aligned_tracks + (unaligned_weight * lhs.unaligned_tracks)
                     > rhs.aligned_tracks + (unaligned_weight * rhs.unaligned_tracks);
            });

  double best_inlier_ratio = 0;
  int best_inliers = 0;
  std::vector<init::Pose> best_poses;
  std::unordered_map<image_t, int> best_image_idx_map;

  // Loop over the possible initialization sets.
  // Just taking the set with the most tracks often leads to a very small
  // baseline which makes it hard to get good constraints.
  // We limit the number of tries, otherwise this runs forever.
  const size_t max_num_init_tries = 10;
  const int num_test_sets = std::min(track_numbers.size(), max_num_init_tries);
  for (int init_set_idx = 0; init_set_idx < num_test_sets; ++init_set_idx) {
    // Assemble the lines for the image set
    const auto init_image_set = track_numbers.at(init_set_idx).image_set;

    std::cout << StringPrintf("Choose initialization set: %u %u %u %u\n",
                              init_image_set.at(0), init_image_set.at(1),
                              init_image_set.at(2), init_image_set.at(3));

    std::vector<FeatureLines> lines(4);
    std::unordered_map<image_t, int> image_idx_map;
    std::vector<Eigen::Vector3d> images_gravity(4);

    int next_idx = 0;
    for (const auto& im_id : init_image_set) {
      image_idx_map[im_id] = next_idx;
      next_idx += 1;
    }

    for (const auto& im_id_idx : image_idx_map) {
      images_gravity.at(im_id_idx.second) =
          aligned_db_cache.Image(im_id_idx.first).GravityDirection();
    }

    // First, collect all the aligned lines
    const auto& aligned_tracks = all_aligned_tracks.at(init_image_set);
    for (const auto& track : aligned_tracks) {
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const int line_idx = track.at(i);
        lines.at(image_idx_map.at(image_id))
            .emplace_back(aligned_db_cache.Image(image_id).Line(line_idx));
        CHECK(aligned_db_cache.Image(image_id).Line(line_idx).IsAligned());
      }
    }

    // Collect all unaligned lines
    const auto& unaligned_tracks = all_unaligned_tracks.at(init_image_set);
    for (const auto& track : unaligned_tracks) {
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const int line_idx = track.at(i);
        lines.at(image_idx_map.at(image_id))
            .emplace_back(aligned_db_cache.Image(image_id).Line(line_idx));
        CHECK(!aligned_db_cache.Image(image_id).Line(line_idx).IsAligned());
      }
    }

    std::vector<init::Pose> poses;
    double inlier_ratio = 0;
    init::InitOptions init_options;
    init_options.min_num_inliers = options.init_min_num_inliers;
    init_options.min_tri_angle = options.init_min_tri_angle;
    init_options.max_error = [&]() {
      double min_error = std::numeric_limits<double>::max();
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const camera_t camera_id = aligned_db_cache.Image(image_id).CameraId();
        const Camera& camera = aligned_db_cache.Camera(camera_id);
        min_error = std::min(min_error, camera.ImageToWorldThreshold(options.init_max_error));
      }
      return min_error;
    }();
    const bool success =
        init::initialize_reconstruction(
            lines, images_gravity, init_options, &poses, &inlier_ratio);

    // We choose the best initialization according to the inlier ratio
    if (success && poses.size() > 0) {
      std::cout << "Got poses with an inlier ratio of " << inlier_ratio << "\n";

      if (inlier_ratio > best_inlier_ratio) {
        std::vector<image_t> prev_best_images;
        for (const auto& im_id_idx : best_image_idx_map) {
          prev_best_images.push_back(im_id_idx.first);
        }
        if (prev_best_images.size() != 4) {
          prev_best_images = {0, 0, 0, 0};
        }
        std::cout << "Update best model "
                  << StringPrintf("(%u %u %u %u)",
                                  prev_best_images.at(0), prev_best_images.at(1),
                                  prev_best_images.at(2), prev_best_images.at(3))
                  << " -> "
                  << StringPrintf("(%u %u %u %u)",
                                  init_image_set.at(0), init_image_set.at(1),
                                  init_image_set.at(2), init_image_set.at(3))
                  << " - (" << best_inlier_ratio << " vs. " << inlier_ratio << ")\n";
        best_inlier_ratio = inlier_ratio;
        best_poses = poses;
        best_image_idx_map = image_idx_map;
        best_inliers = inlier_ratio * lines.at(0).size();
      }
    } else {
      std::cout << "Could not estimate poses\n";
    }
  }

  if (best_poses.size() == 0) {
    std::cerr << "Could not estimate initial image poses\n";
    return false;
  }

  if (best_inliers < options.init_min_num_inliers) {
    std::cerr << "Not enough inliers for initialization (" << best_inliers << " instead of " << options.init_min_num_inliers << ")\n";
    return false;
  }

  for (const auto& im_id_idx : best_image_idx_map) {
    const image_t image_id = im_id_idx.first;
    Image& image = reconstruction_->Image(image_id);
    init::Pose pose = best_poses.at(best_image_idx_map.at(image_id));
    image.SetQvec(RotationMatrixToQuaternion(pose.leftCols<3>()));
    image.SetTvec(pose.rightCols<1>());

    reconstruction_->RegisterImage(image_id);
    RegisterImageEvent(image_id);
  }

  // Triangulate points
  const auto tri_options = IncrementalTriangulator::Options();

  for (const image_t& image_id : reconstruction_->RegImageIds()) {
    TriangulateImage(tri_options, image_id);
  }

  // Merge and retriangulate
  CompleteTracks(tri_options);
  MergeTracks(tri_options);

  std::cout << "Initialized successfully (" << reconstruction_->NumPoints3D() << " points)\n";

  return true;
}

bool IncrementalMapper::RegisterNextImage(const Options& options,
                                          const image_t image_id) {
  CHECK_NOTNULL(reconstruction_);
  CHECK_GE(reconstruction_->NumRegImages(), 2);

  CHECK(options.Check());

  Image& image = reconstruction_->Image(image_id);
  Camera& camera = reconstruction_->Camera(image.CameraId());

  CHECK(!image.IsRegistered()) << "Image cannot be registered multiple times";

  num_reg_trials_[image_id] += 1;

  // Check if enough 2D-3D correspondences.
  if (image.NumVisiblePoints3D() <
      static_cast<size_t>(options.abs_pose_min_num_inliers)) {
    return false;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Search for 2D-3D correspondences
  //////////////////////////////////////////////////////////////////////////////

  const int kCorrTransitivity = 1;

  std::vector<std::pair<point2D_t, point3D_t>> tri_corrs;
  FeatureLines tri_lines2D;
  std::vector<Eigen::Vector3d> tri_lines2D_params;
  std::vector<Eigen::Vector3d> tri_points3D;

  for (point2D_t line_idx = 0; line_idx < image.NumLines();
       ++line_idx) {
    const FeatureLine& line = image.Line(line_idx);
    const CorrespondenceGraph& correspondence_graph =
        database_cache_->CorrespondenceGraph();
    const std::vector<CorrespondenceGraph::Correspondence> corrs =
        correspondence_graph.FindTransitiveCorrespondences(
            image_id, line_idx, kCorrTransitivity);

    std::unordered_set<point3D_t> point3D_ids;

    for (const auto corr : corrs) {
      const Image& corr_image = reconstruction_->Image(corr.image_id);
      if (!corr_image.IsRegistered()) {
        continue;
      }

      const FeatureLine& corr_line = corr_image.Line(corr.line_idx);
      if (!corr_line.HasPoint3D()) {
        continue;
      }

      // Avoid duplicate correspondences.
      if (point3D_ids.count(corr_line.Point3DId()) > 0) {
        continue;
      }

      const Camera& corr_camera =
          reconstruction_->Camera(corr_image.CameraId());

      // Avoid correspondences to images with bogus camera parameters.
      if (corr_camera.HasBogusParams(options.min_focal_length_ratio,
                                     options.max_focal_length_ratio,
                                     options.max_extra_param)) {
        continue;
      }

      const Point3D& point3D =
          reconstruction_->Point3D(corr_line.Point3DId());

      tri_corrs.emplace_back(line_idx, corr_line.Point3DId());
      point3D_ids.insert(corr_line.Point3DId());
      tri_lines2D.push_back(line);
      tri_lines2D_params.push_back(line.Line());
      tri_points3D.push_back(point3D.XYZ());
    }
  }

  // The size of `next_image.num_tri_obs` and `tri_corrs_line_idxs.size()`
  // can only differ, when there are images with bogus camera parameters, and
  // hence we skip some of the 2D-3D correspondences.
  // We need at least 6 correspondences to estimate a pose from lines
  if (tri_lines2D.size() <
          static_cast<size_t>(options.abs_pose_min_num_inliers) ||
          tri_lines2D.size() < 6) {
      return false;
  }

  //////////////////////////////////////////////////////////////////////////////
  // 2D-3D estimation
  //////////////////////////////////////////////////////////////////////////////

  // Only refine / estimate focal length, if no focal length was specified
  // (manually or through EXIF) and if it was not already estimated previously
  // from another image (when multiple images share the same camera
  // parameters)

  AbsolutePoseEstimationOptions abs_pose_options;
  abs_pose_options.num_threads = options.num_threads;
  abs_pose_options.num_focal_length_samples = 30;
  abs_pose_options.min_focal_length_ratio = options.min_focal_length_ratio;
  abs_pose_options.max_focal_length_ratio = options.max_focal_length_ratio;
  abs_pose_options.ransac_options.max_error =
      camera.ImageToWorldThreshold(options.abs_pose_max_error);
  abs_pose_options.ransac_options.min_inlier_ratio =
      options.abs_pose_min_inlier_ratio;
  // Use high confidence to avoid preemptive termination of P3P RANSAC
  // - too early termination may lead to bad registration.
  abs_pose_options.ransac_options.min_num_trials = 100;
  abs_pose_options.ransac_options.max_num_trials = 10000;
  abs_pose_options.ransac_options.confidence = 0.99999;

  AbsolutePoseRefinementOptions abs_pose_refinement_options;
  if (num_reg_images_per_camera_[image.CameraId()] > 0) {
    // Camera already refined from another image with the same camera.
    if (camera.HasBogusParams(options.min_focal_length_ratio,
                              options.max_focal_length_ratio,
                              options.max_extra_param)) {
      // Previously refined camera has bogus parameters,
      // so reset parameters and try to re-refine.
      camera.SetParams(database_cache_->Camera(image.CameraId()).Params());
      abs_pose_options.estimate_focal_length = !camera.HasPriorFocalLength();
      abs_pose_refinement_options.refine_focal_length = false;
      abs_pose_refinement_options.refine_extra_params = false;
    } else {
      abs_pose_options.estimate_focal_length = false;
      abs_pose_refinement_options.refine_focal_length = false;
      abs_pose_refinement_options.refine_extra_params = false;
    }
  } else {
    // Camera not refined before.
    abs_pose_options.estimate_focal_length = false;
    abs_pose_refinement_options.refine_focal_length = false;
    abs_pose_refinement_options.refine_extra_params = false;
  }

  if (!options.abs_pose_refine_focal_length) {
    abs_pose_options.estimate_focal_length = false;
    abs_pose_refinement_options.refine_focal_length = false;
  }

  if (!options.abs_pose_refine_extra_params) {
    abs_pose_refinement_options.refine_extra_params = false;
  }

  size_t num_inliers;
  std::vector<char> inlier_mask;

  if (!EstimateAbsolutePoseFromLines(
        abs_pose_options.ransac_options, tri_lines2D, tri_points3D,
        &image.Qvec(), &image.Tvec(), &num_inliers, &inlier_mask)) {
    return false;
  }

  if (num_inliers < static_cast<size_t>(options.abs_pose_min_num_inliers)) {
    return false;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Pose refinement
  //////////////////////////////////////////////////////////////////////////////

  if (!RefineAbsolutePoseFromLines(abs_pose_refinement_options, inlier_mask,
                                   tri_lines2D_params, tri_points3D, &image.Qvec(),
                                   &image.Tvec(), &camera)) {
    return false;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Continue tracks
  //////////////////////////////////////////////////////////////////////////////

  reconstruction_->RegisterImage(image_id);
  RegisterImageEvent(image_id);

  for (size_t i = 0; i < inlier_mask.size(); ++i) {
    if (inlier_mask[i]) {
      const point2D_t line_idx = tri_corrs[i].first;
      const FeatureLine& line = image.Line(line_idx);
      if (!line.HasPoint3D()) {
        const point3D_t point3D_id = tri_corrs[i].second;
        const TrackElement track_el(image_id, line_idx);
        reconstruction_->AddObservation(point3D_id, track_el);
        triangulator_->AddModifiedPoint3D(point3D_id);
      }
    }
  }

  return true;
}

size_t IncrementalMapper::TriangulateImage(
    const IncrementalTriangulator::Options& tri_options,
    const image_t image_id) {
  CHECK_NOTNULL(reconstruction_);
  return triangulator_->TriangulateImage(tri_options, image_id);
}

size_t IncrementalMapper::CompleteTracks(
    const IncrementalTriangulator::Options& tri_options) {
  CHECK_NOTNULL(reconstruction_);
  return triangulator_->CompleteAllTracks(tri_options);
}

size_t IncrementalMapper::MergeTracks(
    const IncrementalTriangulator::Options& tri_options) {
  CHECK_NOTNULL(reconstruction_);
  return triangulator_->MergeAllTracks(tri_options);
}

IncrementalMapper::LocalBundleAdjustmentReport
IncrementalMapper::AdjustLocalBundle(
    const Options& options, const BundleAdjustmentOptions& ba_options,
    const IncrementalTriangulator::Options& tri_options, const image_t image_id,
    const std::unordered_set<point3D_t>& point3D_ids) {
  CHECK_NOTNULL(reconstruction_);
  CHECK(options.Check());

  LocalBundleAdjustmentReport report;

  // Find images that have most 3D points with given image in common.
  const std::vector<image_t> local_bundle = FindLocalBundle(options, image_id);

  // Do the bundle adjustment only if there is any connected images.
  if (local_bundle.size() > 0) {
    BundleAdjustmentConfig ba_config;
    ba_config.AddImage(image_id);
    for (const image_t local_image_id : local_bundle) {
      ba_config.AddImage(local_image_id);
    }

    // Fix the existing images, if option specified.
    if (options.fix_existing_images) {
      for (const image_t local_image_id : local_bundle) {
        if (existing_image_ids_.count(local_image_id)) {
          ba_config.SetConstantPose(local_image_id);
        }
      }
    }

    // Determine which cameras to fix, when not all the registered images
    // are within the current local bundle.
    std::unordered_map<camera_t, size_t> num_images_per_camera;
    for (const image_t image_id : ba_config.Images()) {
      const Image& image = reconstruction_->Image(image_id);
      num_images_per_camera[image.CameraId()] += 1;
    }

    for (const auto& camera_id_and_num_images_pair : num_images_per_camera) {
      const size_t num_reg_images_for_camera =
          num_reg_images_per_camera_.at(camera_id_and_num_images_pair.first);\
      if (camera_id_and_num_images_pair.second < num_reg_images_for_camera) {
        ba_config.SetConstantCamera(camera_id_and_num_images_pair.first);
      }
    }

    // Fix 7 DOF to avoid scale/rotation/translation drift in bundle adjustment.
    if (local_bundle.size() == 1) {
      ba_config.SetConstantPose(local_bundle[0]);
      ba_config.SetConstantTvec(image_id, {0});
    } else if (local_bundle.size() > 1) {
      const image_t image_id1 = local_bundle[local_bundle.size() - 1];
      const image_t image_id2 = local_bundle[local_bundle.size() - 2];
      ba_config.SetConstantPose(image_id1);
      if (!options.fix_existing_images ||
          !existing_image_ids_.count(image_id2)) {
        ba_config.SetConstantTvec(image_id2, {0});
      }
    }

    // Make sure, we refine all new and short-track 3D points, no matter if
    // they are fully contained in the local image set or not. Do not include
    // long track 3D points as they are usually already very stable and adding
    // to them to bundle adjustment and track merging/completion would slow
    // down the local bundle adjustment significantly.
    std::unordered_set<point3D_t> variable_point3D_ids;
    for (const point3D_t point3D_id : point3D_ids) {
      const Point3D& point3D = reconstruction_->Point3D(point3D_id);
      const size_t kMaxTrackLength = 15;
      if (!point3D.HasError() || point3D.Track().Length() <= kMaxTrackLength) {
        ba_config.AddVariablePoint(point3D_id);
        variable_point3D_ids.insert(point3D_id);
      }
    }

    // Adjust the local bundle.
    BundleAdjuster bundle_adjuster(ba_options, ba_config);
    bundle_adjuster.Solve(reconstruction_);

    report.num_adjusted_observations =
        bundle_adjuster.Summary().num_residuals / 2;

    // Merge refined tracks with other existing points.
    report.num_merged_observations =
        triangulator_->MergeTracks(tri_options, variable_point3D_ids);
    // Complete tracks that may have failed to triangulate before refinement
    // of camera pose and calibration in bundle-adjustment. This may avoid
    // that some points are filtered and it helps for subsequent image
    // registrations.
    report.num_completed_observations =
        triangulator_->CompleteTracks(tri_options, variable_point3D_ids);
    report.num_completed_observations +=
        triangulator_->CompleteImage(tri_options, image_id);
  }

  // Filter both the modified images and all changed 3D points to make sure
  // there are no outlier points in the model. This results in duplicate work as
  // many of the provided 3D points may also be contained in the adjusted
  // images, but the filtering is not a bottleneck at this point.
  std::unordered_set<image_t> filter_image_ids;
  filter_image_ids.insert(image_id);
  filter_image_ids.insert(local_bundle.begin(), local_bundle.end());
  report.num_filtered_observations = reconstruction_->FilterPoints3DInImages(
      options.filter_max_reproj_error, options.filter_min_tri_angle,
      filter_image_ids);
  report.num_filtered_observations += reconstruction_->FilterPoints3D(
      options.filter_max_reproj_error, options.filter_min_tri_angle,
      point3D_ids);

  return report;
}

bool IncrementalMapper::AdjustGlobalBundle(
    const Options& options, const BundleAdjustmentOptions& ba_options) {
  CHECK_NOTNULL(reconstruction_);

  const std::vector<image_t>& reg_image_ids = reconstruction_->RegImageIds();

  CHECK_GE(reg_image_ids.size(), 2) << "At least two images must be "
                                       "registered for global "
                                       "bundle-adjustment";

  // Avoid degeneracies in bundle adjustment.
  reconstruction_->FilterObservationsWithNegativeDepth();

  // Configure bundle adjustment.
  BundleAdjustmentConfig ba_config;
  for (const image_t image_id : reg_image_ids) {
    ba_config.AddImage(image_id);
  }

  // Fix the existing images, if option specified.
  if (options.fix_existing_images) {
    for (const image_t image_id : reg_image_ids) {
      if (existing_image_ids_.count(image_id)) {
        ba_config.SetConstantPose(image_id);
      }
    }
  }

  // Fix 7-DOFs of the bundle adjustment problem.
  ba_config.SetConstantPose(reg_image_ids[0]);
  if (!options.fix_existing_images ||
      !existing_image_ids_.count(reg_image_ids[1])) {
    ba_config.SetConstantTvec(reg_image_ids[1], {0});
  }

  // Run bundle adjustment.
  BundleAdjuster bundle_adjuster(ba_options, ba_config);
  if (!bundle_adjuster.Solve(reconstruction_)) {
    return false;
  }

  // Normalize scene for numerical stability and
  // to avoid large scale changes in viewer.
  reconstruction_->Normalize();

  return true;
}

size_t IncrementalMapper::FilterImages(const Options& options) {
  CHECK_NOTNULL(reconstruction_);
  CHECK(options.Check());

  // Do not filter images in the early stage of the reconstruction, since the
  // calibration is often still refining a lot. Hence, the camera parameters
  // are not stable in the beginning.
  const size_t kMinNumImages = 20;
  if (reconstruction_->NumRegImages() < kMinNumImages) {
    return {};
  }

  const std::vector<image_t> image_ids = reconstruction_->FilterImages(
      options.min_focal_length_ratio, options.max_focal_length_ratio,
      options.max_extra_param);

  for (const image_t image_id : image_ids) {
    DeRegisterImageEvent(image_id);
    filtered_images_.insert(image_id);
  }

  return image_ids.size();
}

size_t IncrementalMapper::FilterPoints(const Options& options) {
  CHECK_NOTNULL(reconstruction_);
  CHECK(options.Check());
  return reconstruction_->FilterAllPoints3D(options.filter_max_reproj_error,
                                            options.filter_min_tri_angle);
}

const Reconstruction& IncrementalMapper::GetReconstruction() const {
  CHECK_NOTNULL(reconstruction_);
  return *reconstruction_;
}

size_t IncrementalMapper::NumTotalRegImages() const {
  return num_total_reg_images_;
}

size_t IncrementalMapper::NumSharedRegImages() const {
  return num_shared_reg_images_;
}

const std::unordered_set<point3D_t>& IncrementalMapper::GetModifiedPoints3D() {
  return triangulator_->GetModifiedPoints3D();
}

void IncrementalMapper::ClearModifiedPoints3D() {
  triangulator_->ClearModifiedPoints3D();
}

std::vector<image_t> IncrementalMapper::FindLocalBundle(
    const Options& options, const image_t image_id) const {
  CHECK(options.Check());

  const Image& image = reconstruction_->Image(image_id);
  CHECK(image.IsRegistered());

  // Extract all images that have at least one 3D point with the query image
  // in common, and simultaneously count the number of common 3D points.

  std::unordered_map<image_t, size_t> shared_observations;

  std::unordered_set<point3D_t> point3D_ids;
  point3D_ids.reserve(image.NumPoints3D());

  for (const FeatureLine& line : image.Lines()) {
    if (line.HasPoint3D()) {
      point3D_ids.insert(line.Point3DId());
      const Point3D& point3D = reconstruction_->Point3D(line.Point3DId());
      for (const TrackElement& track_el : point3D.Track().Elements()) {
        if (track_el.image_id != image_id) {
          shared_observations[track_el.image_id] += 1;
        }
      }
    }
  }

  // Sort overlapping images according to number of shared observations.

  std::vector<std::pair<image_t, size_t>> overlapping_images(
      shared_observations.begin(), shared_observations.end());
  std::sort(overlapping_images.begin(), overlapping_images.end(),
            [](const std::pair<image_t, size_t>& image1,
               const std::pair<image_t, size_t>& image2) {
              return image1.second > image2.second;
            });

  // The local bundle is composed of the given image and its most connected
  // neighbor images, hence the subtraction of 1.

  const size_t num_images =
      static_cast<size_t>(options.local_ba_num_images - 1);
  const size_t num_eff_images = std::min(num_images, overlapping_images.size());

  // Extract most connected images and ensure sufficient triangulation angle.

  std::vector<image_t> local_bundle_image_ids;
  local_bundle_image_ids.reserve(num_eff_images);

  // If the number of overlapping images equals the number of desired images in
  // the local bundle, then simply copy over the image identifiers.
  if (overlapping_images.size() == num_eff_images) {
    for (const auto& overlapping_image : overlapping_images) {
      local_bundle_image_ids.push_back(overlapping_image.first);
    }
    return local_bundle_image_ids;
  }

  // In the following iteration, we start with the most overlapping images and
  // check whether it has sufficient triangulation angle. If none of the
  // overlapping images has sufficient triangulation angle, we relax the
  // triangulation angle threshold and start from the most overlapping image
  // again. In the end, if we still haven't found enough images, we simply use
  // the most overlapping images.

  const double min_tri_angle_rad = DegToRad(options.local_ba_min_tri_angle);

  // The selection thresholds (minimum triangulation angle, minimum number of
  // shared observations), which are successively relaxed.
  const std::array<std::pair<double, double>, 8> selection_thresholds = {{
      std::make_pair(min_tri_angle_rad / 1.0, 0.6 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 1.5, 0.6 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 2.0, 0.5 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 2.5, 0.4 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 3.0, 0.3 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 4.0, 0.2 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 5.0, 0.1 * image.NumPoints3D()),
      std::make_pair(min_tri_angle_rad / 6.0, 0.1 * image.NumPoints3D()),
  }};

  const Eigen::Vector3d proj_center = image.ProjectionCenter();
  std::vector<Eigen::Vector3d> shared_points3D;
  shared_points3D.reserve(image.NumPoints3D());
  std::vector<double> tri_angles(overlapping_images.size(), -1.0);
  std::vector<char> used_overlapping_images(overlapping_images.size(), false);

  for (const auto& selection_threshold : selection_thresholds) {
    for (size_t overlapping_image_idx = 0;
         overlapping_image_idx < overlapping_images.size();
         ++overlapping_image_idx) {
      // Check if the image has sufficient overlap. Since the images are ordered
      // based on the overlap, we can just skip the remaining ones.
      if (overlapping_images[overlapping_image_idx].second <
          selection_threshold.second) {
        break;
      }

      // Check if the image is already in the local bundle.
      if (used_overlapping_images[overlapping_image_idx]) {
        continue;
      }

      const auto& overlapping_image = reconstruction_->Image(
          overlapping_images[overlapping_image_idx].first);
      const Eigen::Vector3d overlapping_proj_center =
          overlapping_image.ProjectionCenter();

      // In the first iteration, compute the triangulation angle. In later
      // iterations, reuse the previously computed value.
      double& tri_angle = tri_angles[overlapping_image_idx];
      if (tri_angle < 0.0) {
        // Collect the commonly observed 3D points.
        shared_points3D.clear();
        for (const FeatureLine& line : image.Lines()) {
          if (line.HasPoint3D() && point3D_ids.count(line.Point3DId())) {
            shared_points3D.push_back(
                reconstruction_->Point3D(line.Point3DId()).XYZ());
          }
        }

        // Calculate the triangulation angle at a certain percentile.
        const double kTriangulationAnglePercentile = 75;
        tri_angle = Percentile(
            CalculateTriangulationAngles(proj_center, overlapping_proj_center,
                                         shared_points3D),
            kTriangulationAnglePercentile);
      }

      // Check that the image has sufficient triangulation angle.
      if (tri_angle >= selection_threshold.first) {
        local_bundle_image_ids.push_back(overlapping_image.ImageId());
        used_overlapping_images[overlapping_image_idx] = true;
        // Check if we already collected enough images.
        if (local_bundle_image_ids.size() >= num_eff_images) {
          break;
        }
      }
    }

    // Check if we already collected enough images.
    if (local_bundle_image_ids.size() >= num_eff_images) {
      break;
    }
  }

  // In case there are not enough images with sufficient triangulation angle,
  // simply fill up the rest with the most overlapping images.

  if (local_bundle_image_ids.size() < num_eff_images) {
    for (size_t overlapping_image_idx = 0;
         overlapping_image_idx < overlapping_images.size();
         ++overlapping_image_idx) {
      // Collect image if it is not yet in the local bundle.
      if (!used_overlapping_images[overlapping_image_idx]) {
        local_bundle_image_ids.push_back(
            overlapping_images[overlapping_image_idx].first);
        used_overlapping_images[overlapping_image_idx] = true;

        // Check if we already collected enough images.
        if (local_bundle_image_ids.size() >= num_eff_images) {
          break;
        }
      }
    }
  }

  return local_bundle_image_ids;
}

void IncrementalMapper::RegisterImageEvent(const image_t image_id) {
  const Image& image = reconstruction_->Image(image_id);
  size_t& num_reg_images_for_camera =
      num_reg_images_per_camera_[image.CameraId()];
  num_reg_images_for_camera += 1;

  size_t& num_regs_for_image = num_registrations_[image_id];
  num_regs_for_image += 1;
  if (num_regs_for_image == 1) {
    num_total_reg_images_ += 1;
  } else if (num_regs_for_image > 1) {
    num_shared_reg_images_ += 1;
  }
}

void IncrementalMapper::DeRegisterImageEvent(const image_t image_id) {
  const Image& image = reconstruction_->Image(image_id);
  size_t& num_reg_images_for_camera =
      num_reg_images_per_camera_.at(image.CameraId());
  CHECK_GT(num_reg_images_for_camera, 0);
  num_reg_images_for_camera -= 1;

  size_t& num_regs_for_image = num_registrations_[image_id];
  num_regs_for_image -= 1;
  if (num_regs_for_image == 0) {
    num_total_reg_images_ -= 1;
  } else if (num_regs_for_image > 0) {
    num_shared_reg_images_ -= 1;
  }
}

}  // namespace colmap
