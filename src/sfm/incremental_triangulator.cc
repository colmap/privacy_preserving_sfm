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

#include "sfm/incremental_triangulator.h"

#include "base/projection.h"
#include "estimators/triangulation.h"
#include "util/misc.h"

namespace colmap {

bool IncrementalTriangulator::Options::Check() const {
  CHECK_OPTION_GE(max_transitivity, 0);
  CHECK_OPTION_GT(create_max_angle_error, 0);
  CHECK_OPTION_GT(continue_max_angle_error, 0);
  CHECK_OPTION_GT(merge_max_reproj_error, 0);
  CHECK_OPTION_GT(complete_max_reproj_error, 0);
  CHECK_OPTION_GE(complete_max_transitivity, 0);
  CHECK_OPTION_GT(re_max_angle_error, 0);
  CHECK_OPTION_GE(re_min_ratio, 0);
  CHECK_OPTION_LE(re_min_ratio, 1);
  CHECK_OPTION_GE(re_max_trials, 0);
  CHECK_OPTION_GT(min_angle, 0);
  return true;
}

IncrementalTriangulator::IncrementalTriangulator(
    const CorrespondenceGraph* correspondence_graph,
    Reconstruction* reconstruction)
    : correspondence_graph_(correspondence_graph),
      reconstruction_(reconstruction) {}

size_t IncrementalTriangulator::TriangulateImage(const Options& options,
                                                 const image_t image_id) {
  CHECK(options.Check());

  size_t num_tris = 0;

  ClearCaches();

  const Image& image = reconstruction_->Image(image_id);
  if (!image.IsRegistered()) {
    return num_tris;
  }

  const Camera& camera = reconstruction_->Camera(image.CameraId());
  if (HasCameraBogusParams(options, camera)) {
    return num_tris;
  }

  // Correspondence data for reference observation in given image. We iterate
  // over all observations of the image and each observation once becomes
  // the reference correspondence.
  CorrData ref_corr_data;
  ref_corr_data.image_id = image_id;
  ref_corr_data.image = &image;
  ref_corr_data.camera = &camera;

  // Container for correspondences from reference observation to other images.
  std::vector<CorrData> corrs_data;

  // Try to triangulate all image observations.
  for (point2D_t line_idx = 0; line_idx < image.NumLines();
       ++line_idx) {
    const size_t num_triangulated =
        Find(options, image_id, line_idx,
             static_cast<size_t>(options.max_transitivity), &corrs_data);
    if (corrs_data.empty()) {
      continue;
    }



    const FeatureLine& line = image.Line(line_idx);
    ref_corr_data.line_idx = line_idx;
    ref_corr_data.line = &line;

    if (num_triangulated == 0) {
      corrs_data.push_back(ref_corr_data);
      num_tris += Create(options, corrs_data);
    } else {
      // Continue correspondences to existing 3D points.
      num_tris += Continue(options, ref_corr_data, corrs_data);
      // Create points from correspondences that are not continued.
      corrs_data.push_back(ref_corr_data);
      num_tris += Create(options, corrs_data);
    }
  }

  return num_tris;
}

size_t IncrementalTriangulator::CompleteImage(const Options& options,
                                              const image_t image_id) {  
  CHECK(options.Check());

  size_t num_tris = 0;

  ClearCaches();

  const Image& image = reconstruction_->Image(image_id);
  if (!image.IsRegistered()) {
    return num_tris;
  }

  const Camera& camera = reconstruction_->Camera(image.CameraId());
  if (HasCameraBogusParams(options, camera)) {
    return num_tris;
  }

  // Setup estimation options.
  EstimateTriangulationOptions tri_options;
  tri_options.min_tri_angle = DegToRad(options.min_angle);
  tri_options.residual_type =
      TriangulationEstimator::ResidualType::REPROJECTION_ERROR;
  tri_options.ransac_options.max_error = options.complete_max_reproj_error;
  tri_options.ransac_options.confidence = 0.9999;
  tri_options.ransac_options.min_inlier_ratio = 0.02;
  tri_options.ransac_options.max_num_trials = 10000;

  // Correspondence data for reference observation in given image. We iterate
  // over all observations of the image and each observation once becomes
  // the reference correspondence.
  CorrData ref_corr_data;
  ref_corr_data.image_id = image_id;
  ref_corr_data.image = &image;
  ref_corr_data.camera = &camera;

  // Container for correspondences from reference observation to other images.
  std::vector<CorrData> corrs_data;

  for (point2D_t line_idx = 0; line_idx < image.NumLines();
       ++line_idx) {
    const FeatureLine& line = image.Line(line_idx);
    if (line.HasPoint3D()) {
      // Complete existing track.
      num_tris += Complete(options, line.Point3DId());
      continue;
    }

    if (options.ignore_two_view_tracks &&
        correspondence_graph_->IsTwoViewObservation(image_id, line_idx)) {
      continue;
    }

    const size_t num_triangulated =
        Find(options, image_id, line_idx,
             static_cast<size_t>(options.max_transitivity), &corrs_data);
    if (num_triangulated || corrs_data.empty()) {
      continue;
    }

    ref_corr_data.line = &line;
    ref_corr_data.line_idx = line_idx;
    corrs_data.push_back(ref_corr_data);

    // Setup data for triangulation estimation.
    std::vector<TriangulationEstimator::PointData> point_data;
    point_data.resize(corrs_data.size());
    std::vector<TriangulationEstimator::PoseData> pose_data;
    pose_data.resize(corrs_data.size());
    for (size_t i = 0; i < corrs_data.size(); ++i) {
      const CorrData& corr_data = corrs_data[i];
      //point_data[i].point = corr_data.point2D->XY();
      //point_data[i].point_normalized =
      //    corr_data.camera->ImageToWorld(point_data[i].point);
      // TODO: I broke it here
      point_data[i].line = corr_data.image->Line(corr_data.line_idx).Line();
      pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
      pose_data[i].proj_center = corr_data.image->ProjectionCenter();
      pose_data[i].camera = corr_data.camera;
    }

    // Enforce exhaustive sampling for small track lengths.
    const size_t kExhaustiveSamplingThreshold = 15;
    if (point_data.size() <= kExhaustiveSamplingThreshold) {
      tri_options.ransac_options.min_num_trials =
          NChooseK(point_data.size(), 2);
    }

    // Estimate triangulation.
    Eigen::Vector3d xyz;
    std::vector<char> inlier_mask;
    if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
                               &xyz)) {
      continue;
    }

    // Add inliers to estimated track.
    Track track;
    track.Reserve(corrs_data.size());
    for (size_t i = 0; i < inlier_mask.size(); ++i) {
      if (inlier_mask[i]) {
        const CorrData& corr_data = corrs_data[i];
        track.AddElement(corr_data.image_id, corr_data.line_idx);
        num_tris += 1;
      }
    }

    const point3D_t point3D_id = reconstruction_->AddPoint3D(xyz, track);
    modified_point3D_ids_.insert(point3D_id);
  }

  return num_tris;
}

size_t IncrementalTriangulator::CompleteTracks(
    const Options& options, const std::unordered_set<point3D_t>& point3D_ids) {
  CHECK(options.Check());

  size_t num_completed = 0;

  ClearCaches();

  for (const point3D_t point3D_id : point3D_ids) {
    num_completed += Complete(options, point3D_id);
  }

  return num_completed;
}

size_t IncrementalTriangulator::CompleteAllTracks(const Options& options) {
  CHECK(options.Check());

  size_t num_completed = 0;

  ClearCaches();

  for (const point3D_t point3D_id : reconstruction_->Point3DIds()) {
    num_completed += Complete(options, point3D_id);
  }

  return num_completed;
}

size_t IncrementalTriangulator::MergeTracks(
    const Options& options, const std::unordered_set<point3D_t>& point3D_ids) {
  CHECK(options.Check());

  size_t num_merged = 0;

  ClearCaches();

  for (const point3D_t point3D_id : point3D_ids) {
    num_merged += Merge(options, point3D_id);
  }

  return num_merged;
}

size_t IncrementalTriangulator::MergeAllTracks(const Options& options) {
  CHECK(options.Check());

  size_t num_merged = 0;

  ClearCaches();

  for (const point3D_t point3D_id : reconstruction_->Point3DIds()) {
    num_merged += Merge(options, point3D_id);
  }

  return num_merged;
}

//size_t IncrementalTriangulator::Retriangulate(const Options& options) {
//
//  CHECK(options.Check());
//
//  size_t num_tris = 0;
//
//  ClearCaches();
//
//  Options re_options = options;
//  re_options.continue_max_angle_error = options.re_max_angle_error;
//
//  for (const auto& image_pair : reconstruction_->ImagePairs()) {
//    // Only perform retriangulation for under-reconstructed image pairs.
//    const double tri_ratio =
//        static_cast<double>(image_pair.second.num_tri_corrs) /
//        static_cast<double>(image_pair.second.num_total_corrs);
//    if (tri_ratio >= options.re_min_ratio) {
//      continue;
//    }
//
//    // Check if images are registered yet.
//
//    image_t image_id1;
//    image_t image_id2;
//    Database::PairIdToImagePair(image_pair.first, &image_id1, &image_id2);
//
//    const Image& image1 = reconstruction_->Image(image_id1);
//    if (!image1.IsRegistered()) {
//      continue;
//    }
//
//    const Image& image2 = reconstruction_->Image(image_id2);
//    if (!image2.IsRegistered()) {
//      continue;
//    }
//
//    // Only perform retriangulation for a maximum number of trials.
//
//    int& num_re_trials = re_num_trials_[image_pair.first];
//    if (num_re_trials >= options.re_max_trials) {
//      continue;
//    }
//    num_re_trials += 1;
//
//    const Camera& camera1 = reconstruction_->Camera(image1.CameraId());
//    const Camera& camera2 = reconstruction_->Camera(image2.CameraId());
//    if (HasCameraBogusParams(options, camera1) ||
//        HasCameraBogusParams(options, camera2)) {
//      continue;
//    }
//
//    // Find correspondences and perform retriangulation.
//
//    const FeatureMatches& corrs =
//        correspondence_graph_->FindCorrespondencesBetweenImages(image_id1,
//                                                                image_id2);
//
//    for (const auto& corr : corrs) {
//      const Point2D& point2D1 = image1.Point2D(corr.point2D_idx1);
//      const Point2D& point2D2 = image2.Point2D(corr.point2D_idx2);
//
//      // Two cases are possible here: both points belong to the same 3D point
//      // or to different 3D points. In the former case, there is nothing
//      // to do. In the latter case, we do not attempt retriangulation,
//      // as retriangulated correspondences are very likely bogus and
//      // would therefore destroy both 3D points if merged.
//      if (point2D1.HasPoint3D() && point2D2.HasPoint3D()) {
//        continue;
//      }
//
//      CorrData corr_data1;
//      corr_data1.image_id = image_id1;
//      corr_data1.point2D_idx = corr.point2D_idx1;
//      corr_data1.image = &image1;
//      corr_data1.camera = &camera1;
//      corr_data1.point2D = &point2D1;
//
//      CorrData corr_data2;
//      corr_data2.image_id = image_id2;
//      corr_data2.point2D_idx = corr.point2D_idx2;
//      corr_data2.image = &image2;
//      corr_data2.camera = &camera2;
//      corr_data2.point2D = &point2D2;
//
//      if (point2D1.HasPoint3D() && !point2D2.HasPoint3D()) {
//        const std::vector<CorrData> corrs_data1 = {corr_data1};
//        num_tris += Continue(re_options, corr_data2, corrs_data1);
//      } else if (!point2D1.HasPoint3D() && point2D2.HasPoint3D()) {
//        const std::vector<CorrData> corrs_data2 = {corr_data2};
//        num_tris += Continue(re_options, corr_data1, corrs_data2);
//      } else if (!point2D1.HasPoint3D() && !point2D2.HasPoint3D()) {
//        const std::vector<CorrData> corrs_data = {corr_data1, corr_data2};
//        // Do not use larger triangulation threshold as this causes
//        // significant drift when creating points (options vs. re_options).
//        num_tris += Create(options, corrs_data);
//      }
//      // Else both points have a 3D point, but we do not want to
//      // merge points in retriangulation.
//    }
//  }
//
//  return num_tris;
//}

void IncrementalTriangulator::AddModifiedPoint3D(const point3D_t point3D_id) {
  modified_point3D_ids_.insert(point3D_id);
}

const std::unordered_set<point3D_t>&
IncrementalTriangulator::GetModifiedPoints3D() {
  // First remove any missing 3D points from the set.
  for (auto it = modified_point3D_ids_.begin();
       it != modified_point3D_ids_.end();) {
    if (reconstruction_->ExistsPoint3D(*it)) {
      ++it;
    } else {
      modified_point3D_ids_.erase(it++);
    }
  }
  return modified_point3D_ids_;
}

void IncrementalTriangulator::ClearModifiedPoints3D() {
  modified_point3D_ids_.clear();
}

void IncrementalTriangulator::ClearCaches() {
  camera_has_bogus_params_.clear();
  merge_trials_.clear();
}

size_t IncrementalTriangulator::Find(const Options& options,
                                     const image_t image_id,
                                     const point2D_t line_idx,
                                     const size_t transitivity,
                                     std::vector<CorrData>* corrs_data) {
  const std::vector<CorrespondenceGraph::Correspondence>& corrs =
      correspondence_graph_->FindTransitiveCorrespondences(
          image_id, line_idx, transitivity);

  corrs_data->clear();
  corrs_data->reserve(corrs.size());

  size_t num_triangulated = 0;

  for (const CorrespondenceGraph::Correspondence corr : corrs) {
    const Image& corr_image = reconstruction_->Image(corr.image_id);
    if (!corr_image.IsRegistered()) {
      continue;
    }

    const Camera& corr_camera = reconstruction_->Camera(corr_image.CameraId());
    if (HasCameraBogusParams(options, corr_camera)) {
      continue;
    }

    CorrData corr_data;
    corr_data.image_id = corr.image_id;
    corr_data.line_idx = corr.line_idx;
    corr_data.image = &corr_image;
    corr_data.camera = &corr_camera;
    corr_data.line = &corr_image.Line(corr.line_idx);

    corrs_data->push_back(corr_data);

    if (corr_data.line->HasPoint3D()) {
      num_triangulated += 1;
    }
  }

  return num_triangulated;
}

size_t IncrementalTriangulator::Create(
    const Options& options, const std::vector<CorrData>& corrs_data) {
  // Extract correspondences without an existing triangulated observation.
  std::vector<CorrData> create_corrs_data;
  create_corrs_data.reserve(corrs_data.size());
  for (const CorrData& corr_data : corrs_data) {
    if (!corr_data.line->HasPoint3D()) {
      create_corrs_data.push_back(corr_data);
    }
  }

  // TODO: Viktor changed things here
  if (create_corrs_data.size() < 3) {
    // Need at least three observations for triangulation.
    return 0;
  } else if (options.ignore_two_view_tracks && create_corrs_data.size() == 2) {
    const CorrData& corr_data1 = create_corrs_data[0];
    if (correspondence_graph_->IsTwoViewObservation(corr_data1.image_id,
                                                    corr_data1.line_idx)) {
      return 0;
    }
  }

  // Setup data for triangulation estimation.
  std::vector<TriangulationEstimator::PointData> point_data;
  point_data.resize(create_corrs_data.size());
  std::vector<TriangulationEstimator::PoseData> pose_data;
  pose_data.resize(create_corrs_data.size());
  int num_random_lines = 0;
  for (size_t i = 0; i < create_corrs_data.size(); ++i) {
    const CorrData& corr_data = create_corrs_data[i];    

    point_data[i].line = corr_data.image->Line(corr_data.line_idx).Line();

    if (!corr_data.image->Line(corr_data.line_idx).IsAligned()) {
      num_random_lines += 1;
    }
        
    pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
    pose_data[i].proj_center = corr_data.image->ProjectionCenter();
    pose_data[i].camera = corr_data.camera;
  }

  // We can not triangulate points from only aligned lines
  if (num_random_lines < 1) {
    return 0;
  }

  // Setup estimation options.
  EstimateTriangulationOptions tri_options;
  tri_options.min_tri_angle = DegToRad(options.min_angle);
  tri_options.residual_type =
      TriangulationEstimator::ResidualType::ANGULAR_ERROR;
  tri_options.ransac_options.max_error =
      DegToRad(options.create_max_angle_error);
  tri_options.ransac_options.confidence = 0.9999;
  tri_options.ransac_options.min_inlier_ratio = 0.02;
  tri_options.ransac_options.max_num_trials = 10000;

  // Enforce exhaustive sampling for small track lengths.
  const size_t kExhaustiveSamplingThreshold = 15;
  if (point_data.size() <= kExhaustiveSamplingThreshold) {
    tri_options.ransac_options.min_num_trials = NChooseK(point_data.size(), 3);
  }

  // Estimate triangulation.
  Eigen::Vector3d xyz;
  std::vector<char> inlier_mask;
  if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
                             &xyz)) {
    return 0;
  }

  // Add inliers to estimated track.
  Track track;
  track.Reserve(create_corrs_data.size());
  for (size_t i = 0; i < inlier_mask.size(); ++i) {
    if (inlier_mask[i]) {
      const CorrData& corr_data = create_corrs_data[i];
      track.AddElement(corr_data.image_id, corr_data.line_idx);
    }
  }

  // Add estimated point to reconstruction.
  const point3D_t point3D_id = reconstruction_->AddPoint3D(xyz, track);
  modified_point3D_ids_.insert(point3D_id);

  const size_t kMinRecursiveTrackLength = 3;
  if (create_corrs_data.size() - track.Length() >= kMinRecursiveTrackLength) {
    return track.Length() + Create(options, create_corrs_data);
  }

  return track.Length();
}

size_t IncrementalTriangulator::Continue(
    const Options& options, const CorrData& ref_corr_data,
    const std::vector<CorrData>& corrs_data) {
  // No need to continue, if the reference observation is triangulated.
  if (ref_corr_data.line->HasPoint3D()) {
    return 0;
  }

  double best_angle_error = std::numeric_limits<double>::max();
  size_t best_idx = std::numeric_limits<size_t>::max();

  for (size_t idx = 0; idx < corrs_data.size(); ++idx) {
    const CorrData& corr_data = corrs_data[idx];
    if (!corr_data.line->HasPoint3D()) {
      continue;
    }

    const Point3D& point3D =
        reconstruction_->Point3D(corr_data.line->Point3DId());

    const double angle_error = CalculateNormalizedLineAngularError(
        ref_corr_data.image->Line(ref_corr_data.line_idx).Line(), point3D.XYZ(), ref_corr_data.image->Qvec(),
        ref_corr_data.image->Tvec(), *ref_corr_data.camera);
    if (angle_error < best_angle_error) {
      best_angle_error = angle_error;
      best_idx = idx;
    }
  }

  const double max_angle_error = DegToRad(options.continue_max_angle_error);
  if (best_angle_error <= max_angle_error &&
      best_idx != std::numeric_limits<size_t>::max()) {
    const CorrData& corr_data = corrs_data[best_idx];
    const TrackElement track_el(ref_corr_data.image_id,
                                ref_corr_data.line_idx);
    reconstruction_->AddObservation(corr_data.line->Point3DId(), track_el);
    modified_point3D_ids_.insert(corr_data.line->Point3DId());
    return 1;
  }

  return 0;
}

size_t IncrementalTriangulator::Merge(const Options& options,
                                      const point3D_t point3D_id) {
  if (!reconstruction_->ExistsPoint3D(point3D_id)) {
    return 0;
  }

  const double max_squared_reproj_error =
      options.merge_max_reproj_error * options.merge_max_reproj_error;

  const auto& point3D = reconstruction_->Point3D(point3D_id);

  for (const auto& track_el : point3D.Track().Elements()) {
    const std::vector<CorrespondenceGraph::Correspondence>& corrs =
        correspondence_graph_->FindCorrespondences(track_el.image_id,
                                                   track_el.line_idx);

    for (const auto corr : corrs) {
      const auto& image = reconstruction_->Image(corr.image_id);
      if (!image.IsRegistered()) {
        continue;
      }

      const FeatureLine& corr_line = image.Line(corr.line_idx);
      if (!corr_line.HasPoint3D() ||
          corr_line.Point3DId() == point3D_id ||
          merge_trials_[point3D_id].count(corr_line.Point3DId()) > 0) {
        continue;
      }

      // Try to merge the two 3D points.

      const Point3D& corr_point3D =
          reconstruction_->Point3D(corr_line.Point3DId());

      merge_trials_[point3D_id].insert(corr_line.Point3DId());
      merge_trials_[corr_line.Point3DId()].insert(point3D_id);

      // Weighted average of point locations, depending on track length.
      const Eigen::Vector3d merged_xyz =
          (point3D.Track().Length() * point3D.XYZ() +
           corr_point3D.Track().Length() * corr_point3D.XYZ()) /
          (point3D.Track().Length() + corr_point3D.Track().Length());

      // Count number of inlier track elements of the merged track.
      bool merge_success = true;
      for (const Track* track : {&point3D.Track(), &corr_point3D.Track()}) {
        for (const auto test_track_el : track->Elements()) {
          const Image& test_image =
              reconstruction_->Image(test_track_el.image_id);
          const Camera& test_camera =
              reconstruction_->Camera(test_image.CameraId());
          const FeatureLine test_line2D =
                  test_image.Line(test_track_el.line_idx);
          if (CalculateSquaredLineReprojectionError(
                  test_line2D.Line(), merged_xyz, test_image.Qvec(),
                  test_image.Tvec(), test_camera) > max_squared_reproj_error) {
            merge_success = false;
            break;
          }
        }
        if (!merge_success) {
          break;
        }
      }

      // Only accept merge if all track elements are inliers.
      if (merge_success) {
        const size_t num_merged =
            point3D.Track().Length() + corr_point3D.Track().Length();

        const point3D_t merged_point3D_id = reconstruction_->MergePoints3D(
            point3D_id, corr_line.Point3DId());

        modified_point3D_ids_.erase(point3D_id);
        modified_point3D_ids_.erase(corr_line.Point3DId());
        modified_point3D_ids_.insert(merged_point3D_id);

        // Merge merged 3D point and return, as the original points are deleted.
        const size_t num_merged_recursive = Merge(options, merged_point3D_id);
        if (num_merged_recursive > 0) {
          return num_merged_recursive;
        } else {
          return num_merged;
        }
      }
    }
  }

  return 0;
}

size_t IncrementalTriangulator::Complete(const Options& options,
                                         const point3D_t point3D_id) {
  size_t num_completed = 0;

  if (!reconstruction_->ExistsPoint3D(point3D_id)) {
    return num_completed;
  }

  const double max_squared_reproj_error =
      options.complete_max_reproj_error * options.complete_max_reproj_error;

  const Point3D& point3D = reconstruction_->Point3D(point3D_id);

  std::vector<TrackElement> queue = point3D.Track().Elements();

  const int max_transitivity = options.complete_max_transitivity;
  for (int transitivity = 0; transitivity < max_transitivity; ++transitivity) {
    if (queue.empty()) {
      break;
    }

    const auto prev_queue = queue;
    queue.clear();

    for (const TrackElement queue_elem : prev_queue) {
      const std::vector<CorrespondenceGraph::Correspondence>& corrs =
          correspondence_graph_->FindCorrespondences(queue_elem.image_id,
                                                     queue_elem.line_idx);

      for (const auto corr : corrs) {
        const Image& image = reconstruction_->Image(corr.image_id);
        if (!image.IsRegistered()) {
          continue;
        }

        const FeatureLine& line = image.Line(corr.line_idx);
        const Eigen::Vector3d& line2D = image.Line(corr.line_idx).Line();
        if (line.HasPoint3D()) {
          continue;
        }

        const Camera& camera = reconstruction_->Camera(image.CameraId());
        if (HasCameraBogusParams(options, camera)) {
          continue;
        }
        // TODO: correct?
        if (CalculateSquaredLineReprojectionError(
                line2D, point3D.XYZ(), image.ProjectionMatrix(), camera)
                > max_squared_reproj_error) {
          continue;
        }

        // Success, add observation to point track.
        const TrackElement track_el(corr.image_id, corr.line_idx);
        reconstruction_->AddObservation(point3D_id, track_el);
        modified_point3D_ids_.insert(point3D_id);

        // Recursively complete track for this new correspondence.
        if (transitivity < max_transitivity - 1) {
          queue.emplace_back(corr.image_id, corr.line_idx);
        }

        num_completed += 1;
      }
    }
  }

  return num_completed;
}

bool IncrementalTriangulator::HasCameraBogusParams(const Options& options,
                                                   const Camera& camera) {
  const auto it = camera_has_bogus_params_.find(camera.CameraId());
  if (it == camera_has_bogus_params_.end()) {
    const bool has_bogus_params = camera.HasBogusParams(
        options.min_focal_length_ratio, options.max_focal_length_ratio,
        options.max_extra_param);
    camera_has_bogus_params_.emplace(camera.CameraId(), has_bogus_params);
    return has_bogus_params;
  } else {
    return it->second;
  }
}

}  // namespace colmap
