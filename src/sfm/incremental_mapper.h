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

#ifndef COLMAP_SRC_SFM_INCREMENTAL_MAPPER_H_
#define COLMAP_SRC_SFM_INCREMENTAL_MAPPER_H_

#include "base/database.h"
#include "base/database_cache.h"
#include "base/reconstruction.h"
#include "optim/bundle_adjustment.h"
#include "sfm/incremental_triangulator.h"
#include "util/alignment.h"

namespace colmap {

// Class that provides all functionality for the incremental reconstruction
// procedure.
class IncrementalMapper {
 public:
  struct Options {
    // Minimum number of inliers for initial image pair.
    int init_min_num_inliers = 20;

    // Maximum error in pixels for relative pose estimation for the initial
    // image set.
    double init_max_error = 5.0;

    // Minimum triangulation angle for initial image pair.
    double init_min_tri_angle = 2.0;

    // Maximum reprojection error in absolute pose estimation.
    double abs_pose_max_error = 12.0;

    // Minimum number of inliers in absolute pose estimation.
    int abs_pose_min_num_inliers = 30;

    // Minimum inlier ratio in absolute pose estimation.
    double abs_pose_min_inlier_ratio = 0.25;

    // Whether to estimate the focal length in absolute pose estimation.
    bool abs_pose_refine_focal_length = false;

    // Whether to estimate the extra parameters in absolute pose estimation.
    bool abs_pose_refine_extra_params = false;

    // Number of images to optimize in local bundle adjustment.
    int local_ba_num_images = 6;

    // Minimum triangulation for images to be chosen in local bundle adjustment.
    double local_ba_min_tri_angle = 6;

    // Thresholds for bogus camera parameters. Images with bogus camera
    // parameters are filtered and ignored in triangulation.
    double min_focal_length_ratio = 0.1;  // Opening angle of ~130deg
    double max_focal_length_ratio = 10;   // Opening angle of ~5deg
    double max_extra_param = 1;

    // Maximum reprojection error in pixels for observations.
    double filter_max_reproj_error = 4.0;

    // Minimum triangulation angle in degrees for stable 3D points.
    double filter_min_tri_angle = 1.5;

    // Maximum number of trials to register an image.
    int max_reg_trials = 3;

    // If reconstruction is provided as input, fix the existing image poses.
    bool fix_existing_images = false;

    // Number of threads.
    int num_threads = -1;

    // Method to find and select next best image to register.
    enum class ImageSelectionMethod {
      MAX_VISIBLE_POINTS_NUM,
      MAX_VISIBLE_POINTS_RATIO,
    };
    // Change the default selection method for the privacy preserving sfm.
    // Min uncertainty will likely not work with line features.
    ImageSelectionMethod image_selection_method =
        ImageSelectionMethod::MAX_VISIBLE_POINTS_RATIO;

    bool Check() const;
  };

  struct LocalBundleAdjustmentReport {
    size_t num_merged_observations = 0;
    size_t num_completed_observations = 0;
    size_t num_filtered_observations = 0;
    size_t num_adjusted_observations = 0;
  };

  // Create incremental mapper. The database cache must live for the entire
  // life-time of the incremental mapper.
  explicit IncrementalMapper(const DatabaseCache* database_cache);

  // Prepare the mapper for a new reconstruction, which might have existing
  // registered images (in which case `RegisterNextImage` must be called) or
  // which is empty (in which case `RegisterInitialImagePair` must be called).
  void BeginReconstruction(Reconstruction* reconstruction);

  // Cleanup the mapper after the current reconstruction is done. If the
  // model is discarded, the number of total and shared registered images will
  // be updated accordingly.
  void EndReconstruction(const bool discard);

  // Find best next image to register in the incremental reconstruction. The
  // images should be passed to `RegisterNextImage`. This function automatically
  // ignores images that failed to registered for `max_reg_trials`.
  std::vector<image_t> FindNextImages(const Options& options);

  bool RegisterInitialLineImages(const Options& options, const DatabaseCache& aligned_db_cache);

  // Attempt to register image to the existing model. This requires that
  // a previous call to `RegisterInitialImagePair` was successful.
  bool RegisterNextImage(const Options& options, const image_t image_id);

  // Triangulate observations of image.
  size_t TriangulateImage(const IncrementalTriangulator::Options& tri_options,
                          const image_t image_id);

  // Complete tracks by transitively following the scene graph correspondences.
  // This is especially effective after bundle adjustment, since many cameras
  // and point locations might have improved. Completion of tracks enables
  // better subsequent registration of new images.
  size_t CompleteTracks(const IncrementalTriangulator::Options& tri_options);

  // Merge tracks by using scene graph correspondences. Similar to
  // `CompleteTracks`, this is effective after bundle adjustment and improves
  // the redundancy in subsequent bundle adjustments.
  size_t MergeTracks(const IncrementalTriangulator::Options& tri_options);

  // Adjust locally connected images and points of a reference image. In
  // addition, refine the provided 3D points. Only images connected to the
  // reference image are optimized. If the provided 3D points are not locally
  // connected to the reference image, their observing images are set as
  // constant in the adjustment.
  LocalBundleAdjustmentReport AdjustLocalBundle(
      const Options& options, const BundleAdjustmentOptions& ba_options,
      const IncrementalTriangulator::Options& tri_options,
      const image_t image_id, const std::unordered_set<point3D_t>& point3D_ids);

  // Global bundle adjustment using Ceres Solver or PBA.
  bool AdjustGlobalBundle(const Options& options,
                          const BundleAdjustmentOptions& ba_options);

  // Filter images and point observations.
  size_t FilterImages(const Options& options);
  size_t FilterPoints(const Options& options);

  const Reconstruction& GetReconstruction() const;

  // Number of images that are registered in at least on reconstruction.
  size_t NumTotalRegImages() const;

  // Number of shared images between current reconstruction and all other
  // previous reconstructions.
  size_t NumSharedRegImages() const;

  // Get changed 3D points, since the last call to `ClearModifiedPoints3D`.
  const std::unordered_set<point3D_t>& GetModifiedPoints3D();

  // Clear the collection of changed 3D points.
  void ClearModifiedPoints3D();

 private:
  // Find local bundle for given image in the reconstruction. The local bundle
  // is defined as the images that are most connected, i.e. maximum number of
  // shared 3D points, to the given image.
  std::vector<image_t> FindLocalBundle(const Options& options,
                                       const image_t image_id) const;

  // Register / De-register image in current reconstruction and update
  // the number of shared images between all reconstructions.
  void RegisterImageEvent(const image_t image_id);
  void DeRegisterImageEvent(const image_t image_id);

  // Class that holds all necessary data from database in memory.
  const DatabaseCache* database_cache_;

  // Class that holds data of the reconstruction.
  Reconstruction* reconstruction_;

  // Class that is responsible for incremental triangulation.
  std::unique_ptr<IncrementalTriangulator> triangulator_;

  // Number of images that are registered in at least on reconstruction.
  size_t num_total_reg_images_;

  // Number of shared images between current reconstruction and all other
  // previous reconstructions.
  size_t num_shared_reg_images_;

  // Images and image pairs that have been used for initialization. Each image
  // and image pair is only tried once for initialization.
  std::unordered_map<image_t, size_t> init_num_reg_trials_;

  // The number of registered images per camera. This information is used
  // to avoid duplicate refinement of camera parameters and degradation of
  // already refined camera parameters in local bundle adjustment when multiple
  // images share intrinsics.
  std::unordered_map<camera_t, size_t> num_reg_images_per_camera_;

  // The number of reconstructions in which images are registered.
  std::unordered_map<image_t, size_t> num_registrations_;

  // Images that have been filtered in current reconstruction.
  std::unordered_set<image_t> filtered_images_;

  // Number of trials to register image in current reconstruction. Used to set
  // an upper bound to the number of trials to register an image.
  std::unordered_map<image_t, size_t> num_reg_trials_;

  // Images that were registered before beginning the reconstruction.
  // This image list will be non-empty, if the reconstruction is continued from
  // an existing reconstruction.
  std::unordered_set<image_t> existing_image_ids_;
};

}  // namespace colmap

#endif  // COLMAP_SRC_SFM_INCREMENTAL_MAPPER_H_
