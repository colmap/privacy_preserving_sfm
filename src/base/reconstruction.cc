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

#include "base/reconstruction.h"

#include <fstream>

#include "base/database_cache.h"
#include "base/pose.h"
#include "base/projection.h"
#include "base/triangulation.h"
#include "util/bitmap.h"
#include "util/misc.h"
#include "util/ply.h"

namespace colmap {

Reconstruction::Reconstruction()
    : correspondence_graph_(nullptr), num_added_points3D_(0) {}

std::unordered_set<point3D_t> Reconstruction::Point3DIds() const {
  std::unordered_set<point3D_t> point3D_ids;
  point3D_ids.reserve(points3D_.size());

  for (const auto& point3D : points3D_) {
    point3D_ids.insert(point3D.first);
  }

  return point3D_ids;
}

void Reconstruction::Load(const DatabaseCache& database_cache) {
  correspondence_graph_ = nullptr;

  // Add cameras.
  cameras_.reserve(database_cache.NumCameras());
  for (const auto& camera : database_cache.Cameras()) {
    if (!ExistsCamera(camera.first)) {
      AddCamera(camera.second);
    }
    // Else: camera was added before, e.g. with `ReadAllCameras`.
  }

  // Add images.
  images_.reserve(database_cache.NumImages());

  for (const auto& image : database_cache.Images()) {
    if (ExistsImage(image.second.ImageId())) {
      class Image& existing_image = Image(image.second.ImageId());
      CHECK_EQ(existing_image.Name(), image.second.Name());
      if (existing_image.NumLines() == 0) {
        existing_image.SetLines(image.second.Lines());
      } else {
        CHECK_EQ(image.second.NumLines(), existing_image.NumLines());
      }
      existing_image.SetNumObservations(image.second.NumObservations());
      existing_image.SetNumCorrespondences(image.second.NumCorrespondences());
    } else {
      AddImage(image.second);
    }
  }

  // Add image pairs.
  for (const auto& image_pair :
       database_cache.CorrespondenceGraph().NumCorrespondencesBetweenImages()) {
    ImagePairStat image_pair_stat;
    image_pair_stat.num_total_corrs = image_pair.second;
    image_pair_stats_.emplace(image_pair.first, image_pair_stat);
  }
}

void Reconstruction::SetUp(const CorrespondenceGraph* correspondence_graph) {
  CHECK_NOTNULL(correspondence_graph);
  correspondence_graph_ = correspondence_graph;

  // If an existing model was loaded from disk and there were already images
  // registered previously, we need to set observations as triangulated.
  for (const auto image_id : reg_image_ids_) {
    const class Image& image = Image(image_id);
    for (point2D_t line_idx = 0; line_idx < image.NumLines();
         ++line_idx) {
      if (image.Line(line_idx).HasPoint3D()) {
        const bool kIsContinuedPoint3D = false;
        SetObservationAsTriangulated(image_id, line_idx,
                                     kIsContinuedPoint3D);
      }
    }
  }
}

void Reconstruction::TearDown() {
  correspondence_graph_ = nullptr;

  // Remove all not yet registered images.
  std::unordered_set<camera_t> keep_camera_ids;
  for (auto it = images_.begin(); it != images_.end();) {
    if (it->second.IsRegistered()) {
      keep_camera_ids.insert(it->second.CameraId());
      ++it;
    } else {
      it = images_.erase(it);
    }
  }

  // Remove all unused cameras.
  for (auto it = cameras_.begin(); it != cameras_.end();) {
    if (keep_camera_ids.count(it->first) == 0) {
      it = cameras_.erase(it);
    } else {
      ++it;
    }
  }

  // Compress tracks.
  for (auto& point3D : points3D_) {
    point3D.second.Track().Compress();
  }
}

void Reconstruction::AddCamera(const class Camera& camera) {
  CHECK(!ExistsCamera(camera.CameraId()));
  CHECK(camera.VerifyParams());
  cameras_.emplace(camera.CameraId(), camera);
}

void Reconstruction::AddImage(const class Image& image) {
  CHECK(!ExistsImage(image.ImageId()));
  images_[image.ImageId()] = image;
}

point3D_t Reconstruction::AddPoint3D(const Eigen::Vector3d& xyz,
                                     const Track& track,
                                     const Eigen::Vector3ub& color) {
  const point3D_t point3D_id = ++num_added_points3D_;
  CHECK(!ExistsPoint3D(point3D_id));

  class Point3D& point3D = points3D_[point3D_id];

  point3D.SetXYZ(xyz);
  point3D.SetTrack(track);
  point3D.SetColor(color);

  for (const auto& track_el : track.Elements()) {
    class Image& image = Image(track_el.image_id);
    CHECK(!image.Line(track_el.line_idx).HasPoint3D());
    image.SetPoint3DForLine(track_el.line_idx, point3D_id);
    CHECK_LE(image.NumPoints3D(), image.NumLines());
  }

  const bool kIsContinuedPoint3D = false;

  for (const auto& track_el : track.Elements()) {
    SetObservationAsTriangulated(track_el.image_id, track_el.line_idx,
                                 kIsContinuedPoint3D);
  }

  return point3D_id;
}

void Reconstruction::AddObservation(const point3D_t point3D_id,
                                    const TrackElement& track_el) {
  class Image& image = Image(track_el.image_id);
  CHECK(!image.Line(track_el.line_idx).HasPoint3D());

  image.SetPoint3DForLine(track_el.line_idx, point3D_id);
  CHECK_LE(image.NumPoints3D(), image.NumLines());

  class Point3D& point3D = Point3D(point3D_id);
  point3D.Track().AddElement(track_el);

  const bool kIsContinuedPoint3D = true;
  SetObservationAsTriangulated(track_el.image_id, track_el.line_idx,
                               kIsContinuedPoint3D);
}

point3D_t Reconstruction::MergePoints3D(const point3D_t point3D_id1,
                                        const point3D_t point3D_id2) {
  const class Point3D& point3D1 = Point3D(point3D_id1);
  const class Point3D& point3D2 = Point3D(point3D_id2);

  const Eigen::Vector3d merged_xyz =
      (point3D1.Track().Length() * point3D1.XYZ() +
       point3D2.Track().Length() * point3D2.XYZ()) /
      (point3D1.Track().Length() + point3D2.Track().Length());
  const Eigen::Vector3d merged_rgb =
      (point3D1.Track().Length() * point3D1.Color().cast<double>() +
       point3D2.Track().Length() * point3D2.Color().cast<double>()) /
      (point3D1.Track().Length() + point3D2.Track().Length());

  Track merged_track;
  merged_track.Reserve(point3D1.Track().Length() + point3D2.Track().Length());
  merged_track.AddElements(point3D1.Track().Elements());
  merged_track.AddElements(point3D2.Track().Elements());

  DeletePoint3D(point3D_id1);
  DeletePoint3D(point3D_id2);

  const point3D_t merged_point3D_id =
      AddPoint3D(merged_xyz, merged_track, merged_rgb.cast<uint8_t>());

  return merged_point3D_id;
}

void Reconstruction::DeletePoint3D(const point3D_t point3D_id) {
  // Note: Do not change order of these instructions, especially with respect to
  // `Reconstruction::ResetTriObservations`

  const class Track& track = Point3D(point3D_id).Track();

  const bool kIsDeletedPoint3D = true;

  for (const auto& track_el : track.Elements()) {
    ResetTriObservations(track_el.image_id, track_el.line_idx,
                         kIsDeletedPoint3D);
  }

  for (const auto& track_el : track.Elements()) {
    class Image& image = Image(track_el.image_id);
    image.ResetPoint3DForLine(track_el.line_idx);
  }

  points3D_.erase(point3D_id);
}

void Reconstruction::DeleteObservation(const image_t image_id,
                                       const point2D_t line_idx) {
  // Note: Do not change order of these instructions, especially with respect to
  // `Reconstruction::ResetTriObservations`

  class Image& image = Image(image_id);
  const point3D_t point3D_id = image.Line(line_idx).Point3DId();
  class Point3D& point3D = Point3D(point3D_id);

  if (point3D.Track().Length() <= 3) {
    DeletePoint3D(point3D_id);
    return;
  }

  point3D.Track().DeleteElement(image_id, line_idx);

  const bool kIsDeletedPoint3D = false;
  ResetTriObservations(image_id, line_idx, kIsDeletedPoint3D);

  image.ResetPoint3DForLine(line_idx);
}

void Reconstruction::RegisterImage(const image_t image_id) {
  class Image& image = Image(image_id);
  if (!image.IsRegistered()) {
    image.SetRegistered(true);
    reg_image_ids_.push_back(image_id);
  }
}

void Reconstruction::DeRegisterImage(const image_t image_id) {
  class Image& image = Image(image_id);

  for (point2D_t line_idx = 0; line_idx < image.NumLines();
       ++line_idx) {
    if (image.Line(line_idx).HasPoint3D()) {
      DeleteObservation(image_id, line_idx);
    }
  }

  image.SetRegistered(false);

  reg_image_ids_.erase(
      std::remove(reg_image_ids_.begin(), reg_image_ids_.end(), image_id),
      reg_image_ids_.end());
}

void Reconstruction::Normalize(const double extent, const double p0,
                               const double p1, const bool use_images) {
  CHECK_GT(extent, 0);
  CHECK_GE(p0, 0);
  CHECK_LE(p0, 1);
  CHECK_GE(p1, 0);
  CHECK_LE(p1, 1);
  CHECK_LE(p0, p1);

  if ((use_images && reg_image_ids_.size() < 2) ||
      (!use_images && points3D_.size() < 2)) {
    return;
  }

  EIGEN_STL_UMAP(class Image*, Eigen::Vector3d) proj_centers;

  for (size_t i = 0; i < reg_image_ids_.size(); ++i) {
    class Image& image = Image(reg_image_ids_[i]);
    const Eigen::Vector3d proj_center = image.ProjectionCenter();
    proj_centers[&image] = proj_center;
  }

  // Coordinates of image centers or point locations.
  std::vector<float> coords_x;
  std::vector<float> coords_y;
  std::vector<float> coords_z;
  if (use_images) {
    coords_x.reserve(proj_centers.size());
    coords_y.reserve(proj_centers.size());
    coords_z.reserve(proj_centers.size());
    for (const auto& proj_center : proj_centers) {
      coords_x.push_back(static_cast<float>(proj_center.second(0)));
      coords_y.push_back(static_cast<float>(proj_center.second(1)));
      coords_z.push_back(static_cast<float>(proj_center.second(2)));
    }
  } else {
    coords_x.reserve(points3D_.size());
    coords_y.reserve(points3D_.size());
    coords_z.reserve(points3D_.size());
    for (const auto& point3D : points3D_) {
      coords_x.push_back(static_cast<float>(point3D.second.X()));
      coords_y.push_back(static_cast<float>(point3D.second.Y()));
      coords_z.push_back(static_cast<float>(point3D.second.Z()));
    }
  }

  // Determine robust bounding box and mean.

  std::sort(coords_x.begin(), coords_x.end());
  std::sort(coords_y.begin(), coords_y.end());
  std::sort(coords_z.begin(), coords_z.end());

  const size_t P0 = static_cast<size_t>(
      (coords_x.size() > 3) ? p0 * (coords_x.size() - 1) : 0);
  const size_t P1 = static_cast<size_t>(
      (coords_x.size() > 3) ? p1 * (coords_x.size() - 1) : coords_x.size() - 1);

  const Eigen::Vector3d bbox_min(coords_x[P0], coords_y[P0], coords_z[P0]);
  const Eigen::Vector3d bbox_max(coords_x[P1], coords_y[P1], coords_z[P1]);

  Eigen::Vector3d mean_coord(0, 0, 0);
  for (size_t i = P0; i <= P1; ++i) {
    mean_coord(0) += coords_x[i];
    mean_coord(1) += coords_y[i];
    mean_coord(2) += coords_z[i];
  }
  mean_coord /= P1 - P0 + 1;

  // Calculate scale and translation, such that
  // translation is applied before scaling.
  const double old_extent = (bbox_max - bbox_min).norm();
  double scale;
  if (old_extent < std::numeric_limits<double>::epsilon()) {
    scale = 1;
  } else {
    scale = extent / old_extent;
  }

  const Eigen::Vector3d translation = mean_coord;

  // Transform images.
  for (auto& image_proj_center : proj_centers) {
    image_proj_center.second -= translation;
    image_proj_center.second *= scale;
    const Eigen::Quaterniond quat(
        image_proj_center.first->Qvec(0), image_proj_center.first->Qvec(1),
        image_proj_center.first->Qvec(2), image_proj_center.first->Qvec(3));
    image_proj_center.first->SetTvec(quat * -image_proj_center.second);
  }

  // Transform points.
  for (auto& point3D : points3D_) {
    point3D.second.XYZ() -= translation;
    point3D.second.XYZ() *= scale;
  }
}

size_t Reconstruction::FilterPoints3D(
    const double max_reproj_error, const double min_tri_angle,
    const std::unordered_set<point3D_t>& point3D_ids) {
    
  size_t num_filtered = 0;
  num_filtered +=
      FilterPoints3DWithLargeReprojectionError(max_reproj_error, point3D_ids);
  
  num_filtered +=
      FilterPoints3DWithSmallTriangulationAngle(min_tri_angle, point3D_ids);
  return num_filtered;
}

size_t Reconstruction::FilterPoints3DInImages(
    const double max_reproj_error, const double min_tri_angle,
    const std::unordered_set<image_t>& image_ids) {

  std::unordered_set<point3D_t> point3D_ids;
  for (const image_t image_id : image_ids) {
    const class Image& image = Image(image_id);
    for (const FeatureLine& line : image.Lines()) {
      if (line.HasPoint3D()) {
        point3D_ids.insert(line.Point3DId());
      }
    }
  }
  return FilterPoints3D(max_reproj_error, min_tri_angle, point3D_ids);
}

size_t Reconstruction::FilterAllPoints3D(const double max_reproj_error,
                                         const double min_tri_angle) {
  // Important: First filter observations and points with large reprojection
  // error, so that observations with large reprojection error do not make
  // a point stable through a large triangulation angle.
  const std::unordered_set<point3D_t>& point3D_ids = Point3DIds();
  size_t num_filtered = 0;
  num_filtered +=
      FilterPoints3DWithLargeReprojectionError(max_reproj_error, point3D_ids);
  num_filtered +=
      FilterPoints3DWithSmallTriangulationAngle(min_tri_angle, point3D_ids);
  return num_filtered;
}

size_t Reconstruction::FilterObservationsWithNegativeDepth() {
  size_t num_filtered = 0;
  for (const auto image_id : reg_image_ids_) {
    const class Image& image = Image(image_id);
    const Eigen::Matrix3x4d proj_matrix = image.ProjectionMatrix();
    for (point2D_t line_idx = 0; line_idx < image.NumLines();
         ++line_idx) {
      const FeatureLine& line = image.Line(line_idx);
      if (line.HasPoint3D()) {
        const class Point3D& point3D = Point3D(line.Point3DId());
        if (!HasPointPositiveDepth(proj_matrix, point3D.XYZ())) {
          DeleteObservation(image_id, line_idx);
          num_filtered += 1;
        }
      }
    }
  }
  return num_filtered;
}

std::vector<image_t> Reconstruction::FilterImages(
    const double min_focal_length_ratio, const double max_focal_length_ratio,
    const double max_extra_param) {
  std::vector<image_t> filtered_image_ids;
  for (const image_t image_id : RegImageIds()) {
    const class Image& image = Image(image_id);
    const class Camera& camera = Camera(image.CameraId());
    if (image.NumPoints3D() == 0) {
      filtered_image_ids.push_back(image_id);
    } else if (camera.HasBogusParams(min_focal_length_ratio,
                                     max_focal_length_ratio, max_extra_param)) {
      filtered_image_ids.push_back(image_id);
    }
  }

  // Only de-register after iterating over reg_image_ids_ to avoid
  // simultaneous iteration and modification of the vector.
  for (const image_t image_id : filtered_image_ids) {
    DeRegisterImage(image_id);
  }

  return filtered_image_ids;
}

size_t Reconstruction::ComputeNumObservations() const {
  size_t num_obs = 0;
  for (const image_t image_id : reg_image_ids_) {
    num_obs += Image(image_id).NumPoints3D();
  }
  return num_obs;
}

double Reconstruction::ComputeMeanTrackLength() const {
  if (points3D_.empty()) {
    return 0.0;
  } else {
    return ComputeNumObservations() / static_cast<double>(points3D_.size());
  }
}

double Reconstruction::ComputeMeanObservationsPerRegImage() const {
  if (reg_image_ids_.empty()) {
    return 0.0;
  } else {
    return ComputeNumObservations() /
           static_cast<double>(reg_image_ids_.size());
  }
}

double Reconstruction::ComputeMeanReprojectionError() const {
  double error_sum = 0.0;
  size_t num_valid_errors = 0;
  for (const auto& point3D : points3D_) {
    if (point3D.second.HasError()) {
      error_sum += point3D.second.Error();
      num_valid_errors += 1;
    }
  }

  if (num_valid_errors == 0) {
    return 0.0;
  } else {
    return error_sum / num_valid_errors;
  }
}

void Reconstruction::Read(const std::string& path) {

  if (ExistsFile(JoinPaths(path, "cameras.txt")) &&
             ExistsFile(JoinPaths(path, "images.txt")) &&
             ExistsFile(JoinPaths(path, "points3D.txt"))) {
    ReadText(path);
  } else {
    LOG(FATAL) << "cameras, images, points3D files do not exist at " << path;
  }
}

void Reconstruction::Write(const std::string& path) const {
  WriteText(path);
}

void Reconstruction::ReadText(const std::string& path) {
  ReadCamerasText(JoinPaths(path, "cameras.txt"));
  ReadImagesText(JoinPaths(path, "images.txt"));
  ReadPoints3DText(JoinPaths(path, "points3D.txt"));
}

void Reconstruction::WriteText(const std::string& path) const {
  WriteCamerasText(JoinPaths(path, "cameras.txt"));
  WriteImagesText(JoinPaths(path, "images.txt"));
  WritePoints3DText(JoinPaths(path, "points3D.txt"));
}

std::vector<PlyPoint> Reconstruction::ConvertToPLY() const {
  std::vector<PlyPoint> ply_points;
  ply_points.reserve(points3D_.size());

  for (const auto& point3D : points3D_) {
    PlyPoint ply_point;
    ply_point.x = point3D.second.X();
    ply_point.y = point3D.second.Y();
    ply_point.z = point3D.second.Z();
    ply_point.r = point3D.second.Color(0);
    ply_point.g = point3D.second.Color(1);
    ply_point.b = point3D.second.Color(2);
    ply_points.push_back(ply_point);
  }

  return ply_points;
}

void Reconstruction::ImportPLY(const std::string& path) {
  points3D_.clear();

  const auto ply_points = ReadPly(path);

  points3D_.reserve(ply_points.size());

  for (const auto& ply_point : ply_points) {
    AddPoint3D(Eigen::Vector3d(ply_point.x, ply_point.y, ply_point.z), Track(),
               Eigen::Vector3ub(ply_point.r, ply_point.g, ply_point.b));
  }
}

void Reconstruction::ExportPLY(const std::string& path) const {
  const auto ply_points = ConvertToPLY();

  const bool kWriteNormal = false;
  const bool kWriteRGB = true;
  WriteBinaryPlyPoints(path, ply_points, kWriteNormal, kWriteRGB);
}

size_t Reconstruction::FilterPoints3DWithSmallTriangulationAngle(
    const double min_tri_angle,
    const std::unordered_set<point3D_t>& point3D_ids) {
  // Number of filtered points.
  size_t num_filtered = 0;

  // Minimum triangulation angle in radians.
  const double min_tri_angle_rad = DegToRad(min_tri_angle);

  // Cache for image projection centers.
  EIGEN_STL_UMAP(image_t, Eigen::Vector3d) proj_centers;

  for (const auto point3D_id : point3D_ids) {
    if (!ExistsPoint3D(point3D_id)) {
      continue;
    }

    const class Point3D& point3D = Point3D(point3D_id);

    // Calculate triangulation angle for all pairwise combinations of image
    // poses in the track. Only delete point if none of the combinations
    // has a sufficient triangulation angle.
    bool keep_point = false;
    for (size_t i1 = 0; i1 < point3D.Track().Length(); ++i1) {
      const image_t image_id1 = point3D.Track().Element(i1).image_id;

      Eigen::Vector3d proj_center1;
      if (proj_centers.count(image_id1) == 0) {
        const class Image& image1 = Image(image_id1);
        proj_center1 = image1.ProjectionCenter();
        proj_centers.emplace(image_id1, proj_center1);
      } else {
        proj_center1 = proj_centers.at(image_id1);
      }

      for (size_t i2 = 0; i2 < i1; ++i2) {
        const image_t image_id2 = point3D.Track().Element(i2).image_id;
        const Eigen::Vector3d proj_center2 = proj_centers.at(image_id2);

        const double tri_angle = CalculateTriangulationAngle(
            proj_center1, proj_center2, point3D.XYZ());

        if (tri_angle >= min_tri_angle_rad) {
          keep_point = true;
          break;
        }
      }

      if (keep_point) {
        break;
      }
    }

    if (!keep_point) {
      num_filtered += 1;
      DeletePoint3D(point3D_id);
    }
  }

  return num_filtered;
}

size_t Reconstruction::FilterPoints3DWithLargeReprojectionError(
    const double max_reproj_error,
    const std::unordered_set<point3D_t>& point3D_ids) {
  const double max_squared_reproj_error = max_reproj_error * max_reproj_error;

  // Number of filtered points.
  size_t num_filtered = 0;

  for (const auto point3D_id : point3D_ids) {
    if (!ExistsPoint3D(point3D_id)) {
      continue;
    }

    class Point3D& point3D = Point3D(point3D_id);

    bool have_non_aligned = false;
    for (const auto& track_el : point3D.Track().Elements()) {
        if(!Image(track_el.image_id).Line(track_el.line_idx).IsAligned())
          have_non_aligned = true;
    }
    if(!have_non_aligned) {
      DeletePoint3D(point3D_id);
      num_filtered += point3D.Track().Length();
      continue;
    }

    if (point3D.Track().Length() < 3) {
      DeletePoint3D(point3D_id);
      num_filtered += point3D.Track().Length();
      continue;
    }

    double reproj_error_sum = 0.0;

    std::vector<TrackElement> track_els_to_delete;

    for (const auto& track_el : point3D.Track().Elements()) {
      const class Image& image = Image(track_el.image_id);
      const class Camera& camera = Camera(image.CameraId());
      const FeatureLine& line2D = image.Line(track_el.line_idx);
      CHECK_NEAR(line2D.Line().head<2>().norm(), 1.0, 1e-6);
      const double squared_reproj_error = CalculateSquaredLineReprojectionError(
          line2D.Line(), point3D.XYZ(), image.Qvec(), image.Tvec(), camera);
      if (squared_reproj_error > max_squared_reproj_error) {
        track_els_to_delete.push_back(track_el);
      } else {
        reproj_error_sum += std::sqrt(squared_reproj_error);
      }
    }

    if (track_els_to_delete.size() >= point3D.Track().Length() - 3) {
      num_filtered += point3D.Track().Length();
      DeletePoint3D(point3D_id);
    } else {
      num_filtered += track_els_to_delete.size();
      for (const auto& track_el : track_els_to_delete) {
        DeleteObservation(track_el.image_id, track_el.line_idx);
      }
      point3D.SetError(reproj_error_sum / point3D.Track().Length());
    }
  }

  return num_filtered;
}

void Reconstruction::ReadCamerasText(const std::string& path) {
  cameras_.clear();

  std::ifstream file(path);
  CHECK(file.is_open()) << path;

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    StringTrim(&line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream(line);

    class Camera camera;

    // ID
    std::getline(line_stream, item, ' ');
    camera.SetCameraId(std::stoul(item));

    // MODEL
    std::getline(line_stream, item, ' ');
    camera.SetModelIdFromName(item);

    // WIDTH
    std::getline(line_stream, item, ' ');
    camera.SetWidth(std::stoll(item));

    // HEIGHT
    std::getline(line_stream, item, ' ');
    camera.SetHeight(std::stoll(item));

    // PARAMS
    camera.Params().clear();
    while (!line_stream.eof()) {
      std::getline(line_stream, item, ' ');
      camera.Params().push_back(std::stold(item));
    }

    CHECK(camera.VerifyParams());

    cameras_.emplace(camera.CameraId(), camera);
  }
}

void Reconstruction::ReadImagesText(const std::string& path) {
  images_.clear();

  std::ifstream file(path);
  CHECK(file.is_open()) << path;

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    StringTrim(&line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream1(line);

    // ID
    std::getline(line_stream1, item, ' ');
    const image_t image_id = std::stoul(item);

    class Image image;
    image.SetImageId(image_id);

    image.SetRegistered(true);
    reg_image_ids_.push_back(image_id);

    // QVEC (qw, qx, qy, qz)
    std::getline(line_stream1, item, ' ');
    image.Qvec(0) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(1) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(2) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(3) = std::stold(item);

    image.NormalizeQvec();

    // TVEC
    std::getline(line_stream1, item, ' ');
    image.Tvec(0) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Tvec(1) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Tvec(2) = std::stold(item);

    // CAMERA_ID
    std::getline(line_stream1, item, ' ');
    image.SetCameraId(std::stoul(item));

    // NAME
    std::getline(line_stream1, item, ' ');
    image.SetName(item);

    // LINES2D
    if (!std::getline(file, line)) {
        break;
    }

    FeatureLines feature_lines;

    StringTrim(&line);
    std::stringstream line_stream2(line);

    if (!line.empty()) {
      while (!line_stream2.eof()) {

          FeatureLine feature_line;

          Eigen::Vector3d line_dir;

          std::getline(line_stream2, item, ' ');

          line_dir(0) = std::stof(item);

          std::getline(line_stream2, item, ' ');
          line_dir(1) = std::stof(item);

          std::getline(line_stream2, item, ' ');
          line_dir(2) = std::stof(item);

          std::getline(line_stream2, item, ' ');
          if (item == "1") {
              feature_line.SetAligned(true);
          } else {
              CHECK(item == "0");
              feature_line.SetAligned(false);
          }

          std::getline(line_stream2, item, ' ');
          if (item == "-1") {
              feature_line.SetPoint3DId(kInvalidPoint3DId);
          } else {
              feature_line.SetPoint3DId(std::stoll(item));
          }

          // Normalize line to simplify distance computations
          const double normalization_factor = line_dir.head<2>().norm();

          feature_line.SetLine(line_dir / normalization_factor);

          feature_lines.push_back(feature_line);
      }

      image.SetLines(feature_lines);
    }

    images_.emplace(image.ImageId(), image);
  }
}

void Reconstruction::ReadPoints3DText(const std::string& path) {
  points3D_.clear();

  std::ifstream file(path);
  CHECK(file.is_open()) << path;

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    StringTrim(&line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream(line);

    // ID
    std::getline(line_stream, item, ' ');
    const point3D_t point3D_id = std::stoll(item);

    // Make sure, that we can add new 3D points after reading 3D points
    // without overwriting existing 3D points.
    num_added_points3D_ = std::max(num_added_points3D_, point3D_id);

    class Point3D point3D;

    // XYZ
    std::getline(line_stream, item, ' ');
    point3D.XYZ(0) = std::stold(item);

    std::getline(line_stream, item, ' ');
    point3D.XYZ(1) = std::stold(item);

    std::getline(line_stream, item, ' ');
    point3D.XYZ(2) = std::stold(item);

    // Color
    std::getline(line_stream, item, ' ');
    point3D.Color(0) = static_cast<uint8_t>(std::stoi(item));

    std::getline(line_stream, item, ' ');
    point3D.Color(1) = static_cast<uint8_t>(std::stoi(item));

    std::getline(line_stream, item, ' ');
    point3D.Color(2) = static_cast<uint8_t>(std::stoi(item));

    // ERROR
    std::getline(line_stream, item, ' ');
    point3D.SetError(std::stold(item));

    // TRACK
    while (!line_stream.eof()) {
      TrackElement track_el;

      std::getline(line_stream, item, ' ');
      StringTrim(&item);
      if (item.empty()) {
        break;
      }
      track_el.image_id = std::stoul(item);

      std::getline(line_stream, item, ' ');
      track_el.line_idx = std::stoul(item);

      point3D.Track().AddElement(track_el);
    }

    point3D.Track().Compress();

    points3D_.emplace(point3D_id, point3D);
  }
}

void Reconstruction::WriteCamerasText(const std::string& path) const {
  std::ofstream file(path, std::ios::trunc);
  CHECK(file.is_open()) << path;

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# Camera list with one line of data per camera:" << std::endl;
  file << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]" << std::endl;
  file << "# Number of cameras: " << cameras_.size() << std::endl;

  for (const auto& camera : cameras_) {
    std::ostringstream line;

    line << camera.first << " ";
    line << camera.second.ModelName() << " ";
    line << camera.second.Width() << " ";
    line << camera.second.Height() << " ";

    for (const double param : camera.second.Params()) {
      line << param << " ";
    }

    std::string line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);

    file << line_string << std::endl;
  }
}

void Reconstruction::WriteImagesText(const std::string& path) const {
  std::ofstream file(path, std::ios::trunc);
  CHECK(file.is_open()) << path;

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# Image list with two lines of data per image:" << std::endl;
  file << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, "
          "NAME"
       << std::endl;
  file << "#   LINES2D[] as (A, B, C, is_aligned, POINT3D_ID)" << std::endl;
  file << "# Number of images: " << reg_image_ids_.size()
       << ", mean observations per image: "
       << ComputeMeanObservationsPerRegImage() << std::endl;

  for (const auto& image : images_) {
    if (!image.second.IsRegistered()) {
      continue;
    }

    std::ostringstream line;
    std::string line_string;

    line << image.first << " ";

    // QVEC (qw, qx, qy, qz)
    const Eigen::Vector4d normalized_qvec =
        NormalizeQuaternion(image.second.Qvec());
    line << normalized_qvec(0) << " ";
    line << normalized_qvec(1) << " ";
    line << normalized_qvec(2) << " ";
    line << normalized_qvec(3) << " ";

    // TVEC
    line << image.second.Tvec(0) << " ";
    line << image.second.Tvec(1) << " ";
    line << image.second.Tvec(2) << " ";

    line << image.second.CameraId() << " ";

    line << image.second.Name();

    file << line.str() << std::endl;

    line.str("");
    line.clear();

    for (const FeatureLine& feature_line : image.second.Lines()) {
      const Eigen::Vector3d& line_dir = feature_line.Line();
      line << line_dir(0) << " ";
      line << line_dir(1) << " ";
      line << line_dir(2) << " ";
      line << (feature_line.IsAligned() ? "1" : "0") << " ";
      if (feature_line.HasPoint3D()) {
        line << feature_line.Point3DId() << " ";
      } else {
        line << -1 << " ";
      }
    }
    line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);
    file << line_string << std::endl;
  }
}

void Reconstruction::WritePoints3DText(const std::string& path) const {
  std::ofstream file(path, std::ios::trunc);
  CHECK(file.is_open()) << path;

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# 3D point list with one line of data per point:" << std::endl;
  file << "#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, "
          "TRACK[] as (IMAGE_ID, line_idx)"
       << std::endl;
  file << "# Number of points: " << points3D_.size()
       << ", mean track length: " << ComputeMeanTrackLength() << std::endl;

  for (const auto& point3D : points3D_) {
    file << point3D.first << " ";
    file << point3D.second.XYZ()(0) << " ";
    file << point3D.second.XYZ()(1) << " ";
    file << point3D.second.XYZ()(2) << " ";
    file << static_cast<int>(point3D.second.Color(0)) << " ";
    file << static_cast<int>(point3D.second.Color(1)) << " ";
    file << static_cast<int>(point3D.second.Color(2)) << " ";
    file << point3D.second.Error() << " ";

    std::ostringstream line;

    for (const auto& track_el : point3D.second.Track().Elements()) {
      line << track_el.image_id << " ";
      line << track_el.line_idx << " ";
    }

    std::string line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);

    file << line_string << std::endl;
  }
}

void Reconstruction::SetObservationAsTriangulated(
    const image_t image_id, const point2D_t line_idx,
    const bool is_continued_point3D) {
  if (correspondence_graph_ == nullptr) {
    return;
  }

  const class Image& image = Image(image_id);
  const FeatureLine& line = image.Line(line_idx);
  const std::vector<CorrespondenceGraph::Correspondence>& corrs =
      correspondence_graph_->FindCorrespondences(image_id, line_idx);

  CHECK(image.IsRegistered());
  CHECK(line.HasPoint3D());

  for (const auto& corr : corrs) {
    class Image& corr_image = Image(corr.image_id);
    const FeatureLine& corr_line = corr_image.Line(corr.line_idx);
    corr_image.IncrementCorrespondenceHasPoint3D(corr.line_idx);
    // Update number of shared 3D points between image pairs and make sure to
    // only count the correspondences once (not twice forward and backward).
    if (line.Point3DId() == corr_line.Point3DId() &&
        (is_continued_point3D || image_id < corr.image_id)) {
      const image_pair_t pair_id =
          Database::ImagePairToPairId(image_id, corr.image_id);
      image_pair_stats_[pair_id].num_tri_corrs += 1;
      CHECK_LE(image_pair_stats_[pair_id].num_tri_corrs,
               image_pair_stats_[pair_id].num_total_corrs)
          << "The correspondence graph graph must not contain duplicate "
             "matches";
    }
  }
}

void Reconstruction::ResetTriObservations(const image_t image_id,
                                          const point2D_t line_idx,
                                          const bool is_deleted_point3D) {
  if (correspondence_graph_ == nullptr) {
    return;
  }

  const class Image& image = Image(image_id);
  const FeatureLine& line = image.Line(line_idx);
  const std::vector<CorrespondenceGraph::Correspondence>& corrs =
      correspondence_graph_->FindCorrespondences(image_id, line_idx);

  CHECK(image.IsRegistered());
  CHECK(line.HasPoint3D());

  for (const auto& corr : corrs) {
    class Image& corr_image = Image(corr.image_id);
    const FeatureLine& corr_line = corr_image.Line(corr.line_idx);
    corr_image.DecrementCorrespondenceHasPoint3D(corr.line_idx);
    // Update number of shared 3D points between image pairs and make sure to
    // only count the correspondences once (not twice forward and backward).
    if (line.Point3DId() == corr_line.Point3DId() &&
        (!is_deleted_point3D || image_id < corr.image_id)) {
      const image_pair_t pair_id =
          Database::ImagePairToPairId(image_id, corr.image_id);
      image_pair_stats_[pair_id].num_tri_corrs -= 1;
      CHECK_GE(image_pair_stats_[pair_id].num_tri_corrs, 0)
          << "The scene graph graph must not contain duplicate matches";
    }
  }
}

}  // namespace colmap
