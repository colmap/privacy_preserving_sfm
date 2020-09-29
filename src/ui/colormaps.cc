// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
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
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "ui/colormaps.h"

#include "util/bitmap.h"

namespace colmap {

PointColormapBase::PointColormapBase()
    : scale(1.0f),
      min(0.0f),
      max(0.0f),
      range(0.0f),
      min_q(0.0f),
      max_q(1.0f) {}

void PointColormapBase::UpdateScale(std::vector<float>* values) {
  if (values->empty()) {
    min = 0.0f;
    max = 0.0f;
    range = 0.0f;
  } else {
    std::sort(values->begin(), values->end());
    min = (*values)[static_cast<size_t>(min_q * (values->size() - 1))];
    max = (*values)[static_cast<size_t>(max_q * (values->size() - 1))];
    range = max - min;
  }
}

float PointColormapBase::AdjustScale(const float gray) {
  if (range == 0.0f) {
    return 0.0f;
  } else {
    const float gray_clipped = std::min(std::max(gray, min), max);
    const float gray_scaled = (gray_clipped - min) / range;
    return std::pow(gray_scaled, scale);
  }
}

void PointColormapPhotometric::Prepare(EIGEN_STL_UMAP(camera_t, Camera) &
                                           cameras,
                                       EIGEN_STL_UMAP(image_t, Image) & images,
                                       EIGEN_STL_UMAP(point3D_t, Point3D) &
                                           points3D,
                                       std::vector<image_t>& reg_image_ids) {}

Eigen::Vector4f PointColormapPhotometric::ComputeColor(
    const point3D_t point3D_id, const Point3D& point3D) {
  return Eigen::Vector4f(point3D.Color(0) / 255.0f, point3D.Color(1) / 255.0f,
                         point3D.Color(2) / 255.0f, 1.0f);
}

void PointColormapError::Prepare(EIGEN_STL_UMAP(camera_t, Camera) & cameras,
                                 EIGEN_STL_UMAP(image_t, Image) & images,
                                 EIGEN_STL_UMAP(point3D_t, Point3D) & points3D,
                                 std::vector<image_t>& reg_image_ids) {
  std::vector<float> errors;
  errors.reserve(points3D.size());

  for (const auto& point3D : points3D) {
    errors.push_back(static_cast<float>(point3D.second.Error()));
  }

  UpdateScale(&errors);
}

Eigen::Vector4f PointColormapError::ComputeColor(const point3D_t point3D_id,
                                                 const Point3D& point3D) {
  const float gray = AdjustScale(static_cast<float>(point3D.Error()));
  return Eigen::Vector4f(JetColormap::Red(gray), JetColormap::Green(gray),
                         JetColormap::Blue(gray), 1.0f);
}

void PointColormapTrackLen::Prepare(EIGEN_STL_UMAP(camera_t, Camera) & cameras,
                                    EIGEN_STL_UMAP(image_t, Image) & images,
                                    EIGEN_STL_UMAP(point3D_t, Point3D) &
                                        points3D,
                                    std::vector<image_t>& reg_image_ids) {
  std::vector<float> track_lengths;
  track_lengths.reserve(points3D.size());

  for (const auto& point3D : points3D) {
    track_lengths.push_back(point3D.second.Track().Length());
  }

  UpdateScale(&track_lengths);
}

Eigen::Vector4f PointColormapTrackLen::ComputeColor(const point3D_t point3D_id,
                                                    const Point3D& point3D) {
  const float gray = AdjustScale(point3D.Track().Length());
  return Eigen::Vector4f(JetColormap::Red(gray), JetColormap::Green(gray),
                         JetColormap::Blue(gray), 1.0f);
}

const Eigen::Vector4f ImageColormapBase::kDefaultPlaneColor = {1.0f, 0.1f, 0.0f,
                                                               0.6f};
const Eigen::Vector4f ImageColormapBase::kDefaultFrameColor = {0.8f, 0.1f, 0.0f,
                                                               1.0f};

ImageColormapBase::ImageColormapBase() {}

void ImageColormapUniform::Prepare(EIGEN_STL_UMAP(camera_t, Camera) & cameras,
                                   EIGEN_STL_UMAP(image_t, Image) & images,
                                   EIGEN_STL_UMAP(point3D_t, Point3D) &
                                       points3D,
                                   std::vector<image_t>& reg_image_ids) {}

void ImageColormapUniform::ComputeColor(const Image& image,
                                        Eigen::Vector4f* plane_color,
                                        Eigen::Vector4f* frame_color) {
  *plane_color = uniform_plane_color;
  *frame_color = uniform_frame_color;
}

void ImageColormapNameFilter::Prepare(EIGEN_STL_UMAP(camera_t, Camera) &
                                          cameras,
                                      EIGEN_STL_UMAP(image_t, Image) & images,
                                      EIGEN_STL_UMAP(point3D_t, Point3D) &
                                          points3D,
                                      std::vector<image_t>& reg_image_ids) {}

void ImageColormapNameFilter::AddColorForWord(
    const std::string& word, const Eigen::Vector4f& plane_color,
    const Eigen::Vector4f& frame_color) {
  image_name_colors_.emplace_back(word,
                                  std::make_pair(plane_color, frame_color));
}

void ImageColormapNameFilter::ComputeColor(const Image& image,
                                           Eigen::Vector4f* plane_color,
                                           Eigen::Vector4f* frame_color) {
  for (const auto& image_name_color : image_name_colors_) {
    if (StringContains(image.Name(), image_name_color.first)) {
      *plane_color = image_name_color.second.first;
      *frame_color = image_name_color.second.second;
      return;
    }
  }

  *plane_color = kDefaultPlaneColor;
  *frame_color = kDefaultFrameColor;
}

}  // namespace colmap
