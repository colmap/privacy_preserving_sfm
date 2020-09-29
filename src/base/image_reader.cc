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

#include "base/image_reader.h"

#include "base/camera_models.h"
#include "util/misc.h"

namespace colmap {

namespace {
bool ReadGravityFile(const std::string& gravity_file_path,
                     Eigen::Vector3d* gravity_direction) {
  std::ifstream gravity_stream(gravity_file_path);
  if (gravity_stream) {
    Eigen::Vector3d gravity;
    gravity_stream >> gravity(0) >> gravity(1) >> gravity(2);
    *gravity_direction = gravity.normalized();

    return true;
  } else {
    return false;
  }
}

bool ReadCameraModelFile(const std::string& camera_model_file_path,
                         Camera* camera) {
  std::ifstream camera_model_stream(camera_model_file_path);
  std::string camera_model_string;
  Camera tmp_cam;

  if (std::getline(camera_model_stream, camera_model_string)) {
    std::stringstream model_string_stream(camera_model_string);
    std::string camera_model;
    model_string_stream >> camera_model;

    // Get the remaining string after the first element.
    std::string model_params;
    std::getline(model_string_stream, model_params);
    tmp_cam.SetModelIdFromName(camera_model);
    tmp_cam.SetParamsFromString(model_params);

    if (tmp_cam.VerifyParams()) {
      *camera = tmp_cam;
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

std::vector<std::string> FilterNonImageFiles(const std::vector<std::string>& file_paths) {

  std::vector<std::string> image_file_paths;
  for (const std::string& path : file_paths) {
    if (Bitmap::IsImageFileType(path)) {
      image_file_paths.emplace_back(path);
    }
  }

  return image_file_paths;
}
} // namespace

bool ImageReaderOptions::Check() const {
  CHECK_OPTION_GT(default_focal_length_factor, 0.0);
  CHECK_OPTION(ExistsCameraModelWithName(camera_model));
  const int model_id = CameraModelNameToId(camera_model);
  if (!camera_params.empty()) {
    CHECK_OPTION(
        CameraModelVerifyParams(model_id, CSVToVector<double>(camera_params)));
  }
  return true;
}

ImageReader::ImageReader(const ImageReaderOptions& options, Database* database)
    : options_(options), database_(database), image_index_(0) {
  CHECK(options_.Check());

  // Ensure trailing slash, so that we can build the correct image name.
  options_.image_path =
      EnsureTrailingSlash(StringReplace(options_.image_path, "\\", "/"));

  // Get a list of all files in the image path, sorted by image name.
  if (options_.image_list.empty()) {
    options_.image_list = FilterNonImageFiles(GetRecursiveFileList(options_.image_path));
    std::sort(options_.image_list.begin(), options_.image_list.end());
  } else {
    for (auto& image_name : options_.image_list) {
      image_name = JoinPaths(options_.image_path, image_name);
    }
  }

  if (static_cast<camera_t>(options_.existing_camera_id) != kInvalidCameraId) {
    CHECK(database->ExistsCamera(options_.existing_camera_id));
    prev_camera_ = database->ReadCamera(options_.existing_camera_id);
  } else {
    // Set the manually specified camera parameters.
    prev_camera_.SetCameraId(kInvalidCameraId);
    prev_camera_.SetModelIdFromName(options_.camera_model);
    if (!options_.camera_params.empty()) {
      CHECK(prev_camera_.SetParamsFromString(options_.camera_params));
      prev_camera_.SetPriorFocalLength(true);
      default_camera_ = prev_camera_;
    }
  }
}

ImageReader::Status ImageReader::Next(Camera* camera, Image* image,
                                      Bitmap* bitmap, Bitmap* mask) {
  CHECK_NOTNULL(camera);
  CHECK_NOTNULL(image);
  CHECK_NOTNULL(bitmap);

  image_index_ += 1;
  CHECK_LE(image_index_, options_.image_list.size());

  const std::string image_path = options_.image_list.at(image_index_ - 1);

  DatabaseTransaction database_transaction(database_);

  //////////////////////////////////////////////////////////////////////////////
  // Set the image name.
  //////////////////////////////////////////////////////////////////////////////

  image->SetName(image_path);
  image->SetName(StringReplace(image->Name(), "\\", "/"));
  image->SetName(
      image->Name().substr(options_.image_path.size(),
                           image->Name().size() - options_.image_path.size()));

  const std::string image_folder = GetParentDir(image->Name());

  //////////////////////////////////////////////////////////////////////////////
  // Check if image already read.
  //////////////////////////////////////////////////////////////////////////////

  const bool exists_image = database_->ExistsImageWithName(image->Name());

  if (exists_image) {
    *image = database_->ReadImageWithName(image->Name());
    const bool exists_lines = database_->ExistsLineFeatures(image->ImageId());
    const bool exists_descriptors =
        database_->ExistsDescriptors(image->ImageId());

    if (exists_lines && exists_descriptors) {
      return Status::IMAGE_EXISTS;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Read image.
  //////////////////////////////////////////////////////////////////////////////

  if (!bitmap->Read(image_path, false)) {
    return Status::BITMAP_ERROR;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Read mask.
  //////////////////////////////////////////////////////////////////////////////

  if (mask && !options_.mask_path.empty()) {
    const std::string mask_path =
        JoinPaths(options_.mask_path,
                  GetRelativePath(options_.image_path, image_path) + ".png");
    if (ExistsFile(mask_path) && !mask->Read(mask_path, false)) {
      // NOTE: Maybe introduce a separate error type MASK_ERROR?
      return Status::BITMAP_ERROR;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Try to read gravity direction
  //////////////////////////////////////////////////////////////////////////////

  {
    const std::string gravity_file_path = image_path + ".gravity.txt";
    Eigen::Vector3d gravity_direction;
    if (!ReadGravityFile(gravity_file_path, &gravity_direction)) {
      gravity_direction.setConstant(std::numeric_limits<double>::quiet_NaN());
    }
    image->SetGravityDirection(gravity_direction);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Check for well-formed data.
  //////////////////////////////////////////////////////////////////////////////

  if (exists_image) {
    const Camera camera = database_->ReadCamera(image->CameraId());

    if (static_cast<size_t>(bitmap->Width()) != camera.Width() ||
        static_cast<size_t>(bitmap->Height()) != camera.Height()) {
      return Status::CAMERA_EXIST_DIM_ERROR;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Read camera model and parameters
  //////////////////////////////////////////////////////////////////////////////

  {
    const std::string camera_model_file_path = image_path + ".camera_model.txt";
    if (!ReadCameraModelFile(camera_model_file_path, &prev_camera_)) {
      if (options_.camera_params.empty()) {
        std::cout << "Could not read camera model file and no default model given\n";
        return Status::FAILURE;
      }
      prev_camera_ = default_camera_;
    }
    prev_camera_.SetWidth(bitmap->Width());
    prev_camera_.SetHeight(bitmap->Height());
  }

  prev_camera_.SetCameraId(database_->WriteCamera(prev_camera_));
  image->SetCameraId(prev_camera_.CameraId());

  //////////////////////////////////////////////////////////////////////////////
  // Extract GPS data.
  //////////////////////////////////////////////////////////////////////////////

  if (!bitmap->ExifLatitude(&image->TvecPrior(0)) ||
      !bitmap->ExifLongitude(&image->TvecPrior(1)) ||
      !bitmap->ExifAltitude(&image->TvecPrior(2))) {
    image->TvecPrior().setConstant(std::numeric_limits<double>::quiet_NaN());
  }

  *camera = prev_camera_;

  image_folders_.insert(image_folder);

  return Status::SUCCESS;
}

size_t ImageReader::NextIndex() const { return image_index_; }

size_t ImageReader::NumImages() const { return options_.image_list.size(); }

}  // namespace colmap
