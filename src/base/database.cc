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

#include "base/database.h"

#include <fstream>

#include "util/sqlite3_utils.h"
#include "util/string.h"

namespace colmap {
namespace {

typedef Eigen::Matrix<float, Eigen::Dynamic, 4, Eigen::RowMajor>
        FeatureLinesBlob;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    FeatureDescriptorsBlob;
typedef Eigen::Matrix<point2D_t, Eigen::Dynamic, 2, Eigen::RowMajor>
    FeatureMatchesBlob;

void SwapFeatureMatchesBlob(FeatureMatchesBlob* matches) {
  matches->col(0).swap(matches->col(1));
}

FeatureLinesBlob FeatureLinesToBlob(const FeatureLines& lines) {
  FeatureLinesBlob blob(lines.size(), 4);
  for (size_t i = 0; i < lines.size(); ++i) {
    blob.row(i).leftCols<3>() = lines.at(i).Line().transpose().cast<float>();
    blob.row(i)(3) = lines.at(i).IsAligned() ? 1 : 0;
  }
  return blob;
}

FeatureLines FeatureLinesFromBlob(const FeatureLinesBlob& blob) {
  FeatureLines lines(static_cast<size_t>(blob.rows()));
  for (int i = 0; i < blob.rows(); ++i) {
    const Eigen::Vector3d line = blob.row(i).leftCols<3>().transpose().cast<double>();
    // Normalize line to simplify distance computations
    const double norm_factor = line.head<2>().norm();
    lines.at(i).SetLine(line / norm_factor);
    lines.at(i).SetAligned(blob.row(i)(3) > 0);
  }
  return lines;
}

FeatureMatchesBlob FeatureMatchesToBlob(const FeatureMatches& matches) {
  const FeatureMatchesBlob::Index kNumCols = 2;
  FeatureMatchesBlob blob(matches.size(), kNumCols);
  for (size_t i = 0; i < matches.size(); ++i) {
    blob(i, 0) = matches[i].line_idx1;
    blob(i, 1) = matches[i].line_idx2;
  }
  return blob;
}

FeatureMatches FeatureMatchesFromBlob(const FeatureMatchesBlob& blob) {
  CHECK_EQ(blob.cols(), 2);
  FeatureMatches matches(static_cast<size_t>(blob.rows()));
  for (FeatureMatchesBlob::Index i = 0; i < blob.rows(); ++i) {
    matches[i].line_idx1 = blob(i, 0);
    matches[i].line_idx2 = blob(i, 1);
  }
  return matches;
}

template <typename MatrixType>
MatrixType ReadStaticMatrixBlob(sqlite3_stmt* sql_stmt, const int rc,
                                const int col) {
  CHECK_GE(col, 0);

  MatrixType matrix;

  if (rc == SQLITE_ROW) {
    const size_t num_bytes =
        static_cast<size_t>(sqlite3_column_bytes(sql_stmt, col));
    if (num_bytes > 0) {
      CHECK_EQ(num_bytes, matrix.size() * sizeof(typename MatrixType::Scalar));
      memcpy(reinterpret_cast<char*>(matrix.data()),
             sqlite3_column_blob(sql_stmt, col), num_bytes);
    } else {
      matrix = MatrixType::Zero();
    }
  } else {
    matrix = MatrixType::Zero();
  }

  return matrix;
}

template <typename MatrixType>
MatrixType ReadDynamicMatrixBlob(sqlite3_stmt* sql_stmt, const int rc,
                                 const int col) {
  CHECK_GE(col, 0);

  MatrixType matrix;

  if (rc == SQLITE_ROW) {
    const size_t rows =
        static_cast<size_t>(sqlite3_column_int64(sql_stmt, col + 0));
    const size_t cols =
        static_cast<size_t>(sqlite3_column_int64(sql_stmt, col + 1));

    CHECK_GE(rows, 0);
    CHECK_GE(cols, 0);
    matrix = MatrixType(rows, cols);

    const size_t num_bytes =
        static_cast<size_t>(sqlite3_column_bytes(sql_stmt, col + 2));
    CHECK_EQ(matrix.size() * sizeof(typename MatrixType::Scalar), num_bytes);

    memcpy(reinterpret_cast<char*>(matrix.data()),
           sqlite3_column_blob(sql_stmt, col + 2), num_bytes);
  } else {
    const typename MatrixType::Index rows =
        (MatrixType::RowsAtCompileTime == Eigen::Dynamic)
            ? 0
            : MatrixType::RowsAtCompileTime;
    const typename MatrixType::Index cols =
        (MatrixType::ColsAtCompileTime == Eigen::Dynamic)
            ? 0
            : MatrixType::ColsAtCompileTime;
    matrix = MatrixType(rows, cols);
  }

  return matrix;
}

template <typename MatrixType>
void WriteStaticMatrixBlob(sqlite3_stmt* sql_stmt, const MatrixType& matrix,
                           const int col) {
  SQLITE3_CALL(sqlite3_bind_blob(
      sql_stmt, col, reinterpret_cast<const char*>(matrix.data()),
      static_cast<int>(matrix.size() * sizeof(typename MatrixType::Scalar)),
      SQLITE_STATIC));
}

template <typename MatrixType>
void WriteDynamicMatrixBlob(sqlite3_stmt* sql_stmt, const MatrixType& matrix,
                            const int col) {
  CHECK_GE(matrix.rows(), 0);
  CHECK_GE(matrix.cols(), 0);
  CHECK_GE(col, 0);

  const size_t num_bytes = matrix.size() * sizeof(typename MatrixType::Scalar);
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt, col + 0, matrix.rows()));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt, col + 1, matrix.cols()));
  SQLITE3_CALL(sqlite3_bind_blob(sql_stmt, col + 2,
                                 reinterpret_cast<const char*>(matrix.data()),
                                 static_cast<int>(num_bytes), SQLITE_STATIC));
}

Camera ReadCameraRow(sqlite3_stmt* sql_stmt) {
  Camera camera;

  camera.SetCameraId(static_cast<camera_t>(sqlite3_column_int64(sql_stmt, 0)));
  camera.SetModelId(sqlite3_column_int64(sql_stmt, 1));
  camera.SetWidth(static_cast<size_t>(sqlite3_column_int64(sql_stmt, 2)));
  camera.SetHeight(static_cast<size_t>(sqlite3_column_int64(sql_stmt, 3)));

  const size_t num_params_bytes =
      static_cast<size_t>(sqlite3_column_bytes(sql_stmt, 4));
  const size_t num_params = num_params_bytes / sizeof(double);
  CHECK_EQ(num_params, camera.NumParams());
  memcpy(camera.ParamsData(), sqlite3_column_blob(sql_stmt, 4),
         num_params_bytes);

  camera.SetPriorFocalLength(sqlite3_column_int64(sql_stmt, 5) != 0);

  return camera;
}

Image ReadImageRow(sqlite3_stmt* sql_stmt) {
  Image image;

  image.SetImageId(static_cast<image_t>(sqlite3_column_int64(sql_stmt, 0)));
  image.SetName(std::string(
      reinterpret_cast<const char*>(sqlite3_column_text(sql_stmt, 1))));
  image.SetCameraId(static_cast<camera_t>(sqlite3_column_int64(sql_stmt, 2)));

  // NaNs are automatically converted to NULLs in SQLite.
  for (size_t i = 0; i < 4; ++i) {
    if (sqlite3_column_type(sql_stmt, i + 3) != SQLITE_NULL) {
      image.QvecPrior(i) = sqlite3_column_double(sql_stmt, i + 3);
    }
  }

  // NaNs are automatically converted to NULLs in SQLite.
  for (size_t i = 0; i < 3; ++i) {
    if (sqlite3_column_type(sql_stmt, i + 7) != SQLITE_NULL) {
      image.TvecPrior(i) = sqlite3_column_double(sql_stmt, i + 7);
    }
  }

  return image;
}

}  // namespace

const size_t Database::kMaxNumImages =
    static_cast<size_t>(std::numeric_limits<int32_t>::max());

std::mutex Database::update_schema_mutex_;

Database::Database() : database_(nullptr) {}

Database::Database(const std::string& path, const bool create_line_table) : Database() { Open(path, create_line_table); }

Database::~Database() { Close(); }

void Database::Open(const std::string& path, const bool create_line_table) {
  Close();

  // SQLITE_OPEN_NOMUTEX specifies that the connection should not have a
  // mutex (so that we don't serialize the connection's operations).
  // Modifications to the database will still be serialized, but multiple
  // connections can read concurrently.
  SQLITE3_CALL(sqlite3_open_v2(
      path.c_str(), &database_,
      SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_NOMUTEX,
      nullptr));

  // Don't wait for the operating system to write the changes to disk
  SQLITE3_EXEC(database_, "PRAGMA synchronous=OFF", nullptr);

  // Use faster journaling mode
  SQLITE3_EXEC(database_, "PRAGMA journal_mode=WAL", nullptr);

  // Store temporary tables and indices in memory
  SQLITE3_EXEC(database_, "PRAGMA temp_store=MEMORY", nullptr);

  // Disabled by default
  SQLITE3_EXEC(database_, "PRAGMA foreign_keys=ON", nullptr);

  CreateTables();

  if (create_line_table) {
    CreateLineFeaturesTable();
    CreateGravityTable();
  }
  PrepareSQLStatements(create_line_table);
}

void Database::Close() {
  if (database_ != nullptr) {
    FinalizeSQLStatements();
    sqlite3_close_v2(database_);
    database_ = nullptr;
  }
}

void Database::CloneTo(Database* target) {
  sqlite3_backup* backup_obj =
      sqlite3_backup_init(target->database_, "main", database_, "main");
  CHECK_NOTNULL(backup_obj);

  int rc = sqlite3_backup_step(backup_obj, -1);

  CHECK_EQ(rc, SQLITE_DONE);

  rc = sqlite3_backup_finish(backup_obj);
  CHECK_EQ(rc, SQLITE_OK);

  rc = sqlite3_errcode(target->database_);
  CHECK_EQ(rc, SQLITE_OK);
}

bool Database::ExistsCamera(const camera_t camera_id) const {
  return ExistsRowId(sql_stmt_exists_camera_, camera_id);
}

bool Database::ExistsImage(const image_t image_id) const {
  return ExistsRowId(sql_stmt_exists_image_id_, image_id);
}

bool Database::ExistsImageWithName(std::string name) const {
  return ExistsRowString(sql_stmt_exists_image_name_, name);
}

bool Database::ExistsDescriptors(const image_t image_id) const {
  return ExistsRowId(sql_stmt_exists_descriptors_, image_id);
}

bool Database::ExistsMatches(const image_t image_id1,
                             const image_t image_id2) const {
  return ExistsRowId(sql_stmt_exists_matches_,
                     ImagePairToPairId(image_id1, image_id2));
}

bool Database::ExistsLineFeaturesTable() const {
    return ExistsTable("line_features");
}

bool Database::ExistsLineFeatures(const colmap::image_t image_id) const {
    return ExistsRowId(sql_stmt_exists_lines_, image_id);
}

bool Database::ExistsImageGravity(const image_t image_id) const {
  return ExistsRowId(sql_stmt_exists_gravity_, image_id);
}

size_t Database::NumCameras() const { return CountRows("cameras"); }

size_t Database::NumImages() const { return CountRows("images"); }

size_t Database::NumLines() const { return SumColumn("rows", "line_features"); }

size_t Database::NumDescriptors() const {
  return SumColumn("rows", "descriptors");
}

size_t Database::MaxNumDescriptors() const {
  return MaxColumn("rows", "descriptors");
}

size_t Database::NumDescriptorsForImage(const image_t image_id) const {
  return CountRowsForEntry(sql_stmt_num_descriptors_, image_id);
}

size_t Database::NumMatches() const { return SumColumn("rows", "matches"); }

size_t Database::NumMatchedImagePairs() const { return CountRows("matches"); }

Camera Database::ReadCamera(const camera_t camera_id) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_camera_, 1, camera_id));

  Camera camera;

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_camera_));
  if (rc == SQLITE_ROW) {
    camera = ReadCameraRow(sql_stmt_read_camera_);
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_camera_));

  return camera;
}

std::vector<Camera> Database::ReadAllCameras() const {
  std::vector<Camera> cameras;

  while (SQLITE3_CALL(sqlite3_step(sql_stmt_read_cameras_)) == SQLITE_ROW) {
    cameras.push_back(ReadCameraRow(sql_stmt_read_cameras_));
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_cameras_));

  return cameras;
}

Image Database::ReadImage(const image_t image_id) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_image_id_, 1, image_id));

  Image image;

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_image_id_));
  if (rc == SQLITE_ROW) {
    image = ReadImageRow(sql_stmt_read_image_id_);
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_image_id_));

  return image;
}

Image Database::ReadImageWithName(const std::string& name) const {
  SQLITE3_CALL(sqlite3_bind_text(sql_stmt_read_image_name_, 1, name.c_str(),
                                 static_cast<int>(name.size()), SQLITE_STATIC));

  Image image;

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_image_name_));
  if (rc == SQLITE_ROW) {
    image = ReadImageRow(sql_stmt_read_image_name_);
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_image_name_));

  return image;
}

std::vector<Image> Database::ReadAllImages() const {
  std::vector<Image> images;
  images.reserve(NumImages());

  while (SQLITE3_CALL(sqlite3_step(sql_stmt_read_images_)) == SQLITE_ROW) {
    images.push_back(ReadImageRow(sql_stmt_read_images_));
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_images_));

  return images;
}

FeatureLines Database::ReadFeatureLines(const image_t image_id) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_lines_, 1, image_id));

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_lines_));
  const FeatureLinesBlob blob = ReadDynamicMatrixBlob<FeatureLinesBlob>(
          sql_stmt_read_lines_, rc, 0);

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_lines_));

  return FeatureLinesFromBlob(blob);
}

Eigen::Vector3d Database::ReadImageGravity(const image_t image_id) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_gravity_, 1, image_id));
  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_gravity_));

  Eigen::Vector3d gravity;
  if (rc == SQLITE_ROW) {
    gravity(0) =
        static_cast<double>(sqlite3_column_double(sql_stmt_read_gravity_, 0));
    gravity(1) =
        static_cast<double>(sqlite3_column_double(sql_stmt_read_gravity_, 1));
    gravity(2) =
        static_cast<double>(sqlite3_column_double(sql_stmt_read_gravity_, 2));
  } else {
    gravity.setConstant(std::numeric_limits<double>::quiet_NaN());
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_gravity_));

  return gravity;
}

FeatureDescriptors Database::ReadDescriptors(const image_t image_id) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_descriptors_, 1, image_id));

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_descriptors_));
  const FeatureDescriptors descriptors =
      ReadDynamicMatrixBlob<FeatureDescriptors>(sql_stmt_read_descriptors_, rc,
                                                0);

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_descriptors_));

  return descriptors;
}

FeatureMatches Database::ReadMatches(image_t image_id1,
                                     image_t image_id2) const {
  const image_pair_t pair_id = ImagePairToPairId(image_id1, image_id2);
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_read_matches_, 1, pair_id));

  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_matches_));
  FeatureMatchesBlob blob =
      ReadDynamicMatrixBlob<FeatureMatchesBlob>(sql_stmt_read_matches_, rc, 0);

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_matches_));

  if (SwapImagePair(image_id1, image_id2)) {
    SwapFeatureMatchesBlob(&blob);
  }

  return FeatureMatchesFromBlob(blob);
}

std::vector<std::pair<image_pair_t, FeatureMatches>> Database::ReadAllMatches()
    const {
  std::vector<std::pair<image_pair_t, FeatureMatches>> all_matches;

  int rc;
  while ((rc = SQLITE3_CALL(sqlite3_step(sql_stmt_read_matches_all_))) ==
         SQLITE_ROW) {
    const image_pair_t pair_id = static_cast<image_pair_t>(
        sqlite3_column_int64(sql_stmt_read_matches_all_, 0));
    const FeatureMatchesBlob blob = ReadDynamicMatrixBlob<FeatureMatchesBlob>(
        sql_stmt_read_matches_all_, rc, 1);
    all_matches.emplace_back(pair_id, FeatureMatchesFromBlob(blob));
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_matches_all_));

  return all_matches;
}

void Database::ReadNumMatches(
    std::vector<std::pair<image_t, image_t> >* image_pairs,
    std::vector<int>* num_inliers) const {
  const auto num_matches = NumMatches();
  image_pairs->reserve(num_matches);
  num_inliers->reserve(num_matches);

  while(SQLITE3_CALL(sqlite3_step(
          sql_stmt_read_num_matches_)) == SQLITE_ROW) {
    image_t image_id1;
    image_t image_id2;
    const image_pair_t  pair_id = static_cast<image_pair_t>(
        sqlite3_column_int64(sql_stmt_read_num_matches_, 0));
    PairIdToImagePair(pair_id, &image_id1, &image_id2);
    image_pairs->emplace_back(image_id1, image_id2);

    const int rows = static_cast<int>(
        sqlite3_column_int64(sql_stmt_read_num_matches_, 1));
    num_inliers->push_back(rows);
  }
  SQLITE3_CALL(sqlite3_reset(sql_stmt_read_num_matches_));
}

camera_t Database::WriteCamera(const Camera& camera,
                               const bool use_camera_id) const {
  if (use_camera_id) {
    CHECK(!ExistsCamera(camera.CameraId())) << "camera_id must be unique";
    SQLITE3_CALL(
        sqlite3_bind_int64(sql_stmt_add_camera_, 1, camera.CameraId()));
  } else {
    SQLITE3_CALL(sqlite3_bind_null(sql_stmt_add_camera_, 1));
  }

  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_camera_, 2, camera.ModelId()));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_camera_, 3,
                                  static_cast<sqlite3_int64>(camera.Width())));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_camera_, 4,
                                  static_cast<sqlite3_int64>(camera.Height())));

  const size_t num_params_bytes = sizeof(double) * camera.NumParams();
  SQLITE3_CALL(sqlite3_bind_blob(sql_stmt_add_camera_, 5, camera.ParamsData(),
                                 static_cast<int>(num_params_bytes),
                                 SQLITE_STATIC));

  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_camera_, 6,
                                  camera.HasPriorFocalLength()));

  SQLITE3_CALL(sqlite3_step(sql_stmt_add_camera_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_add_camera_));

  return static_cast<camera_t>(sqlite3_last_insert_rowid(database_));
}

image_t Database::WriteImage(const Image& image,
                             const bool use_image_id) const {
  if (use_image_id) {
    CHECK(!ExistsImage(image.ImageId())) << "image_id must be unique";
    SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_image_, 1, image.ImageId()));
  } else {
    SQLITE3_CALL(sqlite3_bind_null(sql_stmt_add_image_, 1));
  }

  SQLITE3_CALL(sqlite3_bind_text(sql_stmt_add_image_, 2, image.Name().c_str(),
                                 static_cast<int>(image.Name().size()),
                                 SQLITE_STATIC));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_add_image_, 3, image.CameraId()));

  // NaNs are automatically converted to NULLs in SQLite.
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 4, image.QvecPrior(0)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 5, image.QvecPrior(1)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 6, image.QvecPrior(2)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 7, image.QvecPrior(3)));

  // NaNs are automatically converted to NULLs in SQLite.
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 8, image.TvecPrior(0)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_add_image_, 9, image.TvecPrior(1)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_add_image_, 10, image.TvecPrior(2)));

  SQLITE3_CALL(sqlite3_step(sql_stmt_add_image_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_add_image_));

  return static_cast<image_t>(sqlite3_last_insert_rowid(database_));
}

void Database::WriteFeatureLines(const image_t image_id,
                              const FeatureLines& lines) const {
  const FeatureLinesBlob blob = FeatureLinesToBlob(lines);

  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_write_lines_, 1, image_id));
  WriteDynamicMatrixBlob(sql_stmt_write_lines_, blob, 2);

  SQLITE3_CALL(sqlite3_step(sql_stmt_write_lines_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_write_lines_));
}

void Database::WriteDescriptors(const image_t image_id,
                                const FeatureDescriptors& descriptors) const {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_write_descriptors_, 1, image_id));
  WriteDynamicMatrixBlob(sql_stmt_write_descriptors_, descriptors, 2);

  SQLITE3_CALL(sqlite3_step(sql_stmt_write_descriptors_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_write_descriptors_));
}

void Database::WriteMatches(const image_t image_id1, const image_t image_id2,
                            const FeatureMatches& matches) const {
  const image_pair_t pair_id = ImagePairToPairId(image_id1, image_id2);
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_write_matches_, 1, pair_id));

  // Important: the swapped data must live until the query is executed.
  FeatureMatchesBlob blob = FeatureMatchesToBlob(matches);
  if (SwapImagePair(image_id1, image_id2)) {
    SwapFeatureMatchesBlob(&blob);
    WriteDynamicMatrixBlob(sql_stmt_write_matches_, blob, 2);
  } else {
    WriteDynamicMatrixBlob(sql_stmt_write_matches_, blob, 2);
  }

  SQLITE3_CALL(sqlite3_step(sql_stmt_write_matches_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_write_matches_));
}

void Database::WriteImageGravity(const image_t image_id,
                                 const Eigen::Vector3d& gravity_direction) {
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_write_gravity_, 1, image_id));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_write_gravity_, 2, gravity_direction(0)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_write_gravity_, 3, gravity_direction(1)));
  SQLITE3_CALL(sqlite3_bind_double(sql_stmt_write_gravity_, 4, gravity_direction(2)));

  SQLITE3_CALL(sqlite3_step(sql_stmt_write_gravity_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_write_gravity_));
}

void Database::UpdateCamera(const Camera& camera) const {
  SQLITE3_CALL(
      sqlite3_bind_int64(sql_stmt_update_camera_, 1, camera.ModelId()));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_update_camera_, 2,
                                  static_cast<sqlite3_int64>(camera.Width())));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_update_camera_, 3,
                                  static_cast<sqlite3_int64>(camera.Height())));

  const size_t num_params_bytes = sizeof(double) * camera.NumParams();
  SQLITE3_CALL(
      sqlite3_bind_blob(sql_stmt_update_camera_, 4, camera.ParamsData(),
                        static_cast<int>(num_params_bytes), SQLITE_STATIC));

  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_update_camera_, 5,
                                  camera.HasPriorFocalLength()));

  SQLITE3_CALL(
      sqlite3_bind_int64(sql_stmt_update_camera_, 6, camera.CameraId()));

  SQLITE3_CALL(sqlite3_step(sql_stmt_update_camera_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_update_camera_));
}

void Database::UpdateImage(const Image& image) const {
  SQLITE3_CALL(
      sqlite3_bind_text(sql_stmt_update_image_, 1, image.Name().c_str(),
                        static_cast<int>(image.Name().size()), SQLITE_STATIC));
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_update_image_, 2, image.CameraId()));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 3, image.QvecPrior(0)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 4, image.QvecPrior(1)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 5, image.QvecPrior(2)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 6, image.QvecPrior(3)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 7, image.TvecPrior(0)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 8, image.TvecPrior(1)));
  SQLITE3_CALL(
      sqlite3_bind_double(sql_stmt_update_image_, 9, image.TvecPrior(2)));

  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_update_image_, 10, image.ImageId()));

  SQLITE3_CALL(sqlite3_step(sql_stmt_update_image_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_update_image_));
}

void Database::DeleteMatches(const image_t image_id1,
                             const image_t image_id2) const {
  const image_pair_t pair_id = ImagePairToPairId(image_id1, image_id2);
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt_delete_matches_, 1,
                                  static_cast<sqlite3_int64>(pair_id)));
  SQLITE3_CALL(sqlite3_step(sql_stmt_delete_matches_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_delete_matches_));
}

void Database::ClearMatches() const {
  SQLITE3_CALL(sqlite3_step(sql_stmt_clear_matches_));
  SQLITE3_CALL(sqlite3_reset(sql_stmt_clear_matches_));
}

void Database::BeginTransaction() const {
  SQLITE3_EXEC(database_, "BEGIN TRANSACTION", nullptr);
}

void Database::EndTransaction() const {
  SQLITE3_EXEC(database_, "END TRANSACTION", nullptr);
}

void Database::PrepareSQLStatements(const bool with_line_features) {
  sql_stmts_.clear();

  std::string sql;

  //////////////////////////////////////////////////////////////////////////////
  // num_*
  //////////////////////////////////////////////////////////////////////////////

  sql = "SELECT rows FROM descriptors WHERE image_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_num_descriptors_, 0));
  sql_stmts_.push_back(sql_stmt_num_descriptors_);

  //////////////////////////////////////////////////////////////////////////////
  // exists_*
  //////////////////////////////////////////////////////////////////////////////
  sql = "SELECT 1 FROM cameras WHERE camera_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_exists_camera_, 0));
  sql_stmts_.push_back(sql_stmt_exists_camera_);

  sql = "SELECT 1 FROM images WHERE image_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_exists_image_id_, 0));
  sql_stmts_.push_back(sql_stmt_exists_image_id_);

  sql = "SELECT 1 FROM images WHERE name = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_exists_image_name_, 0));
  sql_stmts_.push_back(sql_stmt_exists_image_name_);

  sql = "SELECT 1 FROM descriptors WHERE image_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_exists_descriptors_, 0));
  sql_stmts_.push_back(sql_stmt_exists_descriptors_);

  sql = "SELECT 1 FROM matches WHERE pair_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_exists_matches_, 0));
  sql_stmts_.push_back(sql_stmt_exists_matches_);

  if (with_line_features) {
      sql = "SELECT 1 FROM line_features WHERE image_id = ?;";
      SQLITE3_CALL(
          sqlite3_prepare_v2(database_, sql.c_str(), -1,
              &sql_stmt_exists_lines_, 0));
      sql_stmts_.push_back(sql_stmt_exists_lines_);

      sql = "SELECT 1 FROM gravity_directions WHERE image_id = ?;";
      SQLITE3_CALL(
          sqlite3_prepare_v2(database_, sql.c_str(), -1,
              &sql_stmt_exists_gravity_, 0));
      sql_stmts_.push_back(sql_stmt_exists_gravity_);
  }

  //////////////////////////////////////////////////////////////////////////////
  // add_*
  //////////////////////////////////////////////////////////////////////////////
  sql =
      "INSERT INTO cameras(camera_id, model, width, height, params, "
      "prior_focal_length) VALUES(?, ?, ?, ?, ?, ?);";
  SQLITE3_CALL(
      sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt_add_camera_, 0));
  sql_stmts_.push_back(sql_stmt_add_camera_);

  sql =
      "INSERT INTO images(image_id, name, camera_id, prior_qw, prior_qx, "
      "prior_qy, prior_qz, prior_tx, prior_ty, prior_tz) VALUES(?, ?, ?, ?, ?, "
      "?, ?, ?, ?, ?);";
  SQLITE3_CALL(
      sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt_add_image_, 0));
  sql_stmts_.push_back(sql_stmt_add_image_);

  //////////////////////////////////////////////////////////////////////////////
  // update_*
  //////////////////////////////////////////////////////////////////////////////
  sql =
      "UPDATE cameras SET model=?, width=?, height=?, params=?, "
      "prior_focal_length=? WHERE camera_id=?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_update_camera_, 0));
  sql_stmts_.push_back(sql_stmt_update_camera_);

  sql =
      "UPDATE images SET name=?, camera_id=?, prior_qw=?, prior_qx=?, "
      "prior_qy=?, prior_qz=?, prior_tx=?, prior_ty=?, prior_tz=? WHERE "
      "image_id=?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_update_image_, 0));
  sql_stmts_.push_back(sql_stmt_update_image_);

  //////////////////////////////////////////////////////////////////////////////
  // read_*
  //////////////////////////////////////////////////////////////////////////////
  sql = "SELECT * FROM cameras WHERE camera_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_camera_, 0));
  sql_stmts_.push_back(sql_stmt_read_camera_);

  sql = "SELECT * FROM cameras;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_cameras_, 0));
  sql_stmts_.push_back(sql_stmt_read_cameras_);

  sql = "SELECT * FROM images WHERE image_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_image_id_, 0));
  sql_stmts_.push_back(sql_stmt_read_image_id_);

  sql = "SELECT * FROM images WHERE name = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_image_name_, 0));
  sql_stmts_.push_back(sql_stmt_read_image_name_);

  sql = "SELECT * FROM images;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_images_, 0));
  sql_stmts_.push_back(sql_stmt_read_images_);

  if (with_line_features) {
      sql = "SELECT rows, cols, data  FROM line_features WHERE image_id = ?;";
      SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                      &sql_stmt_read_lines_, 0));
      sql_stmts_.push_back(sql_stmt_read_lines_);

      sql = "SELECT x, y, z FROM gravity_directions WHERE image_id = ?;";
      SQLITE3_CALL(
          sqlite3_prepare_v2(database_, sql.c_str(),  -1,
              &sql_stmt_read_gravity_, 0));
      sql_stmts_.push_back(sql_stmt_read_gravity_);
  }

  sql = "SELECT rows, cols, data FROM descriptors WHERE image_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_descriptors_, 0));
  sql_stmts_.push_back(sql_stmt_read_descriptors_);

  sql = "SELECT rows, cols, data FROM matches WHERE pair_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_matches_, 0));
  sql_stmts_.push_back(sql_stmt_read_matches_);

  sql = "SELECT * FROM matches WHERE rows > 0;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_read_matches_all_, 0));
  sql_stmts_.push_back(sql_stmt_read_matches_all_);

  sql = "SELECT pair_id, rows FROM matches WHERE rows > 0;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                        &sql_stmt_read_num_matches_,
                                        0));

  //////////////////////////////////////////////////////////////////////////////
  // write_*
  //////////////////////////////////////////////////////////////////////////////

  if (with_line_features) {
      sql = "INSERT INTO line_features(image_id, rows, cols, data) VALUES(?, ?, ?, ?);";
      SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                      &sql_stmt_write_lines_, 0));
      sql_stmts_.push_back(sql_stmt_write_lines_);

      sql = "INSERT INTO gravity_directions(image_id, x, y, z) VALUES(?, ?, ?, ?);";
      SQLITE3_CALL(
          sqlite3_prepare_v2(database_, sql.c_str(), -1,
              &sql_stmt_write_gravity_, 0));
      sql_stmts_.push_back(sql_stmt_write_gravity_);
  }

  sql =
      "INSERT INTO descriptors(image_id, rows, cols, data) VALUES(?, ?, ?, ?);";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_write_descriptors_, 0));
  sql_stmts_.push_back(sql_stmt_write_descriptors_);

  sql = "INSERT INTO matches(pair_id, rows, cols, data) VALUES(?, ?, ?, ?);";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_write_matches_, 0));
  sql_stmts_.push_back(sql_stmt_write_matches_);

  //////////////////////////////////////////////////////////////////////////////
  // delete_*
  //////////////////////////////////////////////////////////////////////////////

  sql = "DELETE FROM matches WHERE pair_id = ?;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_delete_matches_, 0));
  sql_stmts_.push_back(sql_stmt_delete_matches_);

  //////////////////////////////////////////////////////////////////////////////
  // clear_*
  //////////////////////////////////////////////////////////////////////////////
  sql = "DELETE FROM matches;";
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1,
                                  &sql_stmt_clear_matches_, 0));
  sql_stmts_.push_back(sql_stmt_clear_matches_);
}

void Database::FinalizeSQLStatements() {
  for (const auto& sql_stmt : sql_stmts_) {
    SQLITE3_CALL(sqlite3_finalize(sql_stmt));
  }
}

void Database::CreateTables() const {
  CreateCameraTable();
  CreateImageTable();
  CreateDescriptorsTable();
  CreateMatchesTable();
}

void Database::CreateCameraTable() const {
  const std::string sql =
      "CREATE TABLE IF NOT EXISTS cameras"
      "   (camera_id            INTEGER  PRIMARY KEY AUTOINCREMENT  NOT NULL,"
      "    model                INTEGER                             NOT NULL,"
      "    width                INTEGER                             NOT NULL,"
      "    height               INTEGER                             NOT NULL,"
      "    params               BLOB,"
      "    prior_focal_length   INTEGER                             NOT NULL);";

  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

void Database::CreateImageTable() const {
  const std::string sql = StringPrintf(
      "CREATE TABLE IF NOT EXISTS images"
      "   (image_id   INTEGER  PRIMARY KEY AUTOINCREMENT  NOT NULL,"
      "    name       TEXT                                NOT NULL UNIQUE,"
      "    camera_id  INTEGER                             NOT NULL,"
      "    prior_qw   REAL,"
      "    prior_qx   REAL,"
      "    prior_qy   REAL,"
      "    prior_qz   REAL,"
      "    prior_tx   REAL,"
      "    prior_ty   REAL,"
      "    prior_tz   REAL,"
      "CONSTRAINT image_id_check CHECK(image_id >= 0 and image_id < %d),"
      "FOREIGN KEY(camera_id) REFERENCES cameras(camera_id));"
      "CREATE UNIQUE INDEX IF NOT EXISTS index_name ON images(name);",
      kMaxNumImages);

  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

void Database::CreateDescriptorsTable() const {
  const std::string sql =
      "CREATE TABLE IF NOT EXISTS descriptors"
      "   (image_id  INTEGER  PRIMARY KEY  NOT NULL,"
      "    rows      INTEGER               NOT NULL,"
      "    cols      INTEGER               NOT NULL,"
      "    data      BLOB,"
      "FOREIGN KEY(image_id) REFERENCES images(image_id) ON DELETE CASCADE);";

  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

void Database::CreateMatchesTable() const {
  const std::string sql =
      "CREATE TABLE IF NOT EXISTS matches"
      "   (pair_id  INTEGER  PRIMARY KEY  NOT NULL,"
      "    rows     INTEGER               NOT NULL,"
      "    cols     INTEGER               NOT NULL,"
      "    data     BLOB);";

  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

void Database::CreateLineFeaturesTable() const {
  const std::string sql =
      "CREATE TABLE IF NOT EXISTS line_features"
      "   (image_id INTEGER PRIMARY KEY NOT NULL,"
      "    rows INTEGER NOT NULL,"
      "    cols INTEGER NOT NULL,"
      "    data BLOB,"
      "    FOREIGN KEY(image_id) REFERENCES images(image_id) ON DELETE CASCADE)";
  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

void Database::CreateGravityTable() const {
  const std::string sql =
      "CREATE TABLE IF NOT EXISTS gravity_directions"
      "   (image_id INTEGER PRIMARY KEY NOT NULL,"
      "    x REAL,"
      "    y REAL,"
      "    z REAL);";
  SQLITE3_EXEC(database_, sql.c_str(), nullptr);
}

bool Database::ExistsTable(const std::string& table_name) const {
  const std::string sql =
      "SELECT name FROM sqlite_master WHERE type='table' AND name = ?;";

  sqlite3_stmt* sql_stmt;
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt, 0));

  SQLITE3_CALL(sqlite3_bind_text(sql_stmt, 1, table_name.c_str(),
                                 static_cast<int>(table_name.size()),
                                 SQLITE_STATIC));

  const bool exists = SQLITE3_CALL(sqlite3_step(sql_stmt)) == SQLITE_ROW;

  SQLITE3_CALL(sqlite3_finalize(sql_stmt));

  return exists;
}

bool Database::ExistsColumn(const std::string& table_name,
                            const std::string& column_name) const {
  const std::string sql =
      StringPrintf("PRAGMA table_info(%s);", table_name.c_str());

  sqlite3_stmt* sql_stmt;
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt, 0));

  bool exists_column = false;
  while (SQLITE3_CALL(sqlite3_step(sql_stmt)) == SQLITE_ROW) {
    const std::string result =
        reinterpret_cast<const char*>(sqlite3_column_text(sql_stmt, 1));
    if (column_name == result) {
      exists_column = true;
      break;
    }
  }

  SQLITE3_CALL(sqlite3_finalize(sql_stmt));

  return exists_column;
}

bool Database::ExistsRowId(sqlite3_stmt* sql_stmt,
                           const sqlite3_int64 row_id) const {

  CHECK(sql_stmt) << "SQL statement is not initialized";
  SQLITE3_CALL(
      sqlite3_bind_int64(sql_stmt, 1, static_cast<sqlite3_int64>(row_id)));

  const bool exists = SQLITE3_CALL(sqlite3_step(sql_stmt)) == SQLITE_ROW;

  SQLITE3_CALL(sqlite3_reset(sql_stmt));

  return exists;
}

bool Database::ExistsRowString(sqlite3_stmt* sql_stmt,
                               const std::string& row_entry) const {
  CHECK(sql_stmt) << "SQL statement is not initialized";
  SQLITE3_CALL(sqlite3_bind_text(sql_stmt, 1, row_entry.c_str(),
                                 static_cast<int>(row_entry.size()),
                                 SQLITE_STATIC));

  const bool exists = SQLITE3_CALL(sqlite3_step(sql_stmt)) == SQLITE_ROW;

  SQLITE3_CALL(sqlite3_reset(sql_stmt));

  return exists;
}

size_t Database::CountRows(const std::string& table) const {
  const std::string sql =
      StringPrintf("SELECT COUNT(*) FROM %s;", table.c_str());

  sqlite3_stmt* sql_stmt;
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt, 0));

  size_t count = 0;
  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt));
  if (rc == SQLITE_ROW) {
    count = static_cast<size_t>(sqlite3_column_int64(sql_stmt, 0));
  }

  SQLITE3_CALL(sqlite3_finalize(sql_stmt));

  return count;
}

size_t Database::CountRowsForEntry(sqlite3_stmt* sql_stmt,
                                   const sqlite3_int64 row_id) const {
  CHECK(sql_stmt) << "SQL statement is not initialized";
  SQLITE3_CALL(sqlite3_bind_int64(sql_stmt, 1, row_id));

  size_t count = 0;
  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt));
  if (rc == SQLITE_ROW) {
    count = static_cast<size_t>(sqlite3_column_int64(sql_stmt, 0));
  }

  SQLITE3_CALL(sqlite3_reset(sql_stmt));

  return count;
}

size_t Database::SumColumn(const std::string& column,
                           const std::string& table) const {
  const std::string sql =
      StringPrintf("SELECT SUM(%s) FROM %s;", column.c_str(), table.c_str());

  sqlite3_stmt* sql_stmt;
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt, 0));

  size_t sum = 0;
  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt));
  if (rc == SQLITE_ROW) {
    sum = static_cast<size_t>(sqlite3_column_int64(sql_stmt, 0));
  }

  SQLITE3_CALL(sqlite3_finalize(sql_stmt));

  return sum;
}

size_t Database::MaxColumn(const std::string& column,
                           const std::string& table) const {
  const std::string sql =
      StringPrintf("SELECT MAX(%s) FROM %s;", column.c_str(), table.c_str());

  sqlite3_stmt* sql_stmt;
  SQLITE3_CALL(sqlite3_prepare_v2(database_, sql.c_str(), -1, &sql_stmt, 0));

  size_t max = 0;
  const int rc = SQLITE3_CALL(sqlite3_step(sql_stmt));
  if (rc == SQLITE_ROW) {
    max = static_cast<size_t>(sqlite3_column_int64(sql_stmt, 0));
  }

  SQLITE3_CALL(sqlite3_finalize(sql_stmt));

  return max;
}

DatabaseTransaction::DatabaseTransaction(Database* database)
    : database_(database), database_lock_(database->transaction_mutex_) {
  CHECK_NOTNULL(database_);
  database_->BeginTransaction();
}

DatabaseTransaction::~DatabaseTransaction() { database_->EndTransaction(); }

}  // namespace colmap
