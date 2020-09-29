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

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <algorithm>

#include "base/pose.h"
#include "controllers/automatic_reconstruction.h"
#include "controllers/bundle_adjustment.h"
#include "feature/extraction.h"
#include "feature/matching.h"
#include "feature/utils.h"
#include "ui/main_window.h"
#include "util/opengl_utils.h"
#include "util/random.h"
#include "util/string.h"
#include "util/version.h"

#include "init/initializer.h"

using namespace colmap;

#ifdef CUDA_ENABLED
const bool kUseOpenGL = false;
#else
const bool kUseOpenGL = true;
#endif

int RunGraphicalUserInterface(int argc, char** argv) {
  OptionManager options;

  std::string import_path;

  if (argc > 1) {
    options.AddDefaultOption("import_path", &import_path);
    options.AddAllOptions();
    options.Parse(argc, argv);
  }

#if (QT_VERSION >= QT_VERSION_CHECK(5, 6, 0))
  QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
  QApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
#endif

  Q_INIT_RESOURCE(resources);

  QApplication app(argc, argv);

  MainWindow main_window(options);
  main_window.show();

  if (!import_path.empty()) {
    main_window.ImportReconstruction(import_path);
  }

  return app.exec();
}

int RunAutomaticReconstructor(int argc, char** argv) {
  AutomaticReconstructionController::Options reconstruction_options;
  std::string data_type = "individual";
  std::string quality = "high";

  OptionManager options;
  options.AddRequiredOption("workspace_path",
                            &reconstruction_options.workspace_path);
  options.AddRequiredOption("image_path", &reconstruction_options.image_path);
  options.AddDefaultOption("mask_path", &reconstruction_options.mask_path);
  options.AddDefaultOption("data_type", &data_type,
                           "{individual, video }");
  options.AddDefaultOption("quality", &quality, "{low, medium, high, extreme}");
  options.AddDefaultOption("num_threads", &reconstruction_options.num_threads);
  options.AddDefaultOption("use_gpu", &reconstruction_options.use_gpu);
  options.AddDefaultOption("gpu_index", &reconstruction_options.gpu_index);
  options.Parse(argc, argv);

  StringToLower(&data_type);
  if (data_type == "individual") {
    reconstruction_options.data_type =
        AutomaticReconstructionController::DataType::INDIVIDUAL;
  } else if (data_type == "video") {
    reconstruction_options.data_type =
        AutomaticReconstructionController::DataType::VIDEO;
  } else {
    LOG(FATAL) << "Invalid data type provided";
  }

  StringToLower(&quality);
  if (quality == "low") {
    reconstruction_options.quality =
        AutomaticReconstructionController::Quality::LOW;
  } else if (quality == "medium") {
    reconstruction_options.quality =
        AutomaticReconstructionController::Quality::MEDIUM;
  } else if (quality == "high") {
    reconstruction_options.quality =
        AutomaticReconstructionController::Quality::HIGH;
  } else if (quality == "extreme") {
    reconstruction_options.quality =
        AutomaticReconstructionController::Quality::EXTREME;
  } else {
    LOG(FATAL) << "Invalid quality provided";
  }

  ReconstructionManager reconstruction_manager;

  if (reconstruction_options.use_gpu && kUseOpenGL) {
    QApplication app(argc, argv);
    AutomaticReconstructionController controller(reconstruction_options,
                                                 &reconstruction_manager);
    RunThreadWithOpenGLContext(&controller);
  } else {
    AutomaticReconstructionController controller(reconstruction_options,
                                                 &reconstruction_manager);
    controller.Start();
    controller.Wait();
  }

  return EXIT_SUCCESS;
}

int RunBundleAdjuster(int argc, char** argv) {
  std::string input_path;
  std::string output_path;

  OptionManager options;
  options.AddRequiredOption("input_path", &input_path);
  options.AddRequiredOption("output_path", &output_path);
  options.AddBundleAdjustmentOptions();
  options.Parse(argc, argv);

  if (!ExistsDir(input_path)) {
    std::cerr << "ERROR: `input_path` is not a directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!ExistsDir(output_path)) {
    std::cerr << "ERROR: `output_path` is not a directory" << std::endl;
    return EXIT_FAILURE;
  }

  Reconstruction reconstruction;
  reconstruction.Read(input_path);

  BundleAdjustmentController ba_controller(options, &reconstruction);
  ba_controller.Start();
  ba_controller.Wait();

  reconstruction.Write(output_path);

  return EXIT_SUCCESS;
}

int RunDatabaseCreator(int argc, char** argv) {
  OptionManager options;
  options.AddDatabaseOptions();
  options.Parse(argc, argv);

  Database database(*options.database_path);

  return EXIT_SUCCESS;
}

int RunProjectGenerator(int argc, char** argv) {
  std::string output_path;
  std::string quality = "high";

  OptionManager options;
  options.AddRequiredOption("output_path", &output_path);
  options.AddDefaultOption("quality", &quality, "{low, medium, high, extreme}");
  options.Parse(argc, argv);

  OptionManager output_options;
  output_options.AddAllOptions();

  StringToLower(&quality);
  if (quality == "low") {
    output_options.ModifyForLowQuality();
  } else if (quality == "medium") {
    output_options.ModifyForMediumQuality();
  } else if (quality == "high") {
    output_options.ModifyForHighQuality();
  } else if (quality == "extreme") {
    output_options.ModifyForExtremeQuality();
  } else {
    LOG(FATAL) << "Invalid quality provided";
  }

  output_options.Write(output_path);

  return EXIT_SUCCESS;
}

int RunExhaustiveMatcher(int argc, char** argv) {
  OptionManager options;
  options.AddDatabaseOptions();
  options.AddExhaustiveMatchingOptions();
  options.Parse(argc, argv);

  std::unique_ptr<QApplication> app;
  if (options.sift_matching->use_gpu && kUseOpenGL) {
    app.reset(new QApplication(argc, argv));
  }

  ExhaustiveFeatureMatcher feature_matcher(*options.exhaustive_matching,
                                           *options.sift_matching,
                                           *options.database_path);

  if (options.sift_matching->use_gpu && kUseOpenGL) {
    RunThreadWithOpenGLContext(&feature_matcher);
  } else {
    feature_matcher.Start();
    feature_matcher.Wait();
  }

  return EXIT_SUCCESS;
}

bool VerifyCameraParams(const std::string& camera_model,
                        const std::string& params) {
  if (!ExistsCameraModelWithName(camera_model)) {
    std::cerr << "ERROR: Camera model does not exist" << std::endl;
    return false;
  }

  const std::vector<double> camera_params = CSVToVector<double>(params);
  const int camera_model_id = CameraModelNameToId(camera_model);

  if (camera_params.size() > 0 &&
      !CameraModelVerifyParams(camera_model_id, camera_params)) {
    std::cerr << "ERROR: Invalid camera parameters" << std::endl;
    return false;
  }
  return true;
}

int RunFeatureExtractor(int argc, char** argv) {
  std::string image_list_path;

  OptionManager options;
  options.AddDatabaseOptions();
  options.AddImageOptions();
  options.AddDefaultOption("image_list_path", &image_list_path);
  options.AddExtractionOptions();
  options.Parse(argc, argv);

  ImageReaderOptions reader_options = *options.image_reader;
  reader_options.database_path = *options.database_path;
  reader_options.image_path = *options.image_path;

  if (!image_list_path.empty()) {
    reader_options.image_list = ReadTextFileLines(image_list_path);
    if (reader_options.image_list.empty()) {
      return EXIT_SUCCESS;
    }
  }

  if (!ExistsCameraModelWithName(options.image_reader->camera_model)) {
    std::cerr << "ERROR: Camera model does not exist" << std::endl;
  }

  if (!VerifyCameraParams(options.image_reader->camera_model,
                          options.image_reader->camera_params)) {
    return EXIT_FAILURE;
  }

  std::unique_ptr<QApplication> app;
  if (options.sift_extraction->use_gpu && kUseOpenGL) {
    app.reset(new QApplication(argc, argv));
  }

  SiftFeatureExtractor feature_extractor(reader_options,
                                         *options.sift_extraction);

  if (options.sift_extraction->use_gpu && kUseOpenGL) {
    RunThreadWithOpenGLContext(&feature_extractor);
  } else {
    feature_extractor.Start();
    feature_extractor.Wait();
  }

  return EXIT_SUCCESS;
}

int RunImageFilterer(int argc, char** argv) {
  std::string input_path;
  std::string output_path;
  double min_focal_length_ratio = 0.1;
  double max_focal_length_ratio = 10.0;
  double max_extra_param = 100.0;
  size_t min_num_observations = 10;

  OptionManager options;
  options.AddRequiredOption("input_path", &input_path);
  options.AddRequiredOption("output_path", &output_path);
  options.AddDefaultOption("min_focal_length_ratio", &min_focal_length_ratio);
  options.AddDefaultOption("max_focal_length_ratio", &max_focal_length_ratio);
  options.AddDefaultOption("max_extra_param", &max_extra_param);
  options.AddDefaultOption("min_num_observations", &min_num_observations);
  options.Parse(argc, argv);

  Reconstruction reconstruction;
  reconstruction.Read(input_path);

  const size_t num_reg_images = reconstruction.NumRegImages();

  reconstruction.FilterImages(min_focal_length_ratio, max_focal_length_ratio,
                              max_extra_param);

  std::vector<image_t> filtered_image_ids;
  for (const auto& image : reconstruction.Images()) {
    if (image.second.IsRegistered() &&
        image.second.NumPoints3D() < min_num_observations) {
      filtered_image_ids.push_back(image.first);
    }
  }

  for (const auto image_id : filtered_image_ids) {
    reconstruction.DeRegisterImage(image_id);
  }

  const size_t num_filtered_images =
      num_reg_images - reconstruction.NumRegImages();

  std::cout << StringPrintf("Filtered %d images from a total of %d images",
                            num_filtered_images, num_reg_images)
            << std::endl;

  reconstruction.Write(output_path);

  return EXIT_SUCCESS;
}

int RunMapper(int argc, char** argv) {
  std::string input_path;
  std::string output_path;
  std::string image_list_path;

  OptionManager options;
  options.AddDatabaseOptions();
  options.AddImageOptions();
  options.AddDefaultOption("input_path", &input_path);
  options.AddRequiredOption("output_path", &output_path);
  options.AddDefaultOption("image_list_path", &image_list_path);
  options.AddMapperOptions();
  options.Parse(argc, argv);

  if (!ExistsDir(output_path)) {
    std::cerr << "ERROR: `output_path` is not a directory." << std::endl;
    return EXIT_FAILURE;
  }

  if (!image_list_path.empty()) {
    const auto image_names = ReadTextFileLines(image_list_path);
    options.mapper->image_names =
        std::unordered_set<std::string>(image_names.begin(), image_names.end());
  }

  ReconstructionManager reconstruction_manager;
  if (input_path != "") {
    if (!ExistsDir(input_path)) {
      std::cerr << "ERROR: `input_path` is not a directory." << std::endl;
      return EXIT_FAILURE;
    }
    reconstruction_manager.Read(input_path);
  }

  IncrementalMapperController mapper(options.mapper.get(), *options.image_path,
                                     *options.database_path,
                                     &reconstruction_manager);

  // In case a new reconstruction is started, write results of individual sub-
  // models to as their reconstruction finishes instead of writing all results
  // after all reconstructions finished.
  size_t prev_num_reconstructions = 0;
  if (input_path == "") {
    mapper.AddCallback(
        IncrementalMapperController::LAST_IMAGE_REG_CALLBACK, [&]() {
          // If the number of reconstructions has not changed, the last model
          // was discarded for some reason.
          if (reconstruction_manager.Size() > prev_num_reconstructions) {
            const std::string reconstruction_path = JoinPaths(
                output_path, std::to_string(prev_num_reconstructions));
            const auto& reconstruction =
                reconstruction_manager.Get(prev_num_reconstructions);
            CreateDirIfNotExists(reconstruction_path);
            reconstruction.Write(reconstruction_path);
            options.Write(JoinPaths(reconstruction_path, "project.ini"));
            prev_num_reconstructions = reconstruction_manager.Size();
          }
        });
  }


  mapper.Start();
  mapper.Wait();

  // In case the reconstruction is continued from an existing reconstruction, do
  // not create sub-folders but directly write the results.
  if (input_path != "" && reconstruction_manager.Size() > 0) {
    reconstruction_manager.Get(0).Write(output_path);
  }

  return EXIT_SUCCESS;
}

int RunSequentialMatcher(int argc, char** argv) {
  OptionManager options;
  options.AddDatabaseOptions();
  options.AddSequentialMatchingOptions();
  options.Parse(argc, argv);

  std::unique_ptr<QApplication> app;
  if (options.sift_matching->use_gpu && kUseOpenGL) {
    app.reset(new QApplication(argc, argv));
  }

  SequentialFeatureMatcher feature_matcher(*options.sequential_matching,
                                           *options.sift_matching,
                                           *options.database_path);

  if (options.sift_matching->use_gpu && kUseOpenGL) {
    RunThreadWithOpenGLContext(&feature_matcher);
  } else {
    feature_matcher.Start();
    feature_matcher.Wait();
  }

  return EXIT_SUCCESS;
}

typedef std::function<int(int, char**)> command_func_t;

int ShowHelp(
    const std::vector<std::pair<std::string, command_func_t>>& commands) {
  std::cout << StringPrintf(
                   "%s -- Structure-from-Motion and Multi-View Stereo\n"
                   "              (%s)",
                   GetVersionInfo().c_str(), GetBuildInfo().c_str())
            << std::endl
            << std::endl;

  std::cout << "Usage:" << std::endl;
  std::cout << "  colmap [command] [options]" << std::endl << std::endl;

  std::cout << "Documentation:" << std::endl;
  std::cout << "  https://colmap.github.io/" << std::endl << std::endl;

  std::cout << "Example usage:" << std::endl;
  std::cout << "  colmap help [ -h, --help ]" << std::endl;
  std::cout << "  colmap gui" << std::endl;
  std::cout << "  colmap gui -h [ --help ]" << std::endl;
  std::cout << "  colmap automatic_reconstructor -h [ --help ]" << std::endl;
  std::cout << "  colmap automatic_reconstructor --image_path IMAGES "
               "--workspace_path WORKSPACE"
            << std::endl;
  std::cout << "  colmap feature_extractor --image_path IMAGES --database_path "
               "DATABASE"
            << std::endl;
  std::cout << "  colmap exhaustive_matcher --database_path DATABASE"
            << std::endl;
  std::cout << "  colmap mapper --image_path IMAGES --database_path DATABASE "
               "--output_path MODEL"
            << std::endl;
  std::cout << "  ..." << std::endl << std::endl;

  std::cout << "Available commands:" << std::endl;
  std::cout << "  help" << std::endl;
  for (const auto& command : commands) {
    std::cout << "  " << command.first << std::endl;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

int LineInitializer(int argc, char** argv) {

  std::string database_path;
  std::string gravity_path;
  std::string model_output_path;
  IncrementalMapper::Options incremental_mapper_options;
  OptionManager options(false);
  options.AddRequiredOption("database_path", &database_path);
  options.AddRequiredOption("gravity_path", &gravity_path);
  options.AddRequiredOption("model_output_path", &model_output_path);
  options.AddDefaultOption("max_reprojection_error", &incremental_mapper_options.init_max_error);
  options.AddDefaultOption("min_triangulation_angle", &incremental_mapper_options.init_min_tri_angle);
  options.AddDefaultOption("min_num_inliers", &incremental_mapper_options.init_min_num_inliers);
  options.Parse(argc, argv);

  // Read the gravity.txt file
  EIGEN_STL_UMAP(std::string, Eigen::Vector3d) gravity;
  {
    std::ifstream gravity_stream(gravity_path);
    CHECK(gravity_stream.is_open()) << gravity_path;

    std::string line;
    while (std::getline(gravity_stream, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }

      std::stringstream line_stream(line);
      image_t image_id;
      std::string image_name;
      Eigen::Vector3d gravity_dir;
      line_stream >> image_id >> image_name >> gravity_dir(0) >> gravity_dir(1) >> gravity_dir(2);

      gravity.emplace(image_name, gravity_dir);
    }
  }

  // Make sure the database already exists.
  // If not, this would just create a new one and lead to hard to detect
  // errors.
  CHECK(ExistsFile(database_path)) << "Database " << database_path << " does not exist.";

  Database database(database_path, true);

  EIGEN_STL_UMAP(image_t, Image) images;
  {
    // Get the images with lines from the database
    std::vector<Image> tmp_images = database.ReadAllImages();
    for (const auto& image : tmp_images) {
      images.emplace(image.ImageId(), image);
    }
  }

  CHECK_GT(images.size(), 0);

  for (auto& image : images) {
    image.second.SetLines(database.ReadFeatureLines(image.first));
  }

  std::unordered_map<image_t, std::unordered_set<int>> aligned_lines;

  // Find aligned features
  for (const auto& image : images) {
    const FeatureLines& image_lines = image.second.Lines();
    const int num_lines = image_lines.size();
    for (int i = 0; i < num_lines; ++i) {
      if (image_lines.at(i).IsAligned()) {
        aligned_lines[image.first].emplace(i);
      }
    }
  }

  // We only really need the images with aligned lines here. Create a set
  // of the corresponding image names and pass it to the database cache when
  // loading.
  // For random initialization, this will not really help much, but it also
  // shouldn't hurt.
  std::unordered_set<std::string> aligned_image_names;
  for (const auto& im_lines : aligned_lines) {
    const std::string& image_name = images.at(im_lines.first).Name();
    CHECK(gravity.count(image_name) > 0);
    aligned_image_names.emplace(image_name);
  }

  // Get the track from aligned features
  DatabaseCache db_cache;
  db_cache.Load(database, 0, false, aligned_image_names);

  const CorrespondenceGraph& corr_graph = db_cache.CorrespondenceGraph();

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
  for (const auto& im_id : aligned_lines) {
    image_ids.emplace_back(im_id.first);
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

    const Image& image = images.at(image_id);
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
        if (images.at(corr.image_id).Line(corr.line_idx).IsAligned() ==
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
    return EXIT_FAILURE;
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
          gravity.at(images.at(im_id_idx.first).Name());
    }

    // First, collect all the aligned lines
    const auto& aligned_tracks = all_aligned_tracks.at(init_image_set);
    for (const auto& track : aligned_tracks) {
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const int line_idx = track.at(i);
        lines.at(image_idx_map.at(image_id))
            .emplace_back(images.at(image_id).Line(line_idx));
        CHECK(images.at(image_id).Line(line_idx).IsAligned());
      }
    }

    // Collect all unaligned lines
    const auto& unaligned_tracks = all_unaligned_tracks.at(init_image_set);
    for (const auto& track : unaligned_tracks) {
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const int line_idx = track.at(i);
        lines.at(image_idx_map.at(image_id))
            .emplace_back(images.at(image_id).Line(line_idx));
        CHECK(!images.at(image_id).Line(line_idx).IsAligned());
      }
    }



    std::vector<init::Pose> poses;
    double inlier_ratio = 0;
    init::InitOptions init_options;
    init_options.min_num_inliers = incremental_mapper_options.init_min_num_inliers;
    init_options.min_tri_angle = incremental_mapper_options.init_min_tri_angle;
    init_options.max_error = [&]() {
      double min_error = std::numeric_limits<double>::max();
      for (int i = 0; i < 4; ++i) {
        const image_t image_id = init_image_set.at(i);
        const camera_t camera_id = database.ReadImage(image_id).CameraId();
        const Camera& camera = database.ReadCamera(camera_id);
        min_error = std::min(min_error, camera.ImageToWorldThreshold(incremental_mapper_options.init_max_error));
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
      }
    } else {
      std::cout << "Could not estimate poses\n";
    }
  }

  if (best_poses.size() == 0) {
    std::cerr << "Could not estimate initial image poses\n";
    return EXIT_FAILURE;
  }


  Reconstruction reconstruction;
  reconstruction.SetUp(&corr_graph);

  for (const auto& im_id_idx : best_image_idx_map) {
    const image_t image_id = im_id_idx.first;
    Image& image = images.at(image_id);
    init::Pose pose = best_poses.at(best_image_idx_map.at(image_id));
    image.SetQvec(RotationMatrixToQuaternion(pose.leftCols<3>()));
    image.SetTvec(pose.rightCols<1>());

    if (!reconstruction.ExistsCamera(image.CameraId())) {
      reconstruction.AddCamera(database.ReadCamera(image.CameraId()));
    }
    reconstruction.AddImage(image);
    reconstruction.RegisterImage(image_id);
  }

  reconstruction.WriteText(model_output_path);

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  InitializeGlog(argv);

  std::vector<std::pair<std::string, command_func_t>> commands;
  commands.emplace_back("gui", &RunGraphicalUserInterface);
  commands.emplace_back("automatic_reconstructor", &RunAutomaticReconstructor);
  commands.emplace_back("bundle_adjuster", &RunBundleAdjuster);
  commands.emplace_back("database_creator", &RunDatabaseCreator);
  commands.emplace_back("exhaustive_matcher", &RunExhaustiveMatcher);
  commands.emplace_back("feature_extractor", &RunFeatureExtractor);
  commands.emplace_back("image_filterer", &RunImageFilterer);
  commands.emplace_back("mapper", &RunMapper);
  commands.emplace_back("project_generator", &RunProjectGenerator);
  commands.emplace_back("sequential_matcher", &RunSequentialMatcher);
  commands.emplace_back("line_initializer", &LineInitializer);

  if (argc == 1) {
    return ShowHelp(commands);
  }

  const std::string command = argv[1];
  if (command == "help" || command == "-h" || command == "--help") {
    return ShowHelp(commands);
  } else {
    command_func_t matched_command_func = nullptr;
    for (const auto& command_func : commands) {
      if (command == command_func.first) {
        matched_command_func = command_func.second;
        break;
      }
    }
    if (matched_command_func == nullptr) {
      std::cerr << StringPrintf(
                       "ERROR: Command `%s` not recognized. To list the "
                       "available commands, run `colmap help`.",
                       command.c_str())
                << std::endl;
      return EXIT_FAILURE;
    } else {
      int command_argc = argc - 1;
      char** command_argv = &argv[1];
      command_argv[0] = argv[0];
      return matched_command_func(command_argc, command_argv);
    }
  }

  return ShowHelp(commands);
}
