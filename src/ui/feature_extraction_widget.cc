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

#include "ui/feature_extraction_widget.h"

#include "base/camera_models.h"
#include "feature/extraction.h"
#include "ui/options_widget.h"
#include "ui/qt_utils.h"
#include "ui/thread_control_widget.h"

namespace colmap {

class ExtractionWidget : public OptionsWidget {
 public:
  ExtractionWidget(QWidget* parent, OptionManager* options);

  virtual void Run() = 0;

 protected:
  OptionManager* options_;
  ThreadControlWidget* thread_control_widget_;
};

class SIFTExtractionWidget : public ExtractionWidget {
 public:
  SIFTExtractionWidget(QWidget* parent, OptionManager* options);

  void Run() override;
};

ExtractionWidget::ExtractionWidget(QWidget* parent, OptionManager* options)
    : OptionsWidget(parent),
      options_(options),
      thread_control_widget_(new ThreadControlWidget(this)) {}

SIFTExtractionWidget::SIFTExtractionWidget(QWidget* parent,
                                           OptionManager* options)
    : ExtractionWidget(parent, options) {
  AddOptionDirPath(&options->image_reader->mask_path, "mask_path");
  AddOptionFilePath(&options->image_reader->camera_mask_path,
                    "camera_mask_path");

  AddOptionInt(&options->sift_extraction->max_image_size, "max_image_size");
  AddOptionInt(&options->sift_extraction->max_num_features, "max_num_features");
  AddOptionInt(&options->sift_extraction->first_octave, "first_octave", -5);
  AddOptionInt(&options->sift_extraction->num_octaves, "num_octaves");
  AddOptionInt(&options->sift_extraction->octave_resolution,
               "octave_resolution");
  AddOptionDouble(&options->sift_extraction->peak_threshold, "peak_threshold",
                  0.0, 1e7, 0.00001, 5);
  AddOptionDouble(&options->sift_extraction->edge_threshold, "edge_threshold");
  AddOptionBool(&options->sift_extraction->estimate_affine_shape,
                "estimate_affine_shape");
  AddOptionInt(&options->sift_extraction->max_num_orientations,
               "max_num_orientations");
  AddOptionBool(&options->sift_extraction->upright, "upright");
  AddOptionBool(&options->sift_extraction->domain_size_pooling,
                "domain_size_pooling");
  AddOptionDouble(&options->sift_extraction->dsp_min_scale, "dsp_min_scale",
                  0.0, 1e7, 0.00001, 5);
  AddOptionDouble(&options->sift_extraction->dsp_max_scale, "dsp_max_scale",
                  0.0, 1e7, 0.00001, 5);
  AddOptionInt(&options->sift_extraction->dsp_num_scales, "dsp_num_scales", 1);

  AddOptionInt(&options->sift_extraction->num_threads, "num_threads", -1);
  AddOptionBool(&options->sift_extraction->use_gpu, "use_gpu");
  AddOptionText(&options->sift_extraction->gpu_index, "gpu_index");
}

void SIFTExtractionWidget::Run() {
  WriteOptions();

  ImageReaderOptions reader_options = *options_->image_reader;
  reader_options.database_path = *options_->database_path;
  reader_options.image_path = *options_->image_path;

  Thread* extractor =
      new SiftFeatureExtractor(reader_options, *options_->sift_extraction);
  thread_control_widget_->StartThread("Extracting...", true, extractor);
}

FeatureExtractionWidget::FeatureExtractionWidget(QWidget* parent,
                                                 OptionManager* options)
    : parent_(parent), options_(options) {
  // Do not change flag, to make sure feature database is not accessed from
  // multiple threads
  setWindowFlags(Qt::Window);
  setWindowTitle("Feature extraction");

  QGridLayout* grid = new QGridLayout(this);

  grid->addWidget(CreateCameraModelBox(), 0, 0);

  tab_widget_ = new QTabWidget(this);

  QScrollArea* extraction_widget = new QScrollArea(this);
  extraction_widget->setAlignment(Qt::AlignHCenter);
  extraction_widget->setWidget(new SIFTExtractionWidget(this, options));
  tab_widget_->addTab(extraction_widget, tr("Extract"));

  grid->addWidget(tab_widget_);

  QPushButton* extract_button = new QPushButton(tr("Extract"), this);
  connect(extract_button, &QPushButton::released, this,
          &FeatureExtractionWidget::Extract);
  grid->addWidget(extract_button, grid->rowCount(), 0);
}

QGroupBox* FeatureExtractionWidget::CreateCameraModelBox() {
  camera_model_ids_.clear();

  camera_model_cb_ = new QComboBox(this);

#define CAMERA_MODEL_CASE(CameraModel)                                     \
  camera_model_cb_->addItem(                                               \
      QString::fromStdString(CameraModelIdToName(CameraModel::model_id))); \
  camera_model_ids_.push_back(static_cast<int>(CameraModel::model_id));

  CAMERA_MODEL_CASES

#undef CAMERA_MODEL_CASE

  camera_params_heading_ = new QLabel(tr("Default parameters (optional)"), this);

  camera_params_info_ = new QLabel(tr(""), this);
  QPalette pal = QPalette(camera_params_info_->palette());
  pal.setColor(QPalette::WindowText, QColor(130, 130, 130));
  camera_params_info_->setPalette(pal);

  camera_params_text_ = new QLineEdit(this);
  camera_params_text_->setEnabled(true);

  QGroupBox* box = new QGroupBox(tr("Camera model"), this);

  QVBoxLayout* vbox = new QVBoxLayout(box);
  vbox->addWidget(camera_model_cb_);
  vbox->addWidget(camera_params_info_);

  vbox->addWidget(camera_params_heading_);
  vbox->addWidget(camera_params_text_);
  vbox->addStretch(1);

  box->setLayout(vbox);

  SelectCameraModel(camera_model_cb_->currentIndex());

  connect(camera_model_cb_,
          (void (QComboBox::*)(int)) & QComboBox::currentIndexChanged, this,
          &FeatureExtractionWidget::SelectCameraModel);

  return box;
}

void FeatureExtractionWidget::showEvent(QShowEvent* event) {
  parent_->setDisabled(true);
  ReadOptions();
}

void FeatureExtractionWidget::hideEvent(QHideEvent* event) {
  parent_->setEnabled(true);
  WriteOptions();
}

void FeatureExtractionWidget::ReadOptions() {
  const auto camera_code =
      CameraModelNameToId(options_->image_reader->camera_model);
  for (size_t i = 0; i < camera_model_ids_.size(); ++i) {
    if (camera_model_ids_[i] == camera_code) {
      SelectCameraModel(i);
      camera_model_cb_->setCurrentIndex(i);
      break;
    }
  }
  camera_params_text_->setText(
      QString::fromStdString(options_->image_reader->camera_params));
}

void FeatureExtractionWidget::WriteOptions() {
  options_->image_reader->camera_model =
      CameraModelIdToName(camera_model_ids_[camera_model_cb_->currentIndex()]);
  options_->image_reader->camera_params =
      camera_params_text_->text().toUtf8().constData();
}

void FeatureExtractionWidget::SelectCameraModel(const int idx) {
  const int code = camera_model_ids_[idx];
  camera_params_info_->setText(QString::fromStdString(StringPrintf(
      "<small>Parameters: %s</small>", CameraModelParamsInfo(code).c_str())));
}

void FeatureExtractionWidget::Extract() {
  WriteOptions();

  if (!ExistsCameraModelWithName(options_->image_reader->camera_model)) {
    QMessageBox::critical(this, "", tr("Camera model does not exist"));
    return;
  }

  const std::vector<double> camera_params =
      CSVToVector<double>(options_->image_reader->camera_params);
  const auto camera_code =
      CameraModelNameToId(options_->image_reader->camera_model);

  if (camera_params.size() > 0 &&
      !CameraModelVerifyParams(camera_code, camera_params)) {
    QMessageBox::critical(this, "", tr("Invalid camera parameters"));
    return;
  }

  QWidget* widget =
      static_cast<QScrollArea*>(tab_widget_->currentWidget())->widget();
  static_cast<ExtractionWidget*>(widget)->Run();
}

}  // namespace colmap
