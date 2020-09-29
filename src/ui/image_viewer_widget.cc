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

#include "ui/image_viewer_widget.h"

#include "ui/model_viewer_widget.h"
#include "util/misc.h"

namespace colmap {

const double ImageViewerWidget::kZoomFactor = 1.20;

ImageViewerGraphicsScene::ImageViewerGraphicsScene() {
  setSceneRect(0, 0, 0, 0);
  image_pixmap_item_ = addPixmap(QPixmap::fromImage(QImage()));
  image_pixmap_item_->setZValue(-1);
}

QGraphicsPixmapItem* ImageViewerGraphicsScene::ImagePixmapItem() const {
  return image_pixmap_item_;
}

ImageViewerWidget::ImageViewerWidget(QWidget* parent) : QWidget(parent) {
  setWindowFlags(Qt::Window | Qt::WindowTitleHint |
                 Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint |
                 Qt::WindowCloseButtonHint);

  resize(parent->width() - 20, parent->height() - 20);

  QFont font;
  font.setPointSize(10);
  setFont(font);

  grid_layout_ = new QGridLayout(this);
  grid_layout_->setContentsMargins(5, 5, 5, 5);

  graphics_view_ = new QGraphicsView();
  graphics_view_->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

  graphics_view_->setScene(&graphics_scene_);
  graphics_view_->setAlignment(Qt::AlignLeft | Qt::AlignTop);

  grid_layout_->addWidget(graphics_view_, 1, 0);

  button_layout_ = new QHBoxLayout();

  QPushButton* zoom_in_button = new QPushButton("+", this);
  zoom_in_button->setFont(font);
  zoom_in_button->setFixedWidth(50);
  button_layout_->addWidget(zoom_in_button);
  connect(zoom_in_button, &QPushButton::released, this,
          &ImageViewerWidget::ZoomIn);

  QPushButton* zoom_out_button = new QPushButton("-", this);
  zoom_out_button->setFont(font);
  zoom_out_button->setFixedWidth(50);
  button_layout_->addWidget(zoom_out_button);
  connect(zoom_out_button, &QPushButton::released, this,
          &ImageViewerWidget::ZoomOut);

  QPushButton* save_button = new QPushButton("Save image", this);
  save_button->setFont(font);
  button_layout_->addWidget(save_button);
  connect(save_button, &QPushButton::released, this, &ImageViewerWidget::Save);

  grid_layout_->addLayout(button_layout_, 2, 0, Qt::AlignRight);
}

void ImageViewerWidget::resizeEvent(QResizeEvent* event) {
  QWidget::resizeEvent(event);

  graphics_view_->fitInView(graphics_scene_.sceneRect(), Qt::KeepAspectRatio);
}

void ImageViewerWidget::closeEvent(QCloseEvent* event) {
  graphics_scene_.ImagePixmapItem()->setPixmap(QPixmap());
}

void ImageViewerWidget::ShowBitmap(const Bitmap& bitmap) {
  ShowPixmap(QPixmap::fromImage(BitmapToQImageRGB(bitmap)));
}

void ImageViewerWidget::ShowPixmap(const QPixmap& pixmap) {
  graphics_scene_.ImagePixmapItem()->setPixmap(pixmap);
  graphics_scene_.setSceneRect(pixmap.rect());

  show();
  graphics_view_->fitInView(graphics_scene_.sceneRect(), Qt::KeepAspectRatio);

  raise();
}

void ImageViewerWidget::ZoomIn() {
  graphics_view_->scale(kZoomFactor, kZoomFactor);
}

void ImageViewerWidget::ZoomOut() {
  graphics_view_->scale(1.0 / kZoomFactor, 1.0 / kZoomFactor);
}

void ImageViewerWidget::Save() {
  QString filter("PNG (*.png)");
  const QString save_path =
      QFileDialog::getSaveFileName(this, tr("Select destination..."), "",
                                   "PNG (*.png);;JPEG (*.jpg);;BMP (*.bmp)",
                                   &filter)
          .toUtf8()
          .constData();

  // Selection canceled?
  if (save_path == "") {
    return;
  }

  graphics_scene_.ImagePixmapItem()->pixmap().save(save_path);
}

FeatureImageViewerWidget::FeatureImageViewerWidget(
    QWidget* parent, const std::string& switch_text)
    : ImageViewerWidget(parent),
      switch_state_(true),
      switch_text_(switch_text) {
  switch_button_ = new QPushButton(tr(("Hide " + switch_text_).c_str()), this);
  switch_button_->setFont(font());
  button_layout_->addWidget(switch_button_);
  connect(switch_button_, &QPushButton::released, this,
          &FeatureImageViewerWidget::ShowOrHide);
}

void FeatureImageViewerWidget::ReadAndShowWithLines(
    const std::string& path,
    const LineCoords& lines,
    const std::vector<char>& tri_mask) {
  Bitmap bitmap;
  if (!bitmap.Read(path, true)) {
    std::cerr << "ERROR: Cannot read image at path " << path << std::endl;
  }

  image1_ = QPixmap::fromImage(BitmapToQImageRGB(bitmap));
  image2_ = image1_;

  const size_t num_tri_lines = std::count_if(
      tri_mask.begin(), tri_mask.end(), [](const bool tri) { return tri; });

  LineCoords lines_tri(num_tri_lines);
  LineCoords lines_not_tri(lines.size() - num_tri_lines);

  size_t i_tri = 0;
  size_t i_not_tri = 0;
  for (size_t i = 0; i < tri_mask.size(); ++i) {
    if (tri_mask[i]) {
      lines_tri[i_tri] = lines.at(i);
      i_tri += 1;
    } else {
      lines_not_tri[i_not_tri] = lines.at(i);
      i_not_tri += 1;
    }
  }

  DrawLines(&image2_, lines_tri, Qt::magenta);
//  DrawLines(&image2_, lines_not_tri, Qt::green);

  if (switch_state_) {
    ShowPixmap(image2_);
  } else {
    ShowPixmap(image1_);
  }
}

void FeatureImageViewerWidget::ReadAndShowWithLinesAndProjections(
    const std::string& path,
    const LineCoords& lines,
    const std::vector<Eigen::Vector2d>& projected_points,
    const std::vector<char>& tri_mask) {
  Bitmap bitmap;
  if (!bitmap.Read(path, true)) {
    std::cerr << "ERROR: Cannot read image at path " << path << std::endl;
  }

  image1_ = QPixmap::fromImage(BitmapToQImageRGB(bitmap));
  image2_ = image1_;

  const size_t num_tri_keypoints = std::count_if(
      tri_mask.begin(), tri_mask.end(), [](const bool tri) { return tri; });

  FeatureKeypoints tri_projected(num_tri_keypoints);
  LineCoords lines_tri(num_tri_keypoints);
  LineCoords lines_not_tri(lines.size() - num_tri_keypoints);
  size_t i_tri = 0;
  size_t i_not_tri = 0;
  for (size_t i = 0; i < tri_mask.size(); ++i) {
    if (tri_mask[i]) {
      lines_tri[i_tri] = lines.at(i);
      tri_projected[i_tri].x = projected_points.at(i)(0);
      tri_projected[i_tri].y = projected_points.at(i)(1);
      i_tri += 1;
    } else {
      lines_not_tri[i_not_tri] = lines.at(i);
      i_not_tri += 1;
    }
  }

  DrawLines(&image2_, lines_tri, Qt::magenta);

  if (switch_state_) {
    ShowPixmap(image2_);
  } else {
    ShowPixmap(image1_);
  }
}

void FeatureImageViewerWidget::ShowOrHide() {
  if (switch_state_) {
    switch_button_->setText(std::string("Show " + switch_text_).c_str());
    ShowPixmap(image1_);
    switch_state_ = false;
  } else {
    switch_button_->setText(std::string("Hide " + switch_text_).c_str());
    ShowPixmap(image2_);
    switch_state_ = true;
  }
}

DatabaseImageViewerWidget::DatabaseImageViewerWidget(
    QWidget* parent, ModelViewerWidget* model_viewer_widget,
    OptionManager* options)
    : FeatureImageViewerWidget(parent, "keypoints"),
      model_viewer_widget_(model_viewer_widget),
      options_(options) {
  setWindowTitle("Image information");

  table_widget_ = new QTableWidget(this);
  table_widget_->setColumnCount(2);
  table_widget_->setRowCount(11);

  QFont font;
  font.setPointSize(10);
  table_widget_->setFont(font);

  table_widget_->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

  table_widget_->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table_widget_->setSelectionMode(QAbstractItemView::SingleSelection);
  table_widget_->setShowGrid(true);

  table_widget_->horizontalHeader()->setStretchLastSection(true);
  table_widget_->horizontalHeader()->setVisible(false);
  table_widget_->verticalHeader()->setVisible(false);
  table_widget_->verticalHeader()->setDefaultSectionSize(18);

  int table_row = 0;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("image_id"));
  image_id_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, image_id_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("camera_id"));
  camera_id_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, camera_id_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("camera_model"));
  camera_model_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, camera_model_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("camera_params"));
  camera_params_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, camera_params_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("qw, qx, qy, qz"));
  qvec_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, qvec_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("tx, ty, tz"));
  tvec_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, tvec_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("dims"));
  dimensions_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, dimensions_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("num_points2D"));
  num_points2D_item_ = new QTableWidgetItem();
  num_points2D_item_->setForeground(Qt::red);
  table_widget_->setItem(table_row, 1, num_points2D_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("num_points3D"));
  num_points3D_item_ = new QTableWidgetItem();
  num_points3D_item_->setForeground(Qt::magenta);
  table_widget_->setItem(table_row, 1, num_points3D_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0,
                         new QTableWidgetItem("num_observations"));
  num_obs_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, num_obs_item_);
  table_row += 1;

  table_widget_->setItem(table_row, 0, new QTableWidgetItem("name"));
  name_item_ = new QTableWidgetItem();
  table_widget_->setItem(table_row, 1, name_item_);
  table_row += 1;

  grid_layout_->addWidget(table_widget_, 0, 0);

  delete_button_ = new QPushButton(tr("Delete"), this);
  delete_button_->setFont(font);
  button_layout_->addWidget(delete_button_);
  connect(delete_button_, &QPushButton::released, this,
          &DatabaseImageViewerWidget::DeleteImage);
}

void DatabaseImageViewerWidget::ShowImageWithId(const image_t image_id) {
  if (model_viewer_widget_->images.count(image_id) == 0) {
    return;
  }

  image_id_ = image_id;

  const Image& image = model_viewer_widget_->images.at(image_id);
  const Camera& camera = model_viewer_widget_->cameras.at(image.CameraId());

  image_id_item_->setText(QString::number(image_id));
  camera_id_item_->setText(QString::number(image.CameraId()));
  camera_model_item_->setText(QString::fromStdString(camera.ModelName()));
  camera_params_item_->setText(QString::fromStdString(camera.ParamsToString()));
  qvec_item_->setText(QString::number(image.Qvec(0)) + ", " +
                      QString::number(image.Qvec(1)) + ", " +
                      QString::number(image.Qvec(2)) + ", " +
                      QString::number(image.Qvec(3)));
  tvec_item_->setText(QString::number(image.Tvec(0)) + ", " +
                      QString::number(image.Tvec(1)) + ", " +
                      QString::number(image.Tvec(2)));
  dimensions_item_->setText(QString::number(camera.Width()) + "x" +
                            QString::number(camera.Height()));
  num_points2D_item_->setText(QString::number(image.NumLines()));

  num_points3D_item_->setText(QString::number(image.NumPoints3D()));
  num_obs_item_->setText(QString::number(image.NumObservations()));
  name_item_->setText(QString::fromStdString(image.Name()));

  ResizeTable();

  const std::string path = JoinPaths(*options_->image_path, image.Name());

  LineCoords line_coords(image.NumLines());
  std::vector<char> tri_mask(image.NumLines());
  std::vector<Eigen::Vector2d> projected_points(image.NumLines(), Eigen::Vector2d::Zero());
  FeatureKeypoints keypoints(image.NumLines());

  for (point2D_t i = 0; i < image.NumLines(); ++i) {
    const FeatureLine& line = image.Line(i);
    const Eigen::Vector3d& line_dir = line.Line();

    const Eigen::Vector2d p0{-line_dir(2) / line_dir(0), 0};
    const Eigen::Vector2d offset{line_dir(1), -line_dir(0)};

    const Eigen::Vector2d pp0 = camera.WorldToImage(p0);
    const Eigen::Vector2d p_off = camera.WorldToImage(p0 + 0.01 * offset);

    const Eigen::Vector2d offset_dir = p_off - pp0;
    std::array<double,4> s;
    s[0] = -pp0(0) / offset_dir(0);
    s[1] = (camera.Width() - pp0(0))  / offset_dir(0);
    s[2] = -pp0(1) / offset_dir(1);
    s[3] = (camera.Height() - pp0(1)) / offset_dir(1);

    std::sort(s.begin(), s.end());
    const Eigen::Vector2d pp1 = pp0 + s[1] * offset_dir;
    const Eigen::Vector2d pp2 = pp0 + s[2] * offset_dir;

    line_coords.at(i) = LineCoord(pp1, pp2);
    tri_mask.at(i) = line.HasPoint3D();

    if (line.HasPoint3D()) {
      const point3D_t point3D_id = line.Point3DId();
      projected_points.at(i) = camera.WorldToImage(
          (image.ProjectionMatrix() *
           model_viewer_widget_->points3D.at(point3D_id).XYZ().homogeneous())
              .hnormalized());
    }

  }

  ReadAndShowWithLinesAndProjections(path, line_coords, projected_points, tri_mask);
}

void DatabaseImageViewerWidget::ResizeTable() {
  // Set fixed table dimensions.
  table_widget_->resizeColumnsToContents();
  int height = table_widget_->horizontalHeader()->height() +
               2 * table_widget_->frameWidth();
  for (int i = 0; i < table_widget_->rowCount(); i++) {
    height += table_widget_->rowHeight(i);
  }
  table_widget_->setFixedHeight(height);
}

void DatabaseImageViewerWidget::DeleteImage() {
  QMessageBox::StandardButton reply = QMessageBox::question(
      this, "", tr("Do you really want to delete this image?"),
      QMessageBox::Yes | QMessageBox::No);
  if (reply == QMessageBox::Yes) {
    if (model_viewer_widget_->reconstruction->ExistsImage(image_id_)) {
      model_viewer_widget_->reconstruction->DeRegisterImage(image_id_);
    }
    model_viewer_widget_->ReloadReconstruction();
  }
  hide();
}

}  // namespace colmap
