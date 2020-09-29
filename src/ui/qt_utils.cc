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

#include "ui/qt_utils.h"

#include "base/camera_models.h"
#include "util/misc.h"

namespace colmap {

Eigen::Matrix4f QMatrixToEigen(const QMatrix4x4& matrix) {
  Eigen::Matrix4f eigen;
  for (size_t r = 0; r < 4; ++r) {
    for (size_t c = 0; c < 4; ++c) {
      eigen(r, c) = matrix(r, c);
    }
  }
  return eigen;
}

QMatrix4x4 EigenToQMatrix(const Eigen::Matrix4f& matrix) {
  QMatrix4x4 qt;
  for (size_t r = 0; r < 4; ++r) {
    for (size_t c = 0; c < 4; ++c) {
      qt(r, c) = matrix(r, c);
    }
  }
  return qt;
}

QImage BitmapToQImageRGB(const Bitmap& bitmap) {
  QImage image(bitmap.Width(), bitmap.Height(), QImage::Format_RGB32);
  for (int y = 0; y < image.height(); ++y) {
    QRgb* image_line = (QRgb*)image.scanLine(y);
    for (int x = 0; x < image.width(); ++x) {
      BitmapColor<uint8_t> color;
      if (bitmap.GetPixel(x, y, &color)) {
        image_line[x] = qRgba(color.r, color.g, color.b, 255);
      }
    }
  }
  return image;
}

QPixmap ShowImagesSideBySide(const QPixmap& image1, const QPixmap& image2) {
  QPixmap image = QPixmap(QSize(image1.width() + image2.width(),
                                std::max(image1.height(), image2.height())));

  image.fill(Qt::black);

  QPainter painter(&image);
  painter.drawImage(0, 0, image1.toImage());
  painter.drawImage(image1.width(), 0, image2.toImage());

  return image;
}

void DrawKeypoints(QPixmap* pixmap, const FeatureKeypoints& points,
                   const QColor& color, const int radius_scale) {
  if (pixmap->isNull()) {
    return;
  }

  const int pen_width = std::max(pixmap->width(), pixmap->height()) / 2048 + 1;
  const int radius = radius_scale * pen_width + (3 * pen_width) % 2;
  const float radius2 = radius / 2.0f;

  QPainter painter(pixmap);
  painter.setRenderHint(QPainter::Antialiasing);

  QPen pen;
  pen.setWidth(pen_width);
  pen.setColor(color);
  painter.setPen(pen);

  for (const auto& point : points) {
    painter.drawEllipse(point.x - radius2, point.y - radius2, radius, radius);
  }
}

void DrawLines(
        QPixmap* pixmap,
        const LineCoords& lines,
        const QColor& color) {
  if (pixmap->isNull()) {
    return;
  }

  QPainter painter(pixmap);
  painter.setRenderHint(QPainter::Antialiasing);

  // Draw lines
  const int pen_width = std::max(pixmap->width(), pixmap->height()) / 2048 + 1;

  QPen pen;
  pen.setWidth(pen_width);
  pen.setColor(color);
  painter.setPen(pen);

  for (const auto& line : lines) {
    const auto qp1 = QPoint(line.first(0), line.first(1));
    const auto qp2 = QPoint(line.second(0), line.second(1));
    painter.drawLine(qp1, qp2);
  }
}

}  // namespace colmap
