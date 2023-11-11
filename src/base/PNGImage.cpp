///
///	\file    PNGImage.cpp
///	\author  Paul Ullrich
///	\version August 24, 2023
///
///	<remarks>
///		Copyright 2023 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "PNGImage.h"

#include <iostream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

PNGImage::PNGImage(
	const std::string & strFilename
) {
	unsigned int error =
		lodepng::decode(*this, m_width, m_height, strFilename.c_str());

	if (error) {
		_EXCEPTION2("\n  PNG decoder error %i: %s",
			error, lodepng_error_text(error));
	}
}

///////////////////////////////////////////////////////////////////////////////

PNGImage::PNGImage(
	unsigned int width,
	unsigned int height
) {
	m_width = width;
	m_height = height;

	resize(4 * m_width * m_height, 0);
}

///////////////////////////////////////////////////////////////////////////////

bool PNGImage::write(
	const std::string & strFilename
) {
	unsigned error = lodepng::encode(strFilename, *this, m_width, m_height);
	if (error) {
		std::cout << "PNG encoder error (" << error << "): " << lodepng_error_text(error) << std::endl;
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

void PNGImage::from_subset(
	const PNGImage & img,
	int x1,
	int y1,
	int x2,
	int y2,
	unsigned int target_width,
	unsigned int target_height
) {
	// Source image dimensions
	unsigned int imgwidth = img.width();
	unsigned int imgheight = img.height();

	if (imgwidth == 0) {
		_EXCEPTIONT("In PNGImage::from_subset img.width() == 0");
	}
	if (imgheight == 0) {
		_EXCEPTIONT("In PNGImage::from_subset img.height() == 0");
	}

	// Map x1 and x2 to the range [0,imgwidth)
	x1 = (x1<0)?(x1-(x1/imgwidth-1)*imgwidth):(x1-(x1/imgwidth)*imgwidth);
	x2 = (x2<0)?(x2-(x2/imgwidth-1)*imgwidth):(x2-(x2/imgwidth)*imgwidth);
	if (x1 > x2) {
		x2 += imgwidth;
	}
	_ASSERT(x2 >= x1);
	_ASSERT(x1 >= 0);
	_ASSERT(x1 < imgwidth);

	// Map y1 and y2 to the range [0,imgheight)
	y1 = (y1<0)?(y1-(y1/imgheight-1)*imgheight):(y1-(y1/imgheight)*imgheight);
	y2 = (y2<0)?(y2-(y2/imgheight-1)*imgheight):(y2-(y2/imgheight)*imgheight);
	if (y1 > y2) {
		y2 += imgheight;
	}
	_ASSERT(y2 >= y1);
	_ASSERT(y1 >= 0);
	_ASSERT(y1 < imgheight);

	// Target image size
	if (target_width == 0) {
		m_width = static_cast<unsigned int>(x2 - x1 + 1);
	} else {
		m_width = target_width;
	}

	if (target_height == 0) {
		m_height = static_cast<unsigned int>(y2 - y1 + 1);
	} else {
		m_height = target_height;
	}

	// TODO: Special handling for (m_width == 1) or (m_height == 1)
	_ASSERT(m_width > 1);
	_ASSERT(m_height > 1);

	resize(4 * m_width * m_height);

	// Direct subset
	if ((target_width == 0) && (target_height == 0)) {
		for (int j = 0; j < m_height; j++) {
			int jx = (j + y1) % imgheight;
			for (int i = 0; i < m_width; i++) {
				int ix = (i + x1) % imgwidth;

				(*this)[4*j*m_width+4*i+0] = img[4*jx*imgwidth+4*ix+0];
				(*this)[4*j*m_width+4*i+1] = img[4*jx*imgwidth+4*ix+1];
				(*this)[4*j*m_width+4*i+2] = img[4*jx*imgwidth+4*ix+2];
				(*this)[4*j*m_width+4*i+3] = img[4*jx*imgwidth+4*ix+3];
			}
		}

	// Bilinear interpolation
	} else {
		for (int j = 0; j < m_height; j++) {
			double dB = (static_cast<double>(y1) - 0.5) +
				(static_cast<double>(j) + 0.5)
				/ static_cast<double>(m_height)
				* static_cast<double>(y2 - y1);

			int jx;
			if (dB < 0.0) {
				jx = y1;
				dB = 0.0;
			} else if (dB >= static_cast<double>(y2)) {
				jx = y2 - 1;
				dB = 1.0;
			} else {
				jx = static_cast<int>(dB);
				dB -= floor(dB);
			}

			for (int i = 0; i < m_width; i++) {
				double dA = (static_cast<double>(x1) - 0.5) +
					(static_cast<double>(i) + 0.5)
					/ static_cast<double>(m_width)
					* static_cast<double>(x2 - x1);

				int ix;
				if (dA < 0.0) {
					ix = x1;
					dA = 0.0;
				} else if (dA >= static_cast<double>(x2)) {
					ix = x2 - 1;
					dA = 1.0;
				} else {
					ix = static_cast<int>(dA);
					dA -= floor(dA);
				}

				for (int c = 0; c < 4; c++) {
					(*this)[4*j*m_width+4*i+c] = static_cast<unsigned char>(
						  (1.0 - dA) * (1.0 - dB) * static_cast<double>(img[4* jx   *imgwidth+4* ix   +c])
						+ (      dA) * (1.0 - dB) * static_cast<double>(img[4* jx   *imgwidth+4*(ix+1)+c])
						+ (1.0 - dA) * (      dB) * static_cast<double>(img[4*(jx+1)*imgwidth+4* ix   +c])
						+ (      dA) * (      dB) * static_cast<double>(img[4*(jx+1)*imgwidth+4*(ix+1)+c]));
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void PNGImage::overlay(
	const PNGImage & img
) {
	if ((width() != img.width()) || (height() != img.height())) {
		_EXCEPTIONT("Overlay PNGImage must have the same dimensions as the target");
	}

	size_t sTotalSize = width() * height();

	for (size_t i = 0; i < sTotalSize; i++) {
		unsigned int alpha = img[4*i+3] + 1;
		unsigned int inv_alpha = 256 - img[4*i+3];

		(*this)[4*i+0] = (unsigned char)((alpha * img[4*i+0] + inv_alpha * (*this)[4*i+0]) >> 8);
		(*this)[4*i+1] = (unsigned char)((alpha * img[4*i+1] + inv_alpha * (*this)[4*i+1]) >> 8);
		(*this)[4*i+2] = (unsigned char)((alpha * img[4*i+2] + inv_alpha * (*this)[4*i+2]) >> 8);
	}
}

///////////////////////////////////////////////////////////////////////////////

