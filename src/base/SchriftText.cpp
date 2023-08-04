///////////////////////////////////////////////////////////////////////////////
///
///	\file    SchriftText.cpp
///	\author  Paul Ullrich
///	\version July 22, 2023
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

#include "SchriftText.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

SchriftText::SchriftText(
	const std::string & strFontFile,
	double dFontSize
) {
	m_sft.xScale = dFontSize;
	m_sft.yScale = dFontSize;
	m_sft.xOffset = 0;
	m_sft.yOffset = 0;
	m_sft.flags = SFT_DOWNWARD_Y;

	m_sft.font = sft_loadfile(strFontFile.c_str());
	if (m_sft.font == NULL) {
		_EXCEPTION1("Font file \"%s\" not found", strFontFile.c_str());
	}
}

////////////////////////////////////////////////////////////////////////////////

void SchriftText::CalculateStringMinImageBufferSize(
	const std::string & str,
	int & nMinWidth,
	int & nMinHeight,
	int & nBaseline
) {
	nMinWidth = 0;
	nMinHeight = 0;
	nBaseline = 0;

	for (int i = 0; i < str.length(); i++) {
		SFT_Glyph gid;
		if (sft_lookup(&m_sft, str[i], &gid) < 0) {
			_EXCEPTION1("schrift error: character \"%c\" missing", str[i]);
			return;
		}

		SFT_GMetrics mtx;
		if (sft_gmetrics(&m_sft, gid, &mtx) < 0) {
			_EXCEPTION1("schrift error: character \"%c\" bad glyph metrics", str[i]);
			return;
		}

		if (-mtx.yOffset > nBaseline) {
			nBaseline = -mtx.yOffset;
		}
		if (mtx.minHeight > nMinHeight) {
			nMinHeight = mtx.minHeight;
		}
		if (mtx.minHeight + nBaseline > nMinHeight) {
			nMinHeight = mtx.minHeight + nBaseline;
		}
		if (i != str.length()-1) {
			nMinWidth += mtx.advanceWidth;
		} else {
			nMinWidth += mtx.minWidth;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void SchriftText::DrawCharacter(
	unsigned char c,
	int nX,
	int nY,
	int nCanvasWidth,
	int nCanvasHeight,
	unsigned char * imagedata,
	RGBA rgba,
	int * pwidth,
	int * pheight
) {
	SFT_Glyph gid;
	if (sft_lookup(&m_sft, c, &gid) < 0) {
		_EXCEPTION1("schrift error: character \"%c\" missing", c);
		return;
	}

	SFT_GMetrics mtx;
	if (sft_gmetrics(&m_sft, gid, &mtx) < 0) {
		_EXCEPTION1("schrift error: character \"%c\" bad glyph metrics", c);
		return;
	}

	SFT_Image img;
	img.width = (mtx.minWidth + 3) & ~3;
	img.height = mtx.minHeight;
	unsigned char * pixels = new unsigned char[img.width * img.height];
	img.pixels = pixels;
	if (sft_render(&m_sft, gid, img) < 0) {
		_EXCEPTION1("schrift error: character \"%c\" not rendered", c);
		return;
	}

	if (pwidth != NULL) {
		(*pwidth) = mtx.advanceWidth;
	}
	if (pheight != NULL) {
		(*pheight) = mtx.minHeight;
	}

	int s = 0;
	for (int j = 0; j < img.height; j++) {
		int jx = nY + static_cast<size_t>(j) + mtx.yOffset;
		for (int i = 0; i < img.width; i++) {
			int ix = nX + static_cast<size_t>(i) + mtx.leftSideBearing;
			float dShadingFrac = static_cast<float>(255 - pixels[s]) / 255.0;

			imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 0] =
				rgba.r() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 0] * dShadingFrac;
			imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 1] =
				rgba.g() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 1] * dShadingFrac;
			imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 2] =
				rgba.b() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 2] * dShadingFrac;
			s++;
		}
	}

	//_EXCEPTION();

	delete[] pixels;
}

////////////////////////////////////////////////////////////////////////////////

void SchriftText::DrawString(
	const std::string & str,
	int nX,
	int nY,
	TextAlignment eAlign,
	int nCanvasWidth,
	int nCanvasHeight,
	unsigned char * imagedata,
	RGBA rgba,
	int * pwidth,
	int * pheight
) {
	int cumulative_width = 0;
	int cumulative_height = 0;
	int width = 0;
	int height = 0;

	// Render left-aligned text
	if (eAlign == TextAlignment_Left) {
		for (int i = 0; i < str.length(); i++) {

			DrawCharacter(
				str[i],
				nX,
				nY,
				nCanvasWidth,
				nCanvasHeight,
				imagedata,
				rgba,
				&width,
				&height);

			nX += width;
			cumulative_width += width;
			if (height > cumulative_height) {
				cumulative_height = height;
			}
		}

	// Render right-aligned and center-aligned text using an image buffer
	} else if ((eAlign == TextAlignment_Right) || (eAlign == TextAlignment_Center)) {
		int nMinBufferWidth;
		int nMinBufferHeight;
		int nBaseline;

		CalculateStringMinImageBufferSize(str, nMinBufferWidth, nMinBufferHeight, nBaseline);

		if (nMinBufferWidth == 0) {
			if (pwidth != NULL) {
				(*pwidth) = 0;
			}
			if (pheight != NULL) {
				(*pheight) = 0;
			}
			return;
		}

		cumulative_height = nMinBufferHeight;

		unsigned char * stringbuf = new unsigned char[NDIM * nMinBufferWidth * nMinBufferHeight];
		memset(stringbuf, 255, NDIM * nMinBufferWidth * nMinBufferHeight * sizeof(unsigned char));

		int nPenX = 0;

		for (int i = 0; i < str.length(); i++) {
			DrawCharacter(
				str[i],
				nPenX,
				nBaseline,
				nMinBufferWidth,
				nMinBufferHeight,
				stringbuf,
				RGBA(),
				&width,
				&height);

			nPenX += width;
			cumulative_width += width;
		}

		int nPenBeginX = 0;
		if (eAlign == TextAlignment_Right) {
			nPenBeginX = nX - nPenX;
		} else if (eAlign == TextAlignment_Center) {
			nPenBeginX = nX - nPenX / 2;
		}

		for (int i = 0; i < nMinBufferWidth; i++) {
			int ix = nPenBeginX + i;
			for (int j = 0; j < nMinBufferHeight; j++) {
				int jx = nY - nBaseline + j;

				if (nCanvasWidth * jx + ix >= nCanvasWidth * nCanvasHeight) {
					continue;
				}
				if (j * nMinBufferWidth + i >= nMinBufferWidth * nMinBufferHeight) {
					continue;
				}
				float dShadingFrac = static_cast<float>(stringbuf[NDIM * (j * nMinBufferWidth + i) + 0]) / 255.0; 

				imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 0] =
					rgba.r() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 0] * dShadingFrac;
				imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 1] =
					rgba.g() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 1] * dShadingFrac;
				imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 2] =
					rgba.b() * (1.0 - dShadingFrac) + imagedata[NDIM * nCanvasWidth * jx + NDIM * ix + 2] * dShadingFrac;
			}
		}

		delete[] stringbuf;
	}

	if (pwidth != NULL) {
		(*pwidth) = cumulative_width;
	}
	if (pheight != NULL) {
		(*pheight) = cumulative_height;
	}
}

////////////////////////////////////////////////////////////////////////////////

