///////////////////////////////////////////////////////////////////////////////
///
///	\file    PNGImage.h
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

#ifndef _PNGIMAGE_H_
#define _PNGIMAGE_H_

#include "Exception.h"
#include "lodepng.h"

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

class PNGImage : public std::vector<unsigned char> {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PNGImage() :
		m_width(0),
		m_height(0)
	{ }

	///	<summary>
	///		Constructor with file name.
	///	</summary>
	PNGImage(
		const std::string & strFilename
	);

	///	<summary>
	///		Constructor with image size.
	///	</summary>
	PNGImage(
		unsigned int width,
		unsigned int height
	);

	///	<summary>
	///		Write the PNG to a file.
	///	</summary>
	bool write(
		const std::string & strFilename
	);

	///	<summary>
	///		Subset the PNG.
	///	</summary>
	void from_subset(
		const PNGImage & img,
		int x1,
		int y1,
		int x2,
		int y2,
		unsigned int target_width = 0,
		unsigned int target_height = 0
	);

	///	<summary>
	///		Apply another PNG as an overlay.
	///	</summary>
	void overlay(
		const PNGImage & img
	);

	///	<summary>
	///		Get the image width.
	///	</summary>
	unsigned int width() const {
		return m_width;
	}

	///	<summary>
	///		Get the image height.
	///	</summary>
	unsigned int height() const {
		return m_height;
	}

private:
	///	<summary>
	///		Image width.
	///	</summary>
	unsigned int m_width;

	///	<summary>
	///		Image height.
	///	</summary>
	unsigned int m_height;

};

///////////////////////////////////////////////////////////////////////////////

#endif // _PNGIMAGE_H_

