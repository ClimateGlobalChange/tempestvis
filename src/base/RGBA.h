///////////////////////////////////////////////////////////////////////////////
///
///	\file    RGBA.h
///	\author  Paul Ullrich
///	\version August 14, 2023
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

#ifndef _RGBA_H_
#define _RGBA_H_

#include <cstdio>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class holding a RGBA color.
///	</summary>
class RGBA {
	public:
		///	<summary>
		///		Default constructor.
		///	</summary>
		RGBA() {
			color = 0xFFFFFFFFu;
		}

		///	<summary>
		///		Constructor from integers.
		///	</summary>
		RGBA(
			unsigned int iR,
			unsigned int iG,
			unsigned int iB,
			unsigned int iA = 255
		) {
			iR = iR % 256;
			iG = iG % 256;
			iB = iB % 256;
			iA = iA % 256;

			color = (iR << 24) | (iG << 16) | (iB << 8) | iA;
		}

		///	<summary>
		///		Red accessor.
		///	</summary>
		inline unsigned char r() const {
			return static_cast<unsigned char>((color & 0xFF000000u) >> 24);
		}

		///	<summary>
		///		Green accessor.
		///	</summary>
		inline unsigned char g() const {
			return static_cast<unsigned char>((color & 0x00FF0000u) >> 16);
		}

		///	<summary>
		///		Blue accessor.
		///	</summary>
		inline unsigned char b() const {
			return static_cast<unsigned char>((color & 0x0000FF00u) >> 8);
		}

		///	<summary>
		///		Alpha accessor.
		///	</summary>
		inline unsigned char a() const {
			return static_cast<unsigned char>(color & 0x000000FFu);
		}

		///	<summary>
		///		Red assignment.
		///	</summary>
		inline void r(unsigned char x) {
			color = (color - (color & 0xFF000000u)) | (static_cast<unsigned int>(x) << 24);
		}

		///	<summary>
		///		Green assignment.
		///	</summary>
		inline void g(unsigned char x) {
			color = (color - (color & 0x00FF0000u)) | (static_cast<uint32_t>(x) << 16);
		}

		///	<summary>
		///		Blue assignment.
		///	</summary>
		inline void b(unsigned char x) {
			color = (color - (color & 0x0000FF00u)) | (static_cast<uint32_t>(x) << 8);
		}

		///	<summary>
		///		Alpha assignment.
		///	</summary>
		inline void a(unsigned char x) {
			color = (color - (color & 0x000000FFu)) | static_cast<uint32_t>(x);
		}

	private:
		///	<summary>
		///		Structure holding the color.
		///	</summary>
		uint32_t color;
};

///////////////////////////////////////////////////////////////////////////////

#endif // _RGBA_H_

