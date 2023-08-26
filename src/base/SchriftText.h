///////////////////////////////////////////////////////////////////////////////
///
///	\file    SchriftText.h
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

#ifndef _SCHRIFTTEXT_H_
#define _SCHRIFTTEXT_H_

#include "schrift.h"
#include "RGBA.h"

#include <string>

///////////////////////////////////////////////////////////////////////////////

class SchriftText {

public:
	///	<summary>
	///		Number of bytes per pixel in image buffer.
	///	</summary>
	static const int NDIM = 4;

public:
	///	<summary>
	///		Alignment of text.
	///	</summary>
	enum TextAlignment {
		TextAlignment_Left,
		TextAlignment_Right,
		TextAlignment_Center
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SchriftText(
		const std::string & strFontFile,
		double dFontSize
	);

public:
	///	<summary>
	///		Calculate the minimum image buffer size for holding the given string.
	///	</summary>
	void CalculateStringMinImageBufferSize(
		const std::string & str,
		int & nMinWidth,
		int & nMinHeight,
		int & nBaseline
	);

	///	<summary>
	///		Draw the specified character at the specified location on the image,
	///		using font information from m_sft.  The coordinate (x,y) indicates
	///		the top-left corner of the character.
	///	</summary>
	void DrawCharacter(
		unsigned char c,
		int nX,
		int nY,
		int nCanvasWidth,
		int nCanvasHeight,
		unsigned char * imagedata,
		RGBA rgba,
		int * pwidth = NULL,
		int * pheight = NULL);

	///	<summary>
	///		Draw the specified string at the specified location on the image.
	///		using font information from m_sft.  The coordinate (x,y) indicates
	///		the top-left corner of the string.
	///	</summary>
	void DrawString(
		const std::string & str,
		int nX,
		int nY,
		TextAlignment eAlign,
		int nCanvasWidth,
		int nCanvasHeight,
		unsigned char * imagedata,
		RGBA rgba,
		int * pwidth = NULL,
		int * pheight = NULL);

private:
	///	<summary>
	///		Font information.
	///	</summary>
	SFT m_sft;
};

///////////////////////////////////////////////////////////////////////////////

#endif // _SCHRIFTTEXT_H_

