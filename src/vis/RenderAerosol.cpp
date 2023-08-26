///////////////////////////////////////////////////////////////////////////////
///
///	\file    RenderTC.cpp
///	\author  Paul Ullrich
///	\version July 20, 2023
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "SimpleGrid.h"
#include "NetCDFUtilities.h"
#include "FilenameList.h"
#include "NcFileVector.h"
#include "CoordTransforms.h"
#include "SchriftText.h"
#include "ColorMap.h"
#include "ShpFile.h"
#include "PNGImage.h"
#include "RGBA.h"

#include "lodepng.h"
#include "kdtree.h"

#include <libgen.h>
#include <set>
#include <queue>
#include <fstream>
#include <vector>

///////////////////////////////////////////////////////////////////////////////

template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
  return std::max(lower, std::min(n, upper));
}

///////////////////////////////////////////////////////////////////////////////

void RealCoordToImageCoord(
	double dLat,
	double dLon,
	double dLat0,
	double dLon0,
	double dLat1,
	double dLon1,
	int iWidth,
	int iHeight,
	int & iX,
	int & iY
) {
	double dX = (dLon - dLon0) / (dLon1 - dLon0);
	double dY = (dLat - dLat0) / (dLat1 - dLat0);

	iX = static_cast<int>(dX * (iWidth-1));
	iY = static_cast<int>(dY * (iHeight-1));

	iY = iHeight - 1 - iY;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

// Turn off fatal errors in NetCDF
NcError error(NcError::silent_nonfatal);

try {

	// Input data file
	std::string strInputData;

	// List of input data files
	std::string strInputDataList;

	// Output PNG
	std::string strOutputRoot;

	// Begin latitude of render
	double dLatBegin;

	// End latitude of render
	double dLatEnd;

	// Begin longitude of render
	double dLonBegin;

	// End longitude of render
	double dLonEnd;

	// Number of pixels output image in x direction
	int nImageSizeX;

	// Number of pixels output image in y direction
	int nImageSizeY;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Path to data files
	std::string strDataPath = std::string(dirname(argv[0])) + std::string("/../data");

	// Plot type
	std::string strPlotType;

	// Shapefile
	std::string strShapefile;

	// Shapefile 2
	std::string strShapefile2;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strOutputRoot, "out_root", "img");
		CommandLineDouble(dLatBegin, "lat_begin", -90.0);
		CommandLineDouble(dLatEnd, "lat_end", 90.0);
		CommandLineDouble(dLonBegin, "lon_begin", 0.0);
		CommandLineDouble(dLonEnd, "lon_end", 360.0);
		CommandLineInt(nImageSizeX, "imgx", 360);
		CommandLineInt(nImageSizeY, "imgy", 180);
		CommandLineString(strLatitudeName, "latname", "lat");
		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strDataPath, "datapath", strDataPath.c_str());
		CommandLineStringD(strPlotType, "plot", "COD", "(COD|ALLSKY|PRECT)");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_data) or (--in_data_list) specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list) may be specified");
	}
	if (strOutputRoot.length() == 0) {
		_EXCEPTIONT("No output file (--out_png) specified");
	}
	if (nImageSizeX < 1) {
		_EXCEPTIONT("--imgx must be a positive integer");
	}
	if (nImageSizeY < 1) {
		_EXCEPTIONT("--imgy must be a positive integer");
	}
	if ((dLatBegin < -90.0) || (dLatBegin > 90.0)) {
		_EXCEPTIONT("--lat_begin must be in the range [-90.0, 90.0]");
	}
	if ((dLatEnd < -90.0) || (dLatEnd > 90.0)) {
		_EXCEPTIONT("--lat_end must be in the range [-90.0, 90.0]");
	}
	if ((strPlotType != "COD") && (strPlotType != "ALLSKY") && (strPlotType != "PRECT")) {
		_EXCEPTIONT("--plot must be one of \"COD\", \"ALLSKY\" or \"PRECT\"");
	}

	// Load input file list
	FilenameList vecInputFiles;
	if (strInputData.length() != 0) {
		vecInputFiles.push_back(strInputData);
	} else {
		vecInputFiles.FromFile(strInputDataList, true);
	}

	// Create output image
	PNGImage imgOut(nImageSizeX, nImageSizeY);

	// Subset the blue marble
	AnnounceStartBlock("Processing BlueMarble");
	Announce("Reading image from disk");
	PNGImage imgBlueMarble(strDataPath + std::string("/BlueMarble_May2004_3km_notopo.png"));

	Announce("Rescaling image");
	PNGImage imgBlueMarbleSub;

	{
		int x1 = (dLonBegin + 180.0) / 360.0 * static_cast<double>(imgBlueMarble.width());
		int x2 = (dLonEnd + 180.0) / 360.0 * static_cast<double>(imgBlueMarble.width());
		int y1 = imgBlueMarble.height() - (dLatEnd + 90.0) / 180.0 * static_cast<double>(imgBlueMarble.height());
		int y2 = imgBlueMarble.height() - (dLatBegin + 90.0) / 180.0 * static_cast<double>(imgBlueMarble.height());

		if (y1 >= imgBlueMarble.height()) {
			y1 = imgBlueMarble.height()-1;
		}
		if (y2 >= imgBlueMarble.height()) {
			y2 = imgBlueMarble.height()-1;
		}

		imgBlueMarbleSub.from_subset(imgBlueMarble, x1, y1, x2, y2, nImageSizeX, nImageSizeY);
	}
	AnnounceEndBlock("Done");

	// Subset the black marble
	AnnounceStartBlock("Processing BlackMarble");
	Announce("Reading image from disk");
	PNGImage imgBlackMarble(strDataPath + std::string("/BlackMarble_2016_3km.png"));

	Announce("Rescaling image");
	PNGImage imgBlackMarbleSub;

	{
		int x1 = (dLonBegin + 180.0) / 360.0 * static_cast<double>(imgBlackMarble.width());
		int x2 = (dLonEnd + 180.0) / 360.0 * static_cast<double>(imgBlackMarble.width());
		int y1 = imgBlackMarble.height() - (dLatEnd + 90.0) / 180.0 * static_cast<double>(imgBlackMarble.height());
		int y2 = imgBlackMarble.height() - (dLatBegin + 90.0) / 180.0 * static_cast<double>(imgBlackMarble.height());

		if (y1 >= imgBlackMarble.height()) {
			y1 = imgBlackMarble.height()-1;
		}
		if (y2 >= imgBlackMarble.height()) {
			y2 = imgBlackMarble.height()-1;
		}

		imgBlackMarbleSub.from_subset(imgBlackMarble, x1, y1, x2, y2, nImageSizeX, nImageSizeY);
	}
	AnnounceEndBlock("Done");

	// Initialize text renderer
	SchriftText sctext(strDataPath + std::string("/Ubuntu-Regular.ttf"), 18.0);

	// Data buffer
	std::vector<float> dataCOD;
	std::vector<float> dataSOLIN;
	std::vector<float> dataALLSKY;
	std::vector<float> dataPRECT;

	// ColorMap
	ColorMapLibrary cmaplib(strDataPath);
	ColorMap cmapJet;
	cmaplib.GenerateColorMap("jet", cmapJet);

	// Shapefiles
	Mesh meshShp;
	Mesh meshShp2;
	if (strShapefile != "") {
		ReadShpFileAsMesh(strShapefile, meshShp);
	}
	if (strShapefile2 != "") {
		ReadShpFileAsMesh(strShapefile2, meshShp2);
	}

	// Map
	std::vector<size_t> vecMap(nImageSizeX * nImageSizeY);

	int iGuamIx = 0;

	// Load data longitude and latitude information
	{
		AnnounceStartBlock("Building map");

		// Open files
		NcFileVector vecFiles;
		vecFiles.ParseFromString(vecInputFiles[0]);

		// Load longitude and latitude data
		Announce("Loading latitude/longitude data");

		std::vector<double> dDataLon;
		std::vector<double> dDataLat;

		NcVar * varLon = vecFiles[0]->get_var(strLongitudeName.c_str());
		if (varLon == NULL) {
			_EXCEPTION2("No variable \"%s\" found in file \"%s\"",
				strLongitudeName.c_str(), vecFiles.GetFilename(0).c_str());
		}
		NcVar * varLat = vecFiles[0]->get_var(strLatitudeName.c_str());
		if (varLat == NULL) {
			_EXCEPTION2("No variable \"%s\" found in file \"%s\"",
				strLatitudeName.c_str(), vecFiles.GetFilename(0).c_str());
		}

		dDataLon.resize(varLon->get_dim(0)->size());
		varLon->get(&(dDataLon[0]), varLon->get_dim(0)->size());

		dDataLat.resize(varLat->get_dim(0)->size());
		varLat->get(&(dDataLat[0]), varLat->get_dim(0)->size());

		_ASSERT(dDataLon.size() == dDataLat.size());

		// Allocate data
		dataCOD.resize(dDataLon.size());
		dataSOLIN.resize(dDataLon.size());

		if (strPlotType == "ALLSKY") {
			dataALLSKY.resize(dDataLon.size());
		}
		if (strPlotType == "PRECT") {
			dataPRECT.resize(dDataLon.size());
		}
		// Create a new KD Tree containing all nodes
		Announce("Building kdtree");

		kdtree * kdData = kd_create(3);
		if (kdData == NULL) {
			_EXCEPTIONT("kd_create(3) failed");
		}

		for (size_t s = 0; s < dDataLon.size(); s++) {
			double dLatRad = DegToRad(dDataLat[s]);
			double dLonRad = DegToRad(dDataLon[s]);

			double dX = cos(dLonRad) * cos(dLatRad);
			double dY = sin(dLonRad) * cos(dLatRad);
			double dZ = sin(dLatRad);

			kd_insert3(kdData, dX, dY, dZ, (void*)(s));
		}

		// Sample the map
		Announce("Sampling kdtree");

		for (size_t j = 0; j < nImageSizeY; j++) {
			double dLatRad = static_cast<double>(nImageSizeY - 1 - j) / static_cast<double>(nImageSizeY);
			dLatRad = DegToRad(dLatBegin + dLatRad * (dLatEnd - dLatBegin));

			for (size_t i = 0; i < nImageSizeX; i++) {
				double dLonRad = static_cast<double>(i) / static_cast<double>(nImageSizeX);
				dLonRad = DegToRad(dLonBegin + dLonRad * (dLonEnd - dLonBegin));

				double dX = cos(dLonRad) * cos(dLatRad);
				double dY = sin(dLonRad) * cos(dLatRad);
				double dZ = sin(dLatRad);

				kdres * kdresData = kd_nearest3(kdData, dX, dY, dZ);
				_ASSERT(kdresData != NULL);

				vecMap[j*nImageSizeX+i] = (size_t)(kd_res_item_data(kdresData));

				kd_res_free(kdresData);
			}
		}

		kd_free(kdData);

		AnnounceEndBlock("Done");
	}

	// Load data
	size_t sTimeIndex = 0;

	for (size_t f = 0; f < vecInputFiles.size(); f++) {
		AnnounceStartBlock("Processing file %lu/%lu", f, vecInputFiles.size());

		// Open files
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[f]);

		// Get time information
		const NcTimeDimension & vecTimes = vecNcFiles.GetNcTimeDimension(0);
		if (vecTimes.size() == 0) {
			_EXCEPTION1("Input file \"%s\" contains no time information",
				 vecNcFiles.GetFilename(0).c_str());
		}

		if ((strPlotType == "COD") && (vecNcFiles.size() != 1)) {
			_EXCEPTIONT("For --plot \"COD\" one input file must provide COD and SOLIN");
		}
		if ((strPlotType == "ALLSKY") && (vecNcFiles.size() != 1)) {
			_EXCEPTIONT("For --plot \"ALLSKY\" one input file must provide COD, SOLIN and ALLSKY");
		}
		if ((strPlotType == "PRECT") && (vecNcFiles.size() != 3)) {
			_EXCEPTIONT("For --plot \"PRECT\" three input files must provide COD, SOLIN and PRECT");
		}

		for (size_t t = 0; t < vecTimes.size(); t++) {
			AnnounceStartBlock("Processing time %s", vecTimes[t].ToString().c_str());

			// Load COD data
			NcVar * varCOD = vecNcFiles[0]->get_var("cod");
			if (varCOD == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"cod\"",
					vecNcFiles.GetFilename(0).c_str());
			}

			varCOD->set_cur(t,0);
			varCOD->get(&(dataCOD[0]), 1, dataCOD.size());

			// Load SOLIN data
			NcVar * varSOLIN = vecNcFiles[0]->get_var("SOLIN");
			if (varSOLIN == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"SOLIN\"",
					vecNcFiles.GetFilename(0).c_str());
			}

			varSOLIN->set_cur(t,0);
			varSOLIN->get(&(dataSOLIN[0]), 1, dataSOLIN.size());

			// Load AI_ALLSKY data
			if (strPlotType == "ALLSKY") {
				NcVar * varALLSKY = vecNcFiles[0]->get_var("AI_ALLSKY");
				if (varALLSKY == NULL) {
					_EXCEPTION1("File \"%s\" does not contain variable \"AI_ALLSKY\"",
						vecNcFiles.GetFilename(0).c_str());
				}

				varALLSKY->set_cur(t,0);
				varALLSKY->get(&(dataALLSKY[0]), 1, dataALLSKY.size());
			}

			// Load PRECT data
			if (strPlotType == "PRECT") {
				NcVar * varPRECT = vecNcFiles[2]->get_var("PRECT");
				if (varPRECT == NULL) {
					_EXCEPTION1("File \"%s\" does not contain variable \"PRECT\"",
						vecNcFiles.GetFilename(2).c_str());
				}

				varPRECT->set_cur(t,0);
				varPRECT->get(&(dataPRECT[0]), 1, dataPRECT.size());
			}

			// Draw background
			for (int i = 0; i < vecMap.size(); i++) {

				double dIns = 1.0; //sqrt(clamp(dataSOLIN[vecMap[i]] / 1360.0, 0.0, 1.0));

				imgOut[4*i+0] = dIns * imgBlueMarbleSub[4*i+0] + (1.0 - dIns) * imgBlackMarbleSub[4*i+0];
				imgOut[4*i+1] = dIns * imgBlueMarbleSub[4*i+1] + (1.0 - dIns) * imgBlackMarbleSub[4*i+1];
				imgOut[4*i+2] = dIns * imgBlueMarbleSub[4*i+2] + (1.0 - dIns) * imgBlackMarbleSub[4*i+2];
				imgOut[4*i+3] = 255;

				// Draw COD
				double dCOD = clamp(dataCOD[vecMap[i]], 0.0f, 40.0f);

				double dColor;
				double dAlpha;
				if (dCOD > 40.0) {
					dAlpha = 1.0;
					dColor = 1.0;
				} else {
					dAlpha = dCOD / 40.0;
					dColor = 0.7 + 0.3 * dCOD / 40.0;
				}

				dColor -= 0.5 * (1.0 - dIns);

				char cColor = static_cast<char>(255.0 * dColor);

				imgOut[4*i+0] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+0];
				imgOut[4*i+1] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+1];
				imgOut[4*i+2] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+2];
				imgOut[4*i+3] = 255;

			}

			// ALLSKY plot
			if (strPlotType == "ALLSKY") {
				for (int i = 0; i < vecMap.size(); i++) {

					double dALLSKY = clamp(dataALLSKY[vecMap[i]] / 0.8, 0.0, 1.0);

					unsigned char cR = 218;
					unsigned char cG = 96;
					unsigned char cB = 12;

					double dAlpha = dALLSKY;

					imgOut[4*i+0] = dAlpha * cR + (1.0 - dAlpha) * imgOut[4*i+0];
					imgOut[4*i+1] = dAlpha * cG + (1.0 - dAlpha) * imgOut[4*i+1];
					imgOut[4*i+2] = dAlpha * cB + (1.0 - dAlpha) * imgOut[4*i+2];
				}
			}

			// PRECT plot
			if (strPlotType == "PRECT") {
				for (int i = 1; i < vecMap.size()-1; i++) {

					imgOut[4*i+0] *= 0.4;
					imgOut[4*i+1] *= 0.4;
					imgOut[4*i+2] *= 0.4;

					//double dPRECTmmhr = dataPRECT[vecMap[i]] / 2.778e-7;
					double dPRECTmmhr = (0.2 * dataPRECT[vecMap[i-1]] + 0.6 * dataPRECT[vecMap[i]] + 0.2 * dataPRECT[vecMap[i]]) / 2.778e-7;

					if (dPRECTmmhr > 1.0/24.0) {

						unsigned char cR;
						unsigned char cG;
						unsigned char cB;

						double dA;

						if (dPRECTmmhr < 1.0) {
							cR = 0; cG = 64; cB = 64; dA = 0.3;
						} else if (dPRECTmmhr < 2.0) {
							cR = 0; cG = 64; cB = 0; dA = 0.3;
						} else if (dPRECTmmhr < 3.0) {
							cR = 0; cG = 128; cB = 0; dA = 0.5;
						} else if (dPRECTmmhr < 4.0) {
							cR = 0; cG = 192; cB = 0; dA = 0.7;
						} else if (dPRECTmmhr < 5.0) {
							cR = 128; cG = 192; cB = 0; dA = 0.9;
						} else if (dPRECTmmhr < 10.0) {
							cR = 192; cG = 228; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 15.0) {
							cR = 228; cG = 192; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 20.0) {
							cR = 255; cG = 128; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 35.0) {
							cR = 192; cG = 64; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 50.0) {
							cR = 128; cG = 0; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 75.0) {
							cR = 192; cG = 0; cB = 0; dA = 1.0;
						} else if (dPRECTmmhr < 100.0) {
							cR = 255; cG = 0; cB = 255; dA = 1.0;
						} else if (dPRECTmmhr < 150.0) {
							cR = 128; cG = 0; cB = 192; dA = 1.0;
						} else {
							cR = 240; cG = 240; cB = 240; dA = 1.0;
						}

						imgOut[4*i+0] = dA * cR + (1.0 - dA) * imgOut[4*i+0];
						imgOut[4*i+1] = dA * cG + (1.0 - dA) * imgOut[4*i+1];
						imgOut[4*i+2] = dA * cB + (1.0 - dA) * imgOut[4*i+2];
					}
				}
			}

			// Add text
			Time timeGMT = vecTimes[t];

			sctext.DrawString(
				timeGMT.ToString() + std::string(" (GMT)"),
				10,
				nImageSizeY-12,
				SchriftText::TextAlignment_Left,
				nImageSizeX,
				nImageSizeY,
				&(imgOut[0]),
				RGBA());
/*
			// Add Guam
			int iGuamX = nImageSizeX * (144.794 - LonDegToStandardRange(dLonBegin)) / (LonDegToStandardRange(dLonEnd) - LonDegToStandardRange(dLonBegin));
			int iGuamY = nImageSizeY - 1 - nImageSizeY * (13.444 - dLatBegin) / (dLatEnd - dLatBegin);
			int iGuamIx = nImageSizeX * iGuamY + iGuamX;

			sctext.DrawCharacter(
				'+',
				iGuamX-5,
				iGuamY+6,
				nImageSizeX,
				nImageSizeY,
				&(imgOut[0]),
				SchriftText::RGBA(192,0,0));

			imgOut[4*iGuamIx+0] = 255;
			imgOut[4*iGuamIx+1] = 0;
			imgOut[4*iGuamIx+2] = 0;
*/
			// Output file
			char szBuffer[32];
			std::string strOutputFilename = strOutputRoot;
			snprintf(szBuffer, 32, "%05lu.png", sTimeIndex);
			strOutputFilename += szBuffer;

			// Write file
			imgOut.write(strOutputFilename);

			AnnounceEndBlock("Done");

			sTimeIndex++;
		}
		AnnounceEndBlock("Done");
	}

	// Write image
	//

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

