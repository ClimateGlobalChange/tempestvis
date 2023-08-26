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
#include <cfloat>

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

void DrawMesh(
	Mesh & meshShp,
	double dLatBegin,
	double dLatEnd,
	double dLonBegin,
	double dLonEnd,
	PNGImage & imgOut,
	std::vector<RGBA> * pvecOutlGRBA = NULL,
	std::vector<RGBA> * pvecFillRGBA = NULL
) {
	int nImageSizeX = imgOut.width();
	int nImageSizeY = imgOut.height();

	for (size_t f = 0; f < meshShp.faces.size(); f++) {
		const Face & face = meshShp.faces[f];
		if (face.edges.size() == 0) {
			continue;
		}

		// Convert edges to image coords
		std::vector<int> vecX;
		std::vector<int> vecY;
		vecX.reserve(face.edges.size());
		vecY.reserve(face.edges.size());

		for (int i = 0; i < face.edges.size(); i++) {
			double dLatDeg;
			double dLonDeg;

			XYZtoRLL_Deg(
				meshShp.nodes[face[i]].x,
				meshShp.nodes[face[i]].y,
				meshShp.nodes[face[i]].z,
				dLonDeg,
				dLatDeg);

			int iX, iY;

			RealCoordToImageCoord(
				dLatDeg, dLonDeg,
				dLatBegin, dLonBegin,
				dLatEnd, dLonEnd,
				nImageSizeX, nImageSizeY,
				iX, iY);

			if (i != 0) {
				size_t sLast = vecX.size()-1;
				if ((vecX[sLast] != iX) || (vecY[sLast] != iY)) {
					vecX.push_back(iX);
					vecY.push_back(iY);
				}
			} else {
				vecX.push_back(iX);
				vecY.push_back(iY);
			}
		}

		// Draw fill
		if (pvecFillRGBA != NULL) {
			unsigned char cRf = 255;
			unsigned char cGf = 255;
			unsigned char cBf = 255;
			unsigned char cAf = 0;

			if (pvecFillRGBA->size() == meshShp.faces.size()) {
				cRf = (*pvecFillRGBA)[f].r();
				cGf = (*pvecFillRGBA)[f].g();
				cBf = (*pvecFillRGBA)[f].b();
				cAf = (*pvecFillRGBA)[f].a();

			} else if (pvecFillRGBA->size() == 1) {
				cRf = (*pvecFillRGBA)[0].r();
				cGf = (*pvecFillRGBA)[0].g();
				cBf = (*pvecFillRGBA)[0].b();
				cAf = (*pvecFillRGBA)[0].a();

			} else if (pvecFillRGBA->size() == 0) {
				cRf = 255;
				cGf = 255;
				cBf = 255;
				cAf = 255;

			} else {
				_EXCEPTION1("Anomalous pvecFillRGBA array size (%lu)", pvecFillRGBA->size());
			}

			if (cAf != 0) {
				double dAf = static_cast<double>(cAf) / 255.0;

				// Get the polygon bounds
				int iMinX = vecX[0];
				int iMaxX = vecX[0];

				int iMinY = vecY[0];
				int iMaxY = vecY[0];

				for (int i = 0; i < vecX.size(); i++) {
					if (vecX[i] < iMinX) {
						iMinX = vecX[i];
					}
					if (vecX[i] > iMaxX) {
						iMaxX = vecX[i];
					}
					if (vecY[i] < iMinY) {
						iMinY = vecY[i];
					}
					if (vecY[i] > iMaxY) {
						iMaxY = vecY[i];
					}
				}

				// Draw pixels
				if ((iMinX >= 0) && (iMaxX < nImageSizeX) && (iMinY >= 0) && (iMaxY < nImageSizeY)) {
/*
					if (f == 807) {
						for (int v = 0; v < vecX.size(); v++) {
							printf("%i %i\n", vecX[v], vecY[v]);
						}
					}
*/
					//printf("%i %i %i %i\n", iMinX, iMaxX, iMinY, iMaxY);
					for (int j = iMinY; j <= iMaxY; j++) {
/*
						if (f != 807) {
							continue;
						}
*/
						std::vector<int> nPosParity(iMaxX - iMinX + 1, 0);
						std::vector<int> nNegParity(iMaxX - iMinX + 1, 0);
						for (size_t v1 = 0; v1 < vecX.size(); v1++) {
							size_t v2 = (v1 + 1) % vecX.size();
							if (vecY[v1] == vecY[v2]) {
								continue;
							}

							if (vecY[v2] > vecY[v1]) {
								if ((j > vecY[v2]) || (j < vecY[v1])) {
									continue;
								}
								if (j == vecY[v1]) {
									nPosParity[vecX[v1] - iMinX]++;
									continue;
								}
								if (j == vecY[v2]) {
									nPosParity[vecX[v2] - iMinX]++;
									continue;
								}

								int ix = ((vecX[v2] - vecX[v1]) * (j - vecY[v1])) / (vecY[v2] - vecY[v1]) + vecX[v1] - iMinX;
								_ASSERT((ix >= 0) && (ix < nPosParity.size()));
								nPosParity[ix] += 2;

							} else {
								if ((j > vecY[v1]) || (j < vecY[v2])) {
									continue;
								}
								if (j == vecY[v1]) {
									nNegParity[vecX[v1] - iMinX]++;
									continue;
								}
								if (j == vecY[v2]) {
									nNegParity[vecX[v2] - iMinX]++;
									continue;
								}

								int ix = ((vecX[v2] - vecX[v1]) * (j - vecY[v1])) / (vecY[v2] - vecY[v1]) + vecX[v1] - iMinX;
								_ASSERT((ix >= 0) && (ix < nPosParity.size()));
								nNegParity[ix] += 2;
							}
						}
/*
						printf("%i: ", j);
						for (int v = 0; v < nPosParity.size(); v++) {
							printf("(%i,%i) ", nPosParity[v], nNegParity[v]);
						}
						printf("\n");
*/
						int iParity = 0;
						for (int i = iMinX; i <= iMaxX; i++) {
							iParity += nPosParity[i-iMinX];

							if (iParity > 0) {
								int kx = j * nImageSizeX + i;
								imgOut[4*kx+0] = dAf * cRf + (1.0 - dAf) * imgOut[4*kx+0];
								imgOut[4*kx+1] = dAf * cGf + (1.0 - dAf) * imgOut[4*kx+1];
								imgOut[4*kx+2] = dAf * cBf + (1.0 - dAf) * imgOut[4*kx+2];
							}

							iParity -= nNegParity[i-iMinX];
/*
							if (iParity < 0) {
								_EXCEPTION1("Error in face %i", f);
							}
*/
							//_ASSERT(iParity >= 0);
						}
					}
/*
					if (f == 807) {
						_EXCEPTION();
					}
*/
				}
			}
		}

		// Outline color
		unsigned char cRo = 255;
		unsigned char cGo = 255;
		unsigned char cBo = 255;
		unsigned char cAo = 255;

		if (pvecOutlGRBA != NULL) {
			if (pvecOutlGRBA->size() == meshShp.faces.size()) {
				cRo = (*pvecOutlGRBA)[f].r();
				cGo = (*pvecOutlGRBA)[f].g();
				cBo = (*pvecOutlGRBA)[f].b();
				cAo = (*pvecOutlGRBA)[f].a();

			} else if (pvecOutlGRBA->size() == 1) {
				cRo = (*pvecOutlGRBA)[0].r();
				cGo = (*pvecOutlGRBA)[0].g();
				cBo = (*pvecOutlGRBA)[0].b();
				cAo = (*pvecOutlGRBA)[0].a();

			} else if (pvecOutlGRBA->size() == 0) {
				cRo = 255;
				cGo = 255;
				cBo = 255;
				cAo = 255;

			} else {
				_EXCEPTION1("Anomalous pvecOutlGRBA array size (%lu)", pvecOutlGRBA->size());
			}
		}

		// Draw outline
		if (cAo != 0) {
			double dAo = static_cast<double>(cAo) / 255.0;

			int iXprev = 0;
			int iYprev = 0;

			int iXnext = vecX[0];
			int iYnext = vecY[0];

			for (size_t v = 1; v < vecX.size(); v++) {
				iXprev = iXnext;
				iYprev = iYnext;

				iXnext = vecX[v];
				iYnext = vecY[v];

				int nDistX = abs(iXnext - iXprev);
				int nDistY = abs(iYnext - iYprev);

				int nDistMax = nDistX;
				if (nDistY > nDistX) {
					nDistMax = nDistY;
				}
				if (nDistMax > 0.8 * nImageSizeX) {
					continue;
				}

				double dXstep = static_cast<double>(iXnext - iXprev) / static_cast<double>(nDistMax);
				double dYstep = static_cast<double>(iYnext - iYprev) / static_cast<double>(nDistMax);

				for (int i = 0; i < nDistMax; i+=1) {
					int iXcoord = iXprev + static_cast<int>(dXstep * i);
					int iYcoord = iYprev + static_cast<int>(dYstep * i);

					if ((iXcoord >= 0) && (iXcoord < nImageSizeX) && (iYcoord >= 0) && (iYcoord < nImageSizeY)) {
						size_t ix = iXcoord;
						size_t jx = iYcoord;

						int kx = nImageSizeX * jx + ix;
						imgOut[4*kx+0] = dAo * cRo + (1.0 - dAo) * imgOut[4*kx+0];
						imgOut[4*kx+1] = dAo * cGo + (1.0 - dAo) * imgOut[4*kx+1];
						imgOut[4*kx+2] = dAo * cBo + (1.0 - dAo) * imgOut[4*kx+2];

						if (ix != 0) {
							int kx = nImageSizeX * jx + ix-1;
							imgOut[4*kx+0] = dAo * (0.5 * cRo + 0.5 * imgOut[4*kx+0]) + (1.0 - dAo) * imgOut[4*kx+0];
							imgOut[4*kx+1] = dAo * (0.5 * cGo + 0.5 * imgOut[4*kx+1]) + (1.0 - dAo) * imgOut[4*kx+1];
							imgOut[4*kx+2] = dAo * (0.5 * cBo + 0.5 * imgOut[4*kx+2]) + (1.0 - dAo) * imgOut[4*kx+2];
						}

						if (jx != 0) {
							int kx = nImageSizeX * (jx-1) + ix;
							imgOut[4*kx+0] = dAo * (0.5 * cRo + 0.5 * imgOut[4*kx+0]) + (1.0 - dAo) * imgOut[4*kx+0];
							imgOut[4*kx+1] = dAo * (0.5 * cGo + 0.5 * imgOut[4*kx+1]) + (1.0 - dAo) * imgOut[4*kx+1];
							imgOut[4*kx+2] = dAo * (0.5 * cBo + 0.5 * imgOut[4*kx+2]) + (1.0 - dAo) * imgOut[4*kx+2];
						}
					}
				}
			}
		}
	}
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

	// Light angle
	double dAngleDeg;

	// Stencil width
	int nStencilWidth;

	// Maximum scale
	double dMaxScale;

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
		CommandLineStringD(strPlotType, "plot", "", "(|DATA)");
		CommandLineString(strShapefile, "shp", "");
		CommandLineString(strShapefile2, "shp2", "");
		CommandLineDoubleD(dAngleDeg, "lightangle", 135.0, "(0 = east)");
		CommandLineInt(nStencilWidth, "width", 10);
		CommandLineDouble(dMaxScale, "maxscale", 5000.0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list) may be specified");
	}
	if (strOutputRoot.length() == 0) {
		_EXCEPTIONT("No output file (--out_root) specified");
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
	if ((strPlotType != "") && (strPlotType != "DATA")) {
		_EXCEPTIONT("--plot must be one of \"\" or \"DATA\"");
	}

	// Convert light angle to radians
	double dAngleRad = DegToRad(dAngleDeg);

	// Load input file list
	FilenameList vecInputFiles;
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		if (strPlotType == "DATA") {
			_EXCEPTIONT("For --plot \"DATA\" input data must also be specified");
		}
		vecInputFiles.push_back("!");

	} else if (strInputData.length() != 0) {
		vecInputFiles.push_back(strInputData);
	} else {
		vecInputFiles.FromFile(strInputDataList, true);
	}

	// Create output image
	PNGImage imgOut(nImageSizeX, nImageSizeY);

	// Get the topography data
	AnnounceStartBlock("Loading topography data");

	std::string strTopoFile = strDataPath + std::string("/elevation_1km_GMTEDmd.nc");
	NcFile ncfileTopo(strTopoFile.c_str());
	if (!ncfileTopo.is_valid()) {
		_EXCEPTION1("Unable to open data file \"%s\"", strTopoFile.c_str());
	}

	std::vector<float> dZ;
	std::vector<float> dDZ;
	std::vector<double> dZLatDeg;
	std::vector<double> dZLonDeg;

	{
		NcVar * varLon = ncfileTopo.get_var("lon");
		if (varLon == NULL) {
			_EXCEPTION1("Data file \"%s\" does not not contain variable \"lon\"", strTopoFile.c_str());
		}
		if (varLon->num_dims() != 1) {
			_EXCEPTION1("Data file \"%s\" variable \"lon\" must have one dimension", strTopoFile.c_str());
		}
		dZLonDeg.resize(varLon->get_dim(0)->size());
		varLon->get(&(dZLonDeg[0]), dZLonDeg.size());

		NcVar * varLat = ncfileTopo.get_var("lat");
		if (varLat == NULL) {
			_EXCEPTION1("Data file \"%s\" does not not contain variable \"lat\"", strTopoFile.c_str());
		}
		if (varLat->num_dims() != 1) {
			_EXCEPTION1("Data file \"%s\" variable \"lat\" must have one dimension", strTopoFile.c_str());
		}
		dZLatDeg.resize(varLat->get_dim(0)->size());
		varLat->get(&(dZLatDeg[0]), dZLatDeg.size());

		NcVar * varZ = ncfileTopo.get_var("Band1");
		if (varZ == NULL) {
			_EXCEPTION1("Data file \"%s\" does not not contain variable \"Band1\"", strTopoFile.c_str());
		}
		if (varZ->num_dims() != 2) {
			_EXCEPTION1("Data file \"%s\" variable \"Band1\" must have two dimensions", strTopoFile.c_str());
		}
		dZ.resize(dZLatDeg.size() * dZLonDeg.size());
		varZ->get(&(dZ[0]), dZLatDeg.size(), dZLonDeg.size());

		dDZ.resize(dZLatDeg.size() * dZLonDeg.size(), FLT_MAX);
	}

	size_t sZLat = dZLatDeg.size();
	size_t sZLon = dZLonDeg.size();

	double dZLatDegMin = dZLatDeg[0];
	double dZLatDegMax = dZLatDeg[sZLat-1];
	double dZLonDegMin = dZLonDeg[0];
	double dZLonDegMax = dZLonDeg[sZLon-1];

	AnnounceEndBlock("Done");

	// Initialize text renderer
	SchriftText sctext(strDataPath + std::string("/Ubuntu-Regular.ttf"), 18.0);

	// Data buffer
	std::vector<float> dataFLUT;
	std::vector<float> dataSOLIN;

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

	// Load data longitude and latitude information
	if (strPlotType == "DATA") {
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
		dataFLUT.resize(dDataLon.size());

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

	// Loop through all input files
	size_t sTimeIndex = 0;

	for (size_t f = 0; f < vecInputFiles.size(); f++) {
		AnnounceStartBlock("Processing file %lu/%lu", f, vecInputFiles.size());

		// Open files
		NcTimeDimension vecTimes;
		if (vecInputFiles[f] != "!") {
			NcFileVector vecNcFiles;
			vecNcFiles.ParseFromString(vecInputFiles[f]);

			// Get time information
			vecTimes = vecNcFiles.GetNcTimeDimension(0);
			if (vecTimes.size() == 0) {
				_EXCEPTION1("Input file \"%s\" contains no time information",
					 vecNcFiles.GetFilename(0).c_str());
			}

		} else {
			vecTimes.push_back(Time());
		}

		// Loop through all times
		for (size_t t = 0; t < vecTimes.size(); t++) {
			AnnounceStartBlock("Processing time %s", vecTimes[t].ToString().c_str());

			// Calculate gradients and draw background
			for (size_t j = 0; j < nImageSizeY; j++) {
				double dLatDeg = static_cast<double>(nImageSizeY - 1 - j) / static_cast<double>(nImageSizeY);
				dLatDeg = dLatBegin + dLatDeg * (dLatEnd - dLatBegin);

				int jZ = (dLatDeg - dZLatDegMin) / (dZLatDegMax - dZLatDegMin) * static_cast<double>(sZLat);

				if (jZ < 0) {
					jZ = 0;
				}
				if (jZ > sZLat-1) {
					jZ = sZLat-1;
				}

				for (size_t i = 0; i < nImageSizeX; i++) {
					size_t k = j * nImageSizeX + i;

					double dLonDeg = static_cast<double>(i) / static_cast<double>(nImageSizeX);
					dLonDeg = dLonBegin + dLonDeg * (dLonEnd - dLonBegin);

					int iZ = (dLonDeg - dZLonDegMin) / (dZLonDegMax - dZLonDegMin) * static_cast<double>(sZLon);

					if (iZ < 0) {
						iZ += (1 - (iZ + 1) / sZLon) * sZLon;
					}
					if (iZ >= sZLon) {
						iZ -= (iZ / sZLon) * sZLon;
					}

					int kZ = jZ * sZLon + iZ;

					// Calculate derivative along light angle
					if (dDZ[kZ] == FLT_MAX) {
						int jZnext = jZ + nStencilWidth;
						if (jZnext > sZLat-1) {
							jZnext = sZLat-1;
						}

						int jZprev = jZ - nStencilWidth;
						if (jZprev < 0) {
							jZprev = 0;
						}

						int kZnext = jZnext * sZLon + iZ;
						int kZprev = jZprev * sZLon + iZ;

						double dDlatZ = (dZ[kZnext] - dZ[kZprev]) / (dZLatDeg[jZnext] - dZLatDeg[jZprev]);

						int iZnext = (iZ + 1) % sZLon;
						int iZprev = (iZ + sZLon - 1) % sZLon;

						kZnext = jZ * sZLon + iZnext;
						kZprev = jZ * sZLon + iZprev;

						double dDlonZ = (dZ[kZnext] - dZ[kZprev]) / (dZLonDeg[iZnext] - dZLonDeg[iZprev]) / cos(DegToRad(dZLatDeg[jZ]));

						dDZ[kZ] = - cos(dAngleRad) * dDlonZ - sin(dAngleRad) * dDlatZ;
					}

					double dGradZ = clamp(static_cast<double>(dDZ[kZ]), -dMaxScale, dMaxScale);
					dGradZ = (dGradZ + dMaxScale) / (2.0 * dMaxScale);

					//double dGradZ = clamp(dZ[kZ] / 4000.0f, 0.0f, 1.0f);

					imgOut[4*k+0] = dGradZ * 255;
					imgOut[4*k+1] = dGradZ * 255;
					imgOut[4*k+2] = dGradZ * 255;
					imgOut[4*k+3] = 255;
				}
			}
/*
			// Load data
			NcVar * varFLUT = vecNcFiles[0]->get_var("FLUT");
			if (varFLUT == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"FLUT\"",
					vecNcFiles.GetFilename(0).c_str());
			}

			varFLUT->set_cur(t,0);
			varFLUT->get(&(dataFLUT[0]), 1, dataFLUT.size());
*/
/*
			// Draw background
			for (int i = 0; i < vecMap.size(); i++) {

				double dIns = 0.7; //sqrt(clamp(dataSOLIN[vecMap[i]] / 1360.0, 0.0, 1.0));

				imgOut[4*i+0] = dIns * imgBlueMarbleSub[4*i+0] + (1.0 - dIns) * imgBlackMarbleSub[4*i+0];
				imgOut[4*i+1] = dIns * imgBlueMarbleSub[4*i+1] + (1.0 - dIns) * imgBlackMarbleSub[4*i+1];
				imgOut[4*i+2] = dIns * imgBlueMarbleSub[4*i+2] + (1.0 - dIns) * imgBlackMarbleSub[4*i+2];
				imgOut[4*i+3] = 255;

				double dFLUT = clamp(dataFLUT[vecMap[i]], 100.0f, 300.0f);

				double dColor;
				double dAlpha;
				if (dFLUT < 160.0) {
					dAlpha = 1.0;
					dColor = 0.7 + 0.2 * (160.0 - dFLUT) / (160.0 - 100.0);
				} else {
					dAlpha = 1.0 - (dFLUT - 160.0) / (300.0 - 160.0);
					dColor = 0.5 + 0.2 * (300.0 - dFLUT) / (300.0 - 160.0);
				}

				char cColor = static_cast<char>(255.0 * dColor);

				imgOut[4*i+0] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+0];
				imgOut[4*i+1] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+1];
				imgOut[4*i+2] = dAlpha * cColor + (1.0 - dAlpha) * imgOut[4*i+2];
				imgOut[4*i+3] = 255;

			}
*/
			// Draw shapefiles
			if (strShapefile != "") {
				DrawMesh(
					meshShp,
					dLatBegin,
					dLatEnd,
					dLonBegin,
					dLonEnd,
					imgOut);
			}
			if (strShapefile2 != "") {
				std::vector<RGBA> vecOutlGRBA(meshShp2.faces.size());
				std::vector<RGBA> vecFillRGBA(meshShp2.faces.size());
				double dPx, dPy, dPz;
				RLLtoXYZ_Deg(-75.16, 39.95, dPx, dPy, dPz);

				for (int f = 0; f < meshShp2.faces.size(); f++) {
					double dGCD = GreatCircleDistanceXYZ_Deg(
						meshShp2.nodes[meshShp2.faces[f][0]].x,
						meshShp2.nodes[meshShp2.faces[f][0]].y,
						meshShp2.nodes[meshShp2.faces[f][0]].z,
						dPx,
						dPy,
						dPz
					);

					if (dGCD > 3.0) {
						vecOutlGRBA[f] = RGBA(127, 127, 127, 64);
						vecFillRGBA[f] = RGBA(0, 0, 0, 0);

					} else if (dGCD > 2.0) {
						vecOutlGRBA[f] = RGBA(216, 216, 0, 128);
						vecFillRGBA[f] = RGBA(108, 108, 0, 64);

					} else if (dGCD > 1.0) {
						vecOutlGRBA[f] = RGBA(216, 128, 0, 192);
						vecFillRGBA[f] = RGBA(108, 64, 0, 96);

					} else {
						vecOutlGRBA[f] = RGBA(216, 0, 0, 255);
						vecFillRGBA[f] = RGBA(108, 0, 0, 128);
					}
				}

				DrawMesh(
					meshShp2,
					dLatBegin,
					dLatEnd,
					dLonBegin,
					dLonEnd,
					imgOut,
					&vecOutlGRBA,
					&vecFillRGBA);
			}

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

