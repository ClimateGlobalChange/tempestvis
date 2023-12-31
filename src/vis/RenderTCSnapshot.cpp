///////////////////////////////////////////////////////////////////////////////
///
///	\file    RenderTCSnapshot.cpp
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
#include "PNGImage.h"

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

	// Difference plot
	bool fDelta;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strOutputRoot, "out_root", "img");
		CommandLineInt(nImageSizeX, "imgx", 360);
		CommandLineInt(nImageSizeY, "imgy", 360);
		CommandLineString(strDataPath, "datapath", strDataPath.c_str());
		CommandLineStringD(strPlotType, "plot", "FLUT", "(FLUT|U10|PRECT)");
		CommandLineBool(fDelta, "delta");

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
	if ((strPlotType != "FLUT") && (strPlotType != "U10") && (strPlotType != "PRECT")) {
		_EXCEPTIONT("--plot must be one of \"FLUT\", \"U10\" or \"PRECT\"");
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

	// Initialize text renderer
	SchriftText sctext(strDataPath + std::string("/Ubuntu-Regular.ttf"), 18.0);

	// Data buffer
	std::vector<float> dataFLUT;
	std::vector<float> dataU10;
	std::vector<float> dataPRECT;

	// ColorMap
	ColorMapLibrary cmaplib(strDataPath);
	ColorMap cmapJet;
	cmaplib.GenerateColorMap("jet", cmapJet);
	ColorMap cmapBalance;
	cmaplib.GenerateColorMap("bluered", cmapBalance);

	// Load data
	size_t sTimeIndex = 0;

	for (size_t f = 0; f < vecInputFiles.size(); f++) {
		AnnounceStartBlock("Processing file %lu/%lu", f, vecInputFiles.size());

		// Open files
		NcFileVector vecNcFiles;
		vecNcFiles.ParseFromString(vecInputFiles[f]);

		// Get time information
		if ((strPlotType == "FLUT") && (vecNcFiles.size() != 1)) {
			_EXCEPTIONT("For --plot \"FLUT\" one input file must provide FLUT");
		}
		if ((strPlotType == "U10") && (vecNcFiles.size() != 2)) {
			_EXCEPTIONT("For --plot \"U10\" two input files must provide FLUT and U10");
		}
		if ((strPlotType == "PRECT") && (vecNcFiles.size() != 2)) {
			_EXCEPTIONT("For --plot \"PRECT\" two input files must provide FLUT and PRECT");
		}

		NcDim * dimSnapshot = vecNcFiles[0]->get_dim("snapshot");
		if (dimSnapshot == NULL) {
			_EXCEPTION1("File \"%s\" does not contain dimension \"snapshot\"", vecNcFiles.GetFilename(0).c_str());
		}

		NcDim * dimY = vecNcFiles[0]->get_dim("y");
		if (dimY == NULL) {
			_EXCEPTION1("File \"%s\" does not contain dimension \"y\"", vecNcFiles.GetFilename(0).c_str());
		}

		NcDim * dimX = vecNcFiles[0]->get_dim("x");
		if (dimX == NULL) {
			_EXCEPTION1("File \"%s\" does not contain dimension \"x\"", vecNcFiles.GetFilename(0).c_str());
		}

		size_t nSnapshots = dimSnapshot->size();
		size_t nSnapshotY = dimY->size();
		size_t nSnapshotX = dimX->size();

		dataFLUT.resize(nSnapshotY * nSnapshotX);
		if (strPlotType == "U10") {
			dataU10.resize(nSnapshotY * nSnapshotX);
		}
		if (strPlotType == "PRECT") {
			dataPRECT.resize(nSnapshotY * nSnapshotX);
		}

		// Generate map
		std::vector<size_t> vecMap(nImageSizeX * nImageSizeY);
		for (int j = 0; j < nImageSizeY; j++) {
			int jx = static_cast<int>(static_cast<double>(j)
					/ static_cast<double>(nImageSizeY)
					* static_cast<double>(nSnapshotY));

			jx = nSnapshotY - 1 - jx;

			for (int i = 0; i < nImageSizeX; i++) {
				int ix = static_cast<int>(static_cast<double>(i)
						/ static_cast<double>(nImageSizeX)
						* static_cast<double>(nSnapshotX));
				
				vecMap[j * nImageSizeX + i] = jx * nSnapshotX + ix;
			}
		}

		// Draw each snapshot
		for (size_t t = 0; t < nSnapshots; t++) {
			Time timeGuam(2023, 4, 22, 10 * 3600, 0);
			if (t < 24) {
				continue;
			}
			timeGuam.AddSeconds(t * 15 * 60);

			AnnounceStartBlock("Processing snapshot %i", t);

			// Load FLUT data
			NcVar * varFLUT = vecNcFiles[0]->get_var("snap_FLUT");
			if (varFLUT == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"FLUT\"",
					vecNcFiles.GetFilename(0).c_str());
			}

			varFLUT->set_cur(t,0,0);
			varFLUT->get(&(dataFLUT[0]), 1, nSnapshotY, nSnapshotX);

			// Load U10 data
			if (strPlotType == "U10") {
				NcVar * varU10 = vecNcFiles[1]->get_var("snap_U10");
				if (varU10 == NULL) {
					_EXCEPTION1("File \"%s\" does not contain variable \"U10\"",
						vecNcFiles.GetFilename(1).c_str());
				}

				varU10->set_cur(t,0,0);
				varU10->get(&(dataU10[0]), 1, nSnapshotY, nSnapshotX);
			}

			// Load PRECT data
			if (strPlotType == "PRECT") {
				NcVar * varPRECT = vecNcFiles[1]->get_var("snap_PRECT");
				if (varPRECT == NULL) {
					_EXCEPTION1("File \"%s\" does not contain variable \"PRECT\"",
						vecNcFiles.GetFilename(1).c_str());
				}

				varPRECT->set_cur(t,0,0);
				varPRECT->get(&(dataPRECT[0]), 1, nSnapshotY, nSnapshotX);
			}

			// Draw FLUT
			for (int i = 0; i < vecMap.size(); i++) {
				double dFLUT = clamp(dataFLUT[vecMap[i]], 80.0f, 300.0f);

				double dColor;
				if (dFLUT < 120.0) {
					dColor = 0.7 + 0.2 * (120.0 - dFLUT) / (120.0 - 80.0);
				} else {
					dColor = 0.5 + 0.2 * (300.0 - dFLUT) / (300.0 - 120.0);
				}

				char cColor = static_cast<char>(255.0 * dColor);

				imgOut[4*i+0] = cColor;
				imgOut[4*i+1] = cColor;
				imgOut[4*i+2] = cColor;
				imgOut[4*i+3] = 255;
			}

			// U10 plot
			if ((strPlotType == "U10") && !fDelta) {
				for (int i = 0; i < vecMap.size(); i++) {

					imgOut[4*i+0] *= 0.4;
					imgOut[4*i+1] *= 0.4;
					imgOut[4*i+2] *= 0.4;

					if (dataU10[vecMap[i]] > 10.0) {
						double dU10 = clamp(dataU10[vecMap[i]] / 60.0, 0.0, 1.0);

						unsigned char cR;
						unsigned char cG;
						unsigned char cB;

						cmapJet.Sample(dU10, 0.0, 1.0, cR, cG, cB);

						if (dataU10[vecMap[i]] < 20.0) {
							double dA = (dataU10[vecMap[i]] - 10.0) / (20.0 - 10.0);
							imgOut[4*i+0] = dA * cR + (1.0 - dA) * imgOut[4*i+0];
							imgOut[4*i+1] = dA * cG + (1.0 - dA) * imgOut[4*i+1];
							imgOut[4*i+2] = dA * cB + (1.0 - dA) * imgOut[4*i+2];
						} else {
							imgOut[4*i+0] = cR;
							imgOut[4*i+1] = cG;
							imgOut[4*i+2] = cB;
						}
					}
				}
			}

			// PRECT plot
			if ((strPlotType == "U10") && fDelta) {
				for (int i = 1; i < vecMap.size()-1; i++) {

					imgOut[4*i+0] *= 0.2;
					imgOut[4*i+1] *= 0.2;
					imgOut[4*i+2] *= 0.2;

					double dA = 1.0;
					double dScale = clamp(dataU10[vecMap[i]] / 10.0, -1.0, 1.0);

					unsigned char cR;
					unsigned char cG;
					unsigned char cB;

					if (dScale < 0.0) {
						cR = 64 * sqrt(-dScale);
						cG = 128 * sqrt(-dScale);
						cB = 255 * sqrt(-dScale);
						dA = sqrt(-dScale);
					} else {
						cR = 255 * sqrt(dScale);
						cG = 64 * sqrt(dScale);
						cB = 64 * sqrt(dScale);
						dA = sqrt(dScale);

					}

					imgOut[4*i+0] = dA * cR + (1.0 - dA) * imgOut[4*i+0];
					imgOut[4*i+1] = dA * cG + (1.0 - dA) * imgOut[4*i+1];
					imgOut[4*i+2] = dA * cB + (1.0 - dA) * imgOut[4*i+2];
				}
			}

			// PRECT plot
			if ((strPlotType == "PRECT") && !fDelta) {
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

			// PRECT plot
			if ((strPlotType == "PRECT") && fDelta) {
				for (int i = 1; i < vecMap.size()-1; i++) {

					imgOut[4*i+0] *= 0.2;
					imgOut[4*i+1] *= 0.2;
					imgOut[4*i+2] *= 0.2;

					double dPRECTmmhr = (0.2 * dataPRECT[vecMap[i-1]] + 0.6 * dataPRECT[vecMap[i]] + 0.2 * dataPRECT[vecMap[i]]) / 2.778e-7;

					double dA = 1.0;
					double dScale = clamp(dPRECTmmhr / 180.0, -1.0, 1.0);

					unsigned char cR;
					unsigned char cG;
					unsigned char cB;

					if (dScale < 0.0) {
						cR = 255 * sqrt(sqrt(-dScale));
						cG = 64 * sqrt(sqrt(-dScale));
						cB = 0;
						dA = sqrt(sqrt(-dScale));
					} else {
						cR = 0;
						cG = 255 * sqrt(sqrt(dScale));
						cB = 255 * sqrt(sqrt(dScale));
						dA = sqrt(sqrt(dScale));
					}
					//cmapBalance.Sample(dScale, -1.0, 1.0, cR, cG, cB);

					//printf("%1.5f %i %i %i\n", dPRECTmmhr, cR, cG, cB);

					imgOut[4*i+0] = dA * cR + (1.0 - dA) * imgOut[4*i+0];
					imgOut[4*i+1] = dA * cG + (1.0 - dA) * imgOut[4*i+1];
					imgOut[4*i+2] = dA * cB + (1.0 - dA) * imgOut[4*i+2];
				}
			}

			// Add text
			sctext.DrawString(
				timeGuam.ToString() + std::string(" (GMT+10)"),
				10,
				nImageSizeY-12,
				SchriftText::TextAlignment_Left,
				nImageSizeX,
				nImageSizeY,
				&(imgOut[0]),
				RGBA());

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

