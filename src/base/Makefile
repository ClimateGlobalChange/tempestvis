# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTVISDIR= ../..

# Load Makefile framework. 
include $(TEMPESTVISDIR)/mk/framework.make

FILES= Announce.cpp \
       NetCDFUtilities.cpp \
       TimeObj.cpp \
	   NcFileVector.cpp \
       Variable.cpp \
	   DataOp.cpp \
       kdtree.cpp \
	   lodepng.cpp \
	   SimpleGrid.cpp \
	   GaussQuadrature.cpp \
	   LegendrePolynomial.cpp \
	   MeshUtilitiesFuzzy.cpp \
	   GridElements.cpp \
	   SchriftText.cpp \
	   schrift.cpp \
	   ColorMap.cpp \
	   ShpFile.cpp \
	   PNGImage.cpp

LIB_TARGET= libvisbase.a

.PHONY: all clean

# Build rules. 
all: $(LIB_TARGET)

$(LIB_TARGET): $(FILES:%.cpp=$(BUILDDIR)/%.o)
	rm -f $(LIB_TARGET)
	ar -cqs $(LIB_TARGET) build/*.o

# Clean rules.
clean:
	rm -f $(LIB_TARGET) 
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
