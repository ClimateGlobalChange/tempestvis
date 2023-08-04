# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTVISDIR= ../..

# Load Makefile framework. 
include $(TEMPESTVISDIR)/mk/framework.make

TEMPESTVISBASEDIR= $(TEMPESTVISDIR)/src/base
TEMPESTVISBASELIB= $(TEMPESTVISBASEDIR)/libvisbase.a

TEMPESTVISNETCDFDIR= $(TEMPESTVISDIR)/src/netcdf-cxx-4.2
TEMPESTVISNETCDFLIB= $(TEMPESTVISNETCDFDIR)/libnetcdf_c++.a

EXEC_FILES= RenderTC.cpp \
            RenderTCSnapshot.cpp

EXEC_TARGETS= $(EXEC_FILES:%.cpp=%)

FILES= $(UTIL_FILES) $(EXEC_FILES)

.PHONY: all clean

# Build rules. 
all: $(EXEC_TARGETS)

$(EXEC_TARGETS): %: $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/%.o $(TEMPESTVISBASELIB) $(TEMPESTVISNETCDFLIB)
	$(CXX) $(LDFLAGS) -o $@ $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/$*.o $(TEMPESTVISBASELIB) $(LIBRARIES)
	mv $@ $(TEMPESTVISDIR)/bin

$(TEMPESTVISBASELIB): %:
	cd $(TEMPESTVISBASEDIR); $(MAKE)

$(TEMPESTVISNETCDFLIB): %:
	cd $(TEMPESTVISNETCDFDIR); $(MAKE)

# Clean rules.
clean:
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
