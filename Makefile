###################################################################
# This Makefile was created using the bat-project script.
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -O2 -Wall -fPIC -Wno-deprecated 
LD           = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ld
LDFLAGS      = -g -O2 

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += `bat-config --cflags`
LIBS := `bat-config --libs`

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS_all      = \
	efficiency_simulation.cxx

CXXSRCS1      = \
        efficiency_simulation.cxx 

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#
CXXOBJS_all      = $(patsubst %.cxx,%.o,$(CXXSRCS_all))

CXXOBJS1      = $(patsubst %.cxx,%.o,$(CXXSRCS1))

MYPROGS_all     = \
	efficiency_simulation 

MYPROGS1     = \
        efficiency_simulation


GARBAGE = $(CXXOBJS_all) *~ link.d $(MYPROGS_all)

# targets
all : efficiency_simulation

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS_all))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS_all) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm -f $(GARBAGE)

efficiency_simulation : $(CXXOBJS1)
	$(CXX) $(LDFLAGS) $(CXXOBJS1) $(LIBS) -o $@



print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS_all)
	@echo c++ objs  : $(CXXOBJS_all)
	@echo c++ flags : $(CXXFLAGS)
	@echo libs      : $(LIBS)



