# Makefile for c++ files using Pythia and Fastjet libraries

#===========================================
# Include the configuration.
#===========================================
-include Makefile.inc

#===========================================
# Rules for building files
#===========================================
# Telling Make to compile with g++, using pythia and fastjet libraries
# ```make XXX.cc``` produces the executable ```XXX``, called with ```./XXX```
%: $(FASTJET) $(PYTHIA) %.cc
	$(CXX) $@.cc src/ewoc_cmdln.cc src/ewoc_utils.cc src/jet_utils.cc src/general_utils.cc -o $@ $(CXX_COMMON)
