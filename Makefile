# Makefile for c++ files using Pythia and Fastjet libraries

#===========================================
# Include the configuration.
#===========================================
-include Makefile.inc

#===========================================
# Rules for building files
#===========================================
.PHONY : setup write_tools plot_tools scripts ewocs

.DEFAULT_GOAL := ewocs

# Telling Make to compile C++ code with g++, using pythia and fastjet libraries
write_tools: $(FASTJET) $(PYTHIA) write_tools/src/write_ewocs.cc
	# Compiling c++ code for writing EWOC data to the executable ```write_tools/write_ewocs```:
	$(CXX) write_tools/src/write_ewocs.cc write_tools/src/ewoc_cmdln.cc write_tools/src/ewoc_utils.cc write_tools/src/jet_utils.cc write_tools/src/general_utils.cc -o write_tools/write_ewocs $(CXX_COMMON);

plot_tools:
	# Allowing EWOC plotting tools to be run as executables:
	chmod +x plot_tools/plot_ewocs;
	chmod +x plot_tools/scripts/*;

scripts:
	# Allowing EWOC writing and plotting scripts to be run as executables:
	chmod +x scripts/*;

ewocs: write_tools plot_tools scripts
