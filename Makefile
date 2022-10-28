# Makefile for c++ files using Pythia and Fastjet libraries

#===========================================
# Include the configuration.
#===========================================
-include Makefile.inc
# Includes compiler information and library names,
# such as for Pythia8 and FastJet


#===========================================
# Rules for building files
#===========================================

#------------------------------------
# Main Functions
#------------------------------------
# - - - - - - - - - - - - - - -
# Rules:
# - - - - - - - - - - - - - - -
# Possible make targets (to be make with ```make [xxx]```)
.PHONY : write_tools plot_tools scripts ewocs

# Make all main EWOC functions by default
.DEFAULT_GOAL := ewocs

# - - - - - - - - - - - - - - -
# Steps:
# - - - - - - - - - - - - - - -
# Telling Make to compile C++ code with g++, using pythia and fastjet libraries
write_tools: $(FASTJET) $(PYTHIA) write_tools/src/write_ewocs.cc
	# =======================================================
	# Compiling c++ code for writing EWOC data:
	# =======================================================
	# Compiling to the executable ```write_tools/write_ewocs```
	$(CXX) write_tools/src/write_ewocs.cc write_tools/src/ewoc_cmdln.cc write_tools/src/ewoc_utils.cc write_tools/src/jet_utils.cc write_tools/src/general_utils.cc -o write_tools/write_ewocs $(CXX_COMMON);
	@printf "\n"

plot_tools:
	# =======================================================
	# Preparing plotting tools
	# =======================================================
	# Allowing EWOC plotting tools to be run as executables:
	chmod +x plot_tools/plot_ewocs;
	chmod +x plot_tools/event_vis;
	@printf "\n"

scripts:
	# =======================================================
	# Preparing scripts
	# =======================================================
	# Allowing EWOC writing and plotting scripts to be run as executables:
	chmod +x scripts/*;

# - - - - - - - - - - - - - - -
# Make all main functions
# - - - - - - - - - - - - - - -
ewocs: write_tools plot_tools scripts

#------------------------------------
# Plots
#------------------------------------
.PHONY : all_plots all_plots_ee show_plots show_plots_ee

# - - - - - - - - - - - - - - -
# Make all plots from scratch
# - - - - - - - - - - - - - - -
all_plots:
	./scripts/rerun_all.sh --generate_events true --overwrite_hists true --show_plots true

all_plots_ee:
	./scripts/rerun_all.sh --generate_events true --overwrite_hists true --show_plots true --ee

# - - - - - - - - - - - - - - -
# Show all plots (no computation)
# - - - - - - - - - - - - - - -
show_plots:
	./scripts/rerun_all.sh --generate_events false --overwrite_hists false --show_plots true

show_plots_ee:
	./scripts/rerun_all.sh --generate_events false --overwrite_hists false --show_plots true --ee

#------------------------------------
# Functions for Testing
#------------------------------------

# - - - - - - - - - - - - - - -
# Jet energy tests
# - - - - - - - - - - - - - - -
plot_jet_energy:
	# Plotting jet energy and pT spectra for QCD events

# - - - - - - - - - - - - - - -
# Momentum Conservation Tests
# - - - - - - - - - - - - - - -
.PHONY : test_momenta clean_test_momenta

test_momenta:
	# =======================================================
	# Compiling c++ code for testing momentum conservation
	# =======================================================
	$(CXX) write_tools/tests/src/test_momentum_conservation.cc write_tools/src/ewoc_cmdln.cc write_tools/src/ewoc_utils.cc write_tools/src/jet_utils.cc write_tools/src/general_utils.cc -o write_tools/tests/test_momentum_conservation $(CXX_COMMON);
	@printf "\n"
	# =======================================================
	# Testing momentum conservation in QCD events
	# =======================================================
	./write_tools/tests/test_momentum_conservation --process qcd --level parton -n 10
	@printf "\n"

clean_test_momenta:
	# =======================================================
	# Removing momentum conservation test files
	# =======================================================
	rm write_tools/tests/test_momenta
	@printf "\n"


# - - - - - - - - - - - - - - -
# All tests
# - - - - - - - - - - - - - - -
tests: test_momenta clean_test_momenta
