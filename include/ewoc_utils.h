/**
 * @file    ewoc_utils.h
 *
 * @brief   A utility header file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */
#ifndef EWOC_UTILS
#define EWOC_UTILS

// ---------------------------------
// Basic imports
// ---------------------------------
#include <iostream>
#include <string.h>
#include <cmath>
#include <limits>
#include <locale>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>


// ---------------------------------
// HEP imports
// ---------------------------------
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// Local imports:
#include "../include/general_utils.h"
#include "../include/jet_utils.h"

// EWOCs:
#include "../include/ewoc_utils.h"
#include "../include/ewoc_cmdln.h"

using namespace Pythia8;
using namespace fastjet;


// =====================================
// Global Flags
// =====================================

// Recombination schemes for jet finding
extern RecombinationScheme _jet_recomb_scheme;
extern RecombinationScheme _sub_recomb_scheme;


// =====================================
// typedefs
// =====================================
//typedef unsigned long size_t;

// Type definitions for fastjet and pythia
//typedef std::vector<PseudoJet> PseudoJets;

// =====================================
// Utility functions
// =====================================

// ---------------------------------
// Plotting/Labelling Utilities
// ---------------------------------
std::string ewoc_folder(int argc, char* argv[]);

std::string ewoc_file_label(int argc, char* argv[]);


// ---------------------------------
// EWOC Text Storage
// ---------------------------------

void store_event_subpair_info(PseudoJets particles,
                           JetAlgorithm jet_algorithm,
                           double jetR,
                           JetAlgorithm subjet_algorithm,
                           double subjetR,
                           double pt_min, double pt_max,
                           std::ofstream& file);

void store_jet_subpair_info(PseudoJet jet,
                         JetAlgorithm subjet_algorithm,
                         double subjetR,
                         std::ofstream& file);

#endif
