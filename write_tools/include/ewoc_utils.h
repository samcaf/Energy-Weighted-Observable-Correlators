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
// typedefs
// =====================================
//typedef unsigned long size_t;

// Type definitions for fastjet and pythia
//typedef std::vector<PseudoJet> PseudoJets;

// =====================================
// EWOC Storage Utilities
// =====================================

void store_event_subpair_info(const PseudoJets particles,
                           const JetAlgorithm jet_algorithm,
                           const double jetR,
                           const RecombinationScheme jet_recomb,
                           const JetAlgorithm subjet_algorithm,
                           const double subjetR,
                           const RecombinationScheme sub_recomb,
                           const double pt_min, const double pt_max,
                           std::ofstream& file);

void store_jet_subpair_info(const PseudoJet jet,
                         const JetAlgorithm subjet_algorithm,
                         const double subjetR,
                         const RecombinationScheme sub_recomb,
                         std::ofstream& file);

#endif
