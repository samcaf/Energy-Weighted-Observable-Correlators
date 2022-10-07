/**
 * @file    jet_utils.h
 *
 * @brief   A utility header file for jets.
 */
#ifndef JET_UTILS
#define JET_UTILS

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

using namespace fastjet;


// =====================================
// typedefs
// =====================================

// Type definitions for fastjet and pythia
typedef std::vector<PseudoJet> PseudoJets;


// =====================================
// Utility functions
// =====================================

// ---------------------------------
// Error utilities
// ---------------------------------
std::string jetalg_error(JetAlgorithm alg);

double pairwise_error();

// ---------------------------------
// Pseudojet utilities
// ---------------------------------

int is_nu_id(int id);

PseudoJets get_particles_pythia(Pythia8::Event event);

bool equal_pjs(PseudoJet pj1, PseudoJet pj2);

double SumScalarPt(PseudoJets pjs);

double SumEnergy(PseudoJets pjs);

double pair_mass(PseudoJet pj1, PseudoJet pj2);

double pair_cos(PseudoJet pj1, PseudoJet pj2);

double pair_theta(PseudoJet pj1, PseudoJet pj2);

#endif
