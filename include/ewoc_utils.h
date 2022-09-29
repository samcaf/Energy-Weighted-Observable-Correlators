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
typedef unsigned long size_t;

// Type definitions for fastjet and pythia
typedef std::vector<PseudoJet> PseudoJets;

// =====================================
// Utility functions
// =====================================

// ---------------------------------
// General utilities
// ---------------------------------
extern std::string flag_text;

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1);

// ---------------------------------
// Plotting/Labelling Utilities
// ---------------------------------
std::string ewoc_folder(int jet_alg_int, double jet_rad,
                        std::string qcd_level,
                        std::string process_str);

std::string ewoc_file_label(int jet_alg_int, double jet_rad,
                            int sub_alg_int, double sub_rad,
                            int n_events,
                            std::string qcd_level,
                            std::string process_str);

// ---------------------------------
// Error utilities
// ---------------------------------
std::string jetalg_error(JetAlgorithm alg);

double pairwise_error();

// ---------------------------------
// Pseudojet utilities
// ---------------------------------

int is_nu_id(int id);

PseudoJets get_particles(Event event);

bool equal_pjs(PseudoJet pj1, PseudoJet pj2);

double SumScalarPt(PseudoJets pjs);

double SumEnergy(PseudoJets pjs);

double pair_mass(PseudoJet pj1, PseudoJet pj2);

double pair_cos(PseudoJet pj1, PseudoJet pj2);

double pair_theta(PseudoJet pj1, PseudoJet pj2);

// =====================================
// EWOC Storage Utilities
// =====================================

// ---------------------------------
// EWOC Text Storage
// ---------------------------------

void store_event_subpair_info(Event event,
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
