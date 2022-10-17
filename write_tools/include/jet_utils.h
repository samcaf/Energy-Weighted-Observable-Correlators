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
#include <ostream>
#include <cmath>
#include <vector>
#include <limits>
#include <locale>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>

#include <utility>
#include <stdexcept>

#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <stdlib.h>

// ---------------------------------
// HEP imports
// ---------------------------------
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "general_utils.h"

using namespace fastjet;

// =====================================
// typedefs
// =====================================

// Type definitions for fastjet and pythia
typedef std::vector<PseudoJet> PseudoJets;


// =====================================
// Error utilities
// =====================================
std::string jetalg_error(const JetAlgorithm alg);

double pairwise_error();

// =====================================
// Pseudojet utilities
// =====================================

// ---------------------------------
// Particle-level utilities
// ---------------------------------

int is_nu_id(const int id);

// ---------------------------------
// Event-level utilities
// ---------------------------------
//
PseudoJets get_particles_pythia(const Pythia8::Event event);

PseudoJets add_events(const PseudoJets event1, const PseudoJets event2);

double SumScalarPt(const PseudoJets pjs);

double SumEnergy(const PseudoJets pjs);

// ---------------------------------
// Pair-level utilities
// ---------------------------------


double pair_theta(const PseudoJet pj1, const PseudoJet pj2);

double pair_mass(const PseudoJet pj1, const PseudoJet pj2);

double pair_cos(const PseudoJet pj1, const PseudoJet pj2);

// =====================================
// Visualization Utilities
// =====================================

extern std::string event_vis_header;
extern std::string event_vis_footer;

void write_ptyphi_particle(const PseudoJet particle,
                    std::ofstream& file,
                    std::string p_type = "P");

extern bool accept_all_nonempty(const PseudoJet pj);

void write_ptyphis_jets(const PseudoJets particles, const JetDefinition jet_def,
                        std::ofstream& file,
                        bool (*write_jet)(PseudoJet) = &accept_all_nonempty);

// =====================================
// Ghost Particle Utilities
// =====================================
// Base ghost parameters
extern const double _ghost_maxrap;
extern const double _ghost_midrap;
extern const double _ghost_area;
extern const double _mean_ghost_pt;
extern const Selector _no_selector;
// Parameters controlling randomness for the uniform ghosts
extern const double _grid_scatter;
extern const double _pt_scatter;
// Extra information for ghost identification
extern const int _ghost_index;
bool is_ghost(const PseudoJet particle);
bool is_ghost_jet(const PseudoJet jet);
bool not_ghost_jet(const PseudoJet jet);

// ---------------------------------
// Ghost grid setup
// ---------------------------------

PseudoJets uniform_ghosts(const double max_ghost_deltarap = _ghost_maxrap,
                          const double grid_mid_rap = _ghost_midrap,
                          const double ghost_area = _ghost_area,
                          const double mean_ghost_pt = _mean_ghost_pt,
                          const Selector selector = _no_selector,
                          const double grid_scatter = _grid_scatter,
                          const double pt_scatter = _pt_scatter);

// ---------------------------------
// Visualization with ghosts
// ---------------------------------

void write_ptyphis_jets_with_ghosts(const PseudoJets particles,
                                    const JetDefinition jet_def,
                                    std::ofstream& file,
                                    const std::string ghost_type = "active",
                                    bool write_ghost_jets = false);

#endif
