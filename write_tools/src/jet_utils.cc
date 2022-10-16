/**
 * @file    ewoc_utils.cc
 *
 * @brief   A utility file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */

// ---------------------------------
// Basic imports
// ---------------------------------
#include <iostream>
#include <ostream>
#include <string>
#include <string.h>
#include <cmath>
#include <limits>
#include <locale>
#include <string>
#include <vector>
#include <stdexcept>

#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>

// Random numbers for ghost grids
#include <stdlib.h>

// ---------------------------------
// HEP imports
// ---------------------------------
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "../include/general_utils.h"
#include "../include/jet_utils.h"

using namespace fastjet;

// =====================================
// Error utilities
// =====================================
/**
* @brief: Raises an error if an invalid pairwise observable name is specified.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
std::string jetalg_error(const JetAlgorithm alg) {
    throw std::invalid_argument( "Invalid jet algorithm " + std::to_string(alg) );
    return "ERR";
}


/**
* @brief: Raises an error if an invalid pairwise observable name is specified.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
double pairwise_error() {
    throw std::invalid_argument( "Invalid pairwise observable." );
    return 0;
}


// =====================================
// Pseudojet utilities
// =====================================

// ---------------------------------
// Particle-level utilities
// ---------------------------------

/**
* @brief: Function which tests whether the given int is a neutrino
*         particle id.
*
* @param: int id    Given particle id.
*
* @return: bool     True if given int is a neutrino particle id,
*                   false otherwise.
*/
int is_nu_id(const int id) {
    return (id==12 || id==14 || id==16 || id==-12 || id==-14 || id==-16);
}


// ---------------------------------
// Event-level utilities
// ---------------------------------

/**
* @brief: Takes in a Pythia Event and return a vector of PseudoJets
*         containing all particles of a current event.
*
* @param: event     Pythia Event type, as in pythia.event().
*                   The event under consideration.
*
* @return: vector<PseudoJet>    A vector containing all particles in the event.
*/
PseudoJets get_particles_pythia(const Pythia8::Event event) {
    // Storing particles of an event as PseudoJets
    PseudoJets particles;

    for (int ipart = 0; ipart < event.size(); ipart++) {
        // Finding visible final states and converting to pseudojets
        if (event[ipart].isFinal() and !is_nu_id(event[ipart].id())){
        PseudoJet particle(event[ipart].px(), event[ipart].py(),
                           event[ipart].pz(), event[ipart].e());
        particle.set_user_index(ipart);
        particles.push_back(particle);
    }}

    return particles;
}


PseudoJets add_events(const PseudoJets event1, const PseudoJets event2) {
    PseudoJets full_event = event1;
    for (auto particle : event2)
        full_event.push_back(particle);

    return full_event;
}


/**
* @brief: Returns the sum of scalar pT in a vector of pseudojets.
*
* @param: pjs   A vector containing the pseudojets to consider.
*
* @return: double The total energy.
*/
double SumScalarPt(const PseudoJets pjs) {
    // Finding total jet energy
    double pT = 0;
    for (auto pj : pjs) pT += pj.perp();
    return pT;
}


/**
* @brief: Returns the total energy in a vector of pseudojets.
*
* @param: pjs   A vector containing the pseudojets to consider.
*
* @return: double The total energy.
*/
double SumEnergy(const PseudoJets pjs) {
    // Finding total jet energy
    double Q = 0;
    for (auto pj : pjs) Q += pj.E();
    return Q;
}


// ---------------------------------
// Pair-level utilities
// ---------------------------------

/**
* @brief: Returns the angle between a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double   Angle between the pseudojets (in radians).
*/
double pair_theta(const PseudoJet pj1, const PseudoJet pj2) {
    return acos(pair_cos(pj1, pj2));
}


/**
* @brief: Returns mass of a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
double pair_mass(const PseudoJet pj1, const PseudoJet pj2) {
    return sqrt(pj1.m2() + pj2.m2() +
                pj1.E()*pj2.E() - pj1.px()*pj2.px() -
                pj1.py()*pj2.py() - pj1.pz()*pj2.pz());
}


/**
* @brief: Returns the cosine of the angle between a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double   Cosine of the angle between the pseudojets.
*/
double pair_cos(const PseudoJet pj1, const PseudoJet pj2) {
    double cos = (pj1.px()*pj2.px() + pj1.py()*pj2.py() + pj1.pz()*pj2.pz())/(sqrt(pj1.modp2())*sqrt(pj2.modp2()));
    if (cos > 1.0) cos = 1.0; // Need to handle case of numerical overflow
    return cos;
}


// =====================================
// Visualization Utilities
// =====================================

std::string event_vis_header = R"(
# ####################################
# Event Visualization File
# ####################################
# Format:       Header  |   pt  |   y   |   phi |   index
# 1st column:               Header describing information type
#       (J: Jet, P: Particle, SJ: Subjet, G: Ghost)
# 2nd, 3rd, 4th column:     pt, y, phi
# 5th (final) column:       Jet index
#       (-1: no info; N: associated with Nth-hardest jet)
)";


void write_ptyphi_pj(const PseudoJet pj,
                     std::ofstream& file,
                     std::string p_type) {
    file << p_type + " " << pj.pt() << " " << pj.rap()
        << " " << pj.phi();
    if (pj.user_index() != -1) file << " " << pj.user_index();
    file << "\n";
}


void write_ptyphis_event(const PseudoJets particles,
                         std::ofstream& file) {
    file << event_vis_header;
    file << "\n# Storing particles associated with this event.\n\n";
    for (auto part : particles) {
        // Writing particle information for each particle
        write_ptyphi_pj(part, file, (!is_ghost(part) ? "P" : "G"));
    }
}


bool accept_all_nonempty(const PseudoJet pj) {
    if (pj.has_constituents())
        return true;
    return false;
}


void write_ptyphis_jets(const PseudoJets particles,
                        const JetDefinition jet_def,
                        std::ofstream& file,
                        bool (*write_jet)(PseudoJet)) {
    file << event_vis_header;
    file << "\n# Storing particles and jet information associated with this event.\n";
    file << "\n# Using jet definition with description:\n# "
         << jet_def.description() << "\n\n";
    // Finding jets
    ClusterSequence sub_cluster_seq(particles, jet_def);
    PseudoJets jets = sorted_by_pt(sub_cluster_seq.inclusive_jets());

    for (size_t ijet=0; ijet < jets.size(); ijet++) {
        // Looping over jets, hardest first
        PseudoJet jet = jets.at(ijet);

        // If this jet does not pass a selection criterion
        // (Default selection: must have constituents)
        if (!write_jet(jet)) continue;

        jet.set_user_index(ijet+1);
        write_ptyphi_pj(jet, file, "J");
        for (auto part : jet.constituents()) {
            // Storing info about which jet includes this particle
            std::string p_type = (!is_ghost(part) ? "P" : "G");
            part.set_user_index(ijet+1);
            // Writing particle information for each particle
            write_ptyphi_pj(part, file, p_type);
        }
    }
}


// =====================================
// Ghost Particle Utilities
// =====================================

// Base ghost parameters
const double _ghost_maxrap = 5;               // maximum rapidity of the ghosts
const double _ghost_midrap = 0;               // middle of the ghost rapidity grid
const double _mean_ghost_pt = 1e-100;         // mean ghost pt
const double _ghost_area = 0.005;              // changes density of ghosts
const Selector _no_selector = Selector();     // Selector which accepts all particles

// Parameters controlling randomness for the uniform ghosts
const double _grid_scatter = 1.0;             // grid scatter
const double _pt_scatter = 0.1;               // pt scatter

// Extra information for ghost identification
const int _ghost_index = -1000;               // User index for ghost particles
                                              //
bool is_ghost(const PseudoJet particle) {
    return (particle.user_index() == _ghost_index);
}

bool is_ghost_jet(const PseudoJet jet) {
    for (auto part : jet.constituents())
        if (!is_ghost(part)) return false;
    return true;
}

bool not_ghost_jet(const PseudoJet jet) { return !is_ghost_jet(jet); }


// ---------------------------------
// Ghost grid setup
// ---------------------------------

PseudoJets uniform_ghosts(const double max_ghost_deltarap,
                          const double grid_mid_rap,
                          const double ghost_area,
                          const double mean_ghost_pt,
                          const Selector selector,
                          const double grid_scatter,
                          const double pt_scatter) {
    // Initializing ghosts
    PseudoJets ghosts;

    // - - - - - - - - - - - - - - - - -
    // Setting up grid
    // - - - - - - - - - - - - - - - - -
    // Determining grid parameters
    double drap = sqrt(ghost_area);
    double dphi = drap;

    // Using nearest int for number of grid rows/cols
    double nrap = int(max_ghost_deltarap/drap);
    int nphi = int(twopi/dphi);

    // Re-evaluating grid spacing and parameters
    drap = max_ghost_deltarap/nrap;
    dphi = twopi/nphi;
    // ghost_area = dphi * drap;  // not while ghost_area is const

    // int n_ghosts = (2*nrap+1)*nphi;

    // - - - - - - - - - - - - - - - - -
    // Adding to the list of ghosts
    // - - - - - - - - - - - - - - - - -
    // Iterating over grid
    for (int irap = -nrap; irap <= nrap-1; irap++) {
        for (int iphi = 0; iphi < nphi; iphi++) {
            // Grid points and pT, plus random offsets
            double phi = (iphi+0.5) * dphi + dphi*(rand()%1 - 0.5)*grid_scatter;
            double rap = (irap+0.5) * drap + drap*(rand()%1 - 0.5)*grid_scatter
                                                               + grid_mid_rap;
            double pt = mean_ghost_pt*(1 + (rand()%1 - 0.5)*pt_scatter);

            // Initializing ghost particle
            double exprap = exp(+rap);
            double pminus = pt/exprap;
            double pplus  = pt*exprap;
            double px = pt*cos(phi);
            double py = pt*sin(phi);
            PseudoJet mom(px,py,0.5*(pplus-pminus),0.5*(pplus+pminus));
            mom.set_cached_rap_phi(rap,phi);
            mom.set_user_index(_ghost_index);

            // If our particle does not pass the condition placed by an active selector
            if (selector.worker().get() && !selector.pass(mom)) continue;

            // Add the particle at this grid point to our list of ghost particles
            ghosts.push_back(mom);
        }
    }

    return ghosts;
}


// ---------------------------------
// Visualization with ghosts
// ---------------------------------

void write_ptyphis_jets_with_ghosts(const PseudoJets particles,
                                    const JetDefinition jet_def,
                                    std::ofstream& file,
                                    std::string ghost_type,
                                    bool write_ghost_jets) {
    // Ensuring that the ghost type is known
    if (!str_eq(ghost_type, "passive") and !str_eq(ghost_type, "active"))
        throw Error("Invalid ghost_type " + ghost_type + "for function "
                + "```write_ptyphis_jets_with_ghosts```");

    // - - - - - - - - - - - - - - - - -
    // Active Ghosts
    // - - - - - - - - - - - - - - - - -
    if (str_eq(ghost_type, "active")) {
        PseudoJets active_ghost_event = add_events(particles, uniform_ghosts());
        write_ptyphis_jets(active_ghost_event, jet_def, file,
                write_ghost_jets ? &accept_all_nonempty : &not_ghost_jet);
        return;
    }

    // - - - - - - - - - - - - - - - - -
    // Passive Ghosts
    // - - - - - - - - - - - - - - - - -
    if (str_eq(ghost_type, "passive")) {
        // Start by writing full event
        write_ptyphis_jets(particles, jet_def, file);

        // Write ghosts one by one to event
        for (auto ghost : uniform_ghosts()) {
            // Adding single ghost to event (defn of "passive ghosts")
            PseudoJets passive_ghost_event = particles;
            passive_ghost_event.push_back(ghost);

            // Finding jets in passive ghost event
            ClusterSequence sub_cluster_seq(passive_ghost_event, jet_def);
            PseudoJets jets = sorted_by_pt(sub_cluster_seq.inclusive_jets());

            // Looping over jets, hardest first
            size_t ijet = 0;
            bool found_ghost = false;
            while (ijet < jets.size() and !found_ghost) {
                PseudoJet jet = jets.at(ijet);

                // If we only look at non-ghost jets
                if (is_ghost_jet(jet) and !write_ghost_jets) {ijet++; continue;}

                // Finding where the ghost was clustered manually
                // (std::find doesn't play well with fastjet classes)
                for (auto part : jet.constituents()) {
                    if (!is_ghost(part)) continue;
                    // Writing to file
                    ghost.set_user_index(ijet+1);
                    write_ptyphi_pj(ghost, file, "G");
                    break;
                }

                ijet++;
            }
        }
        return;
    }
}
