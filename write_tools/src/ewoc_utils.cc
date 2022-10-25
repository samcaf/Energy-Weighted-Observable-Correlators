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
// No global flags at the moment

// =====================================
// EWOC Storage Utilities
// =====================================

/**
* @brief: Stores event info associated with pairs of subjets
*         for a given jet/subjet algorithm into a file.
*
* @param: event     The Pythia Event under consideration.
* @param: (sub)jet_algorithm, (sub)jetR
*                   Parameters for (sub)jet finding.
* @param: min_pt, max_pt
*                   The minimum and maximum pt for jet finding.
* @param: file      The file to which we will write.
*
* @return: void
*/
void store_event_subpair_info(const PseudoJets particles,
                           const JetAlgorithm jet_algorithm,
                           const double jetR,
                           const RecombinationScheme jet_recomb,
                           const JetAlgorithm subjet_algorithm,
                           const double subjetR,
                           const RecombinationScheme sub_recomb,
                           const double pt_min, const double pt_max,
                           std::ofstream& file) {
    // Finding jets with the given jet radius
    JetDefinition jet_def(jet_algorithm, jetR, jet_recomb);
    ClusterSequence cluster_seq(particles, jet_def);
    // Getting all jets in the event with pt > pt_min
    PseudoJets all_jets = sorted_by_pt(cluster_seq.inclusive_jets(pt_min));

    // Ensuring pt < pt_max
    PseudoJets good_jets;
    for (auto jet : all_jets)
        if (jet.pt() < pt_max) good_jets.push_back(jet);

    for (auto jet : good_jets)
        //for (auto subjetR : get_subjetRs(jetR))
        store_jet_subpair_info(jet, subjet_algorithm, subjetR, sub_recomb, file);

    return;
}


/**
* @brief: Stores info associated with subjet pairs
*         for a given jet/subjet algorithm into a file.
*
* @param: jet   The PseudoJet, with consituents assumed, under
*               consideration.
* @param: subjet_algorithm, subjetR
*               Parameters for subjet finding.
* @param: file  The file to which we will write.
*
* @return: void
*/
void store_jet_subpair_info(const PseudoJet jet,
                         const JetAlgorithm subjet_algorithm,
                         const double subjetR,
                         const RecombinationScheme sub_recomb,
                         std::ofstream& file) {
    // Finding subjets using the given subjet radius
    PseudoJets subjets;
    if (subjetR == 0)
        subjets = sorted_by_pt(jet.constituents());
    else {
        JetDefinition subjet_def(subjet_algorithm, subjetR, sub_recomb);
        ClusterSequence sub_cluster_seq(jet.constituents(), subjet_def);
        subjets = sorted_by_pt(sub_cluster_seq.inclusive_jets());
    }

    // Storing jet energy and number of subjets
    file << "J " << jet.E() << " " << subjets.size() << "\n";

    // Storing information associated with each subjet pair
    for (auto subjet1 : subjets) {
        for (auto subjet2 : subjets) {
        if (subjet1 != subjet2) {
            file << "SP "                                   // Subjet pair
                 << subjet1.E() / jet.E() << " "            // z1
                 << subjet2.E() / jet.E() << " "            // z2
                 << pair_cos(subjet1, subjet2) << "\n\n";   // costheta
                 // << pair_mass(subjet1, subjet2) << "\n";    // mass
        } }
    }
    file << "\n";
}
