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
RecombinationScheme _jet_recomb_scheme = WTA_pt_scheme;
RecombinationScheme _sub_recomb_scheme = WTA_pt_scheme;

// =====================================
// Utility functions
// =====================================

// ---------------------------------
// Plotting/Labelling Utilities
// ---------------------------------
/**
* @brief: Returns a folder name associated with ewoc information
*         for the given params.
*
* @param: argc/argv          Command line input.
*
* @return: string            Folder name.
*/
std::string ewoc_folder(int argc, char* argv[]) {
    // Getting arguments from command line input
    std::string qcd_level     = level_cmdln(argc, argv);
    std::string process_str   = process_cmdln(argc, argv);
    int         jet_alg_int   = jetalg_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);

    // Making directory if it does not exist
    std::string process_folder = "output/" + process_str + "_" + qcd_level + "/";
    // https://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html
    // also, https://stackoverflow.com/a/21589609
    mkdir(process_folder.c_str(),
          S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    std::string jet_rad_folder = "jetR" + str_round(jet_rad, 1) + "/";
    jet_rad_folder = periods_to_hyphens(jet_rad_folder);

    mkdir((process_folder + jet_rad_folder).c_str(),
          S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    std::string jet_alg_folder =
        (jet_alg_int == 2) ? "aktjet/" : (
        (jet_alg_int == 1) ? "cajet/" : (
        (jet_alg_int == 0) ? "ktjet/" :
        "_ERRjet/"));

    std::string folder_name = process_folder + jet_rad_folder + jet_alg_folder;
    int new_dir = mkdir(folder_name.c_str(),
                        S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // If any of the above directories did not exist, let us know
    if (new_dir == 1)
        std::cout << "\n\nDirectory not found; "
                  << "making new directory to store EWOC files.\n";

    return folder_name;
}



/**
* @brief: Returns a label for a file associated with ewoc information
*         for the given params.
*
* @param: argc/argv          Command line input.
*
* @return: string            Label for the file.
*/
std::string ewoc_file_label(int argc, char* argv[]) {
    // Getting directory:
    std::string folder_name = ewoc_folder(argc, argv);

    // Getting arguments from command line input
    // Required basic settings
    int         n_events      = nevents_cmdln(argc, argv);
    int         sub_alg_int   = subalg_cmdln(argc, argv);
    double      sub_rad       = subrad_cmdln(argc, argv);
    // Optional basic settings
    double      pt_min        = ptmin_cmdln(argc, argv);
    double      pt_max        = ptmax_cmdln(argc, argv);

    double      E_cm          = Ecm_cmdln(argc, argv);

    // Advanced Settings
    double      frag_temp     = fragtemp_cmdln(argc, argv);

    // Filename:
    std::string sub_label =
        (sub_alg_int == 2) ? "_aktsub" : (
        (sub_alg_int == 1) ? "_casub" : (
        (sub_alg_int == 0) ? "_ktsub" :
        "_ERRsub"));

    std::string ewoc_file = "subR" + str_round(sub_rad, 2) + sub_label
                            + "_" + std::to_string(n_events) + "evts";

    // Optional basic arguments
    if (pt_min != -1)
        ewoc_file += "_ptmin"+str_round(pt_min, 1);
    if (pt_max != std::numeric_limits<double>::max())
        ewoc_file += "_ptmax"+str_round(pt_max, 1);
    if (E_cm != _Ecm_default)
        ewoc_file += "_Ecm"+str_round(E_cm, 1);

    // Optional advanced arguments
    if (frag_temp != _frag_temp_default) {
        ewoc_file += "_temp"+str_round(frag_temp, 1);
    }

    // Final filename
    ewoc_file = periods_to_hyphens(ewoc_file);
    ewoc_file += ".txt";

    return folder_name + ewoc_file;
}


// =====================================
// EWOC Storage Utilities
// =====================================

// ---------------------------------
// EWOC Text Storage
// ---------------------------------

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
void store_event_subpair_info(PseudoJets particles,
                           JetAlgorithm jet_algorithm,
                           double jetR,
                           JetAlgorithm subjet_algorithm,
                           double subjetR,
                           double pt_min, double pt_max,
                           std::ofstream& file) {
    // Finding jets with the given jet radius
    JetDefinition jet_def(jet_algorithm, jetR, _jet_recomb_scheme);
    ClusterSequence cluster_seq(particles, jet_def);
    // Getting all jets in the event with pt > pt_min
    PseudoJets all_jets = sorted_by_pt(cluster_seq.inclusive_jets(pt_min));

    // Ensuring pt < pt_max
    PseudoJets good_jets;
    for (auto jet : all_jets)
        if (jet.pt() < pt_max) good_jets.push_back(jet);

    for (auto jet : good_jets)
        //for (auto subjetR : get_subjetRs(jetR))
        store_jet_subpair_info(jet, subjet_algorithm, subjetR, file);

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
void store_jet_subpair_info(PseudoJet jet,
                         JetAlgorithm subjet_algorithm,
                         double subjetR,
                         std::ofstream& file) {
    // Finding subjets using the given subjet radius
    PseudoJets subjets;
    if (subjetR == 0)
        subjets = jet.constituents();
    else {
        JetDefinition subjet_def(subjet_algorithm, subjetR, _sub_recomb_scheme);
        ClusterSequence sub_cluster_seq(jet.constituents(), subjet_def);
        subjets = sorted_by_pt(sub_cluster_seq.inclusive_jets());
    }

    // Storing jet energy and number of subjets
    file << "J " << jet.E() << " " << subjets.size() << "\n";

    // Storing information associated with each subjet pair
    for (auto subjet1 : subjets) {
        for (auto subjet2 : subjets) {
        if (!equal_pjs(subjet1, subjet2)) {
            file << "SP "                                   // Subjet pair
                 << subjet1.E() / jet.E() << " "            // z1
                 << subjet2.E() / jet.E() << " "            // z2
                 << pair_cos(subjet1, subjet2) << "\n\n";   // costheta
                 // << pair_mass(subjet1, subjet2) << "\n";    // mass
        } }
    }
    file << "\n";
}
