/**
 * @file    subject_eecs.cc
 *
 * @brief   A file which is used to generate jet and subjet EEC information
 *          with Pythia 8.
 *          Version 0.2
 */

// ---------------------------------
// TO DO:
// ---------------------------------
/*
    * Analytic results
        * Eqn (8) of https://arxiv.org/pdf/1801.03219.pdf
*/

// ---------------------------------
// Basic imports
// ---------------------------------
#include <iostream>
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

#include <chrono>
using namespace std::chrono;

// ---------------------------------
// HEP imports
// ---------------------------------
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "include/ewoc_utils.h"
#include "include/ewoc_cmdln.h"

using namespace Pythia8;
using namespace fastjet;

// =====================================
// Global switches
// =====================================
// No global switches at the moment

// ####################################
// Main
// ####################################
/**
* @brief: Generates e+ e- -> hadrons events with Pythia and analyzes jet and subjet EECs.
*         (Outputs to files for later analysis with python).
*
* @return: int
*/
int main (int argc, char* argv[]){
    // Starting timer
    auto start = high_resolution_clock::now();

    // =====================================
    // Command line setup
    // =====================================
    // Ensure valid command line inputs
    if (checkCmdInputs(argc, argv) == 1) return 1;

    // ---------------------------------
    // Convert command line variables
    // ---------------------------------
    // Basic Settings
    int         n_events      = nevents_cmdln(argc, argv);
    std::string qcd_level     = level_cmdln(argc, argv);
    std::string process_str   = process_cmdln(argc, argv);
    int         jet_alg_int   = jetalg_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);
    int         sub_alg_int   = subalg_cmdln(argc, argv);
    double      sub_rad       = subrad_cmdln(argc, argv);

    // Optional Settings
    double      pt_min        = ptmin_cmdln(argc, argv);
    double      pt_max        = ptmax_cmdln(argc, argv);

    double      E_cm          = Ecm_cmdln(argc, argv);
    std::string s_channel     = schannel_cmdln(argc, argv);

    int         verbose       = verbose_cmdln(argc, argv);

    // Advanced Settings
    double      frag_temp     = fragtemp_cmdln(argc, argv);

    // Set up (sub)jet algorithm
    JetAlgorithm jet_alg = JetAlgorithm(jet_alg_int);
    JetAlgorithm sub_alg = JetAlgorithm(sub_alg_int);

    // Open up an appropriate file
    std::ofstream file;
    std::string filename = ewoc_file_label(argc, argv);
    std::cout << filename;
    file.open(filename);

    // Add header with additional information
    file << "# ==================================\n"
         << "# EWOC information\n"
         << "# ==================================\n"
         << "# Function call ```";
    while( argc-- ) file << *(argv++) << " ";
    file << "```\n\n"
         << "# Process:\n"
         << "process_str = " + process_str << "\n"
         << "# Level (parton/hadron):\n"
         << "level = " + qcd_level << "\n"
         << "# Number of events:\n"
         << "n_events = " + std::to_string(n_events) << "\n"
         << "# Jet information:\n"
         << "jet_alg = " + std::to_string(jet_alg_int) << "\n"
         << "jet_rad = " + std::to_string(jet_rad) << "\n"
         << "# Subjet information:\n"
         << "sub_alg = " + std::to_string(sub_alg_int) << "\n"
         << "sub_rad = " + std::to_string(sub_rad) << "\n\n\n"
         << "# ==================================\n"
         << "# Event output:\n"
         << "# ==================================\n"
         << "# Headers:\n"
         << "#     E:  Event header\n"
         << "#         E   |    event number\n"
         << "#     J:  Jet header\n"
         << "#         J   |    jet num. in event  |    pT\n"
         << "#     SP: Subjet pair header\n"
         << "#         SP  |    pair num. in jet   |   z1*z2  |    angle\n\n\n";

    // ---------------------------------
    // Pythia setup
    // ---------------------------------
    Pythia pythia;

    // Get default options
    pythia.readFile("ewoc_setup.cmnd");

    // Give process information to Pythia
    if (process_str == "quark") {
        std::cout << "pp ->  q qbar, \n";
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");
        pythia.readString("HardQCD:gg2gg = on");
    } else if (process_str == "gluon") {
        std::cout << "pp -> gg, \n";
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");
        pythia.readString("HardQCD:gg2qqbar = on");
    } else if (process_str == "qcd") {
        std::cout << "ee -> hadrons, \n";
        pythia.readString("Beams:idA = 11");
        pythia.readString("Beams:idB = -11");
        // gm, Z, or gmZ in the s-channel
        pythia.readString("WeakSingleBoson:ffbar2ffbar(s:"
                + s_channel + ") = on");
    } else if (process_str == "w") {
        std::cout << "ee -> W W, \n";
        pythia.readString("Beams:idA = 11");
        pythia.readString("Beams:idB = -11");
        pythia.readString("WeakDoubleBoson:ffbar2WW = on");
    }

    // Beam energy
    if (E_cm != _Ecm_default) {
        std::cout << "beam E_cm: "
                  << std::to_string(E_cm/1000.) << " TeV.\n";
        pythia.readString("Beams:eCM = " + std::to_string(E_cm));
    }
    else {
        std::cout << "beam E_cm: 4 TeV.\n";
        pythia.readString("Beams:eCM = 4000");
    }

    // Parton or hadron level
    if (qcd_level == "parton")
        pythia.readString("HadronLevel:all = off");

    // Cut on the phase to produce events with a jet satisfying min pt requirement
    pythia.readString("PhaseSpace:pTHatMin = "+std::to_string(pt_min-20));

    // Misc. options
    if (verbose < 2) {
        pythia.readString("Init:showProcesses = off");
        pythia.readString("Init:showMultipartonInteractions = off");
        pythia.readString("Init:showChangedSettings = off");
        pythia.readString("Init:showChangedParticleData = off");

        pythia.readString("Next:numberCount = 10000");
        pythia.readString("Next:numberShowLHA = 0");
        pythia.readString("Next:numberShowInfo = 0");
        pythia.readString("Next:numberShowProcess = 0");
        pythia.readString("Next:numberShowEvent = 0");

        pythia.readString("Stat:showProcessLevel = off");
        pythia.readString("Stat:showErrors = off");
    }

    // Advanced Options
    if (frag_temp != _frag_temp_default) {
        std::cout << "Turning on thermal model for string pT with temperature "
                  << std::to_string(frag_temp) << "\n";
        pythia.readString("StringPT:thermalModel = on");
        pythia.readString("StringPT:temperature = " + std::to_string(frag_temp));
    }

    pythia.init();

    // ---------------------------------
    // Analyzing events
    // ---------------------------------
    for (int iev = 0; iev < n_events; iev++) {
        // Considering next event, if valid
        if(!pythia.next()) continue;

        // Writing header for this event
        file << "# =========================\n"
             << "E " << std::to_string(iev+1) << "\n";

        // Storing EWOC info for this event in the output file
        store_event_subpair_info(pythia.event,
                              jet_alg, jet_rad,
                              sub_alg, sub_rad,
                              pt_min, pt_max,
                              file);

    } // end event loop

    file.close();

    // Verifying successful run
    if (verbose >= 0)
        std::cout << "\nComplete!\n";
    if (verbose >= 1) {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Analyzed and saved data from " + std::to_string(n_events)
        << " events in " << std::to_string(float(duration.count())/pow(10, 6)) << " seconds.";
    }

    return 0;
}
