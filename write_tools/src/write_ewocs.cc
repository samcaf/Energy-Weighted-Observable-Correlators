/**
 * @file    write_ewocs.cc
 *
 * @brief   A file which is used to generate jet and subjet EEC and EWOC
 *          information.
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
#include <string.h>
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

// Local imports:
#include "../include/general_utils.h"
#include "../include/jet_utils.h"

// EWOCs:
#include "../include/ewoc_utils.h"
#include "../include/ewoc_cmdln.h"

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
int main (int argc, char* argv[]) {
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
    // - - - - - - - - - - - - - - - - -
    // Basic Settings
    // - - - - - - - - - - - - - - - - -
    int         n_events      = nevents_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Jet Settings
    // - - - - - - - - - - - - - - - - -
    int         jet_alg_int   = jetalg_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);
    int         sub_alg_int   = subalg_cmdln(argc, argv);
    double      sub_rad       = subrad_cmdln(argc, argv);

    JetAlgorithm jet_alg = JetAlgorithm(jet_alg_int);
    JetAlgorithm sub_alg = JetAlgorithm(sub_alg_int);

    // - - - - - - - - - - - - - - - - -
    // Optional Settings
    // - - - - - - - - - - - - - - - - -
    double      pt_min        = ptmin_cmdln(argc, argv);
    double      pt_max        = ptmax_cmdln(argc, argv);

    int         verbose       = verbose_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Advanced Settings
    // - - - - - - - - - - - - - - - - -
    std::string event_gen     = eventgen_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Output Settings
    // - - - - - - - - - - - - - - - - -
    // Determines whether EWOC output is written
    bool        write_ewocs   = writeewocs_cmdln(argc, argv);
    // Set up EWOC output file
    std::ofstream ewoc_outfile;
    if (write_ewocs) {
        cout << "\nWriting EWOC data to "
             << ewoc_file_label(argc, argv);
        write_ewoc_header_cmdln(argc, argv);
        ewoc_outfile.open(ewoc_file_label(argc, argv), std::ios_base::app);
    }
    else cout << "\nNot writing EWOC data.\n";

    // Determines whether pT spectrum is written for given PID
    int         write_pid_pt  = writepidpt_cmdln(argc, argv);
    bool        seen_pid_pt   = false;
    // Set up PID pT output file
    std::ofstream pid_pt_outfile;

    // Determines whether event output is written
    bool        write_event   = writeevent_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Set up event generator
    // - - - - - - - - - - - - - - - - -
    // Basic Declarations
    Pythia8::Pythia pythia;

    // Setup
    if (str_eq(event_gen, "pythia")) {
        setup_pythia_ewoc_cmdln(pythia, argc, argv);
    }
    else if (str_eq(event_gen, "herwig"))
        throw Error("Invalid event generator herwig");
    else
        throw Error("Invalid event generator");

    // ---------------------------------
    // Analyzing events
    // ---------------------------------
    for (int iev = 0; iev < n_events; iev++) {
        // Initializing particles in this event
        PseudoJets particles;

        // Considering next event, if valid
        if(!pythia.next()) continue;
        Pythia8::Event event = pythia.event;

        if (str_eq(event_gen, "pythia")) {
            // Getting event and particles in pseudojet form
            particles = get_particles_pythia(event);
        }
        else if (str_eq(event_gen, "herwig"))
            throw Error("Invalid event generator herwig");
        else
            throw Error("Invalid event generator");


        // - - - - - - - - - - - - - - - - -
        // Writing EWOC Data
        // - - - - - - - - - - - - - - - - -
        if (write_ewocs) {
            // Writing header for this event
            ewoc_outfile << "# =========================\n"
                    << "E " << std::to_string(iev+1) << "\n";

            // Storing EWOC info for this event in the output file
            store_event_subpair_info(particles,
                                     jet_alg, jet_rad,
                                     sub_alg, sub_rad,
                                     pt_min, pt_max,
                                     ewoc_outfile);
        }

        // - - - - - - - - - - - - - - - - -
        // Writing pT associated with PID
        // - - - - - - - - - - - - - - - - -
        if (write_pid_pt != 0) {
            // Looking for particle with given PID in this event
            for (int ipart = 0; ipart < event.size(); ipart++) {
                if (event[ipart].id() == write_pid_pt) {
                    // Opening file if this is the first sighting
                    if (not seen_pid_pt) {
                        seen_pid_pt = true;
                        std::string pid_pt_filename = process_folder(argc, argv)
                            + std::to_string(write_pid_pt) + "-pt-spectrum.txt";
                        cout << "Writing pT spectrum for PID " << write_pid_pt
                             << " to " << pid_pt_filename << "\n";
                        pid_pt_outfile.open(pid_pt_filename);
                    }
                    // Storing pT
                    pid_pt_outfile << sqrt(pow(event[ipart].px(), 2.) + pow(event[ipart].py(), 2.))
                                   << "\n";
                }
            }
        }

        // - - - - - - - - - - - - - - - - -
        // Visualizing Subjets
        // - - - - - - - - - - - - - - - - -
        if (iev == 0 and write_event) {
            std::ofstream event_vis_file;
            std::string event_filename =
                remove_extension(ewoc_file_label(argc, argv))
                + "_vis-evt" + std::to_string(iev) + ".txt";
            event_vis_file.open(event_filename);

            if (sub_rad == 0){
                // Visualize the whole event
                JetDefinition jet_def(jet_alg, jet_rad, _jet_recomb_scheme);
                write_ptyphis_jets_with_ghosts(particles, jet_def,
                                               event_vis_file, "passive");
            }
            else {
                // Visualize the subjets of the leading jet in the event
                JetDefinition jet_def(jet_alg, jet_rad, _jet_recomb_scheme);
                ClusterSequence jet_cluster_seq(particles, jet_def);
                PseudoJets jets = sorted_by_pt(jet_cluster_seq.inclusive_jets());

                PseudoJet lead_jet = jets.at(0);
                JetDefinition sub_def(sub_alg, sub_rad, _sub_recomb_scheme);
                write_ptyphis_jets_with_ghosts(lead_jet.constituents(),
                                               sub_def, event_vis_file,
                                               "passive");
            }
            std::cout << "Printed event " << iev << " to " << event_filename << "\n";
            event_vis_file.close();
        }
    } // end event loop

    // Verifying successful run
    if (verbose >= 0)
        std::cout << "\nComplete!\n";
    if (verbose >= 1) {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Analyzed and saved data from " + std::to_string(n_events)
        << " events in " << std::to_string(float(duration.count())/pow(10, 6)) << " seconds.";
    }

    // - - - - - - - - - - - - - - - - -
    // Closing output files
    // - - - - - - - - - - - - - - - - -
    // EWOC Output
    if (write_ewocs) ewoc_outfile.close();
    // pT Output for given particle ID
    if (seen_pid_pt) pid_pt_outfile.close();
    else if (write_pid_pt != 0) cout << "No particles seen with the given PID (for pT analysis).\n";

    return 0;
}
