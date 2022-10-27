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
    std::string jet_alg       = jetalgstr_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);
    std::string sub_alg       = subalgstr_cmdln(argc, argv);
    double      sub_rad       = subrad_cmdln(argc, argv);

    // Optional: recombination schemes
    RecombinationScheme jet_recomb = jetrecomb_cmdln(argc, argv);
    RecombinationScheme sub_recomb = subrecomb_cmdln(argc, argv);

    // Getting jet definitions:
    JetDefinition jet_def = process_JetDef(jet_alg, jet_rad,
                                           jet_recomb);
    JetDefinition sub_def = process_JetDef(sub_alg, sub_rad,
                                           sub_recomb);

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


    // ---------------------------------
    // Output Setup
    // ---------------------------------
    // - - - - - - - - - - - - - - - - -
    // Determines whether EWOC output is written
    // - - - - - - - - - - - - - - - - -
    bool        write_ewocs   = writeewocs_cmdln(argc, argv);
    std::ofstream ewoc_outfile;
    // Set up EWOC output file
    if (write_ewocs) {
        cout << "\nWriting EWOC data to "
             << ewoc_file_label(argc, argv);
        write_ewoc_header_cmdln(argc, argv);
        ewoc_outfile.open(ewoc_file_label(argc, argv), std::ios_base::app);
    }
    else cout << "\nNot writing EWOC data.\n";

    // - - - - - - - - - - - - - - - - -
    // Determines whether pT spectrum is written for given PID
    // - - - - - - - - - - - - - - - - -
    int         write_pid_pt  = writepidpt_cmdln(argc, argv);
    bool        seen_pid_pt   = false;
    // Set up PID pT output file
    std::ofstream pid_pt_outfile;

    // - - - - - - - - - - - - - - - - -
    // Determines whether jet and subjet pT spectrum is written
    // - - - - - - - - - - - - - - - - -
    bool        write_jet_pt  = writejetpt_cmdln(argc, argv);
    std::ofstream jet_pt_outfile;
    std::ofstream subjet_pt_outfile;
    // Set up output files
    if (write_jet_pt) {
        // Jet pT Spectrum
        std::string jet_pt_filename = ewoc_folder(argc, argv) + "jet-pT-spectrum"
            + "_min"+str_round(pt_min, 0) + "_max"+str_round(pt_max, 0) + ".txt";
        cout << "\nWriting jet pt spectrum to " << jet_pt_filename;
        jet_pt_outfile.open(jet_pt_filename);
        jet_pt_outfile << "# Jet pT spectrum for " << jet_alg
                       << " jets with R = " << jet_rad;

        // Subjet pT Spectrum
        std::string subjet_pt_filename = ewoc_folder(argc, argv) + "subjet-pT-spectrum_"
            + alg_label[sub_alg] + "_r"+str_round(sub_rad, 2)
            + "_min" + +"_max";
        subjet_pt_filename = periods_to_hyphens(subjet_pt_filename) + ".txt";
        cout << "\nWriting subjet pt spectrum to " << subjet_pt_filename;
        subjet_pt_outfile.open(subjet_pt_filename);
        subjet_pt_outfile << "# Subjet pT spectrum for " << sub_alg
                          << " subjets with r = " << str_round(sub_rad, 2);
    }

    // - - - - - - - - - - - - - - - - -
    // Determines whether event output is written
    // - - - - - - - - - - - - - - - - -
    std::string write_event   = writeevent_cmdln(argc, argv);
    // Making files which point to event visualization files and
    //  run the visualization. Allows easier command line interface.
    std::ofstream event_vis_pointer;
    std::ofstream event_vis_script;
    if (not str_eq(write_event, "")) {
        event_vis_pointer.open("event_vis_pointer.txt", std::ios_base::app);
        event_vis_script.open("event_vis_script.sh", std::ios_base::app);
    }

    // ---------------------------------
    // Event Generator Setup
    // ---------------------------------
    // Declarations (muting Pythia banner)
    std::streambuf *old = cout.rdbuf();
    stringstream ss; ss.str("");
    if (verbose < 3) cout.rdbuf (ss.rdbuf());  // Redirect output

    Pythia8::Pythia pythia;  // Declaring Pythia8

    cout.rdbuf (old);  // Restore output

    // Setup
    if (str_eq(event_gen, "pythia"))
        setup_pythia_ewoc_cmdln(pythia, argc, argv);
    else if (str_eq(event_gen, "herwig"))
        throw Error("Invalid event generator herwig");
    else
        throw Error("Invalid event generator");

    // ---------------------------------
    // Analyzing events
    // ---------------------------------
    for (int iev = 0; iev < n_events; iev++) {
        // Initializing particles and info for this event
        PseudoJets particles;

        int num_narrow_emissions = 0;

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
        // Writing Jet Data
        // - - - - - - - - - - - - - - - - -
        // Writing headers for this event
        if (ewoc_outfile.is_open()) {
            ewoc_outfile << "# =========================\n"
                    << "E " << std::to_string(iev+1) << "\n";
        }

        // - - - - - - - - - - - - - - - - -
        // Storing EWOC or pT info for this event in the output file
        std::vector<int> narrow_emission_info;
        narrow_emission_info = store_event_subpair_info(particles,
                                                 jet_def, sub_def,
                                                 pt_min, pt_max,
                                                 ewoc_outfile,
                                                 jet_pt_outfile,
                                                 subjet_pt_outfile,
                                                 "num_narrow_emissions");

        // Counting the number of narrow emissions in the leading jet
        if (narrow_emission_info.size() > 0)
            num_narrow_emissions = narrow_emission_info.at(0);
        // - - - - - - - - - - - - - - - - -


        // - - - - - - - - - - - - - - - - -
        // Visualizing Subjets
        bool visualize_event = (str_eq(write_event, "last") and iev == n_events-1)
            or (str_eq(write_event, "narrow_emissions") and num_narrow_emissions > 0);

        if (visualize_event) {
            // - - - - - - - - - - - - - - - - -
            // Saving data for event visualization
            // - - - - - - - - - - - - - - - - -
            std::ofstream event_vis_file;
            std::string event_filename =
                remove_extension(ewoc_file_label(argc, argv))
                + "_vis-evt" + std::to_string(iev) + ".txt";
            event_vis_file.open(event_filename);

            if (sub_rad == 0){
                // Visualize the whole event
                write_ptyphis_jets_with_ghosts(particles, jet_def,
                                               event_vis_file, "passive");
            }
            else {
                // Visualize the subjets of the leading jet in the event
                ClusterSequence jet_cluster_seq(particles, jet_def);
                PseudoJets jets = sorted_by_pt(jet_cluster_seq.inclusive_jets());

                PseudoJet lead_jet = jets.at(0);
                write_ptyphis_jets_with_ghosts(lead_jet.constituents(),
                                               sub_def, event_vis_file,
                                               "passive");
            }
            std::cout << "Printed event " << iev << " to " << event_filename << "\n";
            event_vis_file.close();


            // - - - - - - - - - - - - - - - - -
            // Scripting for plot generation
            // - - - - - - - - - - - - - - - - -
            event_vis_pointer << event_filename << " ";
            event_vis_script << "./plot_tools/event_vis --filename " << event_filename
                             << " scatter_vis\n";
        }
        // - - - - - - - - - - - - - - - - -


        // - - - - - - - - - - - - - - - - -
        // Writing pT associated with PID
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
    } // end event loop

    // Verifying successful run
    if (verbose >= 0)
        std::cout << "\nComplete!\n";
    if (verbose >= 1) {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Analyzed and saved data from " + std::to_string(n_events)
        << " events in " << std::to_string(float(duration.count())/pow(10, 6)) << " seconds.\n";
    }
    if (verbose >= 2) {
        if (str_eq(event_gen, "pythia"))
            pythia.stat();
    }

    // - - - - - - - - - - - - - - - - -
    // Closing output files
    // - - - - - - - - - - - - - - - - -
    // EWOC output
    if (write_ewocs) ewoc_outfile.close();

    // pT output for given particle ID
    if (seen_pid_pt) pid_pt_outfile.close();
    else if (write_pid_pt != 0) cout << "No particles seen with the given PID (for pT analysis).\n";

    // Jet pT spectrum output
    if (write_jet_pt) {
        jet_pt_outfile.close();
        subjet_pt_outfile.close();
    }

    // Event visualization output
    if (not str_eq(write_event, "")) {
        event_vis_pointer.close();
        event_vis_script.close();
    }

    return 0;
}
