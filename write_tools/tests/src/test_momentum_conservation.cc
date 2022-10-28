/**
 * @file    write_ewocs.cc
 *
 * @brief   A file which is used to generate jet and subjet EEC and EWOC
 *          information.
 */

// =====================================
// Global switches
// =====================================

const int print_every = 1;
const int base_verbose = 0;


// =====================================
// Imports
// =====================================

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
#include "../../include/general_utils.h"
#include "../../include/jet_utils.h"

// EWOCs:
#include "../../include/ewoc_utils.h"
#include "../../include/ewoc_cmdln.h"


/**
* @brief: Tests for momentum conservation (p^\mu = 0)
*         in the list of particles `particles`.
*
* @param: particles         A vector<PseudoJet> which
*                           we want to test for momentum
*                           conservation.
* @param: verbose           How talky is function?
*
* @return: double           A measure of the degree of
*                           the degree of violation:
*                             √(Ptot^mu Ptot^mu),
*                           the Euclidean norm.
*/
double test_momentum_conservation(PseudoJets particles,
                                  int verbose = 0) {
    // Initializing total energy and momenta
    double Etot = 0;
    double pxtot = 0;
    double pytot = 0;
    double pztot = 0;

    for (auto part : particles) {
        Etot  += part.E();
        pxtot += part.px();
        pytot += part.py();
        pztot += part.pz();
    }

    if (verbose > 1)
    cout << "\t# == Total 4 Momentum == #\n"
         << "\t# -----------------------------"
         << "--------------------------------- #\n"
         << "\t|\tEtot\t|\tpx\t|\tpy\t|\tpz\t|\n"
         << "\t|\t" << str_round(Etot, 2) << "\t|\t"
         << str_round(pxtot, 2) << "\t|\t"
         << str_round(pytot, 2) << "\t|\t"
         << str_round(pztot, 2) << "\t|\n"
         << "\t# -----------------------------"
         << "--------------------------------- #\n\n";

    if (Etot!=0 or pxtot!=0 or pytot!=0 or pztot!=0) {
        if (verbose > 0)
            cout << "\n\n\n\t# ! # ! # ! # ! # ! # ! # ! # ! # ! #\n"
                 << "\t# ! 4-Mom. Conservation Violation ! #\n"
                 << "\t# ! # ! # ! # ! # ! # ! # ! # ! # ! #\n\n\n\n";
        return sqrt(pow(Etot,2) + pow(pxtot,2) + pow(pytot,2) + pow(pztot,2));
    }


    return 0;
}

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
    double      E_cm          = Ecm_cmdln(argc, argv);
    double      pt_min        = ptmin_cmdln(argc, argv);
    double      pt_max        = ptmax_cmdln(argc, argv);

    int         verbose       = verbose_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Advanced Settings
    // - - - - - - - - - - - - - - - - -
    std::string event_gen     = eventgen_cmdln(argc, argv);

    // - - - - - - - - - - - - - - - - -
    // Set up event generator
    // - - - - - - - - - - - - - - - - -
    // Declarations (muting banners)
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
    // 4-Momentum of the entire event
    PseudoJet event_4mom(0, 0, 0, -E_cm);
    for (int iev = 0; iev < n_events; iev++) {
        // Initializing particles and info for this event
        PseudoJets particles;

        // Considering next event, if valid
        if(!pythia.next()) continue;
        Pythia8::Event event = pythia.event;

        particles = get_particles_pythia(event);

        // - - - - - - - - - - - - - - - - -
        // Full Event
        // - - - - - - - - - - - - - - - - -
        cout << "\n# ##############################"
             << " Event Level "
             << "############################### #\n";

        cout << "\n# ==========================="
             << " All Particles "
             << "============================ #\n";
        // Momentum conservation in the whole event
        /* if ((iev+1)%print_every == 0) { */
        /*     cout << "\nEvent " << iev + 1 << "\n"; */
        /*     cout << "\t# ###   Particles:   ### #\n"; } */
        PseudoJets balanced_particles = particles;
        balanced_particles.push_back(event_4mom);
        test_momentum_conservation(balanced_particles,
               base_verbose + int((iev+1)%print_every == 0));


        cout << "\n# ==========================="
             << " QCD Particles "
             << "============================ #\n";
        PseudoJets qcd_particles = get_particles_pythia(event, qcd_pids);
        qcd_particles.push_back(event_4mom);
        test_momentum_conservation(qcd_particles,
               base_verbose + 2*int((iev+1)%print_every == 0));

        // - - - - - - - - - - - - - - - - -
        // Jets
        // - - - - - - - - - - - - - - - - -
        cout << "\n# ##############################"
             << " Jet Level "
             << "############################### #\n";

        // Muting FastJet banner
        cout.rdbuf (ss.rdbuf());  // Redirect output
        ClusterSequence cluster_seq(particles, jet_def);
        cout.rdbuf (old);  // Restore output

        // Momentum conservation of all jets in the event: ( 1 )
        cout << "\n# ==========================="
             << " All Jets "
             << "============================ #\n";
        PseudoJets all_jets = sorted_by_pt(cluster_seq.inclusive_jets());
        all_jets.push_back(event_4mom);
        test_momentum_conservation(all_jets,
               base_verbose + 2*int((iev+1)%print_every == 0));

        // Momentum conservation of all jets with pt > pt_min ( 2 )
        cout << "\n# ==========================="
             << " Hard Jets "
             << "============================ #\n";
        PseudoJets hard_jets = sorted_by_pt(cluster_seq.inclusive_jets(pt_min));
        hard_jets.push_back(event_4mom);
        test_momentum_conservation(hard_jets,
               base_verbose + 2*int((iev+1)%print_every == 0));

        // =====================================
        // TODO:
        // =====================================

        // Add bins to plotting options -- add to npz file labels

        // Momentum conservation within jets: Particle level: ( 4 )
        // Cluster

        // - - - - - - - - - - - - - - - - -
        // Subjets
        // - - - - - - - - - - - - - - - - -
        // Momentum conservation within jets: Subjet level: ( 3 ):

        // Momentum conservation within subjets: ( 4 )

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

    return 0;
}
