/**
 * @file    ewoc_cmdl.cc
 *
 * @brief   A command line utility file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */
#include <string>
#include <string.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

// Local imports:
#include "../include/general_utils.h"
#include "../include/ewoc_cmdln.h"


// =====================================
// Command Line Basics
// =====================================

// ---------------------------------
// Large text
// ---------------------------------
std::string flag_text = R"(
     ______           _____                         _____                _____
 ___|\     \         |\    \   _____           ____|\    \           ___|\    \
|     \     \        | |    | /    /|         /     /\    \         /    /\    \
|     ,_____/|       \/     / |    ||        /     /  \    \       |    |  |    |
|     \--'\_|/       /     /_  \   \/       |     |    |    |      |    |  |____|
|     /___/|        |     // \  \   \       |     |    |    |      |    |   ____
|     \____|\       |    |/   \ |    |      |\     \  /    /|      |    |  |    |
|____ '     /|      |\ ___/\   \|   /|      | \_____\/____/ |      |\ ___\/    /|
|    /_____/ |      | |   | \______/ |       \ |    ||    | /      | |   /____/ |
|____|     | /       \|___|/\ |    | |        \|____||____|/        \|___|    | /
  \( |_____|/           \(   \|____|/            \(    )/             \( |____|/
   '    )/               '      )/                '    '               '   )/
        '                       '                                          '
)";


std::string help_text = R"(
###################################
# EWOC Parameters:
###################################

  # ================================================
  # Basic Options (Required):
  # ================================================
  # Event Generation
  [<--n_events|-n> <int>] :                     Number of generated events;
  [<--level|-l> <str: parton|hadron>] :         QCD level of generated events (parton or hadron);
  [<--process|-p> <str>] :                      Process for generated events (quark, gluon, qcd, w);

  # Jet and Subjet information
  [<--jet_alg|-j> <int>] :                      Jet algorithm (0 [kt], 1 [ca], 2 [antikt]);
  [<--sub_alg|-s> <int>] :                      Subjet algorithm (0 [kt], 1 [ca], 2 [antikt]);
  [<--jet_rad|-R> <double>] :                   Jet radius (can be 'inf' or 'infty');
  [<--sub_rad|-r> <double>] :                   Subjet radius;
  # If these parameters are not given, the default is [Jet, subjet] : [(AKT, 1), (ca, .1)]

  # ================================================
  # Basic Options (Optional):
  # ================================================
  # Phase Space constraints
  [<--pt_min> <double>] :                       Minimum value of jet p_T;
  [<--pt_max> <double>] :                       Maximum value of jet p_T;

  # Misc. Options
  [<--verbose|-v> <int>] :                      Verbosity;


  # ================================================
  # Advanced Options (Optional):
  # ================================================
  # Event generator
  [<--event_generator|-G> <double>] :           Choice of event generator ("pythia" or "herwig");

  # Fragmentation Temperature
  [<--frag_temp|-T> <double>] :                 Temperature associated with string
                                                fragmentation;

  # Coming soon!
  # Rope related parameters

  # ================================================
  # Output Options (Optional):
  # ================================================
  [<--write_ewocs> <bool>] :                    Determines whether EWOC data is written;
  [<--write_pid_pt> <int>] :                    Determines whether pT spectrum for the particle
                                                with the given particle ID is written;

  # ================================================
  # Help:
  # ================================================
  [<--help|-h> ] :                                      Get this help message!


  # ================================================
  # Example commands:
  # ================================================

  # ---------------------------------
  # Single Subjet Radius
  # ---------------------------------
  # Generate EWOC data for 5000 parton-level qcd events
  # R=1 anti-kt jets, r=.1 kt subjets
  ```
    ./write_ewocs -n 5000 -l parton -p qcd -j 2 -s 0 -R 1.0 -r 0.1 --pt_min 50 --pt_max 3000

    cd python
    python3 plot_ewocs_sub_rads.py -n 5000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r 0.1 --pt_min 50 --pt_max 3000
    cd ..
  ```

  # ---------------------------------
  # Multiple Subjet Radii
  # ---------------------------------
  # Generate EWOC data for 50000 parton-level qcd events
  # R=1 anti-kt jets, r in [.0,.1,.2,.3,.4,.5] kt subjets
  ```
    for rsub in 0.0 0.1 0.2 0.3 0.4 0.5
      ./write_ewocs -n 50000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r $rsub --pt_min 50 --pt_max 3000
    end

    cd python
    python3 plot_ewocs_sub_rads.py -n 50000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r [0.0,0.1,0.2,0.3,0.4,0.5] --pt_min 50 --pt_max 3000
    cd ..
  ```

  // ---------------------------------
  // Plotting options
  // ---------------------------------
  # Notice that the short scripts above call the python function ```plot_ewocs_sub_rads.py``` in the python folder.
  # One can also plot ewocs while varying parameters other than subjet radii.
  # See python folder for event generation and plotting with more options.

  # There are also pre-generated scripts for plotting particular EWOCs and features:
  ```
  ls scripts/
  ls python/scripts/
  ```
)";


// ---------------------------------
// Command Line Parameter Verification
// ---------------------------------

/**
* @brief: checks to ensure that valid command line input was given.
*
* @param: argc  Number of command line inputs;
* @param: argv  Vector of command line inputs.
*
* @return: int  1 if failed, 0 if successful
*/
int checkCmdInputs(int argc, char* argv[]) {
    // ---------------------------------
    // Getting command line variables
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
    std::string event_gen     = eventgen_cmdln(argc, argv);
    double      frag_temp     = fragtemp_cmdln(argc, argv);

    // ---------------------------------
    // Help
    // ---------------------------------
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-h") or
          str_eq(argv[iarg], "--help")) {
            std::cout << flag_text << "\n\n";
            std::cout << help_text;
            return 1;
        }
    }

    // Verbose output
    if (verbose >= 0) {
        std::cout << flag_text << "\n\n";
        std::cout << "Function call:\n    ";
        while( argc-- ) std::cout << *(argv++) << " ";
        std::cout << "\n\n";
    }


    // ---------------------------------
    // Checks on input parameters
    // ---------------------------------
    // Event generator checks
    if (n_events <= 0){
        std::cout << "Number of events must be given as a positive integer.\n";
        return 1;
    }

    if (!str_eq(qcd_level, "parton") and !str_eq(qcd_level, "hadron")) {
        std::cout << "Invalid qcd level (given " << qcd_level
                  << "), must be 'parton' or 'hadron'.\n";
        return 1;
    }

    if (!str_eq(process_str, "quark") and !str_eq(process_str, "gluon")
      and !str_eq(process_str, "qcd") and !str_eq(process_str, "w")
      and !str_eq(process_str, "top")) {
        std::cout << "Invalid process (given " << process_str
                  << "), must be one of 'quark', 'gluon', 'qcd', 'w', 'top'.\n";
        return 1;
    }

    if (jet_alg_int != 0 and jet_alg_int != 1 and jet_alg_int != 2){
        std::cout << "Jet algorithm must be 0 (kt), 1 (ca), or 2 (antikt).\n";
        return 1;
    }

    if (jet_rad <= 0){
        std::cout << "Jet radius must be positive.\n";
        return 1;
    }

    if (sub_alg_int != 0 and sub_alg_int != 1 and sub_alg_int != 2){
        std::cout << "Subjet algorithm must be 0 (ca), 1 (kt), or 2 (antikt).\n";
        return 1;
    }

    if (sub_rad < 0){
        std::cout << "Subjet radius must be positive or zero.\n";
        return 1;
    }

    if (pt_min < 0){
        std::cout << "pt_min must be positive or zero.\n";
        return 1;
    }

    if (pt_max <= pt_min){
        std::cout << "pt_max must be larger than pt_min.\n";
        return 1;
    }

    if (E_cm != _ECM_DEFAULT and E_cm <= 0){
        std::cout << "Center of mass energy must be positive.\n";
        return 1;
    }

    if (!str_eq(s_channel, _SCHANNEL_DEFAULT)) {
        if (!str_eq(process_str, "qcd") and !str_eq(process_str, "top")) {
            std::cout << "No s channel production for desired process ("
               << process_str << ").";
            return 1;
        }
        else if (!str_eq(s_channel, "gm") and !str_eq(s_channel, "Z")) {
            std::cout << "Invalid s channel option: " << s_channel;
            return 1;
        }
    }

    // Advanced Options
    if (!str_eq(event_gen, _PYTHIA_STR) and !str_eq(event_gen, _HERWIG_STR)) {
        std::cout << "Event generator must be " + _PYTHIA_STR + " or " + _HERWIG_STR
                  << " (given: " + event_gen + ").\n";
        return 1;
    }

    if (frag_temp != _FRAG_TEMP_DEFAULT and
     (!str_eq(qcd_level, "hadron") or !str_eq(event_gen, _PYTHIA_STR))) {
        std::cout << "Thermal string fragmentation "
                  << "(given as " << std::to_string(frag_temp)
                  << ") requires hadron level event generation in Pythia.\n";
        return 1;
    }

    // Return 0 if all checks passed
    return 0;
}


// =====================================
// EWOC Setup Utilities
// =====================================

void setup_pythia_ewoc_cmdln(Pythia8::Pythia &pythia, int argc, char* argv[]) {
    // ---------------------------------
    // Convert command line variables
    // ---------------------------------
    // Basic Settings
    std::string qcd_level     = level_cmdln(argc, argv);
    std::string process_str   = process_cmdln(argc, argv);

    // Optional Settings
    double      pt_min        = ptmin_cmdln(argc, argv);

    double      E_cm          = Ecm_cmdln(argc, argv);
    std::string s_channel     = schannel_cmdln(argc, argv);

    int         verbose       = verbose_cmdln(argc, argv);

    // Advanced Settings
    double      frag_temp     = fragtemp_cmdln(argc, argv);


    // ---------------------------------
    // Pythia setup
    // ---------------------------------
    // Get default options
    pythia.readFile("write_tools/src/ewoc_setup.cmnd");

    // Give process information to Pythia
    if (str_eq(process_str, "quark") or str_eq(process_str, "gluon")) {
        // quark or gluon --> p p scattering
        pythia.readString("Beams:idA = 2212");
        pythia.readString("Beams:idB = 2212");

        // Process specific interactions
        if (str_eq(process_str, "quark")) {
            std::cout << "p p ->  q qbar, \n";
            pythia.readString("HardQCD:gg2qqbar = on");
        } else if (str_eq(process_str, "gluon")) {
            std::cout << "p p ->  g g, \n";
            pythia.readString("HardQCD:gg2gg = on");
        }
    }
    else {
        // Otherwise, e+ e- scattering
        pythia.readString("Beams:idA = 11");
        pythia.readString("Beams:idB = -11");

        // Process specific interactions
        if (str_eq(process_str, "qcd")) {
            std::cout << "\n# ---------------------------------"
                      << "\ne+ e- -> hadrons, ";
            // gm, Z, or gmZ in the s-channel
            pythia.readString("WeakSingleBoson:ffbar2ffbar(s:"
                    + s_channel + ") = on");
        } else if (str_eq(process_str, "w")) {
            std::cout << "\n# ---------------------------------"
                      << "\ne+ e- -> W W, ";
            pythia.readString("WeakDoubleBoson:ffbar2WW = on");
        } else if (str_eq(process_str, "top")) {
            std::cout << "\n# ---------------------------------"
                      << "\ne+ e- -> t tbar, ";
            // gm, Z, or gmZ in the s-channel
            pythia.readString("Top:ffbar2ttbar(s:"
                    + s_channel + ") = on");
        }
    }

    // Beam energy
    if (E_cm != _ECM_DEFAULT) {
        std::cout << "beam E_cm: "
                  << std::to_string(E_cm/1000.) << " TeV.\n"
                  << "# ---------------------------------\n";
        pythia.readString("Beams:eCM = " + std::to_string(E_cm));
    }
    else {
        std::cout << "beam E_cm: 4 TeV.\n"
                  << "# ---------------------------------\n";
        pythia.readString("Beams:eCM = 4000");
    }

    // Parton or hadron level
    if (str_eq(qcd_level, "parton"))
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
    }
    if (verbose < 1) {
        pythia.readString("Stat:showProcessLevel = off");
        pythia.readString("Stat:showErrors = off");
    }

    // Advanced Options
    if (frag_temp != _FRAG_TEMP_DEFAULT) {
        std::cout << "Turning on thermal model for string pT with temperature "
                  << std::to_string(frag_temp) << "\n";
        pythia.readString("StringPT:thermalModel = on");
        pythia.readString("StringPT:temperature = " + std::to_string(frag_temp));
    }

    pythia.init();
}



void write_ewoc_header_cmdln(int argc, char* argv[]) {
    // ---------------------------------
    // Getting command line variables
    // ---------------------------------
    // Basic Settings
    int         n_events      = nevents_cmdln(argc, argv);
    std::string qcd_level     = level_cmdln(argc, argv);
    std::string process_str   = process_cmdln(argc, argv);
    int         jet_alg_int   = jetalg_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);
    int         sub_alg_int   = subalg_cmdln(argc, argv);
    double      sub_rad       = subrad_cmdln(argc, argv);

    double      E_cm          = Ecm_cmdln(argc, argv);

    // Advanced Settings
    std::string event_gen     = eventgen_cmdln(argc, argv);

    // Open up the appropriate file
    std::ofstream file;
    std::string filename = ewoc_file_label(argc, argv);
    file.open(filename);

    // Add header with additional information
    file << "# ==================================\n"
         << "# EWOC information\n"
         << "# ==================================\n"
         << "# Function call ```";
    while( argc-- ) file << *(argv++) << " ";
    file << "```\n\n"
         << "# Event Generator:\n"
         << "event_gen = " + event_gen << "\n"
         << "# Process (E_cm = " << E_cm << "\n"
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

    file.close();

    return;
}

// =====================================
// Command Line Reading Utilities
// =====================================
// Format
/**
* @brief: Set of functions that return arguments from command
*         line input.
*
* @param: argc/argv         Command line input.
*
* @return: variable type    The desired argument.
*/

// ---------------------------------
// Required Basic Options
// ---------------------------------

// - - - - - - - - - - - - - - - - -
// Number of events
// - - - - - - - - - - - - - - - - -
int nevents_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-n") or
          str_eq(argv[iarg], "--n_events")) {
            return atoi(argv[iarg+1]);
        }
    }
    return -1;
}

// - - - - - - - - - - - - - - - - -
// QCD Level
// - - - - - - - - - - - - - - - - -
std::string level_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "-l") or
          str_eq(argv[iarg], "--level")) {
            return argv[iarg+1];
        }
    }
    return "no level";
}

// - - - - - - - - - - - - - - - - -
// Process
// - - - - - - - - - - - - - - - - -
std::string process_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-p") or
          str_eq(argv[iarg], "--process")) {
            return argv[iarg+1];
        }
    }
    return "no process";
}


// - - - - - - - - - - - - - - - - -
// Jet and subjet algorithms
// - - - - - - - - - - - - - - - - -
const int _JETALG_DEFAULT = 2;  // anti-kt
const int _SUBALG_DEFAULT = 1;  // cambridge-aachen (ca)

int jetalg_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-j") or
          str_eq(argv[iarg], "--jet_alg")) {
            return jet_alg_str_to_int(argv[iarg+1]);
        }
    }
    return _JETALG_DEFAULT;
}

int subalg_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-s") or
          str_eq(argv[iarg], "--sub_alg")) {
            return jet_alg_str_to_int(argv[iarg+1]);
        }
    }
    return _SUBALG_DEFAULT;
}


// - - - - - - - - - - - - - - - - -
// Jet and subjet radii
// - - - - - - - - - - - - - - - - -
const double _JETRAD_DEFAULT = 1.0;
const double _SUBRAD_DEFAULT = 0.1;

double jetrad_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-R") or
          str_eq(argv[iarg], "--jet_rad")) {
            // Allowing infinite jet radius (whole event)
            if(str_eq(argv[iarg+1], "inf") or
             str_eq(argv[iarg+1], "infty") )
                return 1000;
            // Otherwise, we assume it is finite
            return atof(argv[iarg+1]);
        }
    }
    return _JETRAD_DEFAULT;
}

double subrad_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "-r") or
          str_eq(argv[iarg], "--sub_rad")) {
            return atof(argv[iarg+1]);
        }
    }
    return _SUBRAD_DEFAULT;
}

// - - - - - - - - - - - - - - - - -
// Jet and subjet recombination schemes
// - - - - - - - - - - - - - - - - -
// Recombination schemes for jet finding
const fastjet::RecombinationScheme _JET_RECOMB_DEFAULT = fastjet::WTA_pt_scheme;
const fastjet::RecombinationScheme _SUB_RECOMB_DEFAULT = fastjet::WTA_pt_scheme;

fastjet::RecombinationScheme jetrecomb_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "--jet_recomb")) {
            // Only accepting WTA_pt for now
            if(str_eq(argv[iarg+1], "WTA_pt"))
                return fastjet::WTA_pt_scheme;
            else
                throw std::invalid_argument("Invalid jet recombination scheme");
        }
    }
    return _JET_RECOMB_DEFAULT;
}

fastjet::RecombinationScheme subrecomb_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(str_eq(argv[iarg], "--sub_recomb")) {
            // Only accepting WTA_pt for now
            if(str_eq(argv[iarg+1], "WTA_pt"))
                return fastjet::WTA_pt_scheme;
            else
                throw std::invalid_argument("Invalid subjet recombination scheme");
        }
    }
    return _SUB_RECOMB_DEFAULT;
}


// ---------------------------------
// Optional Basic Options
// ---------------------------------
// - - - - - - - - - - - - - - - - -
// Phase space options
// - - - - - - - - - - - - - - - - -
double ptmin_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "--pt_min"))
            return atof(argv[iarg+1]);
    }
    return 0;
}

double ptmax_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "--pt_max"))
            return atof(argv[iarg+1]);
    }
    return 100000;
}


const double _ECM_DEFAULT = -1;
double Ecm_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "--E_cm")
         or str_eq(argv[iarg], "-E"))
            return atof(argv[iarg+1]);
    }
    return _ECM_DEFAULT;
}

const std::string _SCHANNEL_DEFAULT = "gmZ";
std::string schannel_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "--s_channel"))
            return argv[iarg+1];
    }
    return _SCHANNEL_DEFAULT;
}


// - - - - - - - - - - - - - - - - -
// Misc. Options
// - - - - - - - - - - - - - - - - -
int verbose_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "--verbose")
         or str_eq(argv[iarg], "-v"))
            return atoi(argv[iarg+1]);
    }
    return 0;
}


// - - - - - - - - - - - - - - - - -
// Writing Options (for later plotting)
// - - - - - - - - - - - - - - - - -
bool writeewocs_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Determines whether EWOC output is written
        if(str_eq(argv[iarg], "--write_ewocs"))
            return str_to_bool(argv[iarg+1]);
    }
    return true;
}

int writepidpt_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Determines whether PID pT output is written
        if(str_eq(argv[iarg], "--write_pid_pt"))
            return atoi(argv[iarg+1]);
    }
    return 0;  // Default no PID pt
}


bool writejetpt_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Determines whether jet pT output is written
        if(str_eq(argv[iarg], "--write_jet_pt"))
            return str_to_bool(argv[iarg+1]);
    }
    return false;
}


std::string writeevent_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Determines whether event output is written
        if(str_eq(argv[iarg], "--write_event"))
            return argv[iarg+1];
    }
    return "";  // Default: do not write an event
}


// ---------------------------------
// Advanced Options
// ---------------------------------

// - - - - - - - - - - - - - - - - -
// Choosing event generator
// - - - - - - - - - - - - - - - - -
const std::string _PYTHIA_STR = "pythia";
const std::string _HERWIG_STR = "herwig";
std::string eventgen_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Choice of event generator
        if(str_eq(argv[iarg], "--event_generator")
         or str_eq(argv[iarg], "-G"))
            return argv[iarg+1];
    }
    return _PYTHIA_STR;  // Pythia is the default generator
}


// - - - - - - - - - - - - - - - - -
// Thermodynamical Fragmentation
// - - - - - - - - - - - - - - - - -
const double _FRAG_TEMP_DEFAULT = -10;
double fragtemp_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(str_eq(argv[iarg], "-T")
         or str_eq(argv[iarg], "--frag_temp"))
            return atof(argv[iarg+1]);
    }
    return _FRAG_TEMP_DEFAULT;
}

// - - - - - - - - - - - - - - - - -
// Rope Hadronization
// - - - - - - - - - - - - - - - - -
// Coming soon!
double _ropeshove_amp_default = -10;

// Others
// Coming soon!


// =====================================
// File Labelling Utilities
// =====================================
/**
* @brief: Returns a folder name associated with the given process.
*
* @param: argc/argv          Command line input.
*
* @return: string            Folder name.
*/
std::string process_folder(int argc, char* argv[]) {
    // Getting arguments from command line input
    std::string process_str   = process_cmdln(argc, argv);
    std::string qcd_level     = level_cmdln(argc, argv);
    std::string s_channel     = schannel_cmdln(argc, argv);

    // Making directory if it does not exist
    std::string proc_folder = "output/" + process_str + "_" + qcd_level;
    if (!str_eq(s_channel, _SCHANNEL_DEFAULT))
        proc_folder += "_schan_" + s_channel;
    proc_folder += "/";

    // https://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html
    // also, https://stackoverflow.com/a/21589609
    mkdir(proc_folder.c_str(),
          S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    return proc_folder;
}


std::string jet_alg_int_to_str(int alg_int) {
    assert(alg_int == 2 or alg_int == 1 or alg_int == 0);

    return (alg_int == 2) ? "akt" : (
        (alg_int == 1) ? "ca" : (
        (alg_int == 0) ? "kt" :
        "_ERR"));
}


int jet_alg_str_to_int(std::string alg_str) {
    assert(str_eq(alg_str, "akt") or str_eq(alg_str, "antikt")
           or str_eq(alg_str, "ca") or str_eq(alg_str, "kt")
           or str_eq(alg_str, "2") or str_eq(alg_str, "1")
           or str_eq(alg_str, "0"));

    return (str_eq(alg_str, "akt") or str_eq(alg_str, "2")
            or str_eq(alg_str, "antikt")) ? 2 : (
        (str_eq(alg_str, "ca") or str_eq(alg_str, "1")) ? 1 : (
        (str_eq(alg_str, "kt") or str_eq(alg_str, "0")) ? 0 :
        -1));
}


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
    int         jet_alg_int   = jetalg_cmdln(argc, argv);
    double      jet_rad       = jetrad_cmdln(argc, argv);

    std::string jet_rad_folder = "jetR" + str_round(jet_rad, 1) + "/";
    jet_rad_folder = periods_to_hyphens(jet_rad_folder);

    mkdir((process_folder(argc, argv) + jet_rad_folder).c_str(),
          S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    std::string jet_alg_folder = jet_alg_int_to_str(jet_alg_int) + "jet/";

    std::string folder_name = process_folder(argc, argv) + jet_rad_folder + jet_alg_folder;
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
    std::string sub_label = "_" + jet_alg_int_to_str(sub_alg_int) + "sub";

    std::string ewoc_file = "subR" + str_round(sub_rad, 2) + sub_label
                            + "_" + std::to_string(n_events) + "evts";

    // Optional basic arguments
    if (pt_min != -1)
        ewoc_file += "_ptmin"+str_round(pt_min, 1);
    if (pt_max != std::numeric_limits<double>::max())
        ewoc_file += "_ptmax"+str_round(pt_max, 1);
    if (E_cm != _ECM_DEFAULT)
        ewoc_file += "_Ecm"+str_round(E_cm, 1);

    // Optional advanced arguments
    if (frag_temp != _FRAG_TEMP_DEFAULT) {
        ewoc_file += "_temp"+str_round(frag_temp, 1);
    }

    // Final filename
    ewoc_file = periods_to_hyphens(ewoc_file);
    ewoc_file += ".txt";

    return folder_name + ewoc_file;
}
