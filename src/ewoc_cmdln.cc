/**
 * @file    ewoc_cmdl.cc
 *
 * @brief   A command line utility file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */

#include <string>
#include <iostream>

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
// =====================================
// Old Syntax
// =====================================
The old syntax is

    ```./subjet_ewocs [number of events] [qcd_level] [process] [jet algorithm] [jet radius] [subjet algorithm][subjet radius][pt_min] [pt_max]```.

This is supported for now.


Example command:

    ```./pythia_ewocs 5000 parton qcd 1 1.0 0 0.10 50 3000```.

This produces 5000 parton-level QCD events for e+e- => hadrons, with jets (50 GeV < pT < 3 TeV) clustered with kT and radius 1.0, and subjets clustered with C-A and radius 0.1.


To plot, run pythia_ewocs_plot.py with the same arguments except pt_min/max For example:

    ```python3 plot_pythia_ewocs.py 5000 parton qcd 1 1.0 0 0.10```

plots the output of the binary (```./pythia_ewocs```) above.
The python plotting can also take in a list of subjet radii:
    ```python3 plot_pythia_ewocs.py 5000 parton qcd 1 1.0 0 [0.10,0.20,0.30]```



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

  # ================================================
  # Basic Options (Optional):
  # ================================================
  # Phase Space constraints
  [<--pt_min> <double>] :                       Minimum value of jet p_T;
  [<--pt_max> <double>] :                       Maximum value of jet p_T;

  # Misc. Options
  [<--verbose> <int>] :                         Verbosity (coming soon);


  # ================================================
  # Advanced Options (Optional):
  # ================================================
  # Coming soon!
  # Temperature
  # Rope related parameters


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
  #R=1 anti-kt jets, r=.1 kt subjets
  ```
    ./pythia_ewocs -n 5000 -l parton -p qcd -j 2 -s 0 -R 1.0 -r 0.1 --pt_min 50 --pt_max 3000

    cd python
    python3 plot_pythia_ewocs.py -n 5000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r 0.1 --pt_min 50 --pt_max 300
    cd ..
  ```

  # ---------------------------------
  # Multiple Subjet Radii
  # ---------------------------------
  # Generate EWOC data for 50000 parton-level qcd events
  # R=1 anti-kt jets, r in [.0,.1,.2,.3,.4,.5] kt subjets
  ```
    for rsub in 0.0 0.1 0.2 0.3 0.4 0.5
      ./pythia_ewocs -n 50000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r $rsub --pt_min 50 --pt_max 3000
    end

    cd python
    python3 plot_pythia_ewocs.py -n 50000 -l parton -p qcd -j 2 -s 1 -R 1.0 -r [0.0,0.1,0.2,0.3,0.4,0.5] --pt_min 50 --pt_max 3000
    cd ..
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
    std::cout << flag_text << "\n\n";

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

    // Advanced Settings
    double      frag_temp     = fragtemp_cmdln(argc, argv);

    // ---------------------------------
    // Help
    // ---------------------------------
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-h") == 0 or
          strcmp(argv[iarg], "--help") == 0) {
            std::cout << help_text;
            return 1;
        }
    }

    // ---------------------------------
    // Default/backwards compatible arguments
    // ---------------------------------
    if (argc == 10 and n_events == -1){
        //Convert command line parameters to variables
        n_events    = atoi(argv[1]);
        qcd_level   = argv[2];
        process_str = argv[3];
        jet_alg_int = atoi(argv[4]);
        jet_rad     = atof(argv[5]);
        sub_alg_int = atoi(argv[6]);
        sub_rad     = atof(argv[7]);
        pt_min      = atof(argv[8]);
        pt_max      = atof(argv[9]);
    }

    // If not backwards compatible, no help message, and there aren't enough arguments
    else if (argc < 19) {
        std::cout << "Missing required command line arguments "
        << "(received " << argc-1 << ").\n"
        << "Use ```--help``` for more detail.\n";
        return 1;
    }

    // ---------------------------------
    // Checks on input parameters
    // ---------------------------------
    if (n_events <= 0){
        std::cout << "Number of events must be a positive integer.\n";
        return 1;
    }

    if (qcd_level != "parton" and qcd_level != "hadron") {
        std::cout << "Invalid qcd level " << qcd_level << ", must be 'parton' or 'hadron'.\n";
        return 1;
    }

    if (process_str != "quark" and process_str != "gluon"
      and process_str != "qcd" and process_str != "w") {
        std::cout << "Invalid process " << process_str << ", must be one of (quark, gluon, qcd, w).\n";
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

    if (pt_min <= 0 and pt_min != -1){
        std::cout << "pt_min must be positive.\n";
        return 1;
    }

    if (pt_max <= pt_min){
        std::cout << "pt_max must be larger than pt_min.\n";
        return 1;
    }

    // Advanced Options
    if (frag_temp != _frag_temp_default and qcd_level != "hadron") {
        std::cout << "Thermal string fragmentation requires hadron level event generation.\n";
        return 1;
    }

    // Return 0 if all checks passed
    return 0;
}


/**
* @brief: Set of functions that return arguments from command
*         line input.
*
* @param: argc/argv         Command line input.
*
* @return: variable type    The desired argument.
*/
// =====================================
// Command Line Reading Utilities
// =====================================

// ---------------------------------
// Required Basic Options
// ---------------------------------

// Number of events
int nevents_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-n") == 0 or
          strcmp(argv[iarg], "--n_events") == 0) {
            return atoi(argv[iarg+1]);
        }
    }
    return -1;
}

// QCD Level
std::string level_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(strcmp(argv[iarg], "-l") == 0 or
          strcmp(argv[iarg], "--level") == 0) {
            return argv[iarg+1];
        }
    }
    return "No level";
}

// Process
std::string process_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-p") == 0 or
          strcmp(argv[iarg], "--process") == 0) {
            return argv[iarg+1];
        }
    }
    return "No process";
}


// Jet and subjet algorithm
int jetalg_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-j") == 0 or
          strcmp(argv[iarg], "--jet_alg") == 0) {
            return atoi(argv[iarg+1]);
        }
    }
    return -1;
}

int subalg_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-s") == 0 or
          strcmp(argv[iarg], "--sub_alg") == 0) {
            return atoi(argv[iarg+1]);
        }
    }
    return -1;
}


// Jet and subjet radii
double jetrad_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-R") == 0 or
          strcmp(argv[iarg], "--jet_rad") == 0) {
            // Allowing infinite jet radius (whole event)
            if(strcmp(argv[iarg+1], "inf") == 0 or
             strcmp(argv[iarg+1], "infty") == 0 )
                return 1000;
            // Otherwise, we assume it is finite
            return atof(argv[iarg+1]);
        }
    }
    return -1.0;
}

double subrad_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        if(strcmp(argv[iarg], "-r") == 0 or
          strcmp(argv[iarg], "--sub_rad") == 0) {
            return atof(argv[iarg+1]);
        }
    }
    return -1;
}

// ---------------------------------
// Optional Basic Options
// ---------------------------------
// Phase space options
double ptmin_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(strcmp(argv[iarg], "--pt_min") == 0)
            return atof(argv[iarg+1]);
    }
    return -1;
}

double ptmax_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(strcmp(argv[iarg], "--pt_max") == 0)
            return atof(argv[iarg+1]);
    }
    return std::numeric_limits<double>::max();
}

// Misc. Options
int verbose_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(strcmp(argv[iarg], "--verbose") == 0)
            return atoi(argv[iarg+1]);
    }
    return 0;
}


// ---------------------------------
// Advanced Options
// ---------------------------------
// Thermodynamical Fragmentation
double _frag_temp_default = -10;
double fragtemp_cmdln(int argc, char* argv[]) {
    for(int iarg=0;iarg<argc;iarg++) {
        // Event Generator Parameters
        if(strcmp(argv[iarg], "-T") == 0
         or strcmp(argv[iarg], "--frag_temp"))
            return atof(argv[iarg+1]);
    }
    return _frag_temp_default;
}

// Rope Hadronization
// Coming soon!
double _ropeshove_amp_default = -10;

// Others
// Coming soon!
