/**
 * @file    ewoc_cmdl.h
 *
 * @brief   A command line utility header file for subjet energy-weighted observable correlators.
 */
#ifndef EWOC_CMDLN
#define EWOC_CMDLN

#include <string>
#include <string.h>
#include <iostream>

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"

// =====================================
// Command Line Basics
// =====================================
extern std::string flag_text;
extern std::string help_text;

int checkCmdInputs(int argc, char* argv[]);

// =====================================
// EWOC Setup Utilities
// =====================================
void setup_pythia_ewoc_cmdln(Pythia8::Pythia &pythia, int argc, char* argv[]);
void write_ewoc_header_cmdln(int argc, char* argv[]);

// =====================================
// Command Line Reading Utilities
// =====================================

// ---------------------------------
// Basic Options from Command Line
// ---------------------------------
// Event generation
int nevents_cmdln(int argc, char* argv[]);
std::string level_cmdln(int argc, char* argv[]);
std::string process_cmdln(int argc, char* argv[]);

// Jet and subjet info
extern const int _JETALG_DEFAULT;
extern const int _SUBALG_DEFAULT;
int jetalg_cmdln(int argc, char* argv[]);
int subalg_cmdln(int argc, char* argv[]);

extern const double _JETRAD_DEFAULT;
extern const double _SUBRAD_DEFAULT;
double jetrad_cmdln(int argc, char* argv[]);
double subrad_cmdln(int argc, char* argv[]);

// Recombination scheme info
extern const fastjet::RecombinationScheme _JET_RECOMB_DEFAULT;
extern const fastjet::RecombinationScheme _SUB_RECOMB_DEFAULT;
fastjet::RecombinationScheme jetrecomb_cmdln(int argc, char* argv[]);
fastjet::RecombinationScheme subrecomb_cmdln(int argc, char* argv[]);

// ---------------------------------
// Optional Basic Options
// ---------------------------------

// Phase space options
double ptmin_cmdln(int argc, char* argv[]);
double ptmax_cmdln(int argc, char* argv[]);

extern const double _ECM_DEFAULT;
double Ecm_cmdln(int argc, char* argv[]);

extern const std::string _SCHANNEL_DEFAULT;
std::string schannel_cmdln(int argc, char* argv[]);

// Misc. Options
int verbose_cmdln(int argc, char* argv[]);

// Writing Options (for later plotting)
bool writeewocs_cmdln(int argc, char* argv[]);
int writepidpt_cmdln(int argc, char* argv[]);
std::string writeevent_cmdln(int argc, char* argv[]);

// ---------------------------------
// Advanced Options
// ---------------------------------
extern const std::string _PYTHIA_STR;
extern const std::string _HERWIG_STR;
std::string eventgen_cmdln(int argc, char* argv[]);

extern const double _FRAG_TEMP_DEFAULT;
double fragtemp_cmdln(int argc, char* argv[]);


// =====================================
// File Labelling Utilities
// =====================================

std::string process_folder(int argc, char* argv[]);

std::string jet_alg_int_to_str(int alg_int);
int jet_alg_str_to_int(std::string alg_str);

std::string ewoc_folder(int argc, char* argv[]);
std::string ewoc_file_label(int argc, char* argv[]);


#endif
