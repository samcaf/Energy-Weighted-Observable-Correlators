/**
 * @file    ewoc_cmdl.cc
 *
 * @brief   A command line utility file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */
#ifndef EWOC_CMDLN
#define EWOC_CMDLN

#include <string>
#include <iostream>

// =====================================
// Command Line Basics
// =====================================

int checkCmdInputs(int argc, char* argv[]);

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
int jetalg_cmdln(int argc, char* argv[]);
int subalg_cmdln(int argc, char* argv[]);

double jetrad_cmdln(int argc, char* argv[]);
double subrad_cmdln(int argc, char* argv[]);

// ---------------------------------
// Optional Basic Options
// ---------------------------------
// Phase space options
double ptmin_cmdln(int argc, char* argv[]);
double ptmax_cmdln(int argc, char* argv[]);

// Misc. Options
int verbose_cmdln(int argc, char* argv[]);

// ---------------------------------
// Advanced Options
// ---------------------------------
extern double _frag_temp_default;
double fragtemp_cmdln(int argc, char* argv[]);

#endif
