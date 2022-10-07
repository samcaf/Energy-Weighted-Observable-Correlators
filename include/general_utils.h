/**
 * @file    general_utils.h
 *
 * @brief   A utility header file for miscellaneous utilities in C++.
 * */


#ifndef GENERAL_UTILS
#define GENERAL_UTILS

// ---------------------------------
// Basic imports
// ---------------------------------
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>


// =====================================
// typedefs
// =====================================
typedef unsigned long size_t;

// =====================================
// Utility functions
// =====================================
std::string str_round(double x, int num_dec);
std::string periods_to_hyphens(std::string str);

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1);

#endif
