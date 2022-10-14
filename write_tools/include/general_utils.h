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
bool str_eq(const std::string str_a, const std::string str_b);

bool str_to_bool(const std::string boolstr);

std::string str_round(double x, const int num_dec);
std::string periods_to_hyphens(std::string str);
std::string remove_extension(std::string str);

template<typename T>
std::vector<T> arange(const T start, const T stop, const T step = 1);

#endif
