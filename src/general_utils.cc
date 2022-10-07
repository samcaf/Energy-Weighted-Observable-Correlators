/**
 * @file    general_utils.cc
 *
 * @brief   A set of miscellaneous utilities for C++.
 * */

// ---------------------------------
// Basic imports
// ---------------------------------
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "../include/general_utils.h"

// =====================================
// Utility functions
// =====================================

// ---------------------------------
// General utilities
// ---------------------------------

/**
* @brief:   Rounds a number to the specified number of
*           decimal places, returning a string.
*
* @param: x         A double we wish to round.
* @param: num_dec   The number of decimal places to
*                   which we would like to round.
*
* @return: std::string  The rounded string.
*/
std::string str_round(double x, int num_dec) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(num_dec);
    if (x < pow(10, -1*num_dec)) x = 0;
    ss << x;
    return ss.str();
}


/**
* @brief:   Replaces all periods in a string by hyphens.
*
* @param: str   The string to modify.
*
* @return: std::string      The modified string.
*/
std::string periods_to_hyphens(std::string str) {
    for (size_t i = 0; i < str.length(); ++i)
        if (str[i] == '.')
            str[i] = '-';
    return str;
}

/**
* @brief: C++ equivalent of numpy.arange function.
*         See https://stackoverflow.com/a/21217377
*
* @param: start, stop, step   Start, stop, and step size.
*
* @return: vector             std::vector associated with
*                             the arange.
*/
template<typename T>
std::vector<T> arange(T start, T stop, T step /*= 1*/) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}
