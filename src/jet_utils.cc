/**
 * @file    ewoc_utils.cc
 *
 * @brief   A utility file for subjet energy-weighted observable correlators.
 *          Version 0.2
 */

// ---------------------------------
// Basic imports
// ---------------------------------
#include <iostream>
#include <string.h>
#include <cmath>
#include <limits>
#include <locale>
#include <string>
#include <vector>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>

// ---------------------------------
// HEP imports
// ---------------------------------
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "../include/jet_utils.h"

using namespace fastjet;


// =====================================
// Utility functions
// =====================================

// ---------------------------------
// Error utilities
// ---------------------------------
/**
* @brief: Raises an error if an invalid pairwise observable name is specified.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
std::string jetalg_error(JetAlgorithm alg) {
    throw std::invalid_argument( "Invalid jet algorithm " + std::to_string(alg) );
    return "ERR";
}


/**
* @brief: Raises an error if an invalid pairwise observable name is specified.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
double pairwise_error() {
    throw std::invalid_argument( "Invalid pairwise observable." );
    return 0;
}


// ---------------------------------
// Pseudojet utilities
// ---------------------------------

/**
* @brief: Function which tests whether the given int is a neutrino
*         particle id.
*
* @param: int id    Given particle id.
*
* @return: bool     True if given int is a neutrino particle id,
*                   false otherwise.
*/
int is_nu_id(int id) {
    return (id==12 || id==14 || id==16 || id==-12 || id==-14 || id==-16);
}


/**
* @brief: Takes in a Pythia Event and return a vector of PseudoJets
*         containing all particles of a current event.
*
* @param: event     Pythia Event type, as in pythia.event().
*                   The event under consideration.
*
* @return: vector<PseudoJet>    A vector containing all particles in the event.
*/
PseudoJets get_particles_pythia(Pythia8::Event event){
    // Storing particles of an event as PseudoJets
    PseudoJets particles;

    for (int ipart = 0; ipart < event.size(); ipart++) {
        // Finding visible final states and converting to pseudojets
        if (event[ipart].isFinal() and !is_nu_id(event[ipart].id())){
        PseudoJet particle(event[ipart].px(), event[ipart].py(),
                           event[ipart].pz(), event[ipart].e());
        particle.set_user_index(ipart);
        particles.push_back(particle);
    }}

    return particles;
}


/**
* @brief: Compares two pseudojets, returning true if they have
*         equal 3-momenta.
*
* @param: pj1, pj2   Pseudojets to compare
*
* @return: bool
*/
bool equal_pjs(PseudoJet pj1, PseudoJet pj2) {
    return (pj1.E() == pj2.E() and
            pj1.px() == pj2.px() and
            pj1.py() == pj2.py() and
            pj1.pz() == pj2.pz());
}


/**
* @brief: Returns the sum of scalar pT in a vector of pseudojets.
*
* @param: pjs   A vector containing the pseudojets to consider.
*
* @return: double The total energy.
*/
double SumScalarPt(PseudoJets pjs) {
    // Finding total jet energy
    double pT = 0;
    for (auto pj : pjs) pT += pj.perp();
    return pT;
}


/**
* @brief: Returns the total energy in a vector of pseudojets.
*
* @param: pjs   A vector containing the pseudojets to consider.
*
* @return: double The total energy.
*/
double SumEnergy(PseudoJets pjs) {
    // Finding total jet energy
    double Q = 0;
    for (auto pj : pjs) Q += pj.E();
    return Q;
}


/**
* @brief: Returns mass of a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double Cosine of the angle between the pseudojets.
*/
double pair_mass(PseudoJet pj1, PseudoJet pj2) {
    return sqrt(pj1.m2() + pj2.m2() +
                pj1.E()*pj2.E() - pj1.px()*pj2.px() -
                pj1.py()*pj2.py() - pj1.pz()*pj2.pz());
}


/**
* @brief: Returns the cosine of the angle between a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double   Cosine of the angle between the pseudojets.
*/
double pair_cos(PseudoJet pj1, PseudoJet pj2) {
    double cos = (pj1.px()*pj2.px() + pj1.py()*pj2.py() + pj1.pz()*pj2.pz())/(sqrt(pj1.modp2())*sqrt(pj2.modp2()));
    if (cos > 1.0) cos = 1.0; // Need to handle case of numerical overflow
    return cos;
}


/**
* @brief: Returns the angle between a pair of pseudojets.
*
* @param: pj1, pj2  The given pseudojets.
*
* @return: double   Angle between the pseudojets (in radians).
*/
double pair_theta(PseudoJet pj1, PseudoJet pj2) {
    return acos(pair_cos(pj1, pj2));
}
