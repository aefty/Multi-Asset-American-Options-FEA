/**
 * Parameters.h
 * =============
 * Global variables for simulation
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 */

#include <iostream>
#include <vector>


/** Mesh Parameters
* ===============
* Refine those cells with the largest estimated error that together
* make up 60 (ERR_CONF_INT_REFINE) per cent of the error, and coarsens
* those cells with the smallest error that make up for a combined
* 40 (ERR_CONF_INT_COURSE) per cent of the error.
*/

const int MESH_REFINE_PERIOD = 1;
const double ERR_CONF_INT_REFINE = 0.8;
const double ERR_CONF_INT_COURSE = 0.2;

// Solver Parameters
const double THETA = 0.5;

// Contract Parameters
const double R = 0.1;
std::vector<double> X = {1.0};
std::vector<double> X_SD = {0.05};

const bool STYLE_AMERICAN = false;

// Domain Parameters
const int DIM = 2;
std::vector<double> X1_RANGE = { -10.0, 10.0};
std::vector<double> X2_RANGE = { -10.0, 10.0};
std::vector<double> X3_RANGE = { -10.0, 10.0};

const double DT = 1.0 / 100;
const double T = 1.0 * 25.0 * DT;

const long double EPS = 0.0000001; // 10^-6
