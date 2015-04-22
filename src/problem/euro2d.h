/**
 * Parameters.h
 * =============
 * Global variables for simulation
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Aryan Eftekhari
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
const double RFR = 0.3;
std::vector<double> ASSET_PRICE = {100.0};
std::vector<double> ASSET_SD = {10.0};
double STRIKE_PRICE = .2;

const bool STYLE_AMERICAN = false;

// Domain Parameters
const int DIM = 2;
std::vector<double> S1_RANGE = {0.0, 1.0};
std::vector<double> S2_RANGE = {0.0, 1.0};
std::vector<double> S3_RANGE = {1.0, 10.0};

const double TIME_STEP = 1.0 / 300;
const double EXPIRE_TIME = 2.0 * 25.0 * TIME_STEP;