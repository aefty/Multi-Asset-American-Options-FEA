/**
 * Parameters.h
 * =============
 * Global variables for simulation
 *
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Aryan Eftekhari & Edoardo Vecchi
 */

#include <iostream>
#include <vector>

/** Mesh Paramters
* ===============
* Refine those cells with the largest estimated error that together
* make up 60 (ERR_CONF_INT_REFINE) per cent of the error, and coarsens
* those cells with the smallest error that make up for a combined
* 40 (ERR_CONF_INT_COURSE) per cent of the error.
*/
const int MESH_REFINE_PERIOD = 10;
const double ERR_CONF_INT_REFINE = 0.6;
const double ERR_CONF_INT_COURSE = 0.4;

// Solver Paramters
const double THETA = 0.5;

// Contract Paramters
std::vector<double> ASSET_PRICE = {100};
std::vector<double> ASSET_SD = {10};
std::vector<double> STRIKE_PRICE = {10};
const bool STYLE_AMERICAN = false;

// Domain Paramters
const int DIM = 2;
std::vector<double> S1_RANGE = { -1, 1};
std::vector<double> S2_RANGE = { -1, 1 };
std::vector<double> S3_RANGE = {1, 10};

const double STRIKE_TIME = .1;
const double TIME_STEP = 1.0 / 200;