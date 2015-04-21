/**
 * main.cc
 * =============
 * Entry function
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Edite by : Aryan Eftekhari & Edoardo Vecchi
 *
 * ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2013
 */



// The program starts with the usual include files, all of which you should
// have seen before by now:
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>
#include <algorithm>

// Then the usual placing of all content of this program into a namespace and
// the importation of the deal.II namespace into the one we will work in:
namespace BS_MAA {
using namespace dealii;


#include "lib/kernel.h"

}


// @sect3{The <code>main</code> function}
//
// Having made it this far,  there is, again, nothing
// much to discuss for the main function of this
// program: it looks like all such functions since step-6.
int main() {
	try {
		using namespace dealii;
		using namespace BS_MAA;

		deallog.depth_console(0);

		BlackScholes<DIM> heat_equation_solver;
		heat_equation_solver.run();

		//BlackScholes.output_results();

	} catch (std::exception &exc) {
		std::cerr << std::endl << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
		std::cerr << "Exception on processing: " << std::endl << exc.what()
		          << std::endl << "Aborting!" << std::endl
		          << "----------------------------------------------------"
		          << std::endl;

		return 1;
	} catch (...) {
		std::cerr << std::endl << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
		std::cerr << "Unknown exception!" << std::endl << "Aborting!"
		          << std::endl
		          << "----------------------------------------------------"
		          << std::endl;
		return 1;
	}

	return 0;
}
