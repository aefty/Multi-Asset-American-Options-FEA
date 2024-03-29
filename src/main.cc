/* ---------------------------------------------------------------------
 * Edited By: Aryan Eftekhari & Edoardo Vecchi 2015 (Based on modification done by Patrick Sanan, May 2015 )
 * Università della Svizzera italiana
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
#include <cmath>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <sstream>

using namespace dealii;
#include "lib/IntrinsicValue.cc"
#include "lib/BlackScholes.cc"


int main(int argc, char ** argv) {

  try {
    using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    deallog.depth_console(0);

    int const dim = 3;
    int AMERICAN = false;
    double R;

    Tensor<2, dim>  D;
    Tensor<1, dim>  C;

    if (dim == 1) {
      R = .1;

      D[0][0] = 0.3 * 0.5;

      C[0] = R - D[0][0];

    } else if (dim == 2) {
      R = .1;

      D[0][0] = .3   * 0.5;
      D[0][1] = .0   * 0.5;

      D[1][0] = .0   * 0.5;
      D[1][1] = .3   * 0.5;

      C[0] = R - D[0][0];
      C[1] = R - D[1][1];

    } else if (dim == 3) {
      R = .1;

      D[0][0] = .3   * 0.5;
      D[0][1] = .0   * 0.5;
      D[0][2] = .0   * 0.5;

      D[1][0] = .0   * 0.5;
      D[1][1] = .3   * 0.5;
      D[1][2] = .0   * 0.5;

      D[2][0] = .0   * 0.5;
      D[2][1] = .0   * 0.5;
      D[2][2] = .3   * 0.5;

      C[0] = R - D[0][0];
      C[1] = R - D[1][1];
      C[2] = R - D[2][2];

    } else {
      std::cout << "Problem Dim Exceeded. \n";
      return 0;
    }

    // Descretization Parameters
    double T_DISC = 1.0 / 12.0;
    double T_RANGE = 2.0 / 12.0;

    double X_DISC = 5;
    std::vector <double> X_RANGE = { -10.0, 2.0};
    double THETA = 0.5;

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {

      std::cout << "BlackScholes Solver " << std::endl;
      std::cout << "====================" << std::endl;
      std::cout << "American Option:    " << AMERICAN << std::endl;
      std::cout << "Dimension:          " << dim << std::endl;
      std::cout << "Expiration:         " << T_RANGE << std::endl;
      std::cout << "Time dt:            " << T_DISC << std::endl;
      std::cout << "Time 1/dt:          " << 1.0 / T_DISC << std::endl;
      std::cout << "Time Steps:         " << T_RANGE / T_DISC << std::endl;
      std::cout << "Space Disc:         " << X_DISC << std::endl;
      std::cout << "Theta:              " << THETA << std::endl;
      std::cout << "Risk Free Rate:     " << R << std::endl;
      std::cout << "====================" << std::endl;
    }


    BlackScholes<dim> backscholes_solver(D, C, R, X_DISC, X_RANGE, T_DISC, T_RANGE, THETA , AMERICAN);
    backscholes_solver.run();

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