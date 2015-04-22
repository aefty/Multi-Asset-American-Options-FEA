/**
 * kernel.h
 * =============
 * Primary heat Equation function
 *
 * HPC : Software Atelier 2015
 * Multi-Asset American Options Finite Element Implementation
 * Edited by : Aryan Eftekhari
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


// @sect3{The <code>BlackScholes</code> class}
//
// The next piece is the declaration of the main class of this program. It
// follows the well trodden path of previous examples. If you have looked at
// step-6, for example, the only thing worth noting here is that we need to
// build two matrices (the mass and Laplace matrix) and keep the current and
// previous time step's solution. We then also need to store the current
// time, the size of the time step, and the number of the current time
// step. The last of the member variables denotes the theta parameter
// discussed in the introduction that allows us to treat the explicit and
// implicit Euler methods as well as the Crank-Nicolson method and other
// generalizations all in one program.
//
// As far as member functions are concerned, the only possible surprise is
// that the <code>refine_mesh</code> function takes arguments for the
// minimal and maximal mesh refinement level. The purpose of this is
// discussed in the introduction.

#include "Payoff.h"
#include "Penelty.h"
#include "Boundary.h"

template<int dim>
class BlackScholes {
  public:
	/**
	 * Constructor
	 */
	BlackScholes(): fe(1), dof_handler(triangulation), time_step(TIME_STEP), theta(THETA) {};

	/**
	 * Run Function
	 */
	void run() {
		const unsigned int initial_global_refinement = 2;
		const unsigned int n_adaptive_pre_refinement_steps = 4;

		if (dim == 1) {
			GridGenerator::hyper_rectangle(triangulation,
			                               Point<dim>(S1_RANGE.at(0)),
			                               Point<dim>(S1_RANGE.at(1)));
		} else if ( dim == 2) {
			GridGenerator::hyper_rectangle(triangulation,
			                               Point<dim>(S1_RANGE.at(0), S2_RANGE.at(0)),
			                               Point<dim>(S1_RANGE.at(1), S2_RANGE.at(1)));
		} else {
			GridGenerator::hyper_rectangle(triangulation,
			                               Point<dim>(S1_RANGE.at(0), S2_RANGE.at(0), S3_RANGE.at(0)),
			                               Point<dim>(S1_RANGE.at(1), S2_RANGE.at(1), S3_RANGE.at(2)));
		}

		triangulation.refine_global (initial_global_refinement);

		setup_system();
		unsigned int pre_refinement_step = 0;
		Vector<double> tmp;
		Vector<double> forcing_terms;

		{
		start_time_iteration:

			tmp.reinit (solution.size());
			forcing_terms.reinit (solution.size());


			VectorTools::interpolate(dof_handler, ZeroFunction<dim>(), old_solution);
			solution = old_solution;

			timestep_number = 0;
			time            = 0;

			// Then we start the main loop until the computed time exceeds our
			// end time of 0.5. The first task is to build the right hand
			// side of the linear system we need to solve in each time step.
			// Recall that it contains the term $MU^{n-1}-(1-\theta)k_n AU^{n-1}$.
			// We put these terms into the variable system_rhs, with the
			// help of a temporary vector:
			while (time <= EXPIRE_TIME) {

				time += time_step;
				++timestep_number;

				std::cout << "Time step " << timestep_number << " at t=" << time << std::endl;

				// The second piece is to compute the contributions of the source
				// terms. This corresponds to the term $k_n
				// \left[ (1-\theta)F^{n-1} + \theta F^n \right]$. The following
				// code calls VectorTools::create_right_hand_side to compute the
				// vectors $F$, where we set the time of the right hand side
				// (source) function before we evaluate it. The result of this
				// all ends up in the forcing_terms variable:

				// Next, we add the forcing terms to the ones that
				// come from the time stepping, and also build the matrix
				// $M+k_n\theta A$ that we have to invert in each time step.
				// The final piece of these operations is to eliminate
				// hanging node constrained degrees of freedom from the
				// linear system:

				// 1) Make right-hand side
				mass_matrix.vmult(system_rhs, old_solution); // MU^(n-1)
				laplace_matrix.vmult(tmp, old_solution); // AU^(n-1)
				system_rhs.add(-(1 - theta) * time_step, tmp); // MU^(n-1) + AU^(n-1) * -(1 - theta) * time_step


				// 2) Make Penlety term if American
				if (STYLE_AMERICAN) {
					Penelty<dim> penelty_term;
					penelty_term.set_time(time);
					VectorTools::create_right_hand_side(dof_handler,
					                                    QGauss<dim>(fe.degree + 1),
					                                    penelty_term,
					                                    tmp);
					forcing_terms = tmp;
					forcing_terms *= time_step * theta;

					penelty_term.set_time(time - time_step);
					VectorTools::create_right_hand_side(dof_handler,
					                                    QGauss<dim>(fe.degree + 1),
					                                    penelty_term,
					                                    tmp);

					forcing_terms.add(time_step * (1 - theta), tmp);


					// 3) Combine source with right-hand site
					system_rhs += forcing_terms;
				}

				// 4) Make matrix to be inverted (left-hand side)
				system_matrix.copy_from(mass_matrix);
				system_matrix.add(theta * time_step, laplace_matrix);

				// 5) Eliminate hanging nodes from the right and left hand side
				constraints.condense (system_matrix, system_rhs);

				// There is one more operation we need to do before we
				// can solve it: boundary values. To this end, we create
				// a boundary value object, set the proper time to the one
				// of the current time step, and evaluate it as we have
				// done many times before. The result is used to also
				// set the correct boundary values in the linear system:
				{
					Boundary<dim, Payoff_Functions> boundary_values_function;
					boundary_values_function.set_time(time);

					std::map<types::global_dof_index, double> boundary_values;
					VectorTools::interpolate_boundary_values(dof_handler,
					        0,
					        boundary_values_function,
					        boundary_values);

					MatrixTools::apply_boundary_values(boundary_values,
					                                   system_matrix,
					                                   solution,
					                                   system_rhs);
				}

				// With this out of the way, all we have to do is solve the
				// system, generate graphical data, and...
				solve_time_step();


				// ...take care of mesh refinement. Here, what we want to do is
				// (i) refine the requested number of times at the very beginning
				// of the solution procedure, after which we jump to the top to
				// restart the time iteration, (ii) refine every fifth time
				// step after that.
				//
				// The time loop and, indeed, the main part of the program ends
				// with starting into the next time step by setting old_solution
				// to the solution we have just computed.
				if ((timestep_number == 1) && (pre_refinement_step < n_adaptive_pre_refinement_steps)) {
					refine_mesh (initial_global_refinement,
					             initial_global_refinement + n_adaptive_pre_refinement_steps);
					++pre_refinement_step;

					tmp.reinit (solution.size());
					forcing_terms.reinit (solution.size());

					std::cout << std::endl;

					goto start_time_iteration;
				} else if ((timestep_number > 0) && (timestep_number % MESH_REFINE_PERIOD == 0)) {
					refine_mesh (initial_global_refinement,
					             initial_global_refinement + n_adaptive_pre_refinement_steps);
					tmp.reinit (solution.size());
					forcing_terms.reinit (solution.size());
				}

				output_results();
				old_solution = solution;
			}
		}
	}

  private:
	/**
	 * Variable Definition
	 */

	// Discretization
	Triangulation<dim>   triangulation;
	FE_Q<dim>            fe;
	DoFHandler<dim>      dof_handler;

	ConstraintMatrix     constraints;

	// Matrix containers
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> mass_matrix;
	SparseMatrix<double> laplace_matrix;
	SparseMatrix<double> complication_matrix;
	SparseMatrix<double> system_matrix;

	Vector<double>       solution;
	Vector<double>       old_solution;
	Vector<double>       system_rhs;

	// Time & Other
	double               time;
	double               time_step;
	unsigned int         timestep_number;
	const double         theta;

	/**
	 * 1) setup_system : General system and matrix setup, also deal with hanging nodes
	 */
	void setup_system() {
		dof_handler.distribute_dofs(fe);

		std::cout << std::endl;
		std::cout << "===========================================";
		std::cout << std::endl;
		std::cout << "Number of active cells: " << triangulation.n_active_cells();
		std::cout << std::endl;
		std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs();
		std::cout << std::endl;
		std::cout << std::endl;

		constraints.clear ();
		DoFTools::make_hanging_node_constraints (dof_handler, constraints);
		constraints.close();

		CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern(dof_handler, c_sparsity, constraints, /*keep_constrained_dofs = */ true);
		sparsity_pattern.copy_from(c_sparsity);

		mass_matrix.reinit(sparsity_pattern);
		laplace_matrix.reinit(sparsity_pattern);
		system_matrix.reinit(sparsity_pattern);

		MatrixCreator::create_mass_matrix(dof_handler, QGauss<dim>(fe.degree + 1), mass_matrix, (const Function<dim> *)0, constraints);
		MatrixCreator::create_laplace_matrix(dof_handler, QGauss<dim>(fe.degree + 1), laplace_matrix, (const Function<dim> *)0, constraints);

		solution.reinit(dof_handler.n_dofs());
		old_solution.reinit(dof_handler.n_dofs());
		system_rhs.reinit(dof_handler.n_dofs());
	}

	/**
	 * 2) solve_time_step : Solve current step of PDE
	 */
	void solve_time_step() {
		SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
		SolverCG<> cg(solver_control);

		PreconditionSSOR<> preconditioner;
		preconditioner.initialize(system_matrix, 1.0);

		cg.solve(system_matrix, solution, system_rhs, preconditioner);

		constraints.distribute(solution);

		std::cout << "     " << solver_control.last_step();
		std::cout << " CG iterations." << std::endl;
	}

	/**
	 * [output_results description]
	 */
	void output_results() const {
		DataOut<dim> data_out;

		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(solution, "V");

		data_out.build_patches();

		const std::string filename = "output/solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk";

		std::cout << filename.c_str() ;
		std::ofstream output(filename.c_str());
		data_out.write_vtk(output);
	}

	/**
	 * [refine_mesh  description]
	 * @param min_grid_level [description]
	 * @param max_grid_level [description]
	 */
	void refine_mesh (const unsigned int min_grid_level, const unsigned int max_grid_level) {

		Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

		KellyErrorEstimator<dim>::estimate (dof_handler,
		                                    QGauss < dim - 1 > (fe.degree + 1),
		                                    typename FunctionMap<dim>::type(),
		                                    solution,
		                                    estimated_error_per_cell);

		GridRefinement::refine_and_coarsen_fixed_fraction (triangulation, estimated_error_per_cell, ERR_CONF_INT_REFINE, ERR_CONF_INT_COURSE);

		if (triangulation.n_levels() > max_grid_level)
			for (typename Triangulation<dim>::active_cell_iterator
			        cell = triangulation.begin_active(max_grid_level);
			        cell != triangulation.end(); ++cell)
			{ cell->clear_refine_flag (); }
		for (typename Triangulation<dim>::active_cell_iterator
		        cell = triangulation.begin_active(min_grid_level);
		        cell != triangulation.end_active(min_grid_level); ++cell)
		{ cell->clear_coarsen_flag (); }

		// As part of mesh refinement we need to transfer the solution vectors
		// from the old mesh to the new one. To this end we use the
		// SolutionTransfer class and we have to prepare the solution vectors that
		// should be transferred to the new grid (we will lose the old grid once
		// we have done the refinement so the transfer has to happen concurrently
		// with refinement). At the point where we call this function, we will
		// have just computed the solution, so we no longer need the old_solution
		// variable (it will be overwritten by the solution just after the mesh
		// may have been refined, i.e., at the end of the time step; see below).
		// In other words, we only need the one solution vector, and we copy it
		// to a temporary object where it is safe from being reset when we further
		// down below call <code>setup_system()</code>.
		//
		// Consequently, we initialize a SolutionTransfer object by attaching
		// it to the old DoF handler. We then prepare the triangulation and the
		// data vector for refinement (in this order).
		SolutionTransfer<dim> solution_trans(dof_handler);

		Vector<double> previous_solution;
		previous_solution = solution;
		triangulation.prepare_coarsening_and_refinement();
		solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

		// Now everything is ready, so do the refinement and recreate the dof
		// structure on the new grid, and initialize the matrix structures and the
		// new vectors in the <code>setup_system</code> function. Next, we actually
		// perform the interpolation of the solution from old to new grid.
		triangulation.execute_coarsening_and_refinement ();
		setup_system ();

		solution_trans.interpolate(previous_solution, solution);
	}

};



