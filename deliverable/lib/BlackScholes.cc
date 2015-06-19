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



template<int dim>
class BlackScholes {

   Triangulation<dim>   triangulation;
   FE_Q<dim>            fe;
   DoFHandler<dim>      dof_handler;

   ConstraintMatrix     hanging_node_constraints;

   double               time;
   const double         time_step;
   unsigned int         timestep_number;

   double               t_max;
   const double         theta;
   std::vector <double> range;
   const double         discretization;

   Tensor<2, dim>  diffusion_const;
   Tensor<1, dim>  convection_const;
   double          reaction_const;
   int american;

   PETScWrappers::MPI::SparseMatrix rhs_matrix;
   PETScWrappers::MPI::SparseMatrix system_matrix;

   PETScWrappers::MPI::Vector solution;
   PETScWrappers::MPI::Vector old_solution;
   PETScWrappers::MPI::Vector system_rhs;

   PETScWrappers::MPI::Vector intrinsicValue;

   MPI_Comm mpi_communicator;
   const unsigned int n_mpi_processes;
   const unsigned int this_mpi_process;

   //PDS: a rank-0 output stream
   ConditionalOStream pcout;

 public:

   BlackScholes(Tensor<2, dim> DIFF_CONST, Tensor<1, dim> CONV_CONST , double R, double X_DISC , std::vector<double> X_RANGE, double T_DISC , double T_RANGE , double THETA , int AMERICAN )
      :
      fe(1),
      dof_handler(triangulation),
      time_step(T_DISC),
      t_max(T_RANGE),
      theta(THETA),
      range(X_RANGE),
      discretization(X_DISC),
      diffusion_const(DIFF_CONST),
      convection_const(CONV_CONST),
      reaction_const(R),
      american(AMERICAN),
      mpi_communicator (MPI_COMM_WORLD),
      n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
      this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
      pcout (std::cout) {
      pcout.set_condition(this_mpi_process == 0);
   };

   ~BlackScholes() {
      dof_handler.clear ();
   };

   void run() {

      double check;
      double t_wall = 0.0;
      if (this_mpi_process == 0) {
         t_wall = MPI_Wtime();
      };


      // Set Grid
      {
         if (dim == 1) {
            GridGenerator::hyper_rectangle(triangulation,
                                           Point<dim>(range[0]),
                                           Point<dim>(range[1]));
         } else if (dim == 2) {
            GridGenerator::hyper_rectangle(triangulation,
                                           Point<dim>(range[0], range[0]),
                                           Point<dim>(range[1], range[1]));
         } else if (dim == 3) {
            GridGenerator::hyper_rectangle(triangulation,
                                           Point<dim>(range[0], range[0], range[0]),
                                           Point<dim>(range[1], range[1], range[1]));
         }
      }

      triangulation.refine_global (discretization);

      setup_system();
      assemble_system();

      PETScWrappers::MPI::Vector penTerm;
      PETScWrappers::MPI::Vector v;
      PETScWrappers::MPI::Vector v_last;

      const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);
      intrinsicValue.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

      // Set Initial Conditions
      VectorTools::interpolate(dof_handler, IntrinsicValue<dim>(), intrinsicValue);
      hanging_node_constraints.distribute(intrinsicValue);

      solution        = intrinsicValue;
      old_solution    = intrinsicValue;
      timestep_number = 0;
      time            = 0;


      output_results();

      while (time <= (t_max - 2.0 * time_step)) {
         time += time_step;
         timestep_number = timestep_number + 1;
         rhs_matrix.vmult(system_rhs, old_solution);
         int niter = solve(solution, system_rhs, system_matrix);

         if (american) {

            penTerm.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
            v.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
            v_last.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

            // We dont check norm
            for (int k = 0; k < 15; ++k) {

               v = solution;// old_solution is the old AMER solution, not EURO!

               // Assemble Local Penalty
               for (int i =  penTerm.local_range().first; i < penTerm.local_range().second; ++i) {
                  penTerm[i] = std::fmax(intrinsicValue[i] - v[i], 0);
               };

               v.add(penTerm);
               solution = v;
            };
         };

         output_results();
         old_solution = solution;
      }

      if (this_mpi_process == 0) {
         pcout << "\n+Run (" << n_mpi_processes << " #Procs) :"  << MPI_Wtime() - t_wall  << std::endl;
      };
   };

 private:
   void setup_system() {

      double t_wall = 0.0;
      if (this_mpi_process == 0) {
         t_wall = MPI_Wtime();
      };

      GridTools::partition_triangulation (n_mpi_processes, triangulation);
      dof_handler.distribute_dofs (fe);

      DoFRenumbering::subdomain_wise (dof_handler);
      const types::global_dof_index n_local_dofs =
         DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);

      rhs_matrix.reinit (mpi_communicator,
                         dof_handler.n_dofs(),
                         dof_handler.n_dofs(),
                         n_local_dofs,
                         n_local_dofs,
                         dof_handler.max_couplings_between_dofs());

      system_matrix.reinit (mpi_communicator,
                            dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            n_local_dofs,
                            n_local_dofs,
                            dof_handler.max_couplings_between_dofs());

      solution.reinit       (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      old_solution.reinit   (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      system_rhs.reinit     (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

      hanging_node_constraints.clear ();
      DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
      hanging_node_constraints.close();

      if (this_mpi_process == 0) {
         pcout << "+Setup System : " << MPI_Wtime() - t_wall  << std::endl;
      };

   };

   void assemble_system() {

      double t_wall = 0.0;
      if (this_mpi_process == 0) {
         t_wall = MPI_Wtime();
      };

      QGauss<dim>  quadrature_formula(2);
      FEValues<dim> fe_values (fe, quadrature_formula,
                               update_values   | update_gradients |
                               update_quadrature_points | update_JxW_values);
      const unsigned int   dofs_per_cell = fe.dofs_per_cell;
      const unsigned int   n_q_points    = quadrature_formula.size();

      FullMatrix<double>   cell_matrix_laplace     (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_mass        (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_convection  (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_omega       (dofs_per_cell, dofs_per_cell);

      FullMatrix<double>   cell_matrix_system     (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_rhs        (dofs_per_cell, dofs_per_cell);

      Vector<double>       cell_rhs            (dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell != endc; ++cell) {
         cell_matrix_laplace = 0;
         cell_matrix_mass    = 0;
         cell_matrix_convection    = 0;
         cell_rhs = 0;
         fe_values.reinit (cell);

         if (cell->subdomain_id() == this_mpi_process) {

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
               for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                     cell_matrix_laplace(i, j) += (
                                                     fe_values.shape_grad (i, q_point) *
                                                     diffusion_const *
                                                     fe_values.shape_grad (j, q_point) *
                                                     fe_values.JxW (q_point)
                                                  );

                     cell_matrix_mass(i, j) += (
                                                  fe_values.shape_value(i, q_point) *
                                                  fe_values.shape_value(j, q_point) *
                                                  fe_values.JxW (q_point)
                                               );

                     cell_matrix_convection(i, j) += (
                                                        fe_values.shape_grad (j, q_point) *
                                                        convection_const *
                                                        fe_values.shape_value (i, q_point) *
                                                        fe_values.JxW (q_point)
                                                     );
                  };
               };
            };

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
               for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                  cell_matrix_omega(i, j) =
                     -1.0 * cell_matrix_laplace(i, j) +
                     cell_matrix_convection(i, j) -
                     reaction_const * cell_matrix_mass(i, j);

                  cell_matrix_system(i, j) =
                     cell_matrix_mass(i, j) -
                     time_step * (theta)  * cell_matrix_omega(i, j);

                  cell_matrix_rhs(i, j) =
                     cell_matrix_mass(i, j) +
                     time_step * (1.0 - theta)  *  cell_matrix_omega(i, j);
               };
            };

            cell->get_dof_indices (local_dof_indices);

            // System Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_system, cell_rhs,
               local_dof_indices, system_matrix
               , system_rhs);

            // RHS Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_rhs, cell_rhs,
               local_dof_indices, rhs_matrix
               , system_rhs);
         };
      };

      system_matrix.compress(VectorOperation::add);
      rhs_matrix.compress(VectorOperation::add);

      if (this_mpi_process == 0) {
         pcout << "+Assemble System : " << MPI_Wtime() - t_wall  << std::endl;
      };
   };

   int solve(PETScWrappers::MPI::Vector & solution, PETScWrappers::MPI::Vector & system_rhs, PETScWrappers::MPI::SparseMatrix & system_matrix) {

      pcout << "." << std::flush;

      SolverControl solver_control (solution.size(), 1e-10 * system_rhs.l2_norm());

      PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

      // PDS: We use a very simple parallel preconditioner, block jacobi, which
      // PDS: uses an ILU(0) preconditioner applied to each local block.
      // PDS: The original example used SOR, which is not trivial to parallelize.
      // PDS: Note that this is means the preconditioner *weakens* with core count
      // PDS: (going from 1 to 2 mpi processes will mean the iteration count will change from 4/5 to 8/9)
      PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

      cg.solve (system_matrix, solution, system_rhs, preconditioner);

      PETScWrappers::Vector localized_solution (solution);
      hanging_node_constraints.distribute (localized_solution);
      solution = localized_solution;

      return solver_control.last_step();
   };

   void output_results() const {

      const PETScWrappers::Vector localized_solution (solution);

      if (this_mpi_process == 0) {
         DataOut<dim> data_out;

         data_out.attach_dof_handler(dof_handler);
         data_out.add_data_vector(localized_solution, "V");

         data_out.build_patches();

         const std::string filename = "output/solution-" + Utilities::int_to_string(timestep_number, 4) + ".vtk";
         std::ofstream output(filename.c_str());
         data_out.write_vtk(output);
      }
   };
};


