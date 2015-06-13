
template<int dim>
class BlackScholes {

   Triangulation<dim>   triangulation;
   FE_Q<dim>            fe;
   DoFHandler<dim>      dof_handler;

   //PDS: We change the name of this member, as in the example
   //PDS:  there were no other explicit constraints
   ConstraintMatrix     hanging_node_constraints;

   double               time;
   const double         time_step; // PDS : we make this constant here
   unsigned int         timestep_number;


   double               t_max;
   const double         theta;
   std::vector <double> range;
   const double         discretization;

   Tensor<2, dim>  diffusion_const;
   Tensor<1, dim>  convection_const;
   double          reaction_const;

   // PDS: replace these data structures with PETSc ones
   PETScWrappers::MPI::SparseMatrix mass_matrix;
   PETScWrappers::MPI::SparseMatrix laplace_matrix;
   PETScWrappers::MPI::SparseMatrix system_matrix;
   PETScWrappers::MPI::SparseMatrix convection_matrix;

   PETScWrappers::MPI::SparseMatrix omega_matrix;

   PETScWrappers::MPI::Vector solution;
   PETScWrappers::MPI::Vector old_solution;
   PETScWrappers::MPI::Vector system_rhs;

   //PDS: data concerning an mpi communicator
   MPI_Comm mpi_communicator;
   const unsigned int n_mpi_processes;
   const unsigned int this_mpi_process;

   //PDS: a rank-0 output stream
   ConditionalOStream pcout;

 public:

   BlackScholes(Tensor<2, dim> DIFF_CONST, Tensor<1, dim> CONV_CONST , double R, double X_DISC , std::vector<double> X_RANGE, double T_DISC , double T_RANGE , double THETA  )
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

      // PDS: Initialize additional data members related to MPI
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
      //PDS : choose a set global refinement here
      const unsigned int initial_global_refinement = discretization;
      //const unsigned int n_adaptive_pre_refinement_steps = 4;



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

      triangulation.refine_global (initial_global_refinement);


      setup_system();

      // PDS: new assemble_system routine assembles the full operator to be inverted
      // PDS: this implies that we cannot change the timestep or theta
      // PDS: without calling this routine again
      assemble_system();

      unsigned int pre_refinement_step = 0;

      PETScWrappers::MPI::Vector tmp;
      PETScWrappers::MPI::Vector forcing_terms;

      //PDS: reinit with the extra local/global info required
      const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);
      tmp.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      forcing_terms.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

      //PDS: For now, we ignore adaptive mesh refinement and all the associated complication

      //InitialValue<dim> inital_value;

      //VectorTools::interpolate(dof_handler, inital_value, old_solution);

      VectorTools::interpolate(dof_handler, InitialValue<dim>(),
                               old_solution);

      solution = old_solution;

      timestep_number = 0;
      time            = 0;

      output_results();

      // Then we start the main loop until the computed time exceeds our
      // end time of 0.5. The first task is to build the right hand
      // side of the linear system we need to solve in each time step.
      // Recall that it contains the term $MU^{n-1}-(1-\theta)k_n AU^{n-1}$.
      // We put these terms into the variable system_rhs, with the
      // help of a temporary vector:

      // laplace_matrix.add(-1.0,  convection_matrix);
      //laplace_matrix.add(1.0,  mass_matrix);



      while (time <= t_max) {
         time += time_step;
         ++timestep_number;

         pcout << "Time step " << timestep_number << " at t=" << time << std::endl;

         // system_rhs = M*old_solution;
         //mass_matrix.vmult(system_rhs, old_solution);

         // tmp = A*old_solution;
         // system_rhs += -h*(1-theta)*temp
         omega_matrix.vmult(system_rhs, old_solution);
         // system_rhs.add((1.0 - theta) * time_step, tmp);

         // The second piece is to compute the contributions of the source
         // terms. This corresponds to the term $k_n
         // \left[ (1-\theta)F^{n-1} + \theta F^n \right]$. The following
         // code calls VectorTools::create_right_hand_side to compute the
         // vectors $F$, where we set the time of the right hand side
         // (source) function before we evaluate it. The result of this
         // all ends up in the forcing_terms variable:

         if (false) {
            RightHandSide<dim> rhs_function;
            rhs_function.set_time(time);

            //PDS: TVectorTools::create_right_hand_side is not available, so we must define our
            //PDS   own function
            // F^n
            assemble_F(rhs_function, tmp);
            forcing_terms = tmp;
            forcing_terms *= time_step * theta;


            //F^{n-1}
            rhs_function.set_time(time - time_step);
            assemble_F(rhs_function, tmp);
            forcing_terms.add(time_step * (1 - theta), tmp);
            // Next, we add the forcing terms to the ones that
            // come from the time stepping, and also build the matrix
            // $M+k_n\theta A$ that we have to invert in each time step.
            // The final piece of these operations is to eliminate
            // hanging node constrained degrees of freedom from the
            // linear system:
            system_rhs += forcing_terms;
         };

         // PDS system_matrix should have already been assembled (assuming a constant
         // PDS   time step and theta

         // There is one more operation we need to do before we
         // can solve it: boundary values. To this end, we create
         // a boundary value object, set the proper time to the one
         // of the current time step, and evaluate it as we have
         // done many times before. The result is used to also
         // set the correct boundary values in the linear system:

         if (false) {
            BoundaryValues<dim> boundary_values_function;
            boundary_values_function.set_time(time);

            std::map<types::global_dof_index, double> boundary_values;
            VectorTools::interpolate_boundary_values(dof_handler,
                  0,
                  boundary_values_function,
                  boundary_values);
            MatrixTools::apply_boundary_values(boundary_values,
                                               system_matrix,
                                               solution,
                                               system_rhs,
                                               false); // PDS : note last argument (see step-17 docs)
         };

         // With this out of the way, all we have to do is solve the
         // system, generate graphical data, and...
         //PDS This returns output now, so output number of iterations
         int niter = solve_time_step();

         // Error in converting to 3 digit number to string  ... ?
         //pcout << "Solve_time_step used " << niter << " iterations" << std::endl;

         output_results();

         // ...take care of mesh refinement. Here, what we want to do is
         // (i) refine the requested number of times at the very beginning
         // of the solution procedure, after which we jump to the top to
         // restart the time iteration, (ii) refine every fifth time
         // step after that.
         //
         // The time loop and, indeed, the main part of the program ends
         // with starting into the next time step by setting old_solution
         // to the solution we have just computed.

         //PDS : We ignore the refinement
         /*
         if ((timestep_number == 1) &&
             (pre_refinement_step < n_adaptive_pre_refinement_steps))
           {
             refine_mesh (initial_global_refinement,
                          initial_global_refinement + n_adaptive_pre_refinement_steps);
             ++pre_refinement_step;

             tmp.reinit (solution.size());
             forcing_terms.reinit (solution.size());

             std::cout << std::endl;

             goto start_time_iteration;
           }
         else if ((timestep_number > 0) && (timestep_number % 5 == 0))
           {
             refine_mesh (initial_global_refinement,
                          initial_global_refinement + n_adaptive_pre_refinement_steps);
             tmp.reinit (solution.size());
             forcing_terms.reinit (solution.size());
           }
         */
         old_solution = solution;
      }
   };


 private:
   void setup_system() {

      //PDS:
      GridTools::partition_triangulation (n_mpi_processes, triangulation);
      dof_handler.distribute_dofs (fe);

      //PDS:
      DoFRenumbering::subdomain_wise (dof_handler);
      const types::global_dof_index n_local_dofs
         = DoFTools::count_dofs_with_subdomain_association (dof_handler,
               this_mpi_process);

      //PDS: cout -- > pcout
      pcout << std::endl
            << "==========================================="
            << std::endl
            << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl
            << std::endl;

      //PDS: much changes below
      system_matrix.reinit (mpi_communicator,
                            dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            n_local_dofs,
                            n_local_dofs,
                            dof_handler.max_couplings_between_dofs());

      laplace_matrix.reinit (mpi_communicator,
                             dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             n_local_dofs,
                             n_local_dofs,
                             dof_handler.max_couplings_between_dofs());

      mass_matrix.reinit (mpi_communicator,
                          dof_handler.n_dofs(),
                          dof_handler.n_dofs(),
                          n_local_dofs,
                          n_local_dofs,
                          dof_handler.max_couplings_between_dofs());

      convection_matrix.reinit (mpi_communicator,
                                dof_handler.n_dofs(),
                                dof_handler.n_dofs(),
                                n_local_dofs,
                                n_local_dofs,
                                dof_handler.max_couplings_between_dofs());

      omega_matrix.reinit (mpi_communicator,
                           dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           n_local_dofs,
                           n_local_dofs,
                           dof_handler.max_couplings_between_dofs());


      solution.reinit     (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      old_solution.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      system_rhs.reinit   (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

      //PDS: Note name change
      hanging_node_constraints.clear ();
      DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
      hanging_node_constraints.close();

      // PDS: Note that it remains to assemble the system
   };

   void assemble_system() {
      QGauss<dim>  quadrature_formula(2);
      FEValues<dim> fe_values (fe, quadrature_formula,
                               update_values   | update_gradients |
                               update_quadrature_points | update_JxW_values);
      const unsigned int   dofs_per_cell = fe.dofs_per_cell;
      const unsigned int   n_q_points    = quadrature_formula.size();

      FullMatrix<double>   cell_matrix_laplace    (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_mass       (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_system     (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_convection (dofs_per_cell, dofs_per_cell);
      FullMatrix<double>   cell_matrix_omega      (dofs_per_cell, dofs_per_cell);

      Vector<double>       cell_rhs            (dofs_per_cell); // PDS: this is just a placeholder
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);


      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (; cell != endc; ++cell) {
         cell_matrix_laplace = 0;
         cell_matrix_mass    = 0;
         cell_matrix_convection    = 0;
         cell_rhs = 0; // PDS : dummy
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
                                                        fe_values.shape_grad (i, q_point) *
                                                        convection_const *
                                                        fe_values.shape_value (j, q_point) *
                                                        fe_values.JxW (q_point)
                                                     );
                  }
               }
            }

            //PDS: the system matrix is built from the mass and laplace matrices,
            //PDS   assuming a constant time step and theta
            /*
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
               for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                  cell_matrix_system(i, j) =
                     cell_matrix_mass(i, j) + (time_step * theta * cell_matrix_laplace(i, j));
               }
            }
            */

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
               for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                  cell_matrix_omega(i, j) =
                     -1.0 * cell_matrix_laplace(i, j) +
                     cell_matrix_convection(i, j) -
                     reaction_const * cell_matrix_mass(i, j);

                  cell_matrix_systemL(i, j) =
                     cell_matrix_mass(i, j) -
                     time_step * (theta) * cell_matrix_omega ;


                  cell_matrix_systemR(i, j) =
                     cell_matrix_mass(i, j) +
                     time_step * (1.0 - theta) * cell_matrix_omega ;
               }
            }

            cell->get_dof_indices (local_dof_indices);

            // System Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_system, cell_rhs,
               local_dof_indices, system_matrix
               , system_rhs); //PDS : system_rhs is garbage

            // Mass Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_mass, cell_rhs,
               local_dof_indices,
               mass_matrix, system_rhs); //PDS : system_rhs *garbage*

            // Lapalace Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_laplace, cell_rhs,
               local_dof_indices,
               laplace_matrix, system_rhs); //PDS : system_rhs *garbage*


            // convection Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_convection, cell_rhs,
               local_dof_indices,
               convection_matrix, system_rhs); //PDS : system_rhs *garbage*


            // omega Matrix
            hanging_node_constraints.distribute_local_to_global(
               cell_matrix_omega, cell_rhs,
               local_dof_indices,
               omega_matrix, system_rhs); //PDS : system_rhs *garbage*
         }
      }


      laplace_matrix.compress(VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);
      system_matrix.compress(VectorOperation::add);
      convection_matrix.compress(VectorOperation::add);
      omega_matrix.compress(VectorOperation::add);

      //PDS: Note that we do not enforce the boundary values here (but do so later on once
      //PDS   we have the rhs for the system, which varies with the time step)
   };

   void assemble_F(const RightHandSide<dim> &rhs_function, PETScWrappers::MPI::Vector &f) {
      QGauss<dim>  quadrature_formula(2);
      FEValues<dim> fe_values (fe, quadrature_formula,
                               update_values   | update_gradients |
                               update_quadrature_points | update_JxW_values);
      const unsigned int   dofs_per_cell = fe.dofs_per_cell;
      const unsigned int   n_q_points    = quadrature_formula.size();
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

      {
         const types::global_dof_index n_local_dofs
            = DoFTools::count_dofs_with_subdomain_association (dof_handler,
                  this_mpi_process);
         f.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
      }

      for (; cell != endc; ++cell) {
         if (cell->subdomain_id() == this_mpi_process) {
            cell->get_dof_indices (local_dof_indices);
            fe_values.reinit (cell);
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
               for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  f(local_dof_indices[i]) +=
                     fe_values.shape_value (i, q_point) *
                     rhs_function.value(fe_values.quadrature_point(q_point)) *
                     fe_values.JxW (q_point);
               }
            }
         }
      }
      f.compress(VectorOperation::add);
   };

   int solve_time_step() {

      SolverControl           solver_control (solution.size(),
                                              1e-8 * system_rhs.l2_norm());
      PETScWrappers::SolverCG cg (solver_control,
                                  mpi_communicator);
      // PDS: We use a very simple parallel preconditioner, block jacobi, which
      // PDS:  uses an ILU(0) preconditioner applied to each local block.
      // PDS: The original example used SOR, which is not trivial to parallelize.
      // PDS:Note that this is means the preconditioner *weakens* with core count
      // PDS: (going from 1 to 2 mpi processes will mean the iteration count will change from 4/5 to 8/9)
      PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

      //PDS : you could use a fancier preconditioner (provided there is a wrapper for it)
      //PETScWrappers::PreconditionBoomerAMG preconditioner(system_matrix);

      cg.solve (system_matrix, solution, system_rhs, preconditioner);

      return solver_control.last_step();
   };


   void output_results() const {

      const PETScWrappers::Vector localized_solution (solution);
      //PDS:  on rank 0 only
      if (!this_mpi_process) {
         DataOut<dim> data_out;

         data_out.attach_dof_handler(dof_handler);
         data_out.add_data_vector(localized_solution, "V");

         data_out.build_patches();

         const std::string filename = "output/solution-"
                                      + Utilities::int_to_string(timestep_number, 3) +
                                      ".vtk";
         std::ofstream output(filename.c_str());
         data_out.write_vtk(output);
      }
   };

};


