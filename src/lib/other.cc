void run() {
  const unsigned int initial_global_refinement = discretization;
  PETScWrappers::MPI::Vector tmp;
  PETScWrappers::MPI::Vector forcing_terms;

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
  assemble_system();

  //PDS: reinit with the extra local/global info required
  const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);
  tmp.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
  forcing_terms.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);


  // Set Initial Conditions
  VectorTools::interpolate(dof_handler, InitialValue<dim>(), old_solution);

  solution = old_solution;
  timestep_number = 0;
  time            = 0;

  output_results();


  while (time <= t_max) {
    time += time_step;
    ++timestep_number;

    pcout << "Time step " << timestep_number << " at t=" << time << std::endl;

    // system_rhs = M*old_solution;
    //mass_matrix.vmult(system_rhs, old_solution);

    // tmp = A*old_solution;
    // system_rhs += -h*(1-theta)*temp
    rhs_matrix.vmult(system_rhs, old_solution);
    // system_rhs.add((1.0 - theta) * time_step, tmp);

    // The second piece is to compute the contributions of the source
    // terms. This corresponds to the term $k_n
    // \left[ (1-\theta)F^{n-1} + \theta F^n \right]$. The following
    // code calls VectorTools::create_right_hand_side to compute the
    // vectors $F$, where we set the time of the right hand side
    // (source) function before we evaluate it. The result of this
    // all ends up in the forcing_terms variable:

    /*
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
    */

    int niter = solve_time_step();
    output_results();
    old_solution = solution;
  }
};