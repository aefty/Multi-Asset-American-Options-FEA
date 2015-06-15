
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
  int american;

  // PDS: replace these data structures with PETSc ones
  PETScWrappers::MPI::SparseMatrix rhs_matrix;
  PETScWrappers::MPI::SparseMatrix system_matrix;

  PETScWrappers::MPI::Vector solution;
  PETScWrappers::MPI::Vector old_solution;
  PETScWrappers::MPI::Vector system_rhs;


  PETScWrappers::MPI::SparseMatrix mass_matrix;
  PETScWrappers::MPI::SparseMatrix laplace_matrix;
  PETScWrappers::MPI::SparseMatrix convection_matrix;
  PETScWrappers::MPI::SparseMatrix omega_matrix;


  //PDS: data concerning an mpi communicator
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

    triangulation.refine_global (discretization);

    setup_system();
    assemble_system();

    //PDS: reinit with the extra local/global info required
    const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);
    tmp.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
    forcing_terms.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

    std::cout << dof_handler.n_dofs();

    // Set Initial Conditions
    VectorTools::interpolate(dof_handler, InitialValue<dim>(), old_solution);

    // Create matrixes
    //rhs_matrix.print(std::cout);
    //convection_matrix.print(std::cout);


    hanging_node_constraints.distribute(old_solution);

    solution = old_solution;
    timestep_number = 0;
    time            = 0;

    output_results();

    while (time <= (t_max - 2.0 * time_step)) {
      time += time_step;
      timestep_number = timestep_number + 1;

      pcout << "Time step " << timestep_number << " at t=" << time << std::endl;

      rhs_matrix.vmult(system_rhs, old_solution);

      int niter = solve_time_step();
      output_results();
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
    const types::global_dof_index n_local_dofs =
      DoFTools::count_dofs_with_subdomain_association (dof_handler, this_mpi_process);

    //PDS: cout -- > pcout
    pcout << std::endl
          << "===========================================" << std::endl
          << "Number of active cells: " << triangulation.n_active_cells() << std::endl
          << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl
          << "===========================================" << std::endl
          << "Dimension : " << dim << std::endl
          << "American : " << american << std::endl
          << "Interest Rate : " << reaction_const << std::endl
          << "Expiration Time : " << t_max << std::endl
          << "dt : " << time_step << std::endl;

    //PDS: much changes below
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


    mass_matrix.reinit (mpi_communicator,
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

    FullMatrix<double>   cell_matrix_laplace     (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   cell_matrix_mass        (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   cell_matrix_convection  (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   cell_matrix_omega       (dofs_per_cell, dofs_per_cell);

    FullMatrix<double>   cell_matrix_system     (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   cell_matrix_rhs        (dofs_per_cell, dofs_per_cell);

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
              cell_matrix_laplace(i, j) -
              cell_matrix_convection(i, j) +
              reaction_const * cell_matrix_mass(i, j);

            cell_matrix_system(i, j) =
              cell_matrix_mass(i, j) +
              time_step * (theta)  * cell_matrix_omega(i, j);

            cell_matrix_rhs(i, j) =
              cell_matrix_mass(i, j) -
              time_step * (1.0 - theta)  *  cell_matrix_omega(i, j);
          };
        };



        cell->get_dof_indices (local_dof_indices);

        // System Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_system, cell_rhs,
          local_dof_indices, system_matrix
          , system_rhs); //PDS : system_rhs is garbage

        // RHS Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_rhs, cell_rhs,
          local_dof_indices, rhs_matrix
          , system_rhs); //PDS : system_rhs is garbage


        // mass Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_mass, cell_rhs,
          local_dof_indices, mass_matrix
          , system_rhs); //PDS : system_rhs is garbage


        // lapalce Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_laplace, cell_rhs,
          local_dof_indices, laplace_matrix
          , system_rhs); //PDS : system_rhs is garbage


        // convection Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_convection, cell_rhs,
          local_dof_indices, convection_matrix
          , system_rhs); //PDS : system_rhs is garbage


        // omega Matrix
        hanging_node_constraints.distribute_local_to_global(
          cell_matrix_omega, cell_rhs,
          local_dof_indices, omega_matrix
          , system_rhs); //PDS : system_rhs is garbage

      };
    };




    system_matrix.compress(VectorOperation::add);
    rhs_matrix.compress(VectorOperation::add);


    mass_matrix.compress(VectorOperation::add);
    laplace_matrix.compress(VectorOperation::add);
    convection_matrix.compress(VectorOperation::add);

    omega_matrix.compress(VectorOperation::add);


    //PDS: Note that we do not enforce the boundary values here (but do so later on once
    //PDS   we have the rhs for the system, which varies with the time step)
  };

  void assemble_F(const RightHandSide<dim> &rhs_function, PETScWrappers::MPI::Vector & f) {
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

    SolverControl solver_control (solution.size(), 1e-10 * system_rhs.l2_norm());

    PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

    // PDS: We use a very simple parallel preconditioner, block jacobi, which
    // PDS:  uses an ILU(0) preconditioner applied to each local block.
    // PDS: The original example used SOR, which is not trivial to parallelize.
    // PDS:Note that this is means the preconditioner *weakens* with core count
    // PDS: (going from 1 to 2 mpi processes will mean the iteration count will change from 4/5 to 8/9)
    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

    //PDS : you could use a fancier preconditioner (provided there is a wrapper for it)
    //PETScWrappers::PreconditionBoomerAMG preconditioner(system_matrix);

    cg.solve (system_matrix, solution, system_rhs, preconditioner);

    PETScWrappers::Vector localized_solution (solution);
    hanging_node_constraints.distribute (localized_solution);
    solution = localized_solution;

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

      const std::string filename = "output/solution-" + Utilities::int_to_string(timestep_number, 4) + ".vtk";
      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
    }
  };
};


