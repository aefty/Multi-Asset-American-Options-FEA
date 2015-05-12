// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------
using namespace std;

#ifndef __deal2__matrix_tools_h
#define __deal2__matrix_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/function_map.h>

#include <map>

// forward declarations
template <int dim> class Quadrature;

template<typename number> class Vector;
template<typename number> class FullMatrix;
template<typename number> class SparseMatrix;

template <typename number> class BlockSparseMatrix;
template <typename Number> class BlockVector;

template <int dim, int spacedim> class Mapping;
template <int dim, int spacedim> class DoFHandler;
template <int dim, int spacedim> class MGDoFHandler;
template <int dim, int spacedim> class FEValues;

namespace hp {
    template <int> class QCollection;
    template <int, int> class MappingCollection;
    template <int, int> class DoFHandler;
}


#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers {
    class SparseMatrix;
    class Vector;
    namespace MPI {
        class SparseMatrix;
        class BlockSparseMatrix;
        class Vector;
        class BlockVector;
    }
}
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers {
    class SparseMatrix;
    class Vector;
    class BlockSparseMatrix;
    class BlockVector;
    namespace MPI {
        class Vector;
        class BlockVector;
    }
}
#endif


namespace MatrixCreator_update {

    // create_advection_matrix
    {
        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const Mapping<dim, spacedim>       &mapping,
                                      const DoFHandler<dim, spacedim>    &dof,
                                      const Quadrature<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const DoFHandler<dim, spacedim>    &dof,
                                      const Quadrature<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const Mapping<dim, spacedim>   &mapping,
                                      const DoFHandler<dim, spacedim> &dof,
                                      const Quadrature<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> &rhs,
                                      Vector<double>           &rhs_vector,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const DoFHandler<dim, spacedim> &dof,
                                      const Quadrature<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> &rhs,
                                      Vector<double>           &rhs_vector,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const hp::MappingCollection<dim, spacedim>       &mapping,
                                      const hp::DoFHandler<dim, spacedim>    &dof,
                                      const hp::QCollection<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const hp::DoFHandler<dim, spacedim>    &dof,
                                      const hp::QCollection<dim>    &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const hp::MappingCollection<dim, spacedim> &mapping,
                                      const hp::DoFHandler<dim, spacedim> &dof,
                                      const hp::QCollection<dim> &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> &rhs,
                                      Vector<double>           &rhs_vector,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, typename number, int spacedim>
        void create_advection_matrix (const hp::DoFHandler<dim, spacedim> &dof,
                                      const hp::QCollection<dim> &q,
                                      SparseMatrix<number>     &matrix,
                                      const Function<spacedim> &rhs,
                                      Vector<double>           &rhs_vector,
                                      const Function<spacedim> *const a = 0,
                                      const ConstraintMatrix   &constraints = ConstraintMatrix());
    }

    // create_boundary_advection_matrix
    {
        template <int dim, int spacedim>
        void create_boundary_advection_matrix (
            const Mapping<dim, spacedim>       &mapping,
            const DoFHandler<dim, spacedim>    &dof,
            const Quadrature < dim - 1 >  &q,
            SparseMatrix<double>     &matrix,
            const typename FunctionMap<spacedim>::type & boundary_functions,
            Vector<double>           &rhs_vector,
            std::vector<types::global_dof_index> &dof_to_boundary_mapping,
            const Function<spacedim> *const weight = 0,
            std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

        /**
         * Calls the create_boundary_advection_matrix() function, see above, with
         * <tt>mapping=MappingQ1@<dim@>()</tt>.
         */
        template <int dim, int spacedim>
        void create_boundary_advection_matrix (
            const DoFHandler<dim, spacedim>    &dof,
            const Quadrature < dim - 1 >  &q,
            SparseMatrix<double>     &matrix,
            const typename FunctionMap<spacedim>::type        & boundary_functions,
            Vector<double>           &rhs_vector,
            std::vector<types::global_dof_index> &dof_to_boundary_mapping,
            const Function<spacedim> *const a = 0,
            std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, int spacedim>
        void create_boundary_advection_matrix (
            const hp::MappingCollection<dim, spacedim>       &mapping,
            const hp::DoFHandler<dim, spacedim>    &dof,
            const hp::QCollection < dim - 1 >  &q,
            SparseMatrix<double>     &matrix,
            const typename FunctionMap<spacedim>::type & boundary_functions,
            Vector<double>           &rhs_vector,
            std::vector<types::global_dof_index> &dof_to_boundary_mapping,
            const Function<spacedim> *const a = 0,
            std::vector<unsigned int> component_mapping = std::vector<unsigned int>());


        /**
         * Same function as above, but for hp objects.
         */
        template <int dim, int spacedim>
        void create_boundary_advection_matrix (
            const hp::DoFHandler<dim, spacedim>    &dof,
            const hp::QCollection < dim - 1 >  &q,
            SparseMatrix<double>     &matrix,
            const typename FunctionMap<spacedim>::type        & boundary_functions,
            Vector<double>           &rhs_vector,
            std::vector<types::global_dof_index> &dof_to_boundary_mapping,
            const Function<spacedim> *const a = 0,
            std::vector<unsigned int> component_mapping = std::vector<unsigned int>());
    }


    // create_laplace_matrix
    {
        template <int dim, int spacedim>
        void create_laplace_matrix (
            const Mapping<dim,
            spacedim> &mapping,
            const DoFHandler<dim,
            spacedim> &dof,
            const Quadrature<dim> &q,
            SparseMatrix<double> &matrix,
            const Function<spacedim> *const a = 0,
            const ConstraintMatrix &constraints = ConstraintMatrix());

        /**
         * Calls the create_laplace_matrix() function, see above, with
         * <tt>mapping=MappingQ1@<dim@>()</tt>.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (
            const DoFHandler<dim,
            spacedim> &dof,
            const Quadrature<dim> &q,
            SparseMatrix<double> &matrix,
            const Function<spacedim> *const a = 0,
            const ConstraintMatrix   &constraints = ConstraintMatrix());
    }

    //create_laplace_matrix
    {
        template <int dim, int spacedim>
        void create_laplace_matrix (const Mapping<dim, spacedim>   &mapping,
                                    const DoFHandler<dim, spacedim> &dof,
                                    const Quadrature<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim> &rhs,
                                    Vector<double>           &rhs_vector,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Calls the create_laplace_matrix() function, see above, with
         * <tt>mapping=MappingQ1@<dim@>()</tt>.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (const DoFHandler<dim, spacedim> &dof,
                                    const Quadrature<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim> &rhs,
                                    Vector<double>           &rhs_vector,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Like the functions above, but for hp dof handlers, mappings, and
         * quadrature collections.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (const hp::MappingCollection<dim, spacedim> &mapping,
                                    const hp::DoFHandler<dim, spacedim> &dof,
                                    const hp::QCollection<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Like the functions above, but for hp dof handlers, mappings, and
         * quadrature collections.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (const hp::DoFHandler<dim, spacedim> &dof,
                                    const hp::QCollection<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Like the functions above, but for hp dof handlers, mappings, and
         * quadrature collections.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (const hp::MappingCollection<dim, spacedim> &mapping,
                                    const hp::DoFHandler<dim, spacedim> &dof,
                                    const hp::QCollection<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim>      &rhs,
                                    Vector<double>           &rhs_vector,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());

        /**
         * Like the functions above, but for hp dof handlers, mappings, and
         * quadrature collections.
         */
        template <int dim, int spacedim>
        void create_laplace_matrix (const hp::DoFHandler<dim, spacedim> &dof,
                                    const hp::QCollection<dim>    &q,
                                    SparseMatrix<double>     &matrix,
                                    const Function<spacedim>      &rhs,
                                    Vector<double>           &rhs_vector,
                                    const Function<spacedim> *const a = 0,
                                    const ConstraintMatrix   &constraints = ConstraintMatrix());
    }

    /**
     * Exception
     */
    DeclException0 (ExcComponentMismatch);
}


namespace MatrixTools {
    /**
     * Import namespace MatrixCreator_update for backward compatibility with older
     * versions of deal.II in which these namespaces were classes and class
     * MatrixTools was publicly derived from class MatrixCreator_update.
     */
    using namespace MatrixCreator_update;

    /**
     * Apply Dirichlet boundary conditions to the system matrix and vectors as
     * described in the general documentation.
     */
    template <typename number>
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                SparseMatrix<number>  &matrix,
                                Vector<number>        &solution,
                                Vector<number>        &right_hand_side,
                                const bool             eliminate_columns = true);


    template <typename number>
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                BlockSparseMatrix<number>           &matrix,
                                BlockVector<number>                 &solution,
                                BlockVector<number>                 &right_hand_side,
                                const bool           eliminate_columns = true);

    #ifdef DEAL_II_WITH_PETSC
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values, PETScWrappers::SparseMatrix  &matrix, PETScWrappers::Vector  &solution, PETScWrappers::Vector  &right_hand_side, const bool eliminate_columns = true);
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values, PETScWrappers::MPI::SparseMatrix  &matrix, PETScWrappers::MPI::Vector  &solution, PETScWrappers::MPI::Vector  &right_hand_side, const bool eliminate_columns = true);
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                PETScWrappers::MPI::SparseMatrix  &matrix,
                                PETScWrappers::Vector       &solution,
                                PETScWrappers::MPI::Vector  &right_hand_side,
                                const bool             eliminate_columns = true);
    void apply_boundary_values (const std::map<types::global_dof_index, double>  &boundary_values,
                                PETScWrappers::MPI::BlockSparseMatrix &matrix,
                                PETScWrappers::MPI::BlockVector        &solution,
                                PETScWrappers::MPI::BlockVector        &right_hand_side,
                                const bool       eliminate_columns = true);
    #endif

    #ifdef DEAL_II_WITH_TRILINOS

    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                TrilinosWrappers::SparseMatrix  &matrix,
                                TrilinosWrappers::Vector        &solution,
                                TrilinosWrappers::Vector        &right_hand_side,
                                const bool             eliminate_columns = true);
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                TrilinosWrappers::BlockSparseMatrix  &matrix,
                                TrilinosWrappers::BlockVector        &solution,
                                TrilinosWrappers::BlockVector        &right_hand_side,
                                const bool                eliminate_columns = true);

    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                TrilinosWrappers::SparseMatrix  &matrix,
                                TrilinosWrappers::MPI::Vector   &solution,
                                TrilinosWrappers::MPI::Vector   &right_hand_side,
                                const bool             eliminate_columns = true);

    /**
     * This function does the same as the one above, except now working on block
     * structures.
     */
    void apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                TrilinosWrappers::BlockSparseMatrix  &matrix,
                                TrilinosWrappers::MPI::BlockVector   &solution,
                                TrilinosWrappers::MPI::BlockVector   &right_hand_side,
                                const bool                eliminate_columns = true);
    #endif

    /**
     * Rather than applying boundary values to the global matrix and vector
     * after creating the global matrix, this function does so during assembly,
     * by modifying the local matrix and vector contributions. If you call this
     * function on all local contributions, the resulting matrix will have the
     * same entries, and the final call to apply_boundary_values() on the global
     * system will not be necessary.
     *
     * Since this function does not have to work on the complicated data
     * structures of sparse matrices, it is relatively cheap. It may therefore
     * be a win if you have many fixed degrees of freedom (e.g. boundary nodes),
     * or if access to the sparse matrix is expensive (e.g. for block sparse
     * matrices, or for PETSc or trilinos matrices). However, it doesn't work as
     * expected if there are also hanging nodes to be considered. More caveats
     * are listed in the general documentation of this class.
     */
    void local_apply_boundary_values (const std::map<types::global_dof_index, double> &boundary_values,
                                      const std::vector<types::global_dof_index> &local_dof_indices,
                                      FullMatrix<double> &local_matrix,
                                      Vector<double>     &local_rhs,
                                      const bool          eliminate_columns);

    /**
     * Exception
     */
    DeclException0 (ExcBlocksDontMatch);
}



//DEAL_II_NAMESPACE_CLOSE

#endif
