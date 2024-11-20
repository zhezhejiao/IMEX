

#include "InverseMatrixMG.h"

namespace Main {
namespace DataTypes{

template <int dim>
InverseMatrixMG<dim>::InverseMatrixMG(const SparseMatrix<double>& S,
		double tol,
		DoFHandler<dim>& dof_handler,
		const Triangulation<dim>& triangulation,
		const MappingQ<dim>& mapping,
		FE_Q<dim>& fe)
:
	S(S),
	tol(tol),
	solver_control(global::config.get<double>("Solver_Settings.Max_itersteps"), tol),
	dof_handler(dof_handler),
	triangulation(triangulation),
	mapping(mapping),
	fe(fe)
{
	if(global::config.get<std::string>("Solver_Settings.Stopping_Criterion") == "ScaledEuclideanNorm")
	{
		solver_control.set_tolerance(tol/global::meshwidth);
	}

	setup_multigrid();
	assemble_multigrid();

	mg_transfer.initialize_constraints(mg_constrained_dofs);
	mg_transfer.build_matrices(dof_handler);
	coarse_matrix.copy_from (mg_matrices[0]);
	coarse_grid_solver.initialize (coarse_matrix);
	mg_smoother.initialize(mg_matrices);
	mg_smoother.set_steps(2);
	mg_smoother.set_symmetric(true);
	mg_matrix.initialize(mg_matrices);
	mg_interface_up.initialize(mg_interface_in);
	mg_interface_down.initialize(mg_interface_out);
	mg = new Multigrid<Vector<double> >(mg_matrix,
	                                coarse_grid_solver,
	                                mg_transfer,
	                                mg_smoother,
	                                mg_smoother,
									0,
									numbers::invalid_unsigned_int,
									Multigrid<Vector<double>>::f_cycle);
	mg->set_edge_matrices(mg_interface_down, mg_interface_up);
	preconditioner = new PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
					(dof_handler, *mg, mg_transfer);
}

template <int dim>
InverseMatrixMG<dim>::~InverseMatrixMG() {
	delete preconditioner;
	delete mg;
}

template <int dim>
unsigned int InverseMatrixMG<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src)
{
	SolverGMRES<Vector<double>> solver(solver_control);
	solver.solve(S,dst,src,*preconditioner);
	return solver_control.last_step();
}

template <int dim>
unsigned int InverseMatrixMG<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const SparseMatrix<double>& IP)
{
	LA::SolverGMRES<Vector<double>> solver(solver_control);
	solver.solve(S,dst,src,*preconditioner,IP);

	return solver_control.last_step();
}

template <int dim>
unsigned int InverseMatrixMG<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const double tol)
{
	solver_control.set_tolerance(tol);
	SolverGMRES<Vector<double>> solver(solver_control);
	solver.solve(S,dst,src,*preconditioner);
	return solver_control.last_step();
}

template <int dim>
void InverseMatrixMG<dim>::setup_multigrid()
{
	dof_handler.distribute_mg_dofs();
    mg_constrained_dofs.clear();
    mg_constrained_dofs.initialize(dof_handler);
    const unsigned int n_levels = triangulation.n_levels();
    mg_interface_in.resize(0, n_levels-1);
    mg_interface_in.clear_elements ();
    mg_interface_out.resize(0, n_levels-1);
    mg_interface_out.clear_elements ();
    mg_matrices.resize(0, n_levels-1);
    mg_matrices.clear_elements ();
    mg_sparsity_patterns.resize(0, n_levels-1);
    for (unsigned int level=0; level<n_levels; ++level)
      {
        DynamicSparsityPattern dsp (dof_handler.n_dofs(level),
                                    dof_handler.n_dofs(level));
        MGTools::make_sparsity_pattern(dof_handler, dsp, level);
        mg_sparsity_patterns[level].copy_from (dsp);
        mg_matrices[level].reinit(mg_sparsity_patterns[level]);
        mg_interface_in[level].reinit(mg_sparsity_patterns[level]);
        mg_interface_out[level].reinit(mg_sparsity_patterns[level]);
      }
}

template <int dim>
void InverseMatrixMG<dim>::assemble_multigrid()
{
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients |update_quadrature_points;
    info_box.add_update_flags_all(update_flags);
    update_flags = update_normal_vectors;
    info_box.add_update_flags_boundary(update_flags);
    info_box.initialize(fe, mapping);
    MeshWorker::DoFInfo<dim> dof_info(dof_handler);
    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(mg_constrained_dofs);
    assembler.initialize(mg_matrices);
    assembler.initialize_interfaces(mg_interface_in, mg_interface_out);
    Assembling::SystemIntegrator<dim> matrix_integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_mg(), dof_handler.end_mg(),
      dof_info, info_box, matrix_integrator, assembler);
    const unsigned int nlevels = triangulation.n_levels();
    for (unsigned int level=0; level<nlevels; ++level)
      {
        for (unsigned int i=0; i<dof_handler.n_dofs(level); ++i)
          if (mg_constrained_dofs.is_boundary_index(level,i) ||
              mg_constrained_dofs.at_refinement_edge(level,i))
            mg_matrices[level].set(i, i, 1.);
      }
}

template class InverseMatrixMG<2>;
template class InverseMatrixMG<3>;
} /* namespace DataTypes */
} /* namespace Main */
