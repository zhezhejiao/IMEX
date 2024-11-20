

#include "InverseMatrix.h"

namespace Main {
namespace DataTypes{

template <int dim>
InverseMatrix<dim>::InverseMatrix(const SparseMatrix<double>& M, double tol, const bool symm)
:
		M(M),
		symm(symm),
		tol(tol),
		preconditioner(),
		preconditioner_symm(),
		solver_control(global::config.get<double>("Solver_Settings.Max_itersteps"), tol)
{
	if(global::config.get<std::string>("Solver_Settings.Stopping_Criterion") == "ScaledEuclideanNorm")
	{
		solver_control.set_tolerance(tol/global::meshwidth);
	}
	if (symm == false)
		preconditioner.initialize(M);
	else
		preconditioner_symm.initialize(M);
}


template <int dim>
unsigned int InverseMatrix<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src)
{
	if (symm == false)
	{
		SolverGMRES<Vector<double>> solver(solver_control);
		solver.solve(M,dst,src,preconditioner);
	}
	else
	{
		SolverCG<Vector<double>> solver(solver_control);
		solver.solve(M,dst,src,preconditioner_symm);
	}
	return solver_control.last_step();
}

template <int dim>
unsigned int InverseMatrix<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const SparseMatrix<double>& IP)
{
	LA::SolverGMRES<Vector<double>> solver(solver_control);
	solver.solve(M,dst,src,preconditioner,IP);

	return solver_control.last_step();
}

template <int dim>
unsigned int InverseMatrix<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const double tol)
{
	solver_control.set_tolerance(tol);
	if (symm == false)
	{
		SolverGMRES<Vector<double>> solver(solver_control);
		solver.solve(M,dst,src,preconditioner);
	}
	else
	{
		SolverCG<Vector<double>> solver(solver_control);
		solver.solve(M,dst,src,preconditioner_symm);
	}
	return solver_control.last_step();
}


template class InverseMatrix<2>;
template class InverseMatrix<3>;
} /* namespace DataTypes */
} /* namespace Main */
