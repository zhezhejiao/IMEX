

#ifndef DATA_TYPES_INVERSEMATRIX_H_
#define DATA_TYPES_INVERSEMATRIX_H_

#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_mic.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>

#include "../main/global.h"
#include "../lac/solver_gmres.h"
#include "InverseMatrixBase.h"


using namespace dealii;

namespace Main {
namespace DataTypes{


template <int dim>
class InverseMatrix : public InverseMatrixBase<dim>{
public:
	InverseMatrix(const SparseMatrix<double>& M, double tol, const bool symm = false);
	InverseMatrix();
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const SparseMatrix<double> &IP);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const double tol);

private:
	const SparseMatrix<double> &M;
	const bool symm;
	double tol;
	SparseILU<double> preconditioner;
	
	PreconditionJacobi<SparseMatrix<double>> preconditioner_symm;
	SolverControl solver_control;
};

} /* namespace DataTypes */
} /* namespace Main */

#endif /* DATA_TYPES_INVERSEMATRIX_H_ */
