

#ifndef DATA_TYPES_INVERSEMATRIXMG_H_
#define DATA_TYPES_INVERSEMATRIXMG_H_

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/output.h>
#include <deal.II/meshworker/loop.h>

#include "../main/global.h"
#include "../lac/solver_gmres.h"
#include "../local_integration/LocalIntegrators.h"
#include "InverseMatrixBase.h"

using namespace dealii;

namespace Main {
namespace DataTypes{


template <int dim>
class InverseMatrixMG : public InverseMatrixBase<dim>{
public:
	InverseMatrixMG(const SparseMatrix<double>& S,
					double tol,
					DoFHandler<dim>& dof_handler,
					const Triangulation<dim>& triangulation,
					const MappingQ<dim>& mapping,
					FE_Q<dim>& fe);
	virtual ~InverseMatrixMG();

	unsigned int vmult(Vector<double> &dst, const Vector<double> &src);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const SparseMatrix<double> &IP);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const double tol);


private:
	void setup_multigrid();
	void assemble_multigrid();

	const SparseMatrix<double> &S;
	double tol;
	SolverControl solver_control;

	DoFHandler<dim>& 			dof_handler;
	const Triangulation<dim>&	triangulation;
	const MappingQ<dim>&		mapping;
	FE_Q<dim>&					fe;

    MGConstrainedDoFs 					 mg_constrained_dofs;
    MGLevelObject<SparsityPattern>       mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double> > mg_matrices;
    MGLevelObject<SparseMatrix<double> > mg_interface_in;
    MGLevelObject<SparseMatrix<double> > mg_interface_out;



    mg::Matrix<Vector<double>> mg_matrix;
    mg::Matrix<Vector<double>> mg_interface_up;
    mg::Matrix<Vector<double>> mg_interface_down;

    MGTransferPrebuilt<Vector<double>> mg_transfer;
    FullMatrix<double> coarse_matrix;
    MGCoarseGridHouseholder<> coarse_grid_solver;
    mg::SmootherRelaxation<PreconditionSOR<SparseMatrix<double> >, Vector<double> > mg_smoother;
    Multigrid<Vector<double>>* mg;
	PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double> > >	*preconditioner;
};

} /* namespace DataTypes */
} /* namespace Main */

#endif /* DATA_TYPES_INVERSEMATRIXMG_H_ */
