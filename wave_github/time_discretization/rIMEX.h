/*
 * rIMEX.h
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#ifndef TIME_INTEGRATION_RIMEX_H_
#define TIME_INTEGRATION_RIMEX_H_

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>


#include "../main/global.h"
#include "../exact_solution/Rhs.h"
#include "../InverseMatrix/InverseMatrix.h"
#include "../InverseMatrix/InverseMatrixMG.h"
#include "TimeIntegrator.h"

using namespace dealii;

namespace Main {
namespace TimeIntegration{

template <int dim>
class rIMEX : public TimeIntegrator<dim>{
public:
	rIMEX(const SparseMatrix<double>& M,
			const SparseMatrix<double>& A,
			const SparseMatrix<double>& B,
			Data::Rhs<dim>& rhs,
			const double stepsize);

	rIMEX(const SparseMatrix<double>& M,
			const SparseMatrix<double>& A,
			const SparseMatrix<double>& B,
			Data::Rhs<dim>& rhs,
			const double stepsize,
			DoFHandler<dim>& dof_handler,
			const Triangulation<dim>& triangulation,
			const MappingQ<dim>& mapping,
			FE_Q<dim>& fe);

	~rIMEX();

	void update_matrices(const double time);

	void integrate_step(Vector<double>& u,
						Vector<double>& v,
						const double time);

private:
	const SparseMatrix<double> &M, &A, &B;
	Data::Rhs<dim>& rhs;
	const double stepsize;
	SparseMatrix<double> S, IP;
	DataTypes::InverseMatrixBase<dim>* S_inv;

	Vector<double> system_rhs, vnh, Fn, scaled_Mv;
};

} /* namespace Time_Integration */
} /* namespace Main */

#endif /* TIME_INTEGRATION_RIMEX_H_ */
