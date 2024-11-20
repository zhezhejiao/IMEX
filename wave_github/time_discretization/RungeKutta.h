/*
 * RungeKutta.h
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#ifndef TIME_INTEGRATION_RUNGEKUTTA_H_
#define TIME_INTEGRATION_RUNGEKUTTA_H_

#include <deal.II/lac/diagonal_matrix.h>

#include "../main/global.h"
#include "../exact_solution/Rhs.h"
#include "../InverseMatrix/InverseMatrix.h"
#include "../InverseMatrix/InverseMatrix_ML.h"
#include "TimeIntegrator.h"

using namespace dealii;

namespace Main {
namespace TimeIntegration{


template <int dim>
class RungeKutta  : public TimeIntegrator<dim>{
public:
	RungeKutta(const SparseMatrix<double>& M,
			const SparseMatrix<double>& A,
			const SparseMatrix<double>& B,
			Data::Rhs<dim>& rhs,
			const double stepsize);
	
	RungeKutta(const SparseMatrix<double>& A,
			const SparseMatrix<double>& B,
			Data::Rhs<dim>& rhs,
			const double stepsize,
			DoFHandler<dim>& dof_handler,
			const MappingQ<dim>& mapping,
			FE_Q<dim>& fe,
			double degree);
	~RungeKutta();
    void update_matrices(const double time);
	void integrate_step(Vector<double>& u,
						Vector<double>& v,
						const double time);

private:
	const SparseMatrix<double> &A, &B;
	Data::Rhs<dim>& rhs;
	const double stepsize;
	DataTypes::InverseMatrixBase<dim>* M_inv;
	Vector<double> u_tmp, V2, V3, V4, F1, F2, F3, F4,Bv,Bv2,Bv3,Bv4;
};


} /* namespace Time_Integration */
} /* namespace Main */

#endif /* TIME_INTEGRATION_RUNGE_H_ */
