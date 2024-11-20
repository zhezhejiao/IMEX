/*
 * Rhs.h
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#ifndef DATA_RHS_H_
#define DATA_RHS_H_

#include <deal.II/fe/mapping_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/numbers.h>
#include <deal.II/lac/identity_matrix.h>

#include "rhsfunctions.h"
#include "../main/global.h"


using namespace dealii;

namespace Main {
namespace Data {


template <int dim>
class Rhs {

public:
	Rhs(const Mapping<dim>& mapping,
		const DoFHandler<dim>& dof_handler,
		std::vector< types::global_dof_index>&  dof_to_boundary_mapping,
		const SparseMatrix<double>& M_bnd,
		const SparseMatrix<double>& M_in);

	virtual ~Rhs();


	void get_rhs_vector(Vector<double>& f_u,
						const double time,
						const Vector<double>& u);

	void get_rhs_vector_tx(Vector<double>& f_u,
						const double time);

	void get_rhs_vector_u(Vector<double>& f_u,
						const Vector<double>& u);

	void get_rhs_der(SparseMatrix<double>& Dfu,
						const double time,
						const Vector<double>& u);


private:
	Rhs_func_int<dim> rhs_int;
	Rhs_func_bnd<dim> rhs_bnd;

	const Mapping<dim>& mapping;
	const DoFHandler<dim>& dof_handler;
	std::vector< types::global_dof_index>&  dof_to_boundary_mapping;

	const SparseMatrix<double> &M_bnd, &M_in;

	Vector<double> tmp_in;
	Vector<double> interior_values;
	Vector<double> tmp_bnd,tmp_bnd_2;

	std::map<types::global_dof_index,double>  boundary_values;

	const double PI2;

};

} /* Data */
} /* namespace Main */

#endif /* DATA_RHS_H_ */
