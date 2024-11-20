/*
 * BS_FEM.h
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#ifndef BS_FEM_H_
#define BS_FEM_H_

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <list>

#include "../exact_solution/ExactSolution.h"
#include "../intial_value/InitialValue.h"
#include "../exact_solution/Rhs.h"
#include "../local_integration/LocalIntegrators.h"
#include "../time_discretization/TimeIntegrator.h"
#include "../time_discretization/CrankNicolson.h"
#include "../time_discretization/vIMEX.h"
#include "../time_discretization/rIMEX.h"
#include "../time_discretization/RungeKutta.h"
#include "OutputMatrix.h"
#include "global.h"
#include <sstream>


using namespace dealii;


namespace Main {

using namespace global;

template <int dim>
class BS_FEM {
public:
	BS_FEM				(const Triangulation<dim>&	triangulation,
			 	 	 	 const double				end_time,
						 const unsigned int			deg);
	~BS_FEM();

	void run			(const double total_timesteps);


private:

	void setup_and_assemble_system	();
	void setup_and_assemble_ref		();

	double compute_error			(double time);
	double compute_error_ref		(double time);

	void output_results				(const unsigned int	timestep_number) const;
	void output_ref_sol 			();

	const unsigned int				degree;
	const MappingQ<dim>				mapping;
	const Triangulation<dim>&		triangulation;
	DoFHandler<dim>					dof_handler;
	FE_Q<dim>						fe;
	std::vector< types::global_dof_index>  dof_to_boundary_mapping;
	const QGauss<dim>				quadrature;
	const QGauss<dim-1>				face_quadrature;
	const double					end_time;

	SparsityPattern					sparsity_pattern;
	SparsityPattern					sparsity_pattern_bnd;
	SparseMatrix<double>			M, M_in, M_bnd, A, B;
	Vector<double>					u, v;

	//For reference solution
	Triangulation<dim>				triangulation_ref;
	DoFHandler<dim>					dof_handler_ref;
	SparsityPattern					sparsity_pattern_ref;
	SparseMatrix<double>			M_ref, A_ref;
	Vector<double>					u_ref, v_ref;

};

} /* namespace Main */

#endif /* BS_FEM_H_ */
