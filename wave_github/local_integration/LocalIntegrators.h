

#ifndef ASSEMBLING_LOCALINTEGRATORS_H_
#define ASSEMBLING_LOCALINTEGRATORS_H_

#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>

#include "../main/global.h"

using namespace dealii;

namespace Main {
namespace Assembling{


template <int dim>
class LocalIntegrator
{
public:
	LocalIntegrator(){}
	virtual ~LocalIntegrator(){}

	void assemble_cell_stiffness(const FEValuesBase<dim> &fe_v,
							FullMatrix<double> 	&cell_stiff_matrix,
							double factor = 1) const;

	void assemble_bnd_stiffness(const FEValuesBase<dim> &fe_v,
							FullMatrix<double> 	&cell_bnd_stiff_matrix,
							double factor = 1) const;

	void assemble_cell_damping(const FEValuesBase<dim> &fe_v,
							FullMatrix<double> 	&cell_damp_matrix,
							double factor = 1,
							double alpha = 1,
							double beta = 1) const;

	void assemble_cell_mass(const FEValuesBase<dim> &fe_v,
							FullMatrix<double> 	&cell_mass_matrix) const;

	void assemble_bnd_mass(const FEValuesBase<dim> &fe_v,
							FullMatrix<double> 	&cell_bnd_mass_matrix) const;
};


template <int dim>
class SystemIntegrator  : public MeshWorker::LocalIntegrator<dim> {
public:
	SystemIntegrator();
	virtual ~SystemIntegrator(){}

	void cell(MeshWorker::DoFInfo<dim> &dinfo,
			  MeshWorker::IntegrationInfo<dim> &info) const;
	void boundary(MeshWorker::DoFInfo<dim> &dinfo,
	              MeshWorker::IntegrationInfo<dim> &info) const;

private:
	LocalIntegrator<dim> li;
};




} 
} /* namespace Main */

#endif /* ASSEMBLING_LOCALINTEGRATORS_H_ */
