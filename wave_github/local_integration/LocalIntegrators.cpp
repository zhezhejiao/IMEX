

#include "LocalIntegrators.h"

namespace Main {
namespace Assembling{

template <int dim>
void LocalIntegrator<dim>::assemble_cell_stiffness(const FEValuesBase<dim> &fe_v,
						FullMatrix<double> 	&cell_stiff_matrix,
						double factor) const
{
	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    for (unsigned int q_point=0; q_point<fe_v.n_quadrature_points; ++q_point)
    {
    	for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
    		for (unsigned int j=0; j<dofs_per_cell; ++j)
    		{
    			cell_stiff_matrix(i,j) += factor * (fe_v.shape_grad(i,q_point) *
                                  fe_v.shape_grad(j,q_point) *
                                 fe_v.JxW(q_point));
    		}
        }
    }
}

template <int dim>
void LocalIntegrator<dim>::assemble_bnd_stiffness(const FEValuesBase<dim> &fe_v,
						FullMatrix<double> 	&cell_bnd_stiff_matrix,
						double factor) const
{
    const std::vector<Tensor <1,dim>> &normals = fe_v.get_normal_vectors ();
	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    for (unsigned int q_point=0; q_point<fe_v.n_quadrature_points; ++q_point)
    {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
        	Tensor<1,dim,double> grad_i = fe_v.shape_grad(i,q_point);
        	
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
            	Tensor<1,dim,double> grad_j = fe_v.shape_grad(j,q_point);
            	
            	cell_bnd_stiff_matrix(i,j) +=  factor *(((grad_i -grad_i*normals[q_point] * normals[q_point] ) *
              									(grad_j -grad_j*normals[q_point] * normals[q_point]) //+  value_i * value_j)*
												)*fe_v.JxW(q_point));
            }
        }
    }
}

template <int dim>
void LocalIntegrator<dim>::assemble_cell_damping(const FEValuesBase<dim> &fe_v,
						FullMatrix<double> 	&cell_damp_matrix,
						double factor,
						double alpha,
						double beta) const
{
	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    for (unsigned int q_point=0; q_point<fe_v.n_quadrature_points; ++q_point)
    {
    	for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
    		for (unsigned int j=0; j<dofs_per_cell; ++j)
    		{
    			cell_damp_matrix(i,j) += factor * ( alpha * fe_v.shape_value(i,q_point) * fe_v.shape_value(j,q_point)* fe_v.JxW(q_point));
    		}
        }
    }
}

template <int dim>
void LocalIntegrator<dim>::assemble_cell_mass(const FEValuesBase<dim> &fe_v,
						FullMatrix<double> 	&cell_mass_matrix) const
{
	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
	cell_mass_matrix = 0;

    for (unsigned int q_point=0; q_point<fe_v.n_quadrature_points; ++q_point)
    {
    	for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
    		for (unsigned int j=0; j<dofs_per_cell; ++j)
    		{
                cell_mass_matrix(i,j) += (fe_v.shape_value(i,q_point) *
                                    fe_v.shape_value(j,q_point)*
                                   fe_v.JxW(q_point));
    		}
        }
    }
}

template <int dim>
void LocalIntegrator<dim>::assemble_bnd_mass(const FEValuesBase<dim> &fe_v,
						FullMatrix<double> 	&cell_bnd_mass_matrix) const
{
	const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    for (unsigned int q_point=0; q_point<fe_v.n_quadrature_points; ++q_point)
    {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
            	cell_bnd_mass_matrix(i,j) += (fe_v.shape_value(i,q_point) *
                                    fe_v.shape_value(j,q_point)*
                                   fe_v.JxW(q_point));
            }
        }
    }
}

//--------------------------------------------------------------------------------------

template <int dim>
SystemIntegrator<dim>::SystemIntegrator()
  :
  MeshWorker::LocalIntegrator<dim>(true, true, false)
{}

template <int dim>
void SystemIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo,
		  MeshWorker::IntegrationInfo<dim> &info) const
{
	li.assemble_cell_mass(info.fe_values(0), dinfo.matrix(0,false).matrix);
	li.assemble_cell_stiffness(info.fe_values(0), dinfo.matrix(0,false).matrix,std::pow(global::stepsize,2)/4);
	li.assemble_cell_damping(info.fe_values(0), dinfo.matrix(0,false).matrix,global::stepsize/2);
}

template <int dim>
void SystemIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo,
              MeshWorker::IntegrationInfo<dim> &info) const
{
	li.assemble_bnd_mass(info.fe_values(0), dinfo.matrix(0,false).matrix);
	li.assemble_bnd_stiffness(info.fe_values(0), dinfo.matrix(0,false).matrix,std::pow(global::stepsize,2)/4);
}



template class LocalIntegrator<2>;
template class LocalIntegrator<3>;

template class SystemIntegrator<2>;
template class SystemIntegrator<3>;

} /* namespace Assembling */
} /* namespace Main */
