

#include "InverseMatrix_ML.h"

namespace Main {
namespace DataTypes{

template <int dim>
InverseMatrix_ML<dim>::InverseMatrix_ML(DoFHandler<dim>& dof_handler,
		const MappingQ<dim>& mapping,
		FE_Q<dim>& fe,
		const unsigned int degree)
{
	
	Quadrature<dim>*   	qf_ml = 0;
	Quadrature<dim-1>*	qf_face_ml = 0;
	if(degree == 1)
	{
		qf_ml = new QTrapez<dim>();
		qf_face_ml = new QTrapez<dim-1>();
	}
	else if(degree == 2)
	{
		qf_ml = new QSimpson<dim>();
		qf_face_ml = new QSimpson<dim-1>();
	}
	M_ml_inv.get_vector().reinit(dof_handler.n_dofs());


    Assembling::LocalIntegrator<dim> li;

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    FullMatrix<double> cell_matrix_M_in (dofs_per_cell, dofs_per_cell),
					   cell_matrix_M_bnd (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim>  fe_values (mapping, fe, *qf_ml, update_values|update_quadrature_points|update_JxW_values);
    FEFaceValues<dim> fe_face_values (mapping,fe, *qf_face_ml, update_values|update_quadrature_points|update_JxW_values);


    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
    	cell->get_dof_indices (local_dof_indices);
        cell_matrix_M_in = 0;
        cell_matrix_M_bnd = 0;
        fe_values.reinit (cell);

        li.assemble_cell_mass(fe_values, cell_matrix_M_in);

        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          if (cell->face(face_number)->at_boundary())
            {
              fe_face_values.reinit (cell, face_number);
              li.assemble_bnd_mass(fe_face_values,cell_matrix_M_bnd);
            }
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
        	M_ml_inv.add(local_dof_indices[i], local_dof_indices[i], cell_matrix_M_in(i,i));
        	if(cell_matrix_M_bnd(i,i)!=0)
        	{
        		M_ml_inv.add (local_dof_indices[i], local_dof_indices[i], cell_matrix_M_bnd(i,i));
        	}
          }
      }
    for(unsigned int i=0; i<dof_handler.n_dofs();i++)
    {
    	M_ml_inv(i,i) = 1/M_ml_inv(i,i);
    }
}


template <int dim>
unsigned int InverseMatrix_ML<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src)
{
	M_ml_inv.vmult(dst,src);
	return 1;
}

template <int dim>
unsigned int InverseMatrix_ML<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const SparseMatrix<double>&)
{
	M_ml_inv.vmult(dst,src);
	return 1;
}

template <int dim>
unsigned int InverseMatrix_ML<dim>::vmult(Vector<double>& dst,
						const Vector<double>& src,
						const double)
{
	M_ml_inv.vmult(dst,src);
	return 1;
}


template class InverseMatrix_ML<2>;
template class InverseMatrix_ML<3>;
} /* namespace DataTypes */
} /* namespace Main */
