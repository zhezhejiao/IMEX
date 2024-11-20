/*
 * Rhs.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include "Rhs.h"

namespace Main {
namespace Data {

template <int dim>
Rhs<dim>::Rhs(const Mapping<dim>& mapping,
			const DoFHandler<dim>& dof_handler,
			std::vector< types::global_dof_index>&  dof_to_boundary_mapping,
			const SparseMatrix<double>& M_bnd,
			const SparseMatrix<double>& M_in)
			:
			rhs_int(),
			rhs_bnd(),
			mapping(mapping),
			dof_handler(dof_handler),
			dof_to_boundary_mapping(dof_to_boundary_mapping),
			M_bnd(M_bnd),
			M_in(M_in),
			tmp_in(dof_handler.n_dofs()),
			interior_values(dof_handler.n_dofs()),
			tmp_bnd(dof_handler.n_boundary_dofs()),
			tmp_bnd_2(dof_handler.n_boundary_dofs()),
			PI2(numbers::PI*numbers::PI)
			{

				VectorTools::interpolate(mapping,dof_handler,rhs_int,interior_values);
				VectorTools::interpolate_boundary_values(mapping, dof_handler,0,rhs_bnd,boundary_values);
			}

template <int dim>
Rhs<dim>::~Rhs() {
	
}

template <int dim>
void Rhs<dim>::get_rhs_vector(Vector<double>& f_u,
					const double time,
					const Vector<double>& u)
{
	if (global::config.get<std::string>("Error.Error_against") == "Exact")
	{
		double sin2PIt = std::sin(2*numbers::PI*time);
		double cos2PIt = std::cos(2*numbers::PI*time);
		tmp_in = interior_values;
		tmp_in *= sin2PIt;

		for(unsigned int i = 0; i<dof_handler.n_dofs(); i++)
		{
			tmp_in(i) = (-fabs(tmp_in(i))  - 4*PI2)*tmp_in(i) + 2*(pow(time+1,-1))*numbers::PI*cos2PIt*interior_values(i) + fabs(u(i))*u(i);
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
			{
				tmp_bnd(dof_to_boundary_mapping[i]) = sin2PIt * boundary_values[i] * (-4*PI2 + 6 - pow(sin2PIt * boundary_values[i],2)) + pow(u[i],3);
			}
		}
		M_in.vmult(f_u,tmp_in);
		M_bnd.vmult(tmp_bnd_2,tmp_bnd);
		for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin(); it!=boundary_values.end(); ++it)
		{
			f_u[it->first] += tmp_bnd_2[dof_to_boundary_mapping[it->first]];
		}
	}


	else
	{
		f_u = 0;
		for(unsigned int i = 0; i<dof_handler.n_dofs(); i++)
		{
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
			{
				tmp_bnd(dof_to_boundary_mapping[i]) = -std::sin(u[i])*std::pow(u[i],2);
				
			}
		}
		M_bnd.vmult(tmp_bnd_2,tmp_bnd);
		for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin(); it!=boundary_values.end(); ++it)
		{
			f_u[it->first] += tmp_bnd_2[dof_to_boundary_mapping[it->first]];
		}
	}
}

template <int dim>
void Rhs<dim>::get_rhs_vector_tx(Vector<double>& f_u,
					const double time)
{

	if (global::config.get<std::string>("Error.Error_against") == "Exact")
	{
		double sin2PIt = std::sin(2*numbers::PI*time);
		double cos2PIt = std::cos(2*numbers::PI*time);
		tmp_in = interior_values;
		tmp_in *= sin2PIt;

		for(unsigned int i = 0; i<dof_handler.n_dofs(); i++)
		{
			tmp_in(i) = (-fabs(tmp_in(i))  - 4*PI2)*tmp_in(i) + (2*pow(time+1,-1))*numbers::PI*cos2PIt*interior_values(i);
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
			{
				tmp_bnd(dof_to_boundary_mapping[i]) = sin2PIt * boundary_values[i] * (-4*PI2 + 6 - pow(sin2PIt * boundary_values[i],2));
			}
		}
		M_in.vmult(f_u,tmp_in);
		M_bnd.vmult(tmp_bnd_2,tmp_bnd);
		for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin(); it!=boundary_values.end(); ++it)
		{
			f_u[it->first] += tmp_bnd_2[dof_to_boundary_mapping[it->first]];
		}
	}


	else
		f_u = 0;
}

template <int dim>
void Rhs<dim>::get_rhs_vector_u(Vector<double>& f_u,
					const Vector<double>& u)
{

	if (global::config.get<std::string>("Error.Error_against") == "Exact")
	{
		for(unsigned int i = 0; i<dof_handler.n_dofs(); i++)
		{
			tmp_in(i) =  fabs(u(i))*u(i);
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
			{
				tmp_bnd(dof_to_boundary_mapping[i]) =  pow(u[i],3);
			}
		}
		M_in.vmult(f_u,tmp_in);
		M_bnd.vmult(tmp_bnd_2,tmp_bnd);
		for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin(); it!=boundary_values.end(); ++it)
		{
			f_u[it->first] += tmp_bnd_2[dof_to_boundary_mapping[it->first]];
		}
	}

	else
	{
		f_u = 0;
		for(unsigned int i = 0; i<dof_handler.n_dofs(); i++)
		{
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
			{
				tmp_bnd(dof_to_boundary_mapping[i]) = -std::sin(u[i])*std::pow(u[i],2);
				
			}
		}
		M_bnd.vmult(tmp_bnd_2,tmp_bnd);
		for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin(); it!=boundary_values.end(); ++it)
		{
			f_u[it->first] += tmp_bnd_2[dof_to_boundary_mapping[it->first]];
		}
	}
}

template <int dim>
void Rhs<dim>::get_rhs_der(SparseMatrix<double>& Dfu,
					const double /*time*/,
					const Vector<double>& u)
{
	const SparsityPattern &sp(M_in.get_sparsity_pattern());
	Dfu.reinit(sp);
	int j;
	double temp;
	for (unsigned int i=0; i<Dfu.m(); ++i)
	{
		typename SparsityPattern::iterator entry = sp.begin(i), end = sp.end(i);
		for (; entry!=end; ++entry){
			j = entry->column();
			temp = 2*M_in(i,j)*fabs(u(j));
			if(dof_to_boundary_mapping[i] != numbers::invalid_dof_index && dof_to_boundary_mapping[j] != numbers::invalid_dof_index)
			{
				temp += 3 * M_bnd(dof_to_boundary_mapping[i],dof_to_boundary_mapping[j]) * u(j)*u(j);
			}
			Dfu.add(i,j,temp);
		}
	}

}

template class Rhs<2>;
template class Rhs<3>;

} /* namespace Data */
} /* namespace Main */
