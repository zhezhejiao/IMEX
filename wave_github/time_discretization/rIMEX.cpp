/*
 * 	rIMEX.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include "rIMEX.h"

namespace Main {
namespace TimeIntegration{

template<int dim>
rIMEX<dim>::rIMEX(const SparseMatrix<double> &M,
		const SparseMatrix<double> &A,
		const SparseMatrix<double> &B,
		Data::Rhs<dim>& rhs,
		const double stepsize)
		:
		M(M),
		A(A),
		B(B),
		rhs(rhs),
		stepsize(stepsize),
		S_inv(0)
{
	
	
	IP.reinit(M.get_sparsity_pattern());
	IP.copy_from(M);
	IP.add(1.,A);
	
	vnh.reinit(M.m());
	system_rhs.reinit(M.m());
}




template<int dim>
void rIMEX<dim>::update_matrices(const double time) {
  
   
	double h1=pow(time+stepsize+1,-1);
	double h=pow(time+1,-1);
    S.reinit(M.get_sparsity_pattern());
	S.copy_from(M);
	S.add(std::pow(stepsize,2)/4,A);
    S.add((h1+h)*stepsize / 4, B); 
 
    delete S_inv;  
    S_inv = new DataTypes::InverseMatrix<dim>(S, stepsize * stepsize * global::config.get<double>("Solver_Settings.Lin_Solver_Tolerance"));
}


template <int dim>
void rIMEX<dim>::integrate_step(Vector<double>& u,
					Vector<double>& Mv,
					const double time)
{
	double h1=pow(time+stepsize+1,-1);
	double h=pow(time+1,-1);  
	
	if(Fn.size()==0)
	{
		Fn.reinit(Mv);
		rhs.get_rhs_vector(Fn,time,u);
	}
	
	A.vmult(system_rhs,u);
	system_rhs -= Fn;
	system_rhs *= -stepsize/2;
	system_rhs += Mv;

	
	if(global::config.get<std::string>("Solver_Settings.Stopping_Criterion") == "EnergyNorm")
	{
		std::cout << "   " << S_inv->vmult(vnh, system_rhs, IP);
	}
	else
	{
		std::cout << "   " << S_inv->vmult(vnh, system_rhs);
	}
	
	u.add(stepsize,vnh);
	
	Mv *= -1;
	M.vmult(system_rhs,vnh);
	Mv.add(2,system_rhs,-stepsize/2,Fn);
	rhs.get_rhs_vector(Fn,time + stepsize ,u);
	Mv.add(stepsize/2,Fn);

}

template<int dim>
rIMEX<dim>::~rIMEX() {
	delete S_inv;
}

template class rIMEX<3>;
template class rIMEX<2>;

} /* namespace TimeIntegration */
} /* namespace Main */
