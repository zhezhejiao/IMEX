/*
 * RungeKutta.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include "RungeKutta.h"

namespace Main {
namespace TimeIntegration{

template<int dim>
RungeKutta<dim>::RungeKutta(const SparseMatrix<double>& M,
		const SparseMatrix<double>& A,
		const SparseMatrix<double> &B,
		Data::Rhs<dim>& rhs,
		const double stepsize)
		:
		A(A),
		B(B),
		rhs(rhs),
		stepsize(stepsize),
		M_inv(0),
		u_tmp(A.n()),
		V2(A.n()),
		V3(A.n()),
		V4(A.n()),
		F1(A.n()),
		F2(A.n()),
		F3(A.n()),
		F4(A.n())
{
	M_inv = new DataTypes::InverseMatrix<dim>(M, stepsize * stepsize * stepsize * global::config.get<double>("Solver_Settings.Lin_Solver_Tolerance"),true);

}

template<int dim>
RungeKutta<dim>::RungeKutta(const SparseMatrix<double>& A,
		const SparseMatrix<double>& B,
		Data::Rhs<dim>& rhs,
		const double stepsize,
		DoFHandler<dim>& dof_handler,
		const MappingQ<dim>& mapping,
		FE_Q<dim>& fe,
		double degree)
		:
		A(A),
		B(B),
		rhs(rhs),
		stepsize(stepsize),
		M_inv(0),
		u_tmp(A.n()),
		V2(A.n()),
		V3(A.n()),
		V4(A.n()),
		F1(A.n()),
		F2(A.n()),
		F3(A.n()),
		F4(A.n())
{
	M_inv = new DataTypes::InverseMatrix_ML<dim>(dof_handler, mapping, fe, degree);
}

template<int dim>
void RungeKutta<dim>::update_matrices(const double time) {
  
          if (time==0.8)
		   {
            std::cout << time << std::endl;
          }
		   
   
}



template <int dim>
void RungeKutta<dim>::integrate_step(Vector<double>& u,
					Vector<double>& v,
					const double time)
{


	if(Bv.size()==0)
	{
		Bv.reinit(v);
		Bv2.reinit(v);
		Bv3.reinit(v);
		Bv4.reinit(v);
	}
    
	double h=1/(time+1);
	rhs.get_rhs_vector(V2,time,u);
	V2 *= -1;
	A.vmult_add(V2,u);
	B.vmult(Bv,v);
	V2.add(h,Bv);
	M_inv->vmult(F1,V2);
	V2 = v;
	V2.add(-stepsize/2,F1);

	u_tmp = u;
	u_tmp.add(stepsize/2,v);



	rhs.get_rhs_vector(V3,time + stepsize/2,u_tmp);
	V3 *= -1;
	A.vmult_add(V3,u_tmp);
    B.vmult(Bv2,V2);
	V3.add(h,Bv2);
	M_inv->vmult(F2,V3);
	V3 = v;
	V3.add(-stepsize/2,F2);

	u_tmp = u;
	u_tmp.add(stepsize/2,V2);


	rhs.get_rhs_vector(V4,time + stepsize/2,u_tmp);
	V4 *= -1;
	A.vmult_add(V4,u_tmp);
    B.vmult(Bv3,V3);
    V4.add(h,Bv3);
	M_inv->vmult(F3,V4);
	V4 = v;
	V4.add(-stepsize,F3);

	u_tmp = u;
	u_tmp.add(stepsize,V3);


	u.add(stepsize/6,v,stepsize/3,V2);
	u.add(stepsize/3,V3,stepsize/6,V4);

	v.add(-stepsize/6,F1, -stepsize/3,F2);
	rhs.get_rhs_vector(F1,time + stepsize,u_tmp);
	F1 *= -1;
	A.vmult_add(F1,u_tmp);
	B.vmult(Bv4,V4);
	F1.add(h,Bv4);
	
	std::cout << "   " << M_inv->vmult(F4,F1);
	v.add(-stepsize/3,F3,-stepsize/6,F4);


}

template<int dim>
RungeKutta<dim>::~RungeKutta() {
	delete M_inv;
}

template class RungeKutta<3>;
template class RungeKutta<2>;

} /* namespace TimeIntegration */
} /* namespace Main */
