/*
 * CrankNicolson.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include "CrankNicolson.h"

namespace Main {
namespace TimeIntegration{

template<int dim>
Crank_Nicolson<dim>::Crank_Nicolson(const SparseMatrix<double> &M,
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
	
	tmp.reinit(M.m());
	system_rhs.reinit(M.m());
	solution_update.reinit(M.m());
}


template<int dim>
Crank_Nicolson<dim>::Crank_Nicolson(const SparseMatrix<double>& M,
				const SparseMatrix<double>& A,
				const SparseMatrix<double>& B,
				Data::Rhs<dim>& rhs,
				const double stepsize,
				DoFHandler<dim>& dof_handler,
				const Triangulation<dim>& triangulation,
				const MappingQ<dim>& mapping,
				FE_Q<dim>& fe)
:
	M(M),
	A(A),
	B(B),
	rhs(rhs),
	stepsize(stepsize)
{
	
	
	tmp.reinit(M.m());
	system_rhs.reinit(M.m());
	solution_update.reinit(M.m());
}

template<int dim>
Crank_Nicolson<dim>::~Crank_Nicolson() {
	if(S_inv != 0)
			delete S_inv;
}

template<int dim>
void Crank_Nicolson<dim>::update_matrices(const double time) {
    
	double h1=pow(time+stepsize+1,-1);
	double h=pow(time+1,-1);
    S.reinit(M.get_sparsity_pattern());
	S.copy_from(M);
	S.add(std::pow(stepsize,2)/4,A);
    S.add(h1*stepsize / 2, B); 
	delete S_inv;  
	if(global::config.get<std::string>("Solver_Settings.Newton_Type") == "simple")
		S_inv =  new DataTypes::InverseMatrix<dim>(S,global::config.get<double>("Solver_Settings.Lin_Solver_Tolerance"));
	else
		system_matrix.reinit(M.get_sparsity_pattern());
    
   
}


template <int dim>
void Crank_Nicolson<dim>::integrate_step(Vector<double>& u, Vector<double>& Mv, double time)
{

	
	if(Fn.size()==0)
	{
		Fn.reinit(Mv);
		rhs.get_rhs_vector(Fn, time, u);
		Au.reinit(Mv);
		Bu.reinit(Mv);
		A.vmult(Au,u);
		B.vmult(Bu,u);
	}

    double h1=pow(time+stepsize+1,-1);
	double h=pow(time+1,-1);
	M.vmult(tmp,u);
	tmp.add(-1*std::pow(stepsize,2)/4,Au,h1*stepsize/2,Bu);
	tmp.add((h1-h)*std::pow(stepsize,2)/4,Mv);
	tmp.add(stepsize,Mv,std::pow(stepsize,2)/4,Fn);

	
	Mv.add((h1-h)*stepsize/2,Mv);
	Mv.add(-1.*stepsize/2, Au, stepsize/2, Fn);
	rhs.get_rhs_vector_tx(Fn,time + stepsize);
	tmp.add(std::pow(stepsize,2)/4,Fn);
	Mv.add(h1,Bu,stepsize/2, Fn);
	

	double du_norm = 1., initial_rhs_norm = 0;
	double tol = global::config.get<double>("Solver_Settings.Newton_Tolerance");
	bool first_iteration = true;



	if(global::config.get<std::string>("Solver_Settings.Newton_Type") == "simple")
	{
	
		do
		  {
			assemble_system_simple (u, time + stepsize);
			solution_update=u;
		    if (first_iteration == true)
		    {
		    		std::cout << "   " << S_inv->vmult(u, system_rhs, stepsize * stepsize * stepsize * 10e-3 * global::config.get<double>("Solver_Settings.Lin_Solver_Tolerance") /global::meshwidth * du_norm);
		    }
		    else
		    {
		    		std::cout << '+' << S_inv->vmult(u, system_rhs,  stepsize * stepsize * stepsize * 10e-3 * global::config.get<double>("Solver_Settings.Lin_Solver_Tolerance") /global::meshwidth * du_norm);
		    }
		    solution_update -= u;
		    du_norm = solution_update.l2_norm();
			first_iteration = false;
			//std::cout << std::endl <<  solution_update.l2_norm() * global::meshwidth< << ",  " << stepsize *  stepsize * stepsize *  tol *< std::endl;
		  }

		while (du_norm * global::meshwidth >  stepsize *  stepsize * stepsize  *  tol );
	}
	else
	{
		std::cout << std::endl << "ACHTUNG: Volles Newoten verfahren, funktioniert momentan nicht!" << std::endl;
		do
		  {
			assemble_system (u, time);
		    if (first_iteration == true)
		      initial_rhs_norm = system_rhs.l2_norm();
			SolverControl solver_control (1000, 1e-4*system_rhs.l2_norm());
			SolverGMRES<> gmres (solver_control);
			PreconditionSOR<> preconditioner;
			preconditioner.initialize(system_matrix);
			solution_update = 0;
			gmres.solve (system_matrix, solution_update,
						system_rhs,
						preconditioner);
			u += solution_update;
		    if (first_iteration == true)
		      std::cout << "    " << solver_control.last_step();
		    else
		      std::cout << '+' << solver_control.last_step();
			first_iteration = false;
		  }
		while (system_rhs.l2_norm()* global::meshwidth > stepsize * stepsize * tol * initial_rhs_norm);
	}



	A.vmult(Au,u);
	B.vmult(Bu,u);
	Mv.add(-1.*stepsize/2,Au,-1*h1,Bu);
	rhs.get_rhs_vector_u(tmp,u);

	Mv.add(stepsize/2,tmp);


	Fn += tmp;
}


template <int dim>
void Crank_Nicolson<dim>::assemble_system(Vector<double>& u, double time)
{
	S.vmult(system_rhs,u);
	system_rhs *= -1;
	rhs.get_rhs_vector_u(solution_update,u);
	system_rhs.add(pow(stepsize,2)/4,solution_update);
	system_rhs += tmp;


	rhs.get_rhs_der(system_matrix, time+stepsize, u);
	system_matrix *= -1.*pow(stepsize,2)/4;
	system_matrix.add(1,S);
}

template <int dim>
void Crank_Nicolson<dim>::assemble_system_simple(Vector<double>& u, double /* time */)
{
	rhs.get_rhs_vector_u(system_rhs,u);
	system_rhs *= std::pow(stepsize,2)/4;
	system_rhs += tmp;
}

template class Crank_Nicolson<3>;
template class Crank_Nicolson<2>;

} /* namespace TimeIntegration */
} /* namespace Main */
