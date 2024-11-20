/*
 * BS_FEM.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include "BS_FEM.h"

namespace Main {

template <int dim>
BS_FEM<dim>::BS_FEM(const Triangulation<dim>& triangulation, const double end_time, const unsigned int deg) :
degree(deg),
mapping(deg),
triangulation(triangulation),
dof_handler(triangulation),
fe (degree),
quadrature(degree+1),
face_quadrature(degree+1),
end_time(end_time)
{
	global::meshwidth = GridTools::maximal_cell_diameter(triangulation);
	setup_and_assemble_system();
	if(!global::config.get<bool>("Ref_Solution.Compute_Reference_Solution") && global::config.get<bool>("Error.Error") && global::config.get<std::string>("Error.Error_against") == "Reference")
	{
		
		std::string str_tria_path = "ref_solution/" + global::config.get<std::string>("Ref_Solution.Name") + "/tria.bin";

		std::ifstream triasrc(str_tria_path.c_str());

        boost::archive::binary_iarchive ar{triasrc};

		triangulation_ref.load(ar,0);
		const SphericalManifold<dim> sphere;
		triangulation_ref.set_all_manifold_ids_on_boundary (0);
		triangulation_ref.set_manifold(0, sphere);
		dof_handler_ref.initialize(triangulation_ref, dof_handler.get_fe());
	    DynamicSparsityPattern dsp(dof_handler_ref.n_dofs());
	    DoFTools::make_sparsity_pattern (dof_handler_ref, dsp);
	    sparsity_pattern_ref.copy_from(dsp);

		std::stringstream filename;
		filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_u.vtk";
		std::ifstream in ((filename.str()).c_str());
		u_ref.block_read(in);
		in.close();
		filename.str("");
		filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_v.vtk";
		in.open((filename.str()).c_str());
		v_ref.block_read(in);
		in.close();

		
		filename.str("");
		filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_M.vtk";
		in.open((filename.str()).c_str());
		M_ref.reinit(sparsity_pattern_ref);
		M_ref.block_read(in);
		in.close();
		filename.str("");
		filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_A.vtk";
		in.open((filename.str()).c_str());
		A_ref.reinit(sparsity_pattern_ref);
		A_ref.block_read(in);
		in.close();


	}
}

template <int dim>
BS_FEM<dim>::~BS_FEM ()
{
	dof_handler.clear ();
}


template <int dim>
void BS_FEM<dim>::run	(const double total_timesteps)
{


	double timestep = end_time / (double) total_timesteps;
	double time = 0;
	unsigned int timestep_number = 0;

	Data::Rhs<dim> rhs = Data::Rhs<dim>(mapping,dof_handler,dof_to_boundary_mapping,M_bnd,M_in);

	TimeIntegration::TimeIntegrator<dim> *integrator(0);


	global::ti_timer.restart();


	std::string ti = global::config.get<std::string>("Time_Integration.Time_Integrator");


	if (global::config.get<std::string>("Error.Error_against") == "Exact")
	{
		Data::ExactSolution_v<dim> iv_v = Data::ExactSolution_v<dim>();
		VectorTools::interpolate(mapping,dof_handler,iv_v,u);

	}
	else
	{
		Data::InitialValue_v<dim> iv_v = Data::InitialValue_v<dim>();
		VectorTools::interpolate(mapping,dof_handler,iv_v,u);
	}

	if(ti == "CrankNicolson" || ti == "vIMEX"||ti == "rIMEX")
	{

		if(ti == "CrankNicolson")
			integrator = new TimeIntegration::Crank_Nicolson<dim>(M,A,B,rhs,timestep);
		else if(ti == "IMEX")
			integrator = new TimeIntegration::vIMEX<dim>(M,A,B,rhs,timestep);
		else
			integrator = new TimeIntegration::rIMEX<dim>(M,A,B,rhs,timestep);
		
		M.vmult(v,u);
		if(ti == "CrankNicolson")
			std::cout << std::endl << "Starting time integration (Crank-Nicolson) with "
					<< total_timesteps << " time steps." << std::endl;
		else if(ti == "vIMEX")
			std::cout << std::endl << "Starting time integration (vIMEX) with "
					<< total_timesteps << " time steps." << std::endl;
		else
			std::cout << std::endl << "Starting time integration (rIMEX) with "
					<< total_timesteps << " time steps." << std::endl;		
	}
	else
	{
		v.swap(u);
		
		
		if(ti == "RungeKutta")
		{
			if(global::config.get<bool>("Solver_Settings.Mass_Lumping") == false)
		
				integrator = new TimeIntegration::RungeKutta<dim>(M,A,B,rhs,timestep);
			else
				integrator = new TimeIntegration::RungeKutta<dim>(A,B,rhs,timestep,dof_handler, mapping, fe, degree);
			std::cout << "Starting time integration (Runge-Kutta) with "
											<< total_timesteps << " time steps." << std::endl;
		}
		
	}
	if (global::config.get<std::string>("Error.Error_against") == "Exact")
	{
		Data::ExactSolution_u<dim> iv_u = Data::ExactSolution_u<dim>();
		VectorTools::interpolate(mapping,dof_handler,iv_u,u);

	}
	else
	{
		Data::InitialValue_u<dim> iv_u = Data::InitialValue_u<dim>();
		VectorTools::interpolate(mapping,dof_handler,iv_u,u);
	}

	global::ti_timer.stop();
	if(global::config.get<bool>("Output.Solution"))
		output_results (0);
	global::ti_timer.start();
	
	for(timestep_number=0, time=0; timestep_number < total_timesteps; time += timestep, ++timestep_number)
	{
		if (timestep_number % 10 == 0)
		{
			std::cout << '\r' << timestep_number << "    Time: "<< time << "             ";
			std::cout.flush();
		}
		

		integrator->update_matrices(time);
		integrator->integrate_step(u,v,time);

		global::ti_timer.stop();
		if(global::config.get<bool>("Output.Solution") and (timestep_number + 1) % global::config.get<int>("Output.Output_timestep_skip") == 0)
		{
			output_results (timestep_number+1);
		}
		global::ti_timer.start();
	} 

	std::cout << '\r' << timestep_number << "    Time: "<< time << "             "<< std::endl;
	global::ti_timer.stop();
	std::cout << "Time integration finished after " << global::ti_timer.cpu_time() << "s CPU-time." << std::endl << std::endl ;
	std::cout.flush();

	delete integrator;

	if (global::config.get<bool>("Error.Error") && global::config.get<bool>("Ref_Solution.Compute_Reference_Solution")== false)
	{

		std::ofstream out_err ("error_std.txt", std::ios_base::app);
		double error;
		if (global::config.get<std::string>("Error.Error_against") == "Exact")
			error = compute_error(end_time);
		else
			error = compute_error_ref(end_time);
		out_err << degree << "\t" << GridTools::maximal_cell_diameter(triangulation) << "\t" << timestep << "\t" << error << "\t" << global::ti_timer.cpu_time() << std::endl;
		out_err.close();
	}
	else if(global::config.get<bool>("Ref_Solution.Compute_Reference_Solution"))
		output_ref_sol ();
}



template <int dim>
void BS_FEM<dim>::setup_and_assemble_system	()
{
	std::cout << "Number of active cells:       "
						<< triangulation.n_active_cells()
						<< std::endl;

	
	dof_handler.distribute_dofs(fe);

	std::cout << "Number of degrees of freedom: "
			<< dof_handler.n_dofs()
			<< std::endl << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    DoFTools::map_dof_to_boundary_indices ( dof_handler, dof_to_boundary_mapping);
    dsp.reinit(dof_handler.n_boundary_dofs(),dof_handler.n_boundary_dofs());
    DoFTools::make_boundary_sparsity_pattern (dof_handler, dof_to_boundary_mapping, dsp);
    sparsity_pattern_bnd.copy_from(dsp);
    M.reinit(sparsity_pattern);
    M_in.reinit(sparsity_pattern);
    M_bnd.reinit(sparsity_pattern_bnd);
    A.reinit(sparsity_pattern);
    B.reinit(sparsity_pattern);
    u.reinit(dof_handler.n_dofs());
    v.reinit(dof_handler.n_dofs());

    
    Assembling::LocalIntegrator<dim> li;

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    FullMatrix<double> cell_matrix_A (dofs_per_cell, dofs_per_cell),
    				   cell_matrix_B (dofs_per_cell, dofs_per_cell),
					   cell_matrix_M_in (dofs_per_cell, dofs_per_cell),
					   cell_matrix_M_bnd (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim>  fe_values (mapping, fe, quadrature, update_values|update_gradients|update_quadrature_points|update_JxW_values);
    FEFaceValues<dim> fe_face_values (mapping,fe, face_quadrature, update_values|update_gradients|update_quadrature_points|update_normal_vectors|update_JxW_values);


    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
    	cell->get_dof_indices (local_dof_indices);
        cell_matrix_A = 0;
        cell_matrix_B = 0;
        cell_matrix_M_in = 0;
        cell_matrix_M_bnd = 0;
        fe_values.reinit (cell);

        li.assemble_cell_mass(fe_values, cell_matrix_M_in);
        li.assemble_cell_stiffness(fe_values,cell_matrix_A);
        li.assemble_cell_damping(fe_values, cell_matrix_B);

        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          if (cell->face(face_number)->at_boundary())
            {
              fe_face_values.reinit (cell, face_number);

              li.assemble_bnd_mass(fe_face_values,cell_matrix_M_bnd);
              li.assemble_bnd_stiffness(fe_face_values,cell_matrix_A);
            }
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              A.add(local_dof_indices[i], local_dof_indices[j], cell_matrix_A(i,j));
              B.add(local_dof_indices[i], local_dof_indices[j], cell_matrix_B(i,j));
              M.add(local_dof_indices[i], local_dof_indices[j], cell_matrix_M_in(i,j));
              M_in.add(local_dof_indices[i], local_dof_indices[j], cell_matrix_M_in(i,j));
              if(cell_matrix_M_bnd(i,j)!=0)
              {
            	M.add (local_dof_indices[i], local_dof_indices[j], cell_matrix_M_bnd(i,j));
              	M_bnd.add(dof_to_boundary_mapping[local_dof_indices[i]], dof_to_boundary_mapping[local_dof_indices[j]], cell_matrix_M_bnd(i,j));
              }
            }
          }
      }

	if (global::config.get<bool>("Output.Matrices"))
	{
		OutputMatrix<dim, SparseMatrix<double>> output;
		std::stringstream filename;
		filename << "temp/mass_matrix";
	
		output.output_matrix_matlab(M,filename);

		filename.str("");
		filename << "temp/stiffness_matrix";

		output.output_matrix_matlab(A, filename);

		filename.str("");
		filename << "temp/damping_matrix";

		output.output_matrix_matlab(B, filename);
	}
}

template <int dim>
double BS_FEM<dim>::compute_error(double time)
{
	double total_error_u (0);
	double total_error_v (0);

	Data::ExactSolution_u<dim> exsol_u(time);
	Data::ExactSolution_v<dim> exsol_v(time);


	const QGauss<dim> quadrature(2*degree);
	const QGauss<dim-1> face_quadrature(2*degree);


	Vector<double> v(this->v.size());
	if(global::config.get<std::string>("Time_Integration.Time_Integrator") == "CrankNicolson" ||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "vIMEX"||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "rIMEX")
	{
	
		SolverControl solver_control (1000, 1e-30);
		SolverCG<> cg (solver_control);
		PreconditionSSOR<SparseMatrix<double>> preconditioner;
		preconditioner.initialize(M);
		cg.solve (M, v, this->v,preconditioner);
	}
	else
	{
		v.swap(this->v);
	}




    const unsigned int n_q_points = quadrature.size();
    const unsigned int n_face_q_points = face_quadrature.size();
    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim>  fe_values (mapping, fe, quadrature, update_values|update_gradients|update_quadrature_points|update_JxW_values);
    FEFaceValues<dim> fe_face_values (mapping,fe, face_quadrature, update_values|update_gradients|update_quadrature_points|update_normal_vectors|update_JxW_values);

 	std::vector<double>	values_u(n_q_points), values_u_bnd(n_face_q_points), values_v(n_q_points), values_v_bnd(n_face_q_points);
 	std::vector< Tensor< 1, dim, double > >  	gradients_u(n_q_points), gradients_u_bnd(n_face_q_points);

	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	  {
		cell->get_dof_indices (local_dof_indices);
		fe_values.reinit (cell);
		fe_values.get_function_values(u,values_u);
		fe_values.get_function_gradients(u,gradients_u);
		fe_values.get_function_values(v,values_v);
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{
		total_error_u += ((values_u[q_point] - exsol_u.value(fe_values.quadrature_point(q_point)))*(values_u[q_point] - exsol_u.value(fe_values.quadrature_point(q_point)))*fe_values.JxW(q_point)) +
						+((gradients_u[q_point] - exsol_u.gradient(fe_values.quadrature_point(q_point)))*(gradients_u[q_point] - exsol_u.gradient(fe_values.quadrature_point(q_point)))*fe_values.JxW(q_point));

		total_error_v += ((values_v[q_point] - exsol_v.value(fe_values.quadrature_point(q_point)))*(values_v[q_point] - exsol_v.value(fe_values.quadrature_point(q_point)))*fe_values.JxW(q_point));
		}
		for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
		  if (cell->face(face_number)->at_boundary())
			{
			  fe_face_values.reinit (cell, face_number);
			  fe_face_values.get_function_values(u,values_u_bnd);
			  fe_face_values.get_function_gradients(u,gradients_u_bnd);
			  fe_face_values.get_function_values(v,values_v_bnd);
			  const std::vector<Tensor <1,dim>> &normals = fe_face_values.get_normal_vectors ();
			  for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
				{
				  total_error_u += ( ( (gradients_u_bnd[q_point] - exsol_u.gradient(fe_face_values.quadrature_point(q_point))) -(gradients_u_bnd[q_point] - exsol_u.gradient(fe_face_values.quadrature_point(q_point)))*normals[q_point] * normals[q_point] ) *
						  ( (gradients_u_bnd[q_point] - exsol_u.gradient(fe_face_values.quadrature_point(q_point))) -(gradients_u_bnd[q_point] - exsol_u.gradient(fe_face_values.quadrature_point(q_point)))*normals[q_point] * normals[q_point] ) *
						 fe_face_values.JxW(q_point)) +
						 ( (values_u_bnd[q_point] - exsol_u.value(fe_face_values.quadrature_point(q_point))) *
						   (values_u_bnd[q_point] - exsol_u.value(fe_face_values.quadrature_point(q_point))) *
								 fe_face_values.JxW(q_point));

				  total_error_v += ( (values_v_bnd[q_point] - exsol_v.value(fe_face_values.quadrature_point(q_point))) *
						   (values_v_bnd[q_point] - exsol_v.value(fe_face_values.quadrature_point(q_point))) *
						   fe_face_values.JxW(q_point));
				}
			}
	  }
	return std::sqrt(total_error_u) + std::sqrt(total_error_v);
}



template <int dim>
double BS_FEM<dim>::compute_error_ref(double /* time */)
{
	double total_error_u (0);
	double total_error_v (0);

	Vector<double> v(this->v.size());
	if(global::config.get<std::string>("Time_Integration.Time_Integrator") == "CrankNicolson" ||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "vIMEX"||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "rIMEX")
	{
		
		SolverControl solver_control (1000, 1e-30);
		SolverCG<> cg (solver_control);
		PreconditionSSOR<SparseMatrix<double>> preconditioner;
		preconditioner.initialize(M);
		cg.solve (M, v, this->v,preconditioner);
	}
	else
	{
		v.swap(this->v);
	}

	Vector<double> u_fine(dof_handler_ref.n_dofs()), v_fine(dof_handler_ref.n_dofs());
	VectorTools::interpolate_to_different_mesh(dof_handler,
												u,
												dof_handler_ref,
												u_fine);

	VectorTools::interpolate_to_different_mesh(dof_handler,
												v,
												dof_handler_ref,
												v_fine);



	u_fine -= u_ref;
	v_fine -= v_ref;

	total_error_v = M_ref.matrix_norm_square(v_fine);
	total_error_u = M_ref.matrix_norm_square(u_fine) + A_ref.matrix_norm_square(u_fine);

	return std::sqrt(total_error_u) + std::sqrt(total_error_v);
}


template <int dim>
void BS_FEM<dim>::output_results (const unsigned int timestep_number) const
{
	std::stringstream filename;
	filename << "results/sol-" << timestep_number << ".vtk";



	std::ofstream vtk_output ((filename.str()).c_str());

	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (u, "solution");
	data_out.build_patches ();

	data_out.write_vtk(vtk_output);
}

template <int dim>
void BS_FEM<dim>::output_ref_sol ()
{
	std::string path = "ref_solution/" + global::config.get<std::string>("Ref_Solution.Name");
	const int dir_err = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1 == dir_err)
	{
		printf("Error creating directory! Maybe there is already a directory with this name?\n");
		exit(1);
	}

	
	std::stringstream filename;
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_u.vtk";
	std::ofstream out ((filename.str()).c_str());
	u.block_write(out);

	Vector<double> w(v.size());
	if(global::config.get<std::string>("Time_Integration.Time_Integrator") == "CrankNicolson" ||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "vIMEX"||
			global::config.get<std::string>("Time_Integration.Time_Integrator") == "rIMEX")
	{
		
		SolverControl solver_control (1000, 1e-30);
		SolverCG<> cg (solver_control);
		PreconditionSSOR<SparseMatrix<double>> preconditioner;
		preconditioner.initialize(M);
		cg.solve (M, w, v,preconditioner);
	}
	else
	{
		w.swap(v);
	}
	out.close();
	filename.str("");
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_v.vtk";
	out.open((filename.str()).c_str());
	w.block_write(out);

	out.close();
	filename.str("");
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_M.vtk";
	out.open((filename.str()).c_str());
	M.block_write(out);

	out.close();
	filename.str("");
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/solution_A.vtk";
	out.open((filename.str()).c_str());
	A.block_write(out);
	out.close();


	
	filename.str("");
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/config.ini";
	std::ifstream src("config.ini", std::ios::binary);
	std::ofstream dst((filename.str()).c_str());
	dst << src.rdbuf();


	filename.str("");
	filename << "ref_solution/" << global::config.get<std::string>("Ref_Solution.Name") << "/tria.bin";
	std::ofstream triadst((filename.str()).c_str());
    boost::archive::binary_oarchive ar{triadst};
	triangulation.save(ar,0);

}




template class BS_FEM<2>;
template class BS_FEM<3>;

} /* namespace Main */
