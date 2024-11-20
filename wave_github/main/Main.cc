/*
 * BS_FEM.cpp
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#include <deal.II/base/timer.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>

#include "global.h"
#include "BS_FEM.h"


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/foreach.hpp>
namespace pt = boost::property_tree;

pt::ptree Main::global::config;


dealii::Timer 	Main::global::standard_timer,Main::global::assemble_timer, Main::global::setup_timer,
			Main::global::solve_timer, Main::global::rhs_assemble_timer,
			Main::global::ti_timer, Main::global::output_timer,
			Main::global::rest_timer;

double Main::global::shift;
double Main::global::stepsize;
double Main::global::meshwidth;



static void run();

int main ()
{
	try
	{
		using namespace dealii;
		using namespace Main;

		pt::ini_parser::read_ini("config.ini", global::config);

	
        const unsigned int dim =2;
		Assert (dim>=2, ExcInternalError());

		std::cout << std::endl << "------------------------------------ vIMEX and rIMEX-------------------------------------" << std::endl << std::endl;
		run();
		std::rename("error_std.txt", "results/1vIMEX.txt");

		//global::config.put("Time_Integration.Time_Integrator","rIMEX");
		//run();
		//std::rename("error_std.txt", "results/rIMEX.txt");

		//global::config.put("Time_Integration.Time_Integrator","CrankNicolson");
		//run();
		//std::rename("error_std.txt", "results/CN.txt");

		//global::config.put("Time_Integration.Time_Integrator","RungeKutta");
		//global::config.put("Time_Integration.Timesteps","6400,12800,25600,51200,102400");
		//run();
		//std::rename("error_std.txt", "results/RKM.txt");


	}
	catch (std::exception &exc)
    {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
				<< exc.what() << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
    }
	catch (...)
	{
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	};
	return 0;
}

static void run()
{
	using namespace dealii;
	using namespace Main;

	const unsigned int dim = 2;

	const double end_time = global::config.get<double>("Time_Integration.End_Time");

	const unsigned int starting_lvl = global::config.get<int>("Space_Discretisation.Starting_Lvl");
	const unsigned int n_runs = global::config.get<int>("Space_Discretisation.N_Runs");
	const unsigned int refinement_factor = global::config.get<int>("Space_Discretisation.Refinement_Factor");

	if (global::config.get<bool>("Ref_Solution.Compute_Reference_Solution"))
	{
		std::cout << "Computing reference solution" << std::endl;
		int degree = global::config.get<int>("Ref_Solution.Ref_Degree");

		const SphericalManifold<dim> sphere;
		Triangulation<dim> triangulation(Triangulation<dim>::limit_level_difference_at_vertices);
		GridGenerator::hyper_ball (triangulation);
		triangulation.set_all_manifold_ids_on_boundary (0);
		triangulation.set_manifold(0, sphere);
		triangulation.refine_global(global::config.get<int>("Ref_Solution.Ref_Lvl"));

		BS_FEM<dim> bs_fem(triangulation, end_time, degree);
		bs_fem.run(global::config.get<int>("Ref_Solution.Ref_Timestep"));
	}
	else
	{
		if (global::config.get<bool>("Error.Error"))
		{
			std::ofstream out_err ("error_std.txt");
			out_err << "degree\th_max\ttimestep\terror\tti_time" <<  std::endl;
			out_err.close();
		}


		std::stringstream ss_d(global::config.get<std::string>("Space_Discretisation.Pol_Degree"));
		std::string d;
		
		while (std::getline(ss_d, d, ','))
		{
			unsigned int degree = std::atoi(d.c_str());
			std::cout << std::endl << "---------polynomial_degree: " << degree << "----------------" << std::endl << std::endl;

			const SphericalManifold<dim> sphere;
			Triangulation<dim> triangulation(Triangulation<dim>::limit_level_difference_at_vertices);
			GridGenerator::hyper_ball (triangulation);
			triangulation.set_all_manifold_ids_on_boundary (0);
			triangulation.set_manifold(0, sphere);
			triangulation.refine_global(starting_lvl);

		

			for (unsigned int n_run = 0; n_run <  n_runs; ++n_run)
			{
				std::cout << std::endl << "----------refinement_level: " << starting_lvl + n_run*refinement_factor << "--------" << std::endl << std::endl;


				

				BS_FEM<dim> bs_fem(triangulation, end_time, degree);

				
				std::stringstream ss_t(global::config.get<std::string>("Time_Integration.Timesteps"));
				std::string k;
				while (std::getline(ss_t, k, ','))
				{
					unsigned int total_timesteps = std::atoi(k.c_str());
					bs_fem.run(total_timesteps);
				} 
				triangulation.refine_global(refinement_factor);
			}
		}
	}
}
