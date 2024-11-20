

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include <deal.II/base/timer.h>

namespace Main {


namespace global {

extern pt::ptree config;


extern dealii::Timer 	standard_timer, assemble_timer, setup_timer, solve_timer, rhs_assemble_timer,
						ti_timer, output_timer, rest_timer ;


extern double shift;

extern double stepsize;


extern double meshwidth;

}
}

#endif /* GLOBAL_H_ */
