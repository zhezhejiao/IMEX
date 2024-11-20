

#include "rhsfunctions.h"

namespace Main {
namespace Data{

template<int dim>
Rhs_func_int<dim>::Rhs_func_int(const double current_time) : Function<dim>(1,current_time){


}

template<int dim>
Rhs_func_int<dim>::~Rhs_func_int() {
	
}

template<int dim>
double Rhs_func_int<dim>::value(const Point<dim>   &p,
					const unsigned int /* component */) const
{

	return p[0]*p[1];
	//return 0;
}





template<int dim>
Rhs_func_bnd<dim>::Rhs_func_bnd(const double current_time) : Function<dim>(1,current_time){
	// TODO Auto-generated constructor stub

}

template<int dim>
Rhs_func_bnd<dim>::~Rhs_func_bnd() {
	// TODO Auto-generated destructor stub
}

template<int dim>
double Rhs_func_bnd<dim>::value(const Point<dim>   &p,
					const unsigned int /* component */) const
{
	//return std::sin(2*numbers::PI*this->get_time())*p[1]*p[2]*(-4*numbers::PI*numbers::PI + 6 - abs_square(std::sin(2*numbers::PI*this->get_time())*p[1]*p[2]));
	return p[0]*p[1];
	//return 0;
}

template class Rhs_func_bnd<3>;
template class Rhs_func_bnd<2>;
template class Rhs_func_int<3>;
template class Rhs_func_int<2>;

}/* namespace Data */
} /* namespace Main */
