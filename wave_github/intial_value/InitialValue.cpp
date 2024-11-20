

#include "InitialValue.h"

namespace Main {
namespace Data {

template <int dim>
InitialValue_u<dim>::InitialValue_u(const double /*current_time*/)
	:
	Function<dim>(1){

}

template <int dim>
InitialValue_u<dim>::~InitialValue_u() {
}


template<int dim>
double InitialValue_u<dim>::value(const Point<dim>&   /*p*/,
					const unsigned int /* component */) const
{
	//return 10*std::cos(10*p[0]*p[1]);
        return 0;

}

template <int dim>
InitialValue_v<dim>::InitialValue_v(const double /*current_time*/)
	:
	Function<dim>(1){

}

template <int dim>
InitialValue_v<dim>::~InitialValue_v() {
}


template<int dim>
double InitialValue_v<dim>::value(const Point<dim>   &p,
					const unsigned int /* component */) const
{
	return 100*std::sin(p[0]*p[1]);
	//return 0;
}

template class InitialValue_u<2>;
template class InitialValue_v<2>;

template class InitialValue_u<3>;
template class InitialValue_v<3>;
} /* namespace Data */
} /* namespace Main */
