

#include "ExactSolution.h"

namespace Main {
namespace Data {

template <int dim>
ExactSolution_u<dim>::ExactSolution_u(const double current_time)
	:
	Function<dim>(1,current_time),
	 sin2PIt(std::sin(2*numbers::PI*current_time)){

}

template <int dim>
ExactSolution_u<dim>::~ExactSolution_u() {
}


template<int dim>
double ExactSolution_u<dim>::value(const Point<dim>   &p,
					const unsigned int /* component */) const
{
	return sin2PIt*p[0]*p[1];


}


template <int dim>
Tensor<1,dim> ExactSolution_u<dim>::gradient (const Point<dim>   &p,
                       const unsigned int) const
{
  Tensor<1,dim> return_value;
  return_value[0] = sin2PIt*p[1];
  return_value[1] = sin2PIt*p[0];
  return return_value;
}

template <int dim>
ExactSolution_v<dim>::ExactSolution_v(const double current_time)
	:
	Function<dim>(1,current_time),
	PI2cos2PIt(2*numbers::PI*std::cos(2*numbers::PI*current_time)){

}

template <int dim>
ExactSolution_v<dim>::~ExactSolution_v() {
}


template<int dim>
double ExactSolution_v<dim>::value(const Point<dim>   &p,
					const unsigned int /* component */) const
{
	return PI2cos2PIt*p[0]*p[1];
	//return 0;
}

template class ExactSolution_u<2>;
template class ExactSolution_v<2>;

template class ExactSolution_u<3>;
template class ExactSolution_v<3>;
} /* namespace Data */
} /* namespace Main */
