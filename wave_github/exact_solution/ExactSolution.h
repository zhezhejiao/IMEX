

#ifndef DATA_EXACTSOLUTION_H_
#define DATA_EXACTSOLUTION_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

namespace Main {
namespace Data {

template <int dim>
class ExactSolution_u : public Function<dim> {
public:
	ExactSolution_u(const double current_time = 0);
	virtual ~ExactSolution_u();
	double value (const Point<dim>   &p,
	              const unsigned int  component = 0) const;
	Tensor<1,dim> gradient(const Point<dim> &p,
					const unsigned int component = 0) const;

private:
	const double sin2PIt;
};


template <int dim>
class ExactSolution_v : public Function<dim> {
public:
	ExactSolution_v(const double current_time = 0);
	virtual ~ExactSolution_v();
	double value (const Point<dim>   &p,
	              const unsigned int  component = 0) const;

private:
	const double PI2cos2PIt;
};

} /* namespace Data */
} /* namespace Main */

#endif /* DATA_EXACTSOLUTION_H_ */
