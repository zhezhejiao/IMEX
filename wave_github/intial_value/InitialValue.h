

#ifndef DATA_INITIALVALUE_H_
#define DATA_INITIALVALUE_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

namespace Main {
namespace Data {

template <int dim>
class InitialValue_u : public Function<dim> {
public:
	InitialValue_u(const double current_time = 0);
	virtual ~InitialValue_u();
	double value (const Point<dim>   &p,
	              const unsigned int  component = 0) const;
};


template <int dim>
class InitialValue_v : public Function<dim> {
public:
	InitialValue_v(const double current_time = 0);
	virtual ~InitialValue_v();
	double value (const Point<dim>   &p,
	              const unsigned int  component = 0) const;
};

} /* namespace Data */
} /* namespace Main */

#endif /* DATA_INITIALVALUE_H_ */
