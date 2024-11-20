
#ifndef DATA_RHSFUNCTIONS_H_
#define DATA_RHSFUNCTIONS_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

namespace Main {
namespace Data {

template <int dim>
class Rhs_func_int : public Function<dim>
{
public:
	Rhs_func_int(const double current_time = 0);
	virtual ~Rhs_func_int();
	double value (const Point<dim>   &p,
	              const unsigned int  component = 0) const;
};

template <int dim>
class Rhs_func_bnd : public Function<dim>
{
public:
	Rhs_func_bnd(const double current_time = 0);
	virtual ~Rhs_func_bnd();
	double value (const Point<dim>   &p,
	          	  const unsigned int  component = 0) const;
};

} /* namespace Data */
} /* namespace Main */

#endif /* DATA_RHSFUNCTIONS_H_ */
