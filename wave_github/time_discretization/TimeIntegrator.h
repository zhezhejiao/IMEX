/*
 * TimeIntegrator.h
 *
 *  Created on: 09.07.2024
 *      Author: Li Yaxu
 */

#ifndef TIME_INTEGRATION_TIMEINTEGRATOR_H_
#define TIME_INTEGRATION_TIMEINTEGRATOR_H_



#include <deal.II/lac/vector.h>

using namespace dealii;

namespace Main {


namespace TimeIntegration{




template<int dim>
class TimeIntegrator {
public:


	virtual ~TimeIntegrator() {};


	virtual void integrate_step(Vector<double>& u,
						Vector<double>& v,
						const double time) = 0;

	virtual void update_matrices(const double time) = 0;

};

} /* namespace TimeIntegration */
} /* namespace MaxwellProblem */

#endif /* TIME_INTEGRATION_TIMEINTEGRATOR_H_ */
