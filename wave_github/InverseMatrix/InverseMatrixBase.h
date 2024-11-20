

#ifndef DATA_TYPES_INVERSEMATRIXBASE_H_
#define DATA_TYPES_INVERSEMATRIXBASE_H_

#include <deal.II/lac/vector.h>

using namespace dealii;

namespace Main {
namespace DataTypes{

template<int dim>
class InverseMatrixBase {
public:

	
	virtual ~InverseMatrixBase() {};


	virtual unsigned int vmult(Vector<double> &dst, const Vector<double> &src) = 0;
	virtual unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const SparseMatrix<double> &IP) = 0;
	virtual unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const double tol) =0;

};

} /* namespace DataTypes */
} /* namespace Main */

#endif /* DATA_TYPES_INVERSEMATRIXBASE_H_ */
