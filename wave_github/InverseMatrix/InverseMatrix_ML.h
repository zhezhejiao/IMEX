

#ifndef DATA_TYPES_INVERSEMATRIX_ML_H_
#define DATA_TYPES_INVERSEMATRIX_ML_H_

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/diagonal_matrix.h>

#include "../main/global.h"
#include "../local_integration/LocalIntegrators.h"
#include "InverseMatrixBase.h"


using namespace dealii;

namespace Main {
namespace DataTypes{


template <int dim>
class InverseMatrix_ML : public InverseMatrixBase<dim>{
public:
	InverseMatrix_ML(DoFHandler<dim>& dof_handler,
					const MappingQ<dim>& mapping,
					FE_Q<dim>& fe,
					const unsigned int degree);
	InverseMatrix_ML();
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const SparseMatrix<double>& IP);
	unsigned int vmult(Vector<double> &dst, const Vector<double> &src, const double);

private:
	DiagonalMatrix<Vector<double>> M_ml_inv;
};

} /* namespace DataTypes */
} /* namespace Main */

#endif /* DATA_TYPES_INVERSEMATRIX_ML_H_ */
