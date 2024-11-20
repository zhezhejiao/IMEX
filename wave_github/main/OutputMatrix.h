

#ifndef OUTPUTMATRIX_H_
#define OUTPUTMATRIX_H_


#include <iostream>
#include <fstream>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

using namespace dealii;
namespace Main{

template<int dim, typename Matrixtype = SparseMatrix<double>>
class OutputMatrix{

public:


	void output_matrix_spy (Matrixtype & matrix, std::stringstream & filename);
	
	void output_matrix (Matrixtype & matrix, std::stringstream & filename) ;

	void output_matrix_matlab (Matrixtype & matrix, std::stringstream & filename);

};
}

#endif /* OUTPUTMATRIX_H_ */
