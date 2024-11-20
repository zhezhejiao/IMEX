

#include "OutputMatrix.h"

namespace Main{
template <int dim, typename Matrixtype>
void OutputMatrix<dim, Matrixtype>::output_matrix_spy (Matrixtype & matrix, std::stringstream & filename)
{
	filename << "_spy.vtk";
	std::ofstream out ((filename.str()).c_str());

	for (unsigned int row=0; row<matrix.n(); ++row)
	{
		typename Matrixtype::const_iterator it = matrix.begin(row);
		typename Matrixtype::const_iterator end = matrix.end(row);
		for (; it!=end; ++it)
			if (fabs((*it).value()) > 10e-17)
				out << row << " -" << (*it).column() << "\t" << std::endl;
	}

	out.close();
}


template <int dim, typename Matrixtype>
void OutputMatrix<dim, Matrixtype>::output_matrix (Matrixtype & matrix, std::stringstream & filename)
{
	filename << ".vtk";
	std::ofstream out ((filename.str()).c_str());

	for (unsigned int row=0; row<matrix.n(); ++row)
	{
		typename Matrixtype::const_iterator it = matrix.begin(row);
		typename Matrixtype::const_iterator end = matrix.end(row);
		for (; it!=end; ++it)
			if (fabs((*it).value()) > 10e-17)
				out << row << " -" << (*it).column() << "\t" << (*it).value() << std::endl;
	}

	out.close();
}


template <int dim, typename Matrixtype>
void OutputMatrix<dim, Matrixtype>::output_matrix_matlab (Matrixtype & matrix, std::stringstream & filename)
{
	filename << "_matlab.vtk";
	std::ofstream out ((filename.str()).c_str());

	for (unsigned int row=0; row<matrix.n(); ++row)
	{
		typename Matrixtype::const_iterator it = matrix.begin(row);
		typename Matrixtype::const_iterator end = matrix.end(row);
		for (; it!=end; ++it)
			if (fabs((*it).value()) > 10e-17)
				out << row+1 << "\t" << (*it).column()+1 << "\t" << (*it).value() << std::endl;
	}

	out.close();
}

template class OutputMatrix<2, SparseMatrix<double>>;
template class OutputMatrix<3, SparseMatrix<double>>;
}
