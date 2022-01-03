#include <algorithm>
#include <cassert>
#include <iostream>
#include <thread>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "ProdMatMat.hpp"

namespace {
const int szBlock = 128;
void prodSubBlocks(int iRowBlkA, int iColBlkB, int iColBlkA, int szBlock,
                   const Matrix& A, const Matrix& B, Matrix& C) {

int i, j, k, ib, jb, kb;
#pragma omp parallel private(i, j, k, ib, jb, kb)
{
	#pragma omp for
	for (int jb = iColBlkB; jb < B.nbCols; jb+=szBlock)
	  for (int kb = iColBlkA; kb < A.nbCols; kb+=szBlock)
		 for (int ib = iRowBlkA; ib < A.nbRows; ib+=szBlock)
			for (int j = jb; j < jb + szBlock; j++)
			  for (int k = kb; k < kb + szBlock; k++)
				 for (int i = ib; i < ib + szBlock; i++)
					C(i, j) += A(i, k) * B(k, j);
	int p = omp_get_num_threads();
	std::cout << p << std::endl;
}
}

}  // namespace


Matrix operator*(const Matrix& A, const Matrix& B) {
  Matrix C(A.nbRows, B.nbCols, 0.0);
  prodSubBlocks(0, 0, 0, std::max({A.nbRows, B.nbCols, A.nbCols}), A, B, C);
  return C;
}
