/*
 * CSparseMatrix.h
 *
 *  Created on: May 2, 2013
 *      Author: root
 */

#ifndef CSPARSEMATRIX_H_
#define CSPARSEMATRIX_H_
#include <mkl.h>
#include "CVector.h"
#include "CMatrix.h"

class CSparseMatrix {
public:
	CSparseMatrix();
	CSparseMatrix(MKL_INT r,MKL_INT c);
	MKL_INT KronEyeVecVec(MKL_INT idimeye,CVector* vec,CVector* pVecOp,CVector* pRes,bool bTransMat = false,double alpha = 1,double beta = 0);//Kron Prod between a general eye matrix and dense vector resulting in a sparse matrix
	MKL_INT KronEyeMatVec(MKL_INT idimeye,CMatrix* mat,CVector* pVecOp,CVector* pRes,bool bTransMat = false,double alpha = 1,double beta = 0);//Kron Prod between a general dense matrix and dense vector resulting in a sparse matrix further multiplied with a vector
	MKL_INT MulMatVec(CVector* opB,CVector* vecRes, double alpha = 1,double  beta = 0,bool bTransA = false);
	void DenseMatrix(CMatrix* pMat);
	virtual ~CSparseMatrix();
	void Print();
	MKL_INT Eye(double alpha = 1);
	MKL_INT Inv(bool isDiag = false);
	double* values;
	MKL_INT* columns;
	MKL_INT* rowIndex;
	MKL_INT  Rows;
	MKL_INT  Cols;
	MKL_INT  numNonZero;
	MKL_INT  numZero;
};

#endif /* CSPARSEMATRIX_H_ */
