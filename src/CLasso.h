/*
 * CLasso.h
 *
 *  Created on: May 4, 2013
 *      Author: root
 */

#ifndef CLASSO_H_
#define CLASSO_H_
#include "Kalman1.h"
#include "CVector.h"
#include "CMatrix.h"
#include "CSparseMatrix.h"

class CLasso {
public:
	CLasso();
	CLasso(CMatrix* mata, CVector* vecx, CVector* vecxhat,double lambda,double rho = 1,double alpha = 1,MKL_INT maxiter = 1000,MKL_INT iverbose = 0);
	MKL_INT Shrinkage(CVector* x, double kappa, CVector* res);
	double Objective(CMatrix* matA,CVector* vecB,double lambda,CVector* vecX,CVector* vecZ);
	virtual ~CLasso();
	MKL_INT DoLasso(CVector* vecx, CVector* vecxhat,double lambda,double rho = 1,double alpha = 1,MKL_INT maxiter = 1000,MKL_INT iverbose = 0);
	MKL_INT BestLambda(CMatrix* mata, CVector* vecx, CVector* vecxhat,double* dlambda,double rho = 1,double alpha = 1,MKL_INT maxiter = 1000,MKL_INT iverbose = 0);
	double dLambda;
	MKL_INT iPathLen;
	double dThresh;
	MKL_INT iMaxIter;
	MKL_INT iVerbose;
	MKL_INT iDim;
};

#endif /* CLASSO_H_ */
