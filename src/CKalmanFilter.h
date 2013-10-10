/*
 * CKalmanFilter.h
 *
 *  Created on: May 2, 2013
 *      Author: root
 */

#ifndef CKALMANFILTER_H_
#define CKALMANFILTER_H_
#include "Kalman1.h"
#include "CLasso.h"
#include "CVector.h"
#include "CMatrix.h"
#include "CMatrix3D.h"
#include "CSparseMatrix.h"

class CKalmanFilter {
public:
	CKalmanFilter();
	CKalmanFilter(CVector** ppvecx, CVector** ppvecy,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dLambda);
	MKL_INT SmoothEstimate();
	MKL_INT Estimate(CMatrix* pMatX, CMatrix* pMatY,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dlambda);
	MKL_INT EstimateGene(CMatrix* pMatX, CVector* vecY,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dlambda,MKL_INT curgene);
	virtual ~CKalmanFilter();
	CVector** pEstA;
	CVector** pEstAApr;
	CVector** pEstSmtA;
	CMatrix* matSmthEst;
	CMatrix* matEstA;
	CMatrix* matEstAApr;
	CMatrix** matPkk;
	CMatrix** matPkp1k;
	CMatrix** matPkN;
	MKL_INT nGenes;
	MKL_INT nEffMeas;
	MKL_INT nTimePts;
	CVector* pgEstA;
	CVector* pgEstAApr;
	CVector* pgEstSmtA;
};

#endif /* CKALMANFILTER_H_ */
