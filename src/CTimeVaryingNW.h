/*
 * CTimeVaryingNW.h
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#ifndef CTIMEVARYINGNW_H_
#define CTIMEVARYINGNW_H_
#include "Kalman1.h"
#include "CVector.h"
#include "CMatrix.h"
//#include "CMatrix3D.h"

class CTimeVaryingNW {
public:
	CTimeVaryingNW();
	CTimeVaryingNW(MKL_UINT nGene, MKL_INT ntimepts,MKL_INT nobs, double pchange,MKL_INT seed = 0);
	MKL_INT CreateNW(MKL_UINT nGene,MKL_INT nTimePts,MKL_INT nObs,double pChange,double dSprstyLvl, MKL_INT seed = 0,int streamnum = 0);
	virtual ~CTimeVaryingNW();
	CMatrix* matGeneInteractions;
	MKL_UINT nGene;
	MKL_INT nTimePts;
	MKL_INT nObs;
	double perChange;
	double sparsitylvl;//defines the level of sparsity in the first connectivity matrix
	CVector** pVecX;
	CVector** pVecY;
	CMatrix* pMatX;
	CMatrix* pMatY;
};

#endif /* CTIMEVARYINGNW_H_ */
