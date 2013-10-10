/*
 * CVector.h
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#ifndef CVECTOR_H_
#define CVECTOR_H_
#include "Kalman1.h"
#include "mmio.h"

class CVector {
public:
	CVector();
	CVector(int l,bool bCol);
	CVector(double* ptr,int l, bool bCol);
	void Copy(CVector* pSrcVec);
	void Print(bool bPrintAsRow = false);
	MKL_INT Create(MKL_INT len,bool bcol);
	void Destroy();
	MKL_INT ScaleI(double alpha);//Inplace scale by alpha
	MKL_INT FillVec();
	double DotProd(CVector* opB);
	MKL_INT SubtractI(CVector* opB,double alpha = 1,double beta = -1);
	MKL_INT AddI(CVector* opB,double alpha = 1,double beta = 1);
	MKL_INT ReadMMFile(char* strFileName);
	MKL_INT WriteMMFile(char* strFileName);
	MKL_INT Ones();
	MKL_INT Max(double d,CVector* vecRes);
	MKL_INT IsZero();
	MKL_INT Cardinality();
	double Sum();
	double Norm(MKL_INT n);
	MKL_INT Threshold(double alpha);
	virtual ~CVector();
	//Public Variables
	bool bColVector;
	int length;
	double* data;
};

#endif /* CVECTOR_H_ */
