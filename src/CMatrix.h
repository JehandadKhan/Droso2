/*
 * CMatrix.h
 *
 *  Created on: Apr 19, 2013
 *      Author: root
 */

#ifndef CMATRIX_H_
#define CMATRIX_H_
#include <mkl.h>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include "CVector.h"


class CMatrix {
public:
//


	CMatrix();
	CMatrix(int m,int n);
	virtual ~CMatrix();
	MKL_INT Create(MKL_INT m,MKL_INT n);
	MKL_INT MultiplyV(CVector* res,CVector* opVector,bool TransposeMat = false,double alpha = 1,double beta = 0);//implements the blas level 2 op ?gemv
	void Set(int m,int n,double val);
	void Copy(CMatrix* dest,bool bTrans = false);
	void Eye(void);
	void Zero(void);
	void Print();
	void BinThreshold(double thr);
	void Threshold(double thr);
	void ElemMul(CMatrix* opb, CMatrix* res);
	MKL_INT FillMat();
	MKL_INT Scale(double alpha,bool bTranspose = false);
	MKL_INT Add(CMatrix* opA,CMatrix* opB,double alpha = 1,double beta = 1);
	MKL_INT IsZero(CVector* idx);
	MKL_INT IsNonZero(CVector* idx);
	MKL_INT GetVec(CVector* pres,MKL_INT num,bool row);//returns the num row / col of the matrix depending on var row
	MKL_INT SetVec(CVector* ptrVec,MKL_INT num,bool row);
	MKL_INT SetVec(CMatrix* pvec,MKL_INT num,bool row);
	MKL_INT VecProd(CVector* opA,CVector* opB);
	MKL_INT EyeSubMat(double alpha = -1);
	MKL_INT InvMat(CMatrix* matRes);
	MKL_INT MMMul(CMatrix* matB, CMatrix* res, double alpha = 1,double beta = 0,bool bTransA = false, bool bTransB = false);
	MKL_INT Sprand(double sparsity);
	MKL_INT Vec(CVector* vec);
	MKL_INT UnVec(CVector* pSrcVec);
	MKL_INT Release();
	MKL_INT Transpose();
	static MKL_INT SaveMat(char* filename,int nVars,...);
	static MKL_INT SaveMat3d(char* filename,int nVars, ...);
	static MKL_INT ReadMat3d(char* filename,char* varname,int dim,CMatrix* retmat);
	static MKL_INT ReadMat(char* filename,char* varname,CMatrix* retmat);
	static MKL_INT DelMat(char* filename,char* varname);

	int rows;
	int cols;
	double* mat;
};

#endif /* CMATRIX_H_ */
