/*
 * CMatrix3D.h
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#ifndef CMATRIX3D_H_
#define CMATRIX3D_H_
#include "Kalman1.h"
#include "CVector.h"
#include "CMatrix.h"
#include "CMatrix3D.h"

class CMatrix3D {
public:
	CMatrix3D(MKL_INT r,MKL_INT c, MKL_INT z);
	void Print(void);
	void Set(CMatrix* mat,MKL_INT z);
	void Eye(void);
	virtual ~CMatrix3D();
	CMatrix* pMats;
	MKL_INT depth;
};

#endif /* CMATRIX3D_H_ */
