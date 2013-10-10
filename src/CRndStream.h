/*
 * CRndStream.h
 *
 *  Created on: Apr 22, 2013
 *      Author: root
 */

#ifndef CRNDSTREAM_H_
#define CRNDSTREAM_H_

#define UNIFORM_METHOD VSL_RNG_METHOD_UNIFORM_STD

#include "Kalman1.h"
#include "CVector.h"
#include "CMatrix.h"

class CRndStream {
public:
	CRndStream(MKL_INT seed = 777,int streamnum = 0);
	virtual ~CRndStream();
	void vRngGaussian(CMatrix* mat,double mu, double sigma);
	MKL_INT vRngGaussian(CVector* vec,double mu,double sigma);
	void vRngGaussianMV(CMatrix* mat, CVector* mu, CMatrix* sigma);
	void vRngUniform(CMatrix* mat,MKL_INT a,MKL_INT b);
	MKL_INT vRngUniformDiscrete(CMatrix* mat,MKL_INT a,MKL_INT b);
	VSLStreamStatePtr stream;
};

#endif /* CRNDSTREAM_H_ */
