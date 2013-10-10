/*
 * CRndStream.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: root
 */

#include "CRndStream.h"

CRndStream::CRndStream(MKL_INT seed,int streamnum) {
	// TODO We need to make these streams parallelizable
//	vslNewStream(&this->stream,VSL_BRNG_MT19937,seed);
	vslNewStream(&this->stream,VSL_BRNG_MT2203 + streamnum,seed);

}

CRndStream::~CRndStream() {
	vslDeleteStream(&this->stream);
}

void CRndStream::vRngGaussian(CMatrix* mat,double mu, double sigma)
{
	double numel = mat->cols * mat->rows;
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,this->stream,numel,mat->mat,mu,sigma);
	return;
}
MKL_INT CRndStream::vRngGaussian(CVector* vec,double mu,double sigma)
{
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,vec->length,vec->data,mu,sigma);
	return 1;
}
void CRndStream::vRngGaussianMV(CMatrix* mat, CVector* mu, CMatrix* sigma)
{
	//assert that mu is vector and that sigma is a square matrix
	//~~
	//if mu is a row vector convert to a column vector
	//
	//each row of mat will have one realization with mean mu and variance sigma
	CMatrix chol(sigma->rows,sigma->cols);
	sigma->Copy(&chol);
	LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',chol.cols,chol.mat,chol.cols);
	vdRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_ICDF,stream,mat->rows,mat->mat,mat->cols,VSL_MATRIX_STORAGE_FULL,mu->data,chol.mat);

}

void CRndStream::vRngUniform(CMatrix* mat, MKL_INT a,MKL_INT b)
{
	for(int i =0;i < mat->rows;i++)
	{
		vdRngUniform(UNIFORM_METHOD,stream,mat->cols,(mat->mat + i*mat->cols),a,b);
	}

}

MKL_INT CRndStream::vRngUniformDiscrete(CMatrix* m,MKL_INT a,MKL_INT b)
{
	int* ptr = (int*) mkl_calloc(1,m->rows * m->cols * sizeof(int),64);
	if(ptr == NULL)
	{
		printf("Memory Init Failed CRndStream::vRngUniformDiscrete\n");
		assert(false);
	}
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,m->rows * m->cols,ptr,a,b);
	for(int i =0;i < m->rows;i++)
		for(int j = 0;j < m->cols;j++)
		{
			m->mat[i*m->cols + j] = (double)ptr[i*m->cols + j];
		}

	mkl_free(ptr);
	return 1;
}
