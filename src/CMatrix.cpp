/*
 * CMatrix.cpp
 *
 *  Created on: Apr 19, 2013
 *      Author: root
 */

#include "CMatrix.h"
#include "CRndStream.h"
#include <matio.h>
#include <stdarg.h>
#include <string.h>

CMatrix::CMatrix()
{
	rows = 0;
	cols = 0;
	mat = NULL;
}

CMatrix::CMatrix(int m, int n) {
	mat = (double*) mkl_calloc(1,(m + 1)*n*sizeof(double),64);
	if(mat == NULL)
	{
		printf("Memory Init Failed CMatrix::CMatrix\n");
		assert(false);
	}
	rows = m;
	cols = n;
}
MKL_INT CMatrix::Create(MKL_INT m,MKL_INT n)
{
	if(mat != NULL)
		mkl_free(mat);

	mat = (double*) mkl_calloc(1,m*n*sizeof(double),64);
	if(mat == NULL)
	{
		printf("Memory Init Failed CMatrix::Create\n");
//		assert(false);
		rows = 0;
		cols = 0;
		return -1;
	}
	rows = m;
	cols = n;
	return 0;
}
CMatrix::~CMatrix() {
	if(mat != NULL)
		mkl_free(mat);
	mat = NULL;
	rows = 0;
	cols = 0;
}
MKL_INT CMatrix::Sprand( double sparsity)
{
	CRndStream rnd(time(NULL));
	CMatrix init(rows,cols);
	CMatrix sp(rows,cols);
	rnd.vRngGaussian(&init,0,1);
//	initial.Scale(5);
//	initial.Threshold(0.5);
//	initial.Print();
	rnd.vRngUniform(&sp,-1,1);
//	printf("The spare matrix\n");
//	sparse.Print();
	//now we set all the elements above the threshold zero, this way we can increase or decrease the sparsity

	sp.BinThreshold(sparsity);
//	printf("the threholded sparse matrix\n");
//	sparse.Print();
	//element by elem multiplication
	init.ElemMul(&sp,this);
	return 0;
}
MKL_INT CMatrix::Release()
{
	if(mat != NULL)
	{
		mkl_free(mat);
		mat = NULL;
		cols = 0;
		rows = 0;
	}
	return 0;
}

MKL_INT CMatrix::Vec(CVector* vec)
{
	//transpose the matrix then copy the mat member to the data member of the vector
	assert(mat != NULL);
	CMatrix m(rows,cols);
	mkl_domatcopy('R','T',rows,cols,1,mat,cols,m.mat,cols);
	cblas_dcopy(rows * cols,m.mat,1,vec->data,1);

	return 0;
}
MKL_INT CMatrix::UnVec(CVector* pvec)
{
	cblas_dcopy(pvec->length,pvec->data, 1, mat,1);
	//now we transpose
	mkl_dimatcopy('R','T',rows,cols,1,mat,cols,rows);
	return 0;
}
MKL_INT CMatrix::Transpose()
{

	mkl_dimatcopy('R','T',rows,cols,1,mat,cols,rows);
	MKL_INT tmp = rows;
	rows = cols;
	cols = tmp;
	return 0;
}
void CMatrix::Set(int m,int n,double val)
{
	this->mat[(m*this->cols) + n] = val;
}

void CMatrix::Copy(CMatrix * dest,bool trans)
{
	// TODO if the trans param is true we would need to change the rows and cols in the matrices
	//for the time being
	if(!trans)
	{
		assert(dest->rows == rows);
		assert(dest->cols == cols);

		mkl_domatcopy('R','N',dest->rows,dest->cols,1,this->mat,this->cols,(*dest).mat,dest->cols);
	}
	else
	{
		assert(dest->cols == rows);
		assert(dest->rows == cols);
		mkl_domatcopy('R','T',rows,cols,1,mat,cols,dest->mat,dest->cols);
	}


//	for(int i = 0; i < rows;i++)
//	{
//		cblas_dcopy(cols,mat + (i*cols),1,dest->mat + (i*cols),1);
//	}
}

void CMatrix::Print()
{
	for (int i = 0; i < rows;i++)
	{
		for (int j = 0;j < cols;j++)
		{
		 	printf("%6.2f ",mat[(cols * i + j)]);
		}
		printf("\n");

	}
	printf("\n");

}

void CMatrix::Eye(void)
{
	//assert that the matrix is square
	for(int i = 0; i < cols;i++)
		mat[cols * i + i] = 1.0f;
}
void CMatrix::BinThreshold(double thr)
{
	for(int i =0; i < rows;i++)
	{
		for(int j = 0;j < cols;j++)
		{
			if(fabs(mat[cols*i + j]) > thr)
				mat[cols*i + j] = 0;
			else if(mat[cols*i + j] > 0)
				mat[cols*i + j] = +1;
			else
				mat[cols*i + j] = -1;
		}
	}
}
void CMatrix::Threshold(double thr)
{
	//values below thr will be set to zero
	for(int i =0;i < rows;i++)
	{
		for(int j=0;j < cols;j++)
		{
			if(fabs(mat[cols*i + j]) < thr)
				mat[cols*i + j] = 0;
		}
	}
}

void CMatrix::ElemMul(CMatrix* opb, CMatrix* res)//this matrix itself is the first operand
{
	// all three matrices must have the same dimensions, and so will the result have

	for(int i =0; i < rows;i++)
	{
		vdMul(cols,mat + (i*cols),opb->mat + (i*cols),res->mat + (i*cols));
	}
}

MKL_INT CMatrix::Scale(double alpha,bool btranspose)
{
	if(btranspose)
	{
		mkl_dimatcopy('R','T',rows,cols,alpha,mat,cols,cols);
	}
	else
		mkl_dimatcopy('R','N',rows,cols,alpha,mat,cols,cols);
	return 1;
}
MKL_INT CMatrix::Add(CMatrix* opA,CMatrix* opB,double alpha,double beta)
{
	if(opA->rows != opB->rows && opA->cols != opB->cols)
	{
		return -1;
	}
	MKL_Domatadd('R','N','N',opA->rows,opA->cols,alpha,opA->mat,cols,beta,opB->mat,cols,mat,cols);
	return 1;
}

MKL_INT CMatrix::IsZero(CVector* idx)
{
	//we assume the vector to be of length = this->rows*cols
	MKL_INT cnt = 0;
	if(idx->length < rows*cols)
		return -1;
	for(int i = 0;i < rows;i++)
	{
		for(int j = 0;j < cols;j++)
		{
			if(mat[i*cols + j] == 0)
			{
				idx->data[cnt] = i*cols + j;
				cnt++;
			}
		}
	}
	return cnt;
}

MKL_INT CMatrix::IsNonZero(CVector* idx)
{
	//we assume the vector to be of length = this->rows*cols
	MKL_INT cnt = 0;
	if(idx->length < rows*cols)
		return -1;
	for(int i = 0;i < rows;i++)
	{
		for(int j = 0;j < cols;j++)
		{
			if(mat[i*cols + j] != 0)
			{
				idx->data[cnt] = i*cols + j;
				cnt++;
			}
		}
	}
	return cnt;
}

MKL_INT CMatrix::GetVec(CVector* pRes,MKL_INT num,bool row)
{
	if(pRes == NULL)
		return -1;
	pRes->bColVector = !row;
	if(row)
	{
		cblas_dcopy(cols,mat + (num*cols),1,pRes->data,1);
	}
	else
	{
		CVector vec(cols,1);
		vec.data[num] = 1;
		cblas_dgemv(CblasRowMajor,CblasNoTrans,rows,cols,1,mat,cols,vec.data,1,0,pRes->data,1);
	}
	return 0;
}

MKL_INT CMatrix::SetVec(CVector* pvec,MKL_INT num,bool row)
{
	if(pvec == NULL || pvec->length != cols)
	{
		assert(false);
		return -1;
	}
	if(row)
	{
		cblas_dcopy(cols,pvec->data,1,mat + (num * cols),1);
	}
	else
	{
		assert(false);
		return -1;//not implemented yet
	}
	return 0;
}
MKL_INT CMatrix::SetVec(CMatrix* pvec,MKL_INT num,bool row)
{
	if(pvec == NULL || pvec->cols != cols || pvec->rows != 1)
	{
		assert(pvec);
		assert(pvec->cols == cols);
		assert(pvec->rows == 1);
		assert(false);
		return -1;
	}
	if(row)
	{
		cblas_dcopy(cols,pvec->mat,1,mat + (num * cols),1);
	}
	else
	{
		assert(false);
		return -1;//not implemented yet
	}
	return 0;
}
MKL_INT CMatrix::MultiplyV(CVector* res,CVector* opVector,bool TransposeMat,double alpha ,double beta )
{
	if(TransposeMat)
	{
		assert(rows == opVector->length);

		cblas_dgemv(CblasRowMajor,CblasTrans,rows,cols,alpha,mat,cols,opVector->data,1,beta,res->data,1);
	}
	else
	{
		assert(cols == opVector->length);
		cblas_dgemv(CblasRowMajor,CblasNoTrans,rows,cols,alpha,mat,cols,opVector->data,1,beta,res->data,1);
	}
	return 0;
}
MKL_INT CMatrix::VecProd(CVector* opa,CVector* opb)
{

	if(opa->bColVector == true && opb->bColVector == false && opa->length == opb->length && opa->length == rows && rows == cols) //sanity check
	{
		for(int i = 0;i < opa->length;i++)
		{
			cblas_daxpy(opa->length,opa->data[i],opb->data,1,mat + (i * opa->length),1);
		}
	}
	else
	{
		printf("Error -> CMatrix::VecProd\n");
		return -1;
	}
	return 0;
}
MKL_INT CMatrix::InvMat(CMatrix* ptrres)
{
	assert(rows == cols);

	MKL_INT* i = NULL;
	lapack_int info;

	i = (lapack_int*) mkl_calloc(1,cols * sizeof(lapack_int),64);
	if( i == NULL)
	{
		printf("Memory Init Failed CMatrix::InvMat\n");
		assert(false);
	}
	//get the tri factorization
	// then we get the inverse of the matrix
	cblas_dcopy(rows * cols,mat,1,ptrres->mat,1);
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,rows,cols,ptrres->mat,cols,i);

	info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,rows,ptrres->mat,cols,i);

	mkl_free(i);
	return 0;
}
MKL_INT CMatrix::MMMul(CMatrix* opb,CMatrix* res,double alpha,double beta,bool btransa, bool btransb)
{
//	res = alpha * op(A) * op(B) + beta * C

	if(!btransa && !btransb)
	{
		assert(cols == opb->rows);
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, rows,opb->cols,cols,alpha,mat,cols,opb->mat,opb->cols,beta,res->mat,res->cols);
	}
	else if(!btransa && btransb)
	{
		assert(cols == opb->cols);
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, rows,opb->rows,cols,alpha,mat,cols,opb->mat,opb->cols,beta,res->mat,res->cols);
	}
	else if(btransa && !btransb)
	{
		assert(rows == opb->rows);
		cblas_dgemm(CblasRowMajor, CblasTrans,CblasNoTrans, cols,opb->cols,rows,alpha,mat,cols,opb->mat,opb->cols,beta,res->mat,res->cols);
	}
	else if(btransa && btransb)
	{
		assert(rows == opb->cols);
		cblas_dgemm(CblasRowMajor, CblasTrans,CblasTrans, cols,opb->rows,rows,alpha,mat,cols,opb->mat,opb->cols,beta,res->mat,res->cols);
	}

	return 0;
}
MKL_INT CMatrix::EyeSubMat(double alpha)
{
	if(rows != cols)
	{
		printf("Error -> CMatrix::EyeSubMat \n");
		return -1;
	}
	//inplace scale the matrix by minus 1
	mkl_dimatcopy('R','N',rows,cols,alpha,mat,cols,cols);
	//now add one to the diagonal
	for(int i = 0;i < rows;i++)
	{
		mat[i + i*cols ] = mat[i + i*cols] + 1;
	}
	return 0;
}
MKL_INT CMatrix::FillMat()
{
	for(int i = 0; i < rows;i++)
	{
		for(int j = 0;j < cols;j++)
		{
			mat[cols*i + j] = cols*i+j;
		}
	}

	return 0;
}



MKL_INT CMatrix::DelMat(char*filename, char* varname)
{
	mat_t* matfile;
	matfile = Mat_Open(filename,MAT_ACC_RDWR | MAT_FT_MAT73);
	int ret = Mat_VarDelete(matfile,varname);
	Mat_Close(matfile);
	return ret;
}
