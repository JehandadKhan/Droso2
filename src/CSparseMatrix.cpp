/*
 * CSparseMatrix.cpp
 *
 *  Created on: May 2, 2013
 *      Author: root
 */

#include "CSparseMatrix.h"
#include "string.h"

CSparseMatrix::CSparseMatrix() {
	numZero = 0;
	numNonZero = 0;
	Rows = 0;
	Cols = 0;
	columns = NULL;
	values = NULL;
	rowIndex = NULL;
}
CSparseMatrix::CSparseMatrix(MKL_INT r,MKL_INT c)
{
	Rows = r;
	Cols = c;

	numZero = 0;
	numNonZero = 0;
	columns = NULL;
	values = NULL;
	rowIndex = NULL;
}
MKL_INT CSparseMatrix::MulMatVec(CVector* veb,CVector* vecres, double alpha,double beta,bool transa)
{
	//implements vecRes := alpha * op(Sp Matrix) * (opB) + beta * vecRes

	assert(Cols == veb->length);
	assert(values != NULL);
	assert(columns != NULL);
	assert(rowIndex != NULL);

	char ch[] = "G  C  ";
	char tr;
	if(transa)
		tr = 'T';
	else
		tr = 'N';

	mkl_dcsrmv(&tr,&Rows,&Cols,&alpha,ch,values,columns,rowIndex,rowIndex + 1,veb->data,&beta,vecres->data);

	return 0;
}
MKL_INT CSparseMatrix::KronEyeVecVec(MKL_INT idimeye,CVector* vec,CVector* pvop,CVector* pRes,bool bTransMat,double alpha,double beta)
{
	if(vec->bColVector == false)
	{
		Rows = idimeye;
		Cols = idimeye * vec->length;
		//this is the common case with us
		if(values != NULL)
		{
			mkl_free(values);
			values = NULL;
		}
		values = (double*) mkl_calloc(1,vec->length * idimeye * sizeof(double),64);
		if(values == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec values memory alloc failed");
			assert(false);
			return -1;
		}
		for(int i =0;i < idimeye;i++)
		{
			cblas_dcopy(vec->length,vec->data,1,values + (i* vec->length),1);
		}


		if(columns != NULL)
		{
			mkl_free(columns);
			columns = NULL;
		}
		columns = (MKL_INT*)mkl_calloc(1,vec->length* idimeye* sizeof(int),64);
		if(columns == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec columns memory alloc failed");
			assert(false);
			return -1;
		}

		for(int i =0; i < vec->length * idimeye;i++)
		{
			columns[i] = i;
		}

//		for(int i = 0;i < vec->length * idimeye;i++)
//		{
//			printf("%d ",columns[i]);
//		}
//		printf("\n");
		if(rowIndex != NULL)
		{
			mkl_free(rowIndex);
			rowIndex = NULL;
		}
		rowIndex = (MKL_INT*) mkl_calloc(1,(idimeye + 1) * sizeof(int) ,64);
		if(rowIndex == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec rowIndex memory alloc failed");
			assert(false);
			return -1;
		}
		for(int i =0;i < idimeye;i++)
		{
			rowIndex[i] = i * vec->length ;
		}
		rowIndex[idimeye] = idimeye * vec->length ;
		//now we multiply the sparse matrix with vector
		if(bTransMat)
		{
			char ch[] = "G  C  ";
			char tr = 'T';
			MKL_INT k = idimeye*vec->length;
			mkl_dcsrmv(&tr,&idimeye,&k,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
		}
		else
		{
			char ch[] = "G  C  ";
			char tr = 'N';
			MKL_INT k = idimeye*vec->length;
			mkl_dcsrmv(&tr,&idimeye,&k,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
		}
	}
	else
	{
		Cols = idimeye;
		Rows = idimeye * vec->length;
		if(values != NULL)
		{
			mkl_free(values);
			values = NULL;
		}
		values = (double*) mkl_calloc(1,vec->length * idimeye * sizeof(double),64);
		if(values == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec values memory alloc failed");
			assert(false);
			return -1;
		}
		for(int i =0;i < idimeye;i++)
		{
			cblas_dcopy(vec->length,vec->data,1,values + (i* vec->length),1);
		}

		if(columns != NULL)
		{
			mkl_free(columns);
			columns = NULL;
		}
		columns = (MKL_INT*)mkl_calloc(1,vec->length* idimeye* sizeof(int),64);
		if(columns == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec columns memory alloc failed");
			assert(false);
			return -1;
		}

		for(int i =0; i < vec->length * idimeye;i++)
		{
			columns[i] = (int) i / vec->length;
		}
		if(rowIndex != NULL)
		{
			mkl_free(rowIndex);
			rowIndex = NULL;
		}
		rowIndex = (MKL_INT*) mkl_calloc(1,(idimeye *vec->length + 1) * sizeof(int) ,64);
		if(rowIndex == NULL)
		{
			printf("CSparseMatrix::KronEyeMatVec rowIndex memory alloc failed");
			return -1;
		}
		for(int i =0;i < idimeye * vec->length;i++)
		{
			rowIndex[i] = i ;
		}
		rowIndex[idimeye * vec->length] = idimeye * vec->length ;
		//now we multiply the sparse matrix with vector
		if(bTransMat)
		{
			char ch[] = "G  C  ";
			char tr = 'T';
			MKL_INT k = idimeye*vec->length;
			mkl_dcsrmv(&tr,&k,&idimeye,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
		}
		else
		{
			char ch[] = "G  C  ";
			char tr = 'N';
			MKL_INT k = idimeye*vec->length;
			mkl_dcsrmv(&tr,&k,&idimeye,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
		}

	}
	return 0;
}
MKL_INT CSparseMatrix::KronEyeMatVec(MKL_INT idimeye,CMatrix* mat,CVector* pvop,CVector* pRes,bool bTransMat,double alpha,double beta)
{
	assert(mat != NULL);
	assert(pvop != NULL);
	assert(pRes != NULL);
	assert(pvop->length == idimeye * mat->rows);
	assert(pRes->length == idimeye * mat->rows);


	Rows = idimeye * mat->rows;
	Cols = idimeye * mat->cols;
	//this is the common case with us
	if(values != NULL)
	{
		mkl_free(values);
		values = NULL;
	}
	values = (double*) mkl_calloc(1,Rows * Cols * idimeye * sizeof(double),64);
	if(values == NULL)
	{
		printf("CSparseMatrix::KronEyeMatVec values memory alloc failed");
		return -1;
	}
	for(int i =0;i < idimeye;i++)
	{
		cblas_dcopy(mat->rows * mat->cols,mat->mat,1,values + (i* (mat->rows * mat->cols)),1);
	}

	if(columns != NULL)
	{
		mkl_free(columns);
		columns = NULL;
	}
	columns = (MKL_INT*)mkl_calloc(1,mat->rows * mat->cols * idimeye* sizeof(int),64);
	if(columns == NULL)
	{
		printf("CSparseMatrix::KronEyeMatVec columns memory alloc failed");
		return -1;
	}

	for(int i =0; i < idimeye;i++)
	{
		for(int j = 0;j < mat->cols * mat->rows;j++)
		{
			columns[(i * mat->cols * mat->rows ) + j] = (j % mat->cols) +  ( i * mat->cols) ;
		}
	}

	if(rowIndex != NULL)
	{
		mkl_free(rowIndex);
		rowIndex = NULL;
	}
	rowIndex = (MKL_INT*) mkl_calloc(1,((idimeye  * mat->rows)+ 1) * sizeof(int) ,64);
	if(rowIndex == NULL)
	{
		printf("CSparseMatrix::KronEyeMatVec rowIndex memory alloc failed");
		return -1;
	}
	for(int i =0;i < idimeye * mat->rows;i++)
	{
		rowIndex[i] = i * mat->cols;
	}
	rowIndex[idimeye * mat->rows] = idimeye * mat->cols * mat->rows;
	//now we multiply the sparse matrix with vector
	if(bTransMat)
	{
		char ch[] = "G  C  ";
		char tr = 'T';
		MKL_INT k = idimeye* mat->cols;
		MKL_INT m = idimeye* mat->rows;
		mkl_dcsrmv(&tr,&m,&k,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
	}
	else
	{
		char ch[] = "G  C  ";
		char tr = 'N';
		MKL_INT k = idimeye*mat->cols;
		MKL_INT m = idimeye* mat->rows;
		mkl_dcsrmv(&tr,&m,&k,&alpha,ch,values,columns,rowIndex,rowIndex + 1,pvop->data,&beta,pRes->data);
	}

	return 0;
}
void CSparseMatrix::DenseMatrix(CMatrix* pMat)
{
	int job[] = {1,0,0,2,8,1};
	MKL_INT info;
	mkl_ddnscsr(job,&(pMat->rows),&(pMat->cols),pMat->mat,&(pMat->cols),values,columns,rowIndex,&info);
}
MKL_INT CSparseMatrix::Eye(double alpha)
{
	//create a sparse identity matrix with diagonal entries alpha
	assert(Rows == Cols);
	values = (double*) mkl_calloc(1,Rows * sizeof(double),64);
	columns = (MKL_INT*) mkl_calloc(1,Rows * sizeof(MKL_INT),64);
	rowIndex = (MKL_INT*) mkl_calloc(1,(Rows + 1) * sizeof(MKL_INT),64);

	for(int i = 0;i < Rows;i++)
	{
		values[i] = alpha;
		columns[i] = i;
		rowIndex[i] = i;
	}
	rowIndex[Rows] = Rows ;
	return 0;
}
MKL_INT CSparseMatrix::Inv(bool isdiag)
{
	assert(isdiag == true); //only diag matrix
	assert(Rows == Cols);// only square matrices

	for(int i = 0;i < Rows;i++)
	{
		values[i] = 1 / values[i];
	}
	return 0;
}
void CSparseMatrix::Print()
{
	printf("Printing Sparse Matrix\n");
	CMatrix densmat(Rows,Cols);
	DenseMatrix(&densmat);
	densmat.Print();
}
CSparseMatrix::~CSparseMatrix() {

	if(values != NULL)
		mkl_free(values);
	if(columns != NULL)
		mkl_free(columns);
	if(rowIndex != NULL)
		mkl_free(rowIndex);
}

