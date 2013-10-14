/*
 * CVector.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#include "CVector.h"
#include <math.h>
CVector::CVector()
{
	data = NULL;
	bColVector = -1;
	length = 0;
}
CVector::CVector(int l,bool bCol) {
	data = NULL;
	bColVector = -1;
	length = 0;
	Create(l,bCol);
}
MKL_INT CVector::Create(MKL_INT l,bool bcol)
{
	if(data != NULL)
	{
		mkl_free(data);
		data = NULL;
	}
	data = (double*) mkl_calloc(1,l*sizeof(double),64);
	if(data == NULL)
	{
		eprintf("Memory Init Failed CVector::Create \n");
//		assert(false);
		return -1;
	}
	bColVector = bcol;
	length = l;
	return 0;
}
void CVector::Destroy()
{
	if(data!= NULL)
		mkl_free(data);
	data = NULL;
	length = 0;
	bColVector = -1;
}
CVector::CVector(double* ptr,int l, bool bCol)
{
	data = (double*) mkl_calloc(1,l*sizeof(double),64);
	if(data == NULL)
	{
		eprintf("Memory Init Failed CVector::CVector \n");
//		assert(false);
	}
	bColVector = bCol;
	length = l;
	cblas_dcopy(l,ptr,1,data,1);
}
void CVector::Copy(CVector* pSrc)
{
//	if(data != NULL)
//		mkl_free(data);
	if(data == NULL)
		data = (double*)mkl_calloc(1,pSrc->length * sizeof(double),64);
	if(data == NULL)
	{
		eprintf("Memory Init Failed CVector::Copy \n");
//		assert(false);
	}
	bColVector = pSrc->bColVector;
	length = pSrc->length;
	cblas_dcopy(length,pSrc->data,1,data,1);

}
void CVector::Print(bool flag)
{
	for(int i=0;i < length;i++)
	{
		printf("%6.2f",data[i]);
		if(!flag)
		{
			if(bColVector)
				printf("\n");
			else
				printf(",");
		}
		else
			printf(",");
	}
	printf("\n");
}

MKL_INT CVector::ScaleI(double a)
{
	cblas_dscal(length,a,data,1);
	return 0;
}
MKL_INT CVector::Threshold(double alpha)
{
	//all values below alpha will be set to zero
	for(int i =0 ; i < length;i++)
	{
		if(abs(data[i]) < alpha)
			data[i] = 0;
	}
	return 0;
}
MKL_INT CVector::FillVec()
{
	for(int i =0;i < length;i++)
		data[i] = i;
	return 0;
}
double CVector::DotProd(CVector* opb)
{
	if(length != opb->length)
		return -1;
	return cblas_ddot(length,data,1,opb->data,1);
}
MKL_INT CVector::SubtractI(CVector* opB,double alpha,double beta)
{
	//same length and same orientation
	if(length == opB->length && bColVector == opB->bColVector)
		cblas_daxpby(length,beta,opB->data,1,alpha,data,1);
	else
		return -1;
	return 0;
}
MKL_INT CVector::AddI(CVector* opB, double alpha, double beta)
{
	//Implements the equation  y := alpha * x + beta * y
	//this :=  alpha * (this) + beta * ( opB)
	//same length and same orientation
	if(length == opB->length && bColVector == opB->bColVector)
		cblas_daxpby(length,beta,opB->data,1,alpha,data,1);
	else
	{
//		assert(false);
		eprintf("CVector:AddI Invalid operands");
		return -1;
	}
	return 0;
}
MKL_INT CVector::ReadMMFile(char* strfilename)
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i;

    if ((f = fopen(strfilename, "r")) == NULL) {
           printf("Could not file vector input file : %s.\n", strfilename);
           return -1;
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner in input file %s.\n", strfilename);
        return -1;
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return -1;
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
        printf("Could not process Matrix Market size in input file: %s.\n", strfilename);
        return -1;
    }


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    this->Create(nz,0);
//    prob->y.reserve(nz);
    int I,J;
    double val;

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I, &J, &val);
        I--;  /* adjust from 1-based to 0-based */
        J--;
        assert(J==0);
//        prob->y.push_back(val);
        data[I] = val;
    }

    if (f !=stdin) fclose(f);
	return 0;
}
MKL_INT CVector::Ones()
{
	for(int i =0;i < length;i++)
		data[i] = 1.0;
	return 0;
}
MKL_INT CVector::Max(double d,CVector* vecRes)
{
	assert(length == vecRes->length);

	for(int i =0; i < length;i++)
	{
		if(d > data[i])
			vecRes->data[i] = d;
		else
			vecRes->data[i] = data[i];
	}
	return 0;
}
double CVector::Norm(MKL_INT n)
{
	double s = 0;
	if(n == 1)
	{
		for(int i =0;i < length;i ++)
			s = s + abs(data[i]);
	}
	else if(n == 2)
	{
		for(int i =0;i < length;i ++)
			s = s + (data[i] * data[i]);
		s = sqrt(s);
	}
	return s;
}
double CVector::Sum()
{
	//return the sum of the elements of the vector
	double s = 0;
	for(int i =0; i < length;i++)
		s = s + data[i];
	return s;
}

MKL_INT CVector::WriteMMFile(char* strfilename)
{
	MM_typecode matcode;
	int i;
	MKL_INT M;
	MKL_INT N;
	MKL_INT nz;//number of nonzero elements in the array
	if(bColVector == false)
	{
		M = 1;
		N = length;
		nz = length;
	}
	else
	{
		M = length;
		N = 1;
		nz = length;
	}

	mm_initialize_typecode(&matcode);
	mm_set_matrix(&matcode);
	mm_set_coordinate(&matcode);
	mm_set_real(&matcode);
	FILE * fp = fopen(strfilename,"w");
	if(fp == NULL)
	{
		char str[512] = {0};
		sprintf(str,"CVector:WriteMMFile. Could not open file: %s ",strfilename);
		eprintf(str);
	}


	mm_write_banner(fp, matcode);
	mm_write_mtx_crd_size(fp, M, N, nz);

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	  of a vector has index 1, not 0.  */

	for (i=0; i<nz; i++)
	{
		if(bColVector)
			fprintf(fp, "%d %d %10.16g\n", i+ 1,N, data[i]);
		else
			fprintf(fp, "%d %d %10.16g\n", M,i+ 1,data[i]);
	}
	fclose(fp);

	return 0;
}

MKL_INT CVector::IsZero()
{
	int flag = 1;
	for(int i =0;i < length;i++)
	{
		if(data[i] != 0)
			flag = 0;
	}
	return flag;
}

MKL_INT CVector::Cardinality()
{
	MKL_INT card = 0;
	for(int i = 0;i < length;i++)
	{
		if(data[i] != 0)
		{
			card++;
		}
	}
	return card;
}
CVector::~CVector()
{
	if(data != NULL)
		mkl_free(data);
	data = NULL;
	bColVector = -1;
	length = 0;
}


