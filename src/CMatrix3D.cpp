/*
 * CMatrix3D.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#include "CMatrix3D.h"

CMatrix3D::CMatrix3D(MKL_INT r,MKL_INT c, MKL_INT z)
{
//Since we are planning on very large non-sparse arrays therefore we will make this a list of CMatrix pointers
//this will enable the creation of large matrix arrays in different memory regions thus avoid wasting memory at the cost of performance

	pMats = (CMatrix*) mkl_malloc(z*sizeof(CMatrix),64);

	depth = z;

	for(int i =0; i < z;i++)
	{

		pMats[i].Create(r,c);
	}
}

CMatrix3D::~CMatrix3D()
{
	//first free each matrix and then release the pointer array

	for(int i = 0;i < depth;i++)
	{
		pMats[i].~CMatrix();//force all the arrays to dealloc
	}
	mkl_free(pMats);

}
void CMatrix3D::Print(void)
{
	for(int i =0;i < depth;i++)
	{
		printf("Depth = %d\n",i);
		pMats[i].Print();
		printf("\n");
	}
}

void CMatrix3D::Set(CMatrix* mat,MKL_INT z)
{
	if(z > 0 && z <= depth)
		mat->Copy(&(pMats[z]));
}

void CMatrix3D::Eye(void)
{
	//This would set all the matrices in the depths to identity
	for(int i=0;i < depth;i++)
	{
		pMats[i].Eye();
	}
}
