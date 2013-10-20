/*
 * CTimeVaryingNW.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: root
 */

#include "CTimeVaryingNW.h"
#include "CRndStream.h"
CTimeVaryingNW::CTimeVaryingNW()
{
	nGene = 0;
	perChange = 0;
	matGeneInteractions = NULL;
	sparsitylvl = 0;
	nTimePts = 0;
	nObs = 0;
	pVecX = NULL;
	pVecY = NULL;
	pMatX = NULL;
	pMatY = NULL;
}
MKL_INT CTimeVaryingNW::CreateNW(MKL_UINT ngene,MKL_INT ntimepts,MKL_INT nobs,double pchange,double sparsitylvl, MKL_INT seed,int streamnum)
{
	//creates network with gene expression in matrix form rather than arrays
	nGene = ngene;
	perChange = pchange;
	nTimePts = ntimepts;
	nObs = nobs;
	pVecX = NULL;
	pVecY = NULL;


	//pMatX and pMatY would be matrices with nGene rows and nObs columns, the total number of the matrices would be nTimePts
	pMatX = new CMatrix[nTimePts];
	pMatY = new CMatrix[nTimePts];

	for(int i =0;i < nTimePts;i++)
	{
		pMatX[i].Create(nGene,nObs);
		pMatY[i].Create(nGene,nObs);
	}

	//now we create the conn mat for time zero
	CMatrix initial(nGene,nGene);
	CMatrix sparse(nGene,nGene);

	//create a 3d mat to store the time varying info
	matGeneInteractions = new CMatrix[nTimePts];

	for(int i =0;i < nTimePts;i++)
	{
		matGeneInteractions[i].Create(nGene,nGene);
	}
//	pVecX = new
	//We generate the connectivity matrix and then multiply it with a uniform random binary matrix that has been thresholded

	CRndStream rnd(seed,streamnum);
	//we create the first connectivity matrix
	rnd.vRngGaussian(&initial,0,1);
	initial.Scale(3);
	initial.Threshold(0.5);
	initial.Print();
	rnd.vRngUniform(&sparse,-1,1);
//	printf("The spare matrix\n");
//	sparse.Print();
	//now we set all the elements above the threshold zero, this way we can increase or decrease the sparsity

//	printf("Sparsity Level is %f\n",sparsitylvl);
	sparse.BinThreshold(sparsitylvl);
//	printf("the threholded sparse matrix\n");
//	sparse.Print();
	//element by elem multiplication
	initial.ElemMul(&sparse,matGeneInteractions);
	//set elements below 0.5 to zero
//	matGeneInteractions->pMats[0]->Print();
	//create the gene expression for time zero
//	printf("Interaction matrix time %d\n",0);
//	matGeneInteractions[0].Print();
//	matGeneInteractions->pMats[0].Print();

	//generate the matrix X with nObs columns and nGene rows
	rnd.vRngUniform(&(pMatX[0]),0,1);
	pMatX[0].Scale(10);
	matGeneInteractions[0].MMMul(&pMatX[0],&pMatY[0]);

//	printf("The gene expression observations for time %d are \n",0);
//	printf("Matrix X \n");
//	pMatX[0].Print();
//	printf("Matrix Y \n");
//	pMatY[0].Print();

	//now that the base matrix is created we can change it for each time step using the pchange evolution

	for(int i = 1;i < nTimePts;i++)
	{
		matGeneInteractions[i-1].Copy(&(matGeneInteractions[i]));
		//recreate the gaussian matrix, scale and add it to the previous matrix
		rnd.vRngGaussian(&initial,0,1);
		initial.ElemMul(&sparse,&initial);
		matGeneInteractions[i].Add(&initial,&(matGeneInteractions[i-1]),0.1,1);
		matGeneInteractions[i].Threshold(0.5);
//
//		printf("the network after perturbation\n\n");
//		matGeneInteractions->pMats[i]->Print();
		//now we need to choose which edges to kill and what new to create

		//find the indices of all the nonzero elements in the connectivity matrix
		CVector idxZero(nGene*nGene,0);
		CVector idxNonZero(nGene*nGene,0);
		MKL_INT cntZero  = matGeneInteractions[i].IsZero(&idxZero);
		MKL_INT cntNonZero = matGeneInteractions[i].IsNonZero(&idxNonZero);

		//here we can check and ensure the number of zero / nonzero entries in the matrix


		//we generate a random sequence of numbers from 0 to cntZero, length being floor(perChange*cntZero)
		//these are the edges that will be born, then set them to gaussian nonzero numbers

		//first we create the edges
		double numNewEdges = floor(perChange*cntNonZero);
		if(numNewEdges > 0 )
		{
			CMatrix newEdge(1,numNewEdges);
			rnd.vRngUniformDiscrete(&newEdge,0,cntZero);


			//now we pick the edge numbers and create an edge there
			for(int j =0;j < numNewEdges;j++)
			{
				int curEdgeNum = (int) idxZero.data[(int)newEdge.mat[j]];
				CMatrix newVal(1,1);
				rnd.vRngGaussian(&newVal,0,1);
				matGeneInteractions[i].mat[curEdgeNum] = newVal.mat[0];
			}
		}
//		printf("The network after new edges\n");
//		matGeneInteractions->pMats[i]->Print();

		double numDelEdges = floor(perChange * cntNonZero);
		if(numDelEdges > 0)
		{
			CMatrix delEdges(1,numDelEdges);
			rnd.vRngUniformDiscrete(&delEdges,0,cntNonZero);
			for(int j =0; j < numDelEdges;j++)
			{
				int curEdgeNum = (int) idxNonZero.data[(int)delEdges.mat[j]];
				matGeneInteractions[i].mat[curEdgeNum] = 0;
			}
		}
//		printf("The final network at time %d \n",i);
//		matGeneInteractions->pMats[i]->Print();

//		printf("Interaction matrix time %d\n",i);
//		matGeneInteractions->pMats[i].Print();

		//generate the matrix X with nObs columns and nGene rows
		rnd.vRngUniform(&(pMatX[i]),0,1);
		pMatX[i].Scale(10);
		matGeneInteractions[i].MMMul(&pMatX[i],&pMatY[i]);

//		printf("The gene expression observations for time %d are \n",i);
//		printf("Matrix X \n");
//		pMatX[i].Print();
//		printf("Matrix Y \n");
//		pMatY[i].Print();
	}


	return 0;
}
CTimeVaryingNW::CTimeVaryingNW(MKL_UINT nG, MKL_INT ntpts, MKL_INT nO,double perCh,MKL_INT seed)
{
	nGene = nG;
	perChange = perCh;
	nTimePts = ntpts;
	nObs = nO;
	pMatX = NULL;
	pMatY = NULL;

	pVecX = NULL;
	pVecY = NULL;

	//we initialize the x and y vectors
	pVecX = new CVector*[nTimePts];
	pVecY = new CVector*[nTimePts];
	for(int i =0;i < nTimePts;i++)
	{
		pVecX[i] = new CVector[nObs];//(nGene,1);
		pVecY[i] = new CVector[nObs];//(nGene,1);
		for(int j = 0;j < nObs;j++)
		{
			pVecX[i][j].Create(nGene,1);
			pVecY[i][j].Create(nGene,1);
		}

	}

	CMatrix initial(nGene,nGene);
	CMatrix sparse(nGene,nGene);
	sparsitylvl = 0.08; //lesser means more zeros

	//create a 3d mat to store the time varying info
	matGeneInteractions = new CMatrix[nTimePts];
	for(int i = 0;i < nTimePts;i++)
	{
		matGeneInteractions[i].Create(nGene,nGene);
	}
//	pVecX = new
	//We generate the connectivity matrix and then multiply it with a uniform random binary matrix that has been thresholded

	CRndStream rnd(seed);
	//we create the first connectivity matrix
	rnd.vRngGaussian(&initial,0,1);
	initial.Scale(3);
	initial.Threshold(0.5);
//	initial.Print();
	rnd.vRngUniform(&sparse,-1,1);
//	printf("The spare matrix\n");
//	sparse.Print();
	//now we set all the elements above the threshold zero, this way we can increase or decrease the sparsity

	sparse.BinThreshold(sparsitylvl);
//	printf("the threholded sparse matrix\n");
//	sparse.Print();
	//element by elem multiplication
	initial.ElemMul(&sparse,matGeneInteractions);
	//set elements below 0.5 to zero
//	matGeneInteractions->pMats[0]->Print();

	//now we update all the time instances for the gene network

	for(int i = 1;i < nTimePts;i++)
	{
		matGeneInteractions[i-1].Copy(&(matGeneInteractions[i]));
		//recreate the gaussian matrix, scale and add it to the previous matrix
		rnd.vRngGaussian(&initial,0,1);
		initial.ElemMul(&sparse,&initial);
		matGeneInteractions[i].Add(&initial,&(matGeneInteractions[i-1]),0.1,1);
//
//		printf("the network after perturbation\n\n");
//		matGeneInteractions->pMats[i]->Print();
		//now we need to choose which edges to kill and what new to create

		//find the indices of all the nonzero elements in the connectivity matrix
		CVector idxZero(nGene*nGene,0);
		CVector idxNonZero(nGene*nGene,0);
		MKL_INT cntZero  = matGeneInteractions[i].IsZero(&idxZero);
		MKL_INT cntNonZero = matGeneInteractions[i].IsNonZero(&idxNonZero);

		//we generate a random sequence of numbers from 0 to cntZero, length being floor(perChange*cntZero)
		//these are the edges that will be born, then set them to gaussian nonzero numbers

		//first we create the edges
		double numNewEdges = floor(perChange*cntNonZero);
		CMatrix newEdge(1,numNewEdges);
		rnd.vRngUniformDiscrete(&newEdge,0,cntZero);
		//now we pick the edge numbers and create an edge there
		for(int j =0;j < numNewEdges;j++)
		{
			int curEdgeNum = (int) idxZero.data[(int)newEdge.mat[j]];
			CMatrix newVal(1,1);
			rnd.vRngGaussian(&newVal,0,1);
			matGeneInteractions[i].mat[curEdgeNum] = newVal.mat[0];
		}

//		printf("The network after new edges\n");
//		matGeneInteractions->pMats[i]->Print();

		double numDelEdges = floor(perChange * cntNonZero);
		CMatrix delEdges(1,numDelEdges);
		rnd.vRngUniformDiscrete(&delEdges,0,cntNonZero);
		for(int j =0; j < numDelEdges;j++)
		{
			int curEdgeNum = (int) idxNonZero.data[(int)delEdges.mat[j]];
			matGeneInteractions[i].mat[curEdgeNum] = 0;
		}
//		printf("The final network at time %d \n",i);
//		matGeneInteractions->pMats[i]->Print();
	}

	//now we create the x and y for each matrix
	for(int i = 0;i < nTimePts;i++)
	{
		//for each time point that we have we generate nObs time X and Y matrices

//		printf("Interaction matrix time %d\n",i);
//		matGeneInteractions->pMats[i].Print();
		for(int j=0; j < nObs;j++)
		{
			rnd.vRngGaussian(&pVecX[i][j],0,1);
			pVecX[i][j].ScaleI(10);
			matGeneInteractions[i].MultiplyV(&pVecY[i][j],&pVecX[i][j],true);
			//both pVecX and pVecY are rows vecs
			pVecX[i][j].bColVector = false;
			pVecY[i][j].bColVector = false;
//			printf("Input %d\n",j);
//			pVecX[i][j].Print();
//			printf("Output \n");
//			pVecY[i][j].Print();
		}
	}

}

CTimeVaryingNW::~CTimeVaryingNW()
{

	for(int i = 0;i < nTimePts;i++)
		matGeneInteractions[i].Release();

	delete [] matGeneInteractions;


	for(int i =0;i < nTimePts;i++)
	{
		for(int j = 0;j < nObs;j++)
		{
			if(pVecX != NULL)
				pVecX[i][j].Destroy();
			if(pVecY != NULL)
				pVecY[i][j].Destroy();
		}
//		delete pVecX[i];
//		delete pVecY[i];
		if(pMatX != NULL)
			pMatX[i].Release();
		if(pMatY != NULL)
			pMatY[i].Release();

	}
	if(pVecX != NULL)
		delete pVecX;
	if(pVecY != NULL)
		delete pVecY;
	if(pMatX != NULL)
		delete [] pMatX;
	if(pMatY != NULL)
		delete [] pMatY;
}

