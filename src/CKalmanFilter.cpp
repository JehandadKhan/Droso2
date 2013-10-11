/*
 * CKalmanFilter.cpp
 *
 *  Created on: May 2, 2013
 *      Author: root
 */

#include "CKalmanFilter.h"
#include <mpi.h>
#define DELTALAMBDA 0.1
#define DELTAFINELAMBDA 0.01

CKalmanFilter::CKalmanFilter()
{
	pEstA = NULL;
	pEstAApr = NULL;
	pEstSmtA = NULL;
	nGenes = 0;
	nEffMeas = 0;
	matSmthEst = NULL;
	matEstA = NULL;
	matEstAApr = NULL;
	nTimePts = 0;

	matPkN = NULL;
	matPkk = NULL;
	matPkp1k = NULL;

	pgEstA = NULL;
	pgEstAApr = NULL;
	pgEstSmtA = NULL;
}

CKalmanFilter::CKalmanFilter(CVector** ppvecx, CVector** ppvecy,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dlambda)
{
	if(ppvecx == NULL)
		return;
	if(ppvecy == NULL)
		return;

	pEstA = NULL;
	pEstAApr = NULL;
	pEstSmtA = NULL;
	nGenes = 0;
	nEffMeas = 0;
	matSmthEst = NULL;
	matEstA = NULL;
	matEstAApr = NULL;
	nTimePts = 0;


	matPkN = NULL;
	matPkk = NULL;
	matPkp1k = NULL;

	pgEstA = NULL;
	pgEstAApr = NULL;
	pgEstSmtA = NULL;
//	nGenes = ngenes;
//	nEffMeas = ntimepts*nobs + 1;
//	CMatrix matSigmaQ(ngenes,ngenes);
//	matSigmaQ.Eye();
//	matSigmaQ.Scale(7);
//
//	double dsigmaR = 0.2;
//	CMatrix3D mat3dbarPkk(ngenes,ngenes,nEffMeas);
//	CMatrix3D mat3dbarPkp1k(ngenes,ngenes,nEffMeas);
//	CMatrix3D mat3dbarPkN(nGenes,nGenes,nEffMeas);
//	pEstA = NULL;
//	pEstAApr = NULL;
//	pEstSmtA = NULL;
//
//	pEstA = new CVector*[nEffMeas];
//	pEstAApr = new CVector*[nEffMeas];
//	pEstSmtA = new CVector*[nEffMeas];
//	for(int i=0; i < nEffMeas;i++)
//	{
//		pEstA[i] = new CVector(ngenes*ngenes,1);
//		pEstAApr[i] = new CVector(ngenes*ngenes,1);
//		pEstSmtA[i] = new CVector(ngenes*ngenes,1);
//	}
//	//CVector vecXnn(ngenes*ngenes,1);
//	//CVector vecXnn1(ngenes*ngenes,1);
//	printf("Kalman Filter Constructor\n");
//
//	pEstA[0]->Copy(ic); //initialize the state variable with the initial condition
//	CVector vecKalmanGain(ngenes,1);
//	int efftime= 0;
//	int effobs = 0;
//	CVector temp(nGenes,1);
//	double* arrayLambda = new double[nEffMeas];//store the final lambda to apply at smoothing time
//	for(int j =1; j < nEffMeas;j++) //this is the estimator time
//	{
//		printf("Kalman filter estimation at j = %d\n",j);
//		CSparseMatrix spMatErr;
//		CVector res1(ngenes,1);
//
//		efftime = (int) ((j - 1) / nobs);
//		effobs = (j - 1) - (efftime * nobs);
//		printf("Effective time for Observations is %d\n",efftime);
//		printf("Effective Observation %d at time %d\n",effobs,efftime);
//
//		ppvecx[efftime][effobs].bColVector = true;
//		ppvecy[efftime][effobs].bColVector = true;
//
//
//		printf("Posterior covariance from j = %d \n",j - 1);
//		mat3dbarPkk.pMats[j-1].Print();
//
//		//propagate the state to the apriori estimate
//		pEstAApr[j]->Copy(pEstA[j-1]);
//		//propagate the covariance
//		mat3dbarPkp1k.pMats[j].Add(mat3dbarPkk.pMats[j-1],&matSigmaQ);
//		printf("Apriori Covariance \n");
//		mat3dbarPkp1k.pMats[j].Print();
//		//calc the kalman gain
//
//		mat3dbarPkp1k.pMats[j].MultiplyV(&temp,&ppvecx[efftime][effobs]);
//		double d = temp.DotProd(&ppvecx[efftime][effobs]) + dsigmaR;
//		mat3dbarPkp1k.pMats[j].MultiplyV(&vecKalmanGain,&ppvecx[efftime][effobs]);
//		d = 1 / d;
//		vecKalmanGain.ScaleI(d);
//		//calc the aposteriori state estimate
//
//		ppvecx[efftime][effobs].bColVector = false;//the next call requires it to be transposed
//		spMatErr.KronEyeVecVec(ngenes,&ppvecx[efftime][effobs],pEstAApr[j],&res1);
//		ppvecx[efftime][effobs].bColVector = true;
//		MKL_INT r = res1.AddI(&ppvecy[efftime][effobs],-1,1);
//		//the first sparse product is kron product of an eye mat and a column vector and multiplied with the res above
//		res1.bColVector = true;
//		spMatErr.KronEyeVecVec(ngenes,&vecKalmanGain,&res1,pEstA[j]);
//
//		//next we add inplace to the aposteriori estimate the apriori estimate
//		pEstA[j]->AddI(pEstAApr[j]);
//
//		//finally we update the error covariance
//		//calc the cross product
//		ppvecx[efftime][effobs].bColVector = false;//transpose
//		CMatrix mat(ngenes,ngenes);
//		mat.VecProd(&vecKalmanGain,&ppvecx[efftime][effobs]);
//		//subract from identity
//		mat.EyeSubMat();
//		//now multiply with apriori cov to form the aposteriori cov
//		mat3dbarPkk.pMats[j].MMMul(&mat,&(mat3dbarPkp1k.pMats[j])); //calc complete
//
//		printf("Aposteriori Covariance\n");
//		mat3dbarPkk.pMats[j].Print();
//		//now we apply the lasso algorithm to constrain the estimate
//
//		printf("Begining search for best lambda\n");
//		double delta = 0.1;
//		int numLambda = (dlambda[1] - dlambda[0]) / delta;
//		double* pMse = new double[numLambda];
//		//find the minimum Mse and select that lambda for fine search
//		double min = pMse[0];
//		int minidx = 0;
//		pMse[0] = INFINITY;
//		for(int k =0; k < numLambda;k++)
//		{
//			printf("Lasso for Lambda = %f\n",k * delta + dlambda[0]);
//
//			CVector vecConstrEst(nGenes * nGenes,1);
//			CMatrix matConn(nGenes,nGenes);
//			CVector vecYhat(nGenes,1);
//
//
//			CLasso lasso;
//			//(pEstA[j],&vecConstrEst,k * delta + dlambda[0]);
//			lasso.DoLasso(pEstA[j],&vecConstrEst,k * delta + dlambda[0],1/*rho*/,1/*alpha*/,100/*Max iteration*/,1/*verbose*/);
//
//			//calculate and store the error and lambda of this lasso estimate
//
//			//to calc the error we need to calculate the product estmatrix * vecxx[][] = vecyy_hat
//			//unvec the connectivity matrix
//			printf("calc the mean sq error for the estimated obsn and the true obsn\n");///
//			matConn.UnVec(&vecConstrEst);
//			matConn.Transpose();
//			printf("The estimated Observation matrix for lambda = %f\n",k * delta + dlambda[0]);
//			matConn.Print();
//			matConn.MultiplyV(&vecYhat,&ppvecx[efftime][effobs]);
//			printf("the estimated observation from the above matrix\n");
//			vecYhat.bColVector = false;
//			vecYhat.Print();
//			vecYhat.bColVector = true;
//			printf("the original observation \n");
//			ppvecy[efftime][effobs].bColVector = false;
//			ppvecy[efftime][effobs].Print();
//			ppvecy[efftime][effobs].bColVector = true;
//			//then we calc the difference between vecyy[][] and vecyy_hat and take its two norm
//			vecYhat.AddI(&ppvecy[efftime][effobs], 1, -1);
//			pMse[k] = vecYhat.Norm(2); //this is the metric for decision
//			printf("The estimated vector for lambda = %f and mse = %f is \n", k * delta, pMse[k]);
//			vecConstrEst.bColVector = false;
//			vecConstrEst.Print();
//			vecConstrEst.bColVector = true;
////			if(k != 0)
////			{
////				if(pMse[k] > pMse[k - 1])
////				{
////					minidx = k-1;
////					printf("Minima found \n");
////					break;
////				}
////			}
//		}
//		if(minidx == 0)
//		{
//			//
//
//			min = pMse[0];
//			minidx = 0;
//			for(int k = 1; k < numLambda;k++)
//			{
//				if(pMse[k] < min)
//				{
//					min = pMse[k];
//					minidx = k;
//				}
//			}
//		}
//		printf("Best coarse Lambda %f, corresponding MSE %f\n",minidx* delta + dlambda[0],pMse[minidx]);
//		delete pMse;
//
//		double newdelta = 0.01;
//		double minlambda = (minidx - 1) * delta + dlambda[0];
//		double maxlambda = (minidx + 1) * delta + dlambda[0];
//		numLambda = (maxlambda - minlambda) / newdelta;
//		pMse = new double[numLambda];
//		//the second search loop for fine grained search
//		printf("search for a finer lambda\n");
//		for(int k =0; k < numLambda;k++)
//		{
//			CVector vecConstrEst(nGenes * nGenes,1);
//			CMatrix matConn(nGenes,nGenes);
//			CVector vecYhat(nGenes,1);
//
//
//			CLasso lasso;
//			//(pEstA[j],&vecConstrEst,k * delta + dlambda[0]);
//			lasso.DoLasso(pEstA[j],&vecConstrEst,k * newdelta + minlambda);
//
//			//calculate and store the error and lambda of this lasso estimate
//
//			//to calc the error we need to calculate the product estmatrix * vecxx[][] = vecyy_hat
//			//unvec the connectivity matrix
//			matConn.UnVec(&vecConstrEst);
//			matConn.Transpose();
//			printf("The estimated Observation matrix for lambda = %f\n",k * newdelta + minlambda);
//			matConn.Print();
//			matConn.MultiplyV(&vecYhat,&ppvecx[efftime][effobs]);
//			printf("the estimated observation from the above matrix\n");
//			vecYhat.bColVector = false;
//			vecYhat.Print();
//			vecYhat.bColVector = true;
//			printf("the original observation \n");
//			ppvecy[efftime][effobs].bColVector = false;
//			ppvecy[efftime][effobs].Print();
//			ppvecy[efftime][effobs].bColVector = true;
//			//then we calc the difference between vecyy[][] and vecyy_hat and take its two norm
//			vecYhat.AddI(&ppvecy[efftime][effobs], 1, -1);
//			pMse[k] = vecYhat.Norm(2); //this is the metric for decision
//			printf("The estimated matrix for lambda = %f and mse = %f is \n", k * newdelta + minlambda, pMse[k]);
//
//
//		}
//		min = pMse[0];
//		minidx = 0;
//
//		for(int k = 1; k < numLambda;k++)
//		{
//			if(pMse[k] < min)
//			{
//				min = pMse[k];
//				minidx = k;
//			}
//		}
//
//		printf("Best fine lambda %f, corresponding Mse %f \n",minidx * newdelta + minlambda,pMse[minidx]);
//		//find the optimum estimate (again) and overwrite the unconstrained one
//		CLasso lasso;
//		CMatrix matEst(nGenes,nGenes);
//		arrayLambda[j] = minidx * newdelta + minlambda;
//		lasso.DoLasso(pEstA[j],pEstA[j],minidx * newdelta + minlambda);
//		printf("Unsmoothed Estimate at time index %d \n",j);
//		matEst.UnVec(pEstA[j]);
//		matEst.Print();
//	}
//	//initialize the smooth estimates at the final times
//
//	//the estimated covariance at the final time is the smoothed covariance
//	pEstSmtA[nEffMeas - 1]->Copy(pEstA[nEffMeas - 1]);
//	// the estimated state at the final time is the smoothed estimate
//	mat3dbarPkk.pMats[nEffMeas - 1].Copy(&(mat3dbarPkN.pMats[nEffMeas - 1]));
//
//
//	for(int j = nEffMeas - 2;j >= 1;j--)
//	{
//		CMatrix temp(nGenes,nGenes);
//		CMatrix temp2(nGenes,nGenes);
//
//		CMatrix barAk(nGenes,nGenes);
//		printf("Pkp1k at time %d \n",j);
//		mat3dbarPkp1k.pMats[j].Print();
//
//		printf("Pkk at time %d \n",j);
//		mat3dbarPkk.pMats[j].Print();
//
//
//
//		mat3dbarPkp1k.pMats[j].InvMat(&temp);
//
//		printf("Pkp1k inverse \n");
//		temp.Print();
//		mat3dbarPkk.pMats[j].MMMul(&temp,&barAk);
//
//		printf("The Ak matrix at time %d \n",j);
//		barAk.Print();
//
//		CSparseMatrix mat;
//		CVector vTmp(nGenes * nGenes,1);
//		vTmp.Copy(pEstSmtA[j+1]);
//		vTmp.SubtractI(pEstAApr[j+1]);
//		mat.KronEyeMatVec(nGenes,&barAk,&vTmp,pEstSmtA[j]);
//		// TODO dealloc vTmp here
//		pEstSmtA[j]->AddI(pEstA[j]);
//		CMatrix tmp(nGenes,nGenes);
//		tmp.UnVec(pEstSmtA[j]);
//		printf("The smoothed estimate at time %d \n",j);
//		tmp.Print();
//
//
//		//update the covariance
//		temp.Add(mat3dbarPkN.pMats[j+1],mat3dbarPkp1k.pMats[j],1,-1);
//
//		barAk.MMMul(&temp,&temp2);
////		temp2.MMMul(&barAk,&temp,1,1,false,true);
//		assert(false);//the above multiplication function has been changed
//		mat3dbarPkN.pMats[j]->Add(mat3dbarPkk.pMats[j],&temp);
//
//		printf("THe smoothed covariance at time %d \n",j);
//		mat3dbarPkN.pMats[j]->Print();
//
////		CLasso lasso;
////		//(pEstSmtA[j],pEstSmtA[j],dlambda);
////		lasso.DoLasso(pEstSmtA[j],pEstSmtA[j],arrayLambda[j]);
//
//		//make the small values in the estimate equal to zero
//		pEstSmtA[j]->Threshold(0.5);
//
//		tmp.UnVec(pEstSmtA[j]);
//		printf("The Reconstrained smoothed estimate at time %d \n",j);
//		tmp.Print();
//	}
//
}
MKL_INT CKalmanFilter::Estimate(CMatrix* pMatX, CMatrix* pMatY,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dLambda)
{
	//Computes the Kalman estimate based on each row of the connectivity matrix (9th May 2013)
//form the input (X) at one time into a matrix, with dimension nObs x nGenes

	int numprocs, rank, namelen, id,mpirank;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Get_processor_name(processor_name,&namelen);


	nGenes = ngenes;
	nEffMeas = ntimepts*nobs + 1;
	ntimepts++; // TODO this is dangerous

	nTimePts = ntimepts;


// TODO should all the genes have the same evolution noise ?
	CMatrix matSigmaQ(nGenes,nGenes); // Noise covariance of the state dynamics
	matSigmaQ.Eye();
	matSigmaQ.Scale(5.5);

	// TODO should all the genes have the same observation noise
	CMatrix matSigmaR(nobs,nobs);
	matSigmaR.Eye();
	matSigmaR.Scale(0.2);


	matPkN = NULL;
	matPkk = NULL;
	matPkp1k = NULL;

	matPkk = new CMatrix*[ntimepts];
	matPkp1k = new CMatrix*[ntimepts];
	matPkN = new CMatrix*[ntimepts];

	for(int i =0;i < ntimepts;i++)
	{
		matPkk[i] = new CMatrix[nGenes];
		matPkp1k[i] = new CMatrix[nGenes];
		matPkN[i] = new CMatrix[nGenes];
		for(int j = 0; j < nGenes;j++)
		{
			matPkk[i][j].Create(nGenes,nGenes);
			matPkp1k[i][j].Create(nGenes,nGenes);
			matPkN[i][j].Create(nGenes,nGenes);
		}
	}



	pEstA = NULL;
	pEstAApr = NULL;
	pEstSmtA = NULL;

	pEstA = new CVector*[ntimepts];
	pEstAApr = new CVector*[ntimepts];
	pEstSmtA = new CVector*[ntimepts];
	for(int i=0; i < ntimepts;i++)
	{
		pEstA[i] = new CVector[nGenes];
		pEstAApr[i] = new CVector[nGenes];
		pEstSmtA[i] = new CVector[nGenes];
		for(int j =0; j < nGenes;j++)
		{
			pEstA[i][j].Create(nGenes,1);
			pEstAApr[i][j].Create(nGenes,1);
			pEstSmtA[i][j].Create(nGenes,1);
		}
	}

	printf("Realization %d Successfully Allocated memory\n",mpirank);
	//initialize the estimate at time zero from the inital condition
	for(int j = 0;j < nGenes;j++)
	{
		pEstA[0][j].Copy(&(ic[j])); //initialize the state variable with the initial condition
	}


	CMatrix* matKalmanGain = new CMatrix[nGenes];
	for(int j =0; j < nGenes;j++)
	{
		matKalmanGain[j].Create(nGenes,nobs);
	}

	int efftime= 0;
	int effobs = 0;

	//for each time k
	for(int i = 1;i < ntimepts;i++)
	{
//		printf("Kalman Filter at time %d\n",i);
		CMatrix H(nobs,nGenes);
//		pMatX[i-1]->Print();
		pMatX[i-1].Copy(&H,true);
//		H.Print();
		CVector vecA(nGenes,1);


		//for each row of the connectivity matrix (that is for each gene) //
#pragma omp parallel for
		for(int j =0;j < nGenes;j++)
		{
			CVector curObs(nobs,1);
			CMatrix temp1(nGenes,nobs);
			CMatrix temp2(nobs,nobs);
			CMatrix temp3(nobs,nobs);
			CMatrix temp4(nGenes,nGenes);
			CMatrix temp5(nGenes,nGenes);
			CVector vectemp1(nobs,1);
			CVector vectemp2(nobs,1);
			CVector vectemp3(nGenes,1);



			// we update the state dynamics   //update the row dynamics
			pEstAApr[i][j].Copy(&pEstA[i-1][j]);
//			printf("The apriori state estimate time:%d and row:%d\n",i,j);
//			pEstAApr[i][j].Print(true);

			//the apriori error covariance
			matPkp1k[i][j].Add(&matPkk[i - 1][j],&matSigmaQ); // P^{-1}_k = P^+_{k-1} + Q_{k-1}
//			printf("The apriori covariance at time:%d and row:%d\n",i,j);
//			matPkp1k[i][j].Print();
			//calculate the Kalman Gain vector
//			 K_k = P^-_k H^T_k ( H_k * P^-_k * H^t_k + R_k )^{-1}
			//first calc the inverse term
//			printf("calculating the Kalman Gain at time:%d for row:%d\n",i,j);
			//temp1 = P^-_k * H^t_k
			matPkp1k[i][j].MMMul(&(pMatX[i-1]),&temp1);//temp1 is used twice
//			printf("P^-_k\n");
//			matPkp1k[i][j].Print();
//			printf("H^t_k\n");
//			pMatX[i-1]->Print();
//			printf("temp1 = P^-_k * H^t_k\n");
//			temp1.Print();
			//temp2 = H_k * P^-_k * H^t_k
//			printf("H \n");
//			H.Print();
			H.MMMul(&temp1,&temp2);
//			printf("temp2 = H_k * P^-_k * H^t_k\n");
//			temp2.Print();
			//temp3 = H_k * P^-_k * H^t_k + R_k
			temp3.Add(&temp2,&matSigmaR);
//			printf("temp3 = H_k * P^-_k * H^t_k + R_k\n");
//			temp3.Print();
			//temp2 = ( H_k * P^-_k * H^t_k + R_k )^{-1}
			temp3.InvMat(&temp2);
//			printf("temp2 = ( H_k * P^-_k * H^t_k + R_k )^{-1}\n");
//			temp2.Print();

			//K_k = P^-_k H^T_k ( H_k * P^-_k * H^t_k + R_k )^{-1}
			temp1.MMMul(&temp2,&matKalmanGain[j]);
//			printf("The kalman gain \n");
//			matKalmanGain[j].Print();

			//aposteriori estimate
			//\hat{x}^+_k = \hat{x}^-_k + k_k * ( y_k - H_k * \hat{x}^-_k )
//			vectemp1 = H_k * \hat{x}^-_k
			pMatX[i-1].MultiplyV(&vectemp1,&pEstAApr[i][j],true);
//			printf("vectemp1 = H_k * \\hat{x}^-_k\n");
//			vectemp1.Print(true);
			//vectemp2 = y_k = Y(j,:)
			pMatY[i-1].GetVec(&curObs,j,true);//extract the row frm the expression matrix
			curObs.bColVector = true;//transpose it
			vectemp2.Copy(&curObs);

//			printf("vectemp2 = y_k = Y(j,:)\n");
//			vectemp2.Print(true);
			//vectemp2 = y_k - H_k * \hat{x}^-_k
			vectemp2.AddI(&vectemp1,1,-1);
//			printf("vectemp2 = y_k - H_k * \\hat{x}^-_k\n");
//			vectemp2.Print(true);
			//vectemp3 = k_k * ( y_k - H_k * \hat{x}^-_k )
			matKalmanGain[j].MultiplyV(&vectemp3,&vectemp2);
//			printf("vectemp3 = k_k * ( y_k - H_k * \\hat{x}^-_k )\n");
//			vectemp3.Print();
			//hat{x}^+_k = \hat{x}^-_k
			pEstA[i][j].Copy(&pEstAApr[i][j]);
			pEstA[i][j].AddI(&vectemp3);
//			printf("\\hat{x}^+_k = \\hat{x}^-_k + vectemp3\n");
//			pEstA[i][j].Print(true);

			//we update the covariance using the Joseph Form
			// P^+_k = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k
			//temp4 = K_k * H_k
			matKalmanGain[j].MMMul(&H,&temp4);
//			printf("temp4 = K_k * H_k\n");
//			temp4.Print();
			//temp4 = ( I - K_k * H_k)
			temp4.EyeSubMat();
//			printf("( I - K_k * H_k)\n");
//			temp4.Print();
			//P^+_k = ( I - K_k * H_k) * P^-_k
			temp4.MMMul(&matPkp1k[i][j],&matPkk[i][j]);
//			printf("P^+_k = ( I - K_k * H_k) * P^-_k \n");
//			matPkk[i][j].Print();
			//temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T
			matPkk[i][j].MMMul(&temp4,&temp5,1,1,false,true);
//			printf("temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T\n");
//			temp5.Print();
			//temp1 = K_k * R_k
			matKalmanGain[j].MMMul(&matSigmaR,&temp1);
//			printf("temp1 = K_k * R_k \n");
//			temp1.Print();
			//temp5 = temp1 * K^T_k + temp5
			//temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k
			temp1.MMMul(&matKalmanGain[j],&temp5,1,1,false,true);
//			printf("temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k\n");
//			temp5.Print();
			temp5.Copy(&matPkk[i][j]);
			//constraint the solution using lasso

			//update the row estimate with the constrained one
			CMatrix mat(nGenes,nGenes);
			mat.Eye();
//			lasso.DoLasso(&pEstA[i][j],&pEstA[i][j],arrLambda[j]); //faulty implementation of the Lasso

			//now we search for the best lambda
			dLambda[1] = 5;
			MKL_INT numLambda = (dLambda[1] - dLambda[0]) / DELTALAMBDA;
			double dLargestLambda = 0;
			double* gcv = new double[numLambda];
			gcv[0] = INFINITY;

#pragma omp parallel for schedule(dynamic)
			for(int m=1; m < numLambda;m++)
			{
				double lambda = dLambda[0] + (m * DELTALAMBDA);
				double pLambda = 0;
				CVector beta(nGenes,1);
				CVector ynew(nobs,1);	CMatrix** matPkk = NULL;
				CMatrix** matPkp1k = NULL;
				CMatrix** matPkN = NULL;
				CLasso l(&mat,&(pEstA[i][j]),&beta,lambda,1,1,1000,0);
//				printf("beta: ");beta.Print(true);
				if(beta.IsZero() && dLargestLambda == 0)
				{
					dLargestLambda = lambda;
					gcv[m] = INFINITY;
				}
				else if( dLargestLambda == 0)
				{
					pLambda = 	beta.Cardinality();
					H.MultiplyV(&ynew,&beta);
					ynew.AddI(&curObs,1,-1);

					gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
				}
				else
					gcv[m] = INFINITY;

			}

			//find the minimum and then finetune the search
			MKL_INT minidx = 0;
			double mingcv = gcv[0];

			for(int m =1; m < numLambda;m++)
			{
				if(gcv[m] < mingcv)
				{
					minidx = m;
					mingcv = gcv[m];
				}
			}

			double minFineLambda = (minidx - 1) * DELTALAMBDA + dLambda[0];
			double maxFineLambda = (minidx + 1) * DELTALAMBDA + dLambda[0];
			MKL_INT numFineLambda = (maxFineLambda - minFineLambda) / DELTAFINELAMBDA;
			dLargestLambda = 0;

			delete [] gcv;

			gcv = new double[numFineLambda];

#pragma omp parallel for schedule(dynamic)
			for(int m = 0; m < numFineLambda;m++)
			{
				double lambda = minFineLambda + (m * DELTAFINELAMBDA);
				double pLambda = 0;
				CVector beta(nGenes,1);
				CVector ynew(nobs,1);
				CLasso l(&mat,&(pEstA[i][j]),&beta,lambda,1,1,1000,0);
				if(beta.IsZero() && dLargestLambda == 0)
				{
					dLargestLambda = lambda;
					gcv[m] = INFINITY;
				}
				else if( dLargestLambda == 0)
				{
					pLambda = 	beta.Cardinality();
					H.MultiplyV(&ynew,&beta);
					ynew.AddI(&curObs,1,-1);

					gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
				}
				else
					gcv[m] = INFINITY;
			}

			minidx = 0;
			mingcv = gcv[0];

			for(int m = 1; m < numFineLambda;m++)
			{
				if(gcv[m] < mingcv)
				{
					minidx = m;
					mingcv = gcv[m];
				}
			}


			delete [] gcv;
			CLasso lasso(&mat,&(pEstA[i][j]),&(pEstA[i][j]),minFineLambda + (minidx * DELTAFINELAMBDA),1,1,1000,0);
//			pEstA[i][j].Print(true);
		}
	}


	for(int j =0; j < nGenes;j++)
	{
		matKalmanGain[j].Release();
	}

	delete [] matKalmanGain;

	/* Smooth the estimate
	 *
	 *
	 */

	//initialize the smoother with the final result of the forward run
	/*
	 *
	 * PkN(:,:,:,end) = Pkk(:,:,:,end);
	 * aHatSmth(:,:,end) = aHatApos(:,:,end);
	 */
#pragma omp parallel for
	for(int j = 0;j < nGenes;j++)
	{

		matPkk[ntimepts - 1][j].Copy(&(matPkN[ntimepts - 1][j]));
		pEstSmtA[ntimepts - 1][j].Copy(&(pEstA[ntimepts - 1][j]));

	}
	//smoothing takes place backwards
	for(int k=(ntimepts - 2);k >= 1;k--)
	{
//		printf("Smoothing at time: %d\n",k);
		CMatrix H(nobs,nGenes);
//		pMatX[i-1]->Print();
		pMatX[k-1].Copy(&H,true);

		//the gene elements should / can be processed in parallel
#pragma omp parallel for
		for(int j = 0;j < nGenes;j++)
		{
			/*
			 * Ak = Pkk(:,:,k,j) / Pkp1k(:,:,k,j);
			 */
			CMatrix mat1(nGenes,nGenes);
			CMatrix mat2(nGenes,nGenes);
			CVector vecTemp1(nGenes,1);

			CMatrix Ak(nGenes,nGenes);
			CVector curObs(nobs,1);
			pMatY[k-1].GetVec(&curObs,j,true);//extract the row frm the expression matrix
			curObs.bColVector = true;//transpose it

			matPkp1k[k][j].InvMat(&mat1);
			matPkk[k][j].MMMul(&mat1,&Ak);

			/*
			 *
			 * aHatSmth(:,k,j) = aHatApos(:,k,j) + Ak * ( aHatSmth(:,k,j+1) - aHatApr(:,k,j));
			 */
			vecTemp1.Copy(&(pEstSmtA[k+1][j]));
			vecTemp1.AddI(&(pEstAApr[k][j]),1,-1);


			Ak.MultiplyV(&(pEstSmtA[k][j]),&vecTemp1);
			pEstSmtA[k][j].AddI(&(pEstA[k][j]));

			/*
			 *  PkN(:,:,k,j) = Pkk(:,:,k,j) + Ak * ( PkN(:,:,k,j+1) - Pkp1k(:,:,k,j) ) * Ak';
			 */

			mat1.Add(&(matPkN[k+1][j]),&(matPkp1k[k][j]),1,-1);
			Ak.MMMul(&mat1,&mat2);

			matPkk[k][j].Copy(&(matPkN[k][j]));

			mat2.MMMul(&Ak,&(matPkN[k][j]),1,1,false,true);


			MKL_INT numLambda = (dLambda[1] - dLambda[0]) / DELTALAMBDA;
			double dLargestLambda = 0;
			double* gcv = new double[numLambda];
			CMatrix mat(nGenes,nGenes);
			mat.Eye();

			gcv[0] = INFINITY;
#pragma omp parallel for schedule(dynamic)
			for(int m=1; m < numLambda;m++)
			{
				double lambda = dLambda[0] + (m * DELTALAMBDA);
				double pLambda = 0;
				CVector beta(nGenes,1);
				CVector ynew(nobs,1);
				CLasso l(&mat,&(pEstSmtA[k][j]),&beta,lambda,1,1,1000,0);
//				printf("beta: ");beta.Print(true);
				if(beta.IsZero() && dLargestLambda == 0)
				{
					dLargestLambda = lambda;
					gcv[m] = INFINITY;
				}
				else if( dLargestLambda == 0)
				{
					pLambda = 	beta.Cardinality();
					H.MultiplyV(&ynew,&beta);
					ynew.AddI(&curObs,1,-1);

					gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
				}
				else
					gcv[m] = INFINITY;

			}

			//find the minimum and then finetune the search
			MKL_INT minidx = 0;
			double mingcv = gcv[0];

			for(int m =1; m < numLambda;m++)
			{
				if(gcv[m] < mingcv)
				{
					minidx = m;
					mingcv = gcv[m];
				}
			}

			double minFineLambda = (minidx - 1) * DELTALAMBDA + dLambda[0];
			double maxFineLambda = (minidx + 1) * DELTALAMBDA + dLambda[0];
			MKL_INT numFineLambda = (maxFineLambda - minFineLambda) / DELTAFINELAMBDA;
			dLargestLambda = 0;

			delete [] gcv;

			gcv = new double[numFineLambda];

#pragma omp parallel for schedule(dynamic)
			for(int m = 0; m < numFineLambda;m++)
			{
				double lambda = minFineLambda + (m * DELTAFINELAMBDA);
				double pLambda = 0;
				CVector beta(nGenes,1);
				CVector ynew(nobs,1);
				CLasso l(&mat,&(pEstSmtA[k][j]),&beta,lambda,1,1,1000,0);
				if(beta.IsZero() && dLargestLambda == 0)
				{
					dLargestLambda = lambda;
					gcv[m] = INFINITY;
				}
				else if( dLargestLambda == 0)
				{
					pLambda = 	beta.Cardinality();
					H.MultiplyV(&ynew,&beta);
					ynew.AddI(&curObs,1,-1);

					gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
				}
				else
					gcv[m] = INFINITY;
			}

			minidx = 0;
			mingcv = gcv[0];

			for(int m = 1; m < numFineLambda;m++)
			{
				if(gcv[m] < mingcv)
				{
					minidx = m;
					mingcv = gcv[m];
				}
			}
			delete [] gcv;

			CLasso lasso(&mat,&(pEstSmtA[k][j]),&(pEstSmtA[k][j]),minFineLambda + (minidx * DELTAFINELAMBDA),1,1,1000,0);
//			pEstSmtA[k][j].Print(true);


		}
	}

	//conver the vector arrays to matrices for further processing / saving

	if(matSmthEst == NULL)
	{
		matSmthEst = new CMatrix[ntimepts - 1];
		matEstA = new CMatrix[ntimepts - 1];
		matEstAApr = new CMatrix[ntimepts - 1];

		for(int m = 0;m < nGenes;m++)
		{
			pEstSmtA[0][m].Destroy();
			pEstA[0][m].Destroy();
			pEstAApr[0][m].Destroy();
		}
		for(int k =  0;k < (ntimepts - 1);k++)
		{
			matSmthEst[k].Create(nGenes,nGenes);
			matEstA[k].Create(nGenes,nGenes);
			matEstAApr[k].Create(nGenes,nGenes);
			for(int j =0; j < nGenes;j++)
			{
				matSmthEst[k].SetVec(&(pEstSmtA[k + 1][j]),j,true);
				pEstSmtA[k + 1][j].Destroy();

				matEstA[k].SetVec(&(pEstA[k+1][j]),j,true);
				pEstA[k+1][j].Destroy();

				matEstAApr[k].SetVec(&(pEstAApr[k+1][j]),j,true);
				pEstAApr[k+1][j].Destroy();
			}
//			printf("Estimated Smooth matrix at time %d \n",k + 1);
//			matSmthEst[k].Print();
//
//			printf("Estimated Apos matrix at time %d \n",k + 1);
//			matEstA[k].Print();
//
//			printf("Estimated Apr matrix at time %d \n",k + 1);
//			matEstAApr[k].Print();

		}
	}


	return 0;
}
MKL_INT CKalmanFilter::EstimateGene(CMatrix* pMatX, CVector* vecY,MKL_INT ngenes,MKL_INT ntimepts,MKL_INT nobs, CVector* ic,double* dLambda,MKL_INT curgene)
{
	//Computes the Kalman estimate based on each row of the connectivity matrix (9th May 2013)
//form the input (X) at one time into a matrix, with dimension nObs x nGenes


	nGenes = ngenes; //the total number of genes
//	nEffMeas = ntimepts*nobs + 1;
	ntimepts++; // TODO this is dangerous

	nTimePts = ntimepts;


// TODO should all the genes have the same evolution noise ?
	CMatrix matSigmaQ(nGenes,nGenes); // Noise covariance of the state dynamics
	matSigmaQ.Eye();
	matSigmaQ.Scale(5.5);

	// TODO should all the genes have the same observation noise
	CMatrix matSigmaR(nobs,nobs);
	matSigmaR.Eye();
	matSigmaR.Scale(0.2);
	matPkN = NULL;
	matPkk = NULL;
	matPkp1k = NULL;
	pEstA = NULL;
	pEstAApr = NULL;
	pEstSmtA = NULL;


	CMatrix* mPkN = NULL;
	CMatrix* mPkk = NULL;
	CMatrix* mPkp1k = NULL;

	mPkk = new CMatrix[ntimepts];
	mPkp1k = new CMatrix[ntimepts];
	mPkN = new CMatrix[ntimepts];

	for(int i =0;i < ntimepts;i++)
	{
		if(mPkk[i].Create(nGenes,nGenes) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}
		if(mPkp1k[i].Create(nGenes,nGenes) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}
		if(mPkN[i].Create(nGenes,nGenes) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}
	}



	pgEstA = NULL;
	pgEstAApr = NULL;
	pgEstSmtA = NULL;

	pgEstA = new CVector[ntimepts];
	pgEstAApr = new CVector[ntimepts];
	pgEstSmtA = new CVector[ntimepts];
	for(int i=0; i < ntimepts;i++)
	{
		if(pgEstA[i].Create(nGenes,1) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}
		if(pgEstAApr[i].Create(nGenes,1) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}
		if(pgEstSmtA[i].Create(nGenes,1) == -1)
		{
			eprintf("Create Failed");
			return -1;
		}

	}

	//initialize the estimate at time zero from the inital condition
	pgEstA[0].Copy(ic); //initialize the state variable with the initial condition


	CMatrix matKalmanGain;
	if(matKalmanGain.Create(nGenes,nobs) == -1)
	{
		eprintf("Create Failed");
		return -1;
	}

	int efftime= 0;
	int effobs = 0;
	char str[256] = {0};

	//for each time k
	for(int i = 1;i < ntimepts;i++)
	{
		sprintf(str,"Kalman Filter at time %d\n",i);
		eprintf(str);
		CMatrix H(nobs,nGenes);
		eprintf("declared H");
//		pMatX[i-1]->Print();
		pMatX[i-1].Copy(&H,true);
//		H.Print();
		CVector vecA(nGenes,1);


		//for each row of the connectivity matrix (that is for each gene) //

		CVector curObs(nobs,1);
		CMatrix temp1(nGenes,nobs);
		CMatrix temp2(nobs,nobs);
		CMatrix temp3(nobs,nobs);
		CMatrix temp4(nGenes,nGenes);
		CMatrix temp5(nGenes,nGenes);
		CVector vectemp1(nobs,1);
		CVector vectemp2(nobs,1);
		CVector vectemp3(nGenes,1);



		// we update the state dynamics   //update the row dynamics
		pgEstAApr[i].Copy(&pgEstA[i-1]);
//			printf("The apriori state estimate time:%d and row:%d\n",i,j);
//			pEstAApr[i][j].Print(true);
		eprintf("Dynamics Updated");
		//the apriori error covariance
		mPkp1k[i].Add(&mPkk[i - 1],&matSigmaQ); // P^{-1}_k = P^+_{k-1} + Q_{k-1}
//			printf("The apriori covariance at time:%d and row:%d\n",i,j);
//			matPkp1k[i][j].Print();
		//calculate the Kalman Gain vector
//			 K_k = P^-_k H^T_k ( H_k * P^-_k * H^t_k + R_k )^{-1}
		//first calc the inverse term
//			printf("calculating the Kalman Gain at time:%d for row:%d\n",i,j);
		//temp1 = P^-_k * H^t_k
		mPkp1k[i].MMMul(&(pMatX[i-1]),&temp1);//temp1 is used twice
//			printf("P^-_k\n");
//			matPkp1k[i][j].Print();
//			printf("H^t_k\n");
//			pMatX[i-1]->Print();
//			printf("temp1 = P^-_k * H^t_k\n");
//			temp1.Print();
		//temp2 = H_k * P^-_k * H^t_k
//			printf("H \n");
//			H.Print();
		H.MMMul(&temp1,&temp2);
//			printf("temp2 = H_k * P^-_k * H^t_k\n");
//			temp2.Print();
		//temp3 = H_k * P^-_k * H^t_k + R_k
		temp3.Add(&temp2,&matSigmaR);
//			printf("temp3 = H_k * P^-_k * H^t_k + R_k\n");
//			temp3.Print();
		//temp2 = ( H_k * P^-_k * H^t_k + R_k )^{-1}
		temp3.InvMat(&temp2);
//			printf("temp2 = ( H_k * P^-_k * H^t_k + R_k )^{-1}\n");
//			temp2.Print();

		//K_k = P^-_k H^T_k ( H_k * P^-_k * H^t_k + R_k )^{-1}
		temp1.MMMul(&temp2,&matKalmanGain);
//			printf("The kalman gain \n");
//			matKalmanGain[j].Print();

		//aposteriori estimate
		//\hat{x}^+_k = \hat{x}^-_k + k_k * ( y_k - H_k * \hat{x}^-_k )
//			vectemp1 = H_k * \hat{x}^-_k
		pMatX[i-1].MultiplyV(&vectemp1,&pgEstAApr[i],true);
//			printf("vectemp1 = H_k * \\hat{x}^-_k\n");
//			vectemp1.Print(true);
		//vectemp2 = y_k = Y(j,:)
//		pMatY[i-1].GetVec(&curObs,j,true);//extract the row frm the expression matrix
		curObs.Copy(&(vecY[i-1]));
		curObs.bColVector = true;//transpose it
		vectemp2.Copy(&curObs);

//			printf("vectemp2 = y_k = Y(j,:)\n");
//			vectemp2.Print(true);
		//vectemp2 = y_k - H_k * \hat{x}^-_k
		vectemp2.AddI(&vectemp1,1,-1);
//			printf("vectemp2 = y_k - H_k * \\hat{x}^-_k\n");
//			vectemp2.Print(true);
		//vectemp3 = k_k * ( y_k - H_k * \hat{x}^-_k )
		matKalmanGain.MultiplyV(&vectemp3,&vectemp2);
//			printf("vectemp3 = k_k * ( y_k - H_k * \\hat{x}^-_k )\n");
//			vectemp3.Print();
		//hat{x}^+_k = \hat{x}^-_k
		pgEstA[i].Copy(&pgEstAApr[i]);
		pgEstA[i].AddI(&vectemp3);
//			printf("\\hat{x}^+_k = \\hat{x}^-_k + vectemp3\n");
//			pEstA[i][j].Print(true);

		//we update the covariance using the Joseph Form
		// P^+_k = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k
		//temp4 = K_k * H_k
		matKalmanGain.MMMul(&H,&temp4);
//			printf("temp4 = K_k * H_k\n");
//			temp4.Print();
		//temp4 = ( I - K_k * H_k)
		temp4.EyeSubMat();
//			printf("( I - K_k * H_k)\n");
//			temp4.Print();
		//P^+_k = ( I - K_k * H_k) * P^-_k
		temp4.MMMul(&mPkp1k[i],&mPkk[i]);
//			printf("P^+_k = ( I - K_k * H_k) * P^-_k \n");
//			matPkk[i][j].Print();
		//temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T
		mPkk[i].MMMul(&temp4,&temp5,1,1,false,true);
//			printf("temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T\n");
//			temp5.Print();
		//temp1 = K_k * R_k
		matKalmanGain.MMMul(&matSigmaR,&temp1);
//			printf("temp1 = K_k * R_k \n");
//			temp1.Print();
		//temp5 = temp1 * K^T_k + temp5
		//temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k
		temp1.MMMul(&matKalmanGain,&temp5,1,1,false,true);
//			printf("temp5 = ( I - K_k * H_k) * P^-_k * (I - K_k * H_k)^T + K_k * R_k * K^T_k\n");
//			temp5.Print();
		temp5.Copy(&mPkk[i]);
		//constraint the solution using lasso

		//update the row estimate with the constrained one
		CMatrix mat(nGenes,nGenes);
		mat.Eye();
//			lasso.DoLasso(&pEstA[i][j],&pEstA[i][j],arrLambda[j]); //faulty implementation of the Lasso

		//now we search for the best lambda
		dLambda[1] = 5;
		MKL_INT numLambda = (dLambda[1] - dLambda[0]) / DELTALAMBDA;
		double dLargestLambda = 0;
		double* gcv = new double[numLambda];
		gcv[0] = INFINITY;

#pragma omp parallel for schedule(dynamic)
		for(int m=1; m < numLambda;m++)
		{
			double lambda = dLambda[0] + (m * DELTALAMBDA);
			double pLambda = 0;
			CVector beta(nGenes,1);
			CVector ynew(nobs,1);
			CLasso l(&mat,&(pgEstA[i]),&beta,lambda,1,1,100,0);
//				printf("beta: ");beta.Print(true);
			if(beta.IsZero() && dLargestLambda == 0)
			{
				dLargestLambda = lambda;
				gcv[m] = INFINITY;
			}
			else if( dLargestLambda == 0)
			{
				pLambda = 	beta.Cardinality();
				H.MultiplyV(&ynew,&beta);
				ynew.AddI(&curObs,1,-1);

				gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
			}
			else
				gcv[m] = INFINITY;

		}

		//find the minimum and then finetune the search
		MKL_INT minidx = 0;
		double mingcv = gcv[0];

		for(int m =1; m < numLambda;m++)
		{
			if(gcv[m] < mingcv)
			{
				minidx = m;
				mingcv = gcv[m];
			}
		}

		double minFineLambda = (minidx - 1) * DELTALAMBDA + dLambda[0];
		double maxFineLambda = (minidx + 1) * DELTALAMBDA + dLambda[0];
		MKL_INT numFineLambda = (maxFineLambda - minFineLambda) / DELTAFINELAMBDA;
		dLargestLambda = 0;

		delete [] gcv;

		gcv = new double[numFineLambda];

#pragma omp parallel for schedule(dynamic)
		for(int m = 0; m < numFineLambda;m++)
		{
			double lambda = minFineLambda + (m * DELTAFINELAMBDA);
			double pLambda = 0;
			CVector beta(nGenes,1);
			CVector ynew(nobs,1);
			CLasso l(&mat,&(pgEstA[i]),&beta,lambda,1,1,100,0);
			if(beta.IsZero() && dLargestLambda == 0)
			{
				dLargestLambda = lambda;
				gcv[m] = INFINITY;
			}
			else if( dLargestLambda == 0)
			{
				pLambda = 	beta.Cardinality();
				H.MultiplyV(&ynew,&beta);
				ynew.AddI(&curObs,1,-1);

				gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
			}
			else
				gcv[m] = INFINITY;
		}

		minidx = 0;
		mingcv = gcv[0];

		for(int m = 1; m < numFineLambda;m++)
		{
			if(gcv[m] < mingcv)
			{
				minidx = m;
				mingcv = gcv[m];
			}
		}


		delete [] gcv;
		CLasso lasso(&mat,&(pgEstA[i]),&(pgEstA[i]),minFineLambda + (minidx * DELTAFINELAMBDA),1,1,100,0);
//		printf("Unsmoothes Estimate at time %d\n",i);
//		pgEstA[i].Print(true);

	}


	matKalmanGain.Release();

	/* Smooth the estimate
	 *
	 *
	 */

	//initialize the smoother with the final result of the forward run
	/*
	 *
	 * PkN(:,:,:,end) = Pkk(:,:,:,end);
	 * aHatSmth(:,:,end) = aHatApos(:,:,end);
	 */

	mPkk[ntimepts - 1].Copy(&(mPkN[ntimepts - 1]));
	pgEstSmtA[ntimepts - 1].Copy(&(pgEstA[ntimepts - 1]));

	//smoothing takes place backwards
	for(int k=(ntimepts - 2);k >= 1;k--)
	{
		printf("Smoothing at time: %d\n",k);
		CMatrix H(nobs,nGenes);
//		pMatX[i-1]->Print();
		pMatX[k-1].Copy(&H,true);

		//the gene elements should / can be processed in parallel
		/*
		 * Ak = Pkk(:,:,k,j) / Pkp1k(:,:,k,j);
		 */
		CMatrix mat1(nGenes,nGenes);
		CMatrix mat2(nGenes,nGenes);
		CVector vecTemp1(nGenes,1);

		CMatrix Ak(nGenes,nGenes);
		CVector curObs(nobs,1);
//		pMatY[k-1].GetVec(&curObs,j,true);//extract the row frm the expression matrix
		curObs.Copy(&(vecY[k]));
		curObs.bColVector = true;//transpose it

		mPkp1k[k].InvMat(&mat1);
		mPkk[k].MMMul(&mat1,&Ak);

		/*
		 *
		 * aHatSmth(:,k,j) = aHatApos(:,k,j) + Ak * ( aHatSmth(:,k,j+1) - aHatApr(:,k,j));
		 */
		vecTemp1.Copy(&(pgEstSmtA[k+1]));
		vecTemp1.AddI(&(pgEstAApr[k]),1,-1);


		Ak.MultiplyV(&(pgEstSmtA[k]),&vecTemp1);
		pgEstSmtA[k].AddI(&(pgEstA[k]));

		/*
		 *  PkN(:,:,k,j) = Pkk(:,:,k,j) + Ak * ( PkN(:,:,k,j+1) - Pkp1k(:,:,k,j) ) * Ak';
		 */

		mat1.Add(&(mPkN[k+1]),&(mPkp1k[k]),1,-1);
		Ak.MMMul(&mat1,&mat2);

		mPkk[k].Copy(&(mPkN[k]));

		mat2.MMMul(&Ak,&(mPkN[k]),1,1,false,true);


		MKL_INT numLambda = (dLambda[1] - dLambda[0]) / DELTALAMBDA;
		double dLargestLambda = 0;
		double* gcv = new double[numLambda];
		CMatrix mat(nGenes,nGenes);
		mat.Eye();

		gcv[0] = INFINITY;
#pragma omp parallel for schedule(dynamic)
		for(int m=1; m < numLambda;m++)
		{
			double lambda = dLambda[0] + (m * DELTALAMBDA);
			double pLambda = 0;
			CVector beta(nGenes,1);
			CVector ynew(nobs,1);
			CLasso l(&mat,&(pgEstSmtA[k]),&beta,lambda,1,1,100,0);
//				printf("beta: ");beta.Print(true);
			if(beta.IsZero() && dLargestLambda == 0)
			{
				dLargestLambda = lambda;
				gcv[m] = INFINITY;
			}
			else if( dLargestLambda == 0)
			{
				pLambda = 	beta.Cardinality();
				H.MultiplyV(&ynew,&beta);
				ynew.AddI(&curObs,1,-1);

				gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
			}
			else
				gcv[m] = INFINITY;

		}

		//find the minimum and then finetune the search
		MKL_INT minidx = 0;
		double mingcv = gcv[0];

		for(int m =1; m < numLambda;m++)
		{
			if(gcv[m] < mingcv)
			{
				minidx = m;
				mingcv = gcv[m];
			}
		}

		double minFineLambda = (minidx - 1) * DELTALAMBDA + dLambda[0];
		double maxFineLambda = (minidx + 1) * DELTALAMBDA + dLambda[0];
		MKL_INT numFineLambda = (maxFineLambda - minFineLambda) / DELTAFINELAMBDA;
		dLargestLambda = 0;

		delete [] gcv;

		gcv = new double[numFineLambda];

#pragma omp parallel for schedule(dynamic)
		for(int m = 0; m < numFineLambda;m++)
		{
			double lambda = minFineLambda + (m * DELTAFINELAMBDA);
			double pLambda = 0;
			CVector beta(nGenes,1);
			CVector ynew(nobs,1);
			CLasso l(&mat,&(pgEstSmtA[k]),&beta,lambda,1,1,100,0);
			if(beta.IsZero() && dLargestLambda == 0)
			{
				dLargestLambda = lambda;
				gcv[m] = INFINITY;
			}
			else if( dLargestLambda == 0)
			{
				pLambda = 	beta.Cardinality();
				H.MultiplyV(&ynew,&beta);
				ynew.AddI(&curObs,1,-1);

				gcv[m] = (ynew.Norm(2) / pow(( 1 - (pLambda / nGenes) ),2)) / nGenes;
			}
			else
				gcv[m] = INFINITY;
		}

		minidx = 0;
		mingcv = gcv[0];

		for(int m = 1; m < numFineLambda;m++)
		{
			if(gcv[m] < mingcv)
			{
				minidx = m;
				mingcv = gcv[m];
			}
		}
		delete [] gcv;

		CLasso lasso(&mat,&(pgEstSmtA[k]),&(pgEstSmtA[k]),minFineLambda + (minidx * DELTAFINELAMBDA),1,1,100,0);
//		printf("Final Estimate at time %d\n",k);
//		pgEstSmtA[k].Print(true);

	}

	//release the covariance matrices


	for(int i =0;i < ntimepts;i++)
	{
		mPkk[i].Release();
		mPkp1k[i].Release();
		mPkN[i].Release();
	}
	delete [] mPkk;
	delete [] mPkp1k;
	delete [] mPkN;
	mPkk = NULL;
	mPkp1k = NULL;
	mPkN = NULL;

	return 0;
}
MKL_INT CKalmanFilter::SmoothEstimate()
{
	/*
	 *
	 *
barPkN = zeros(n_Gene,n_Gene,n_Measurements_Effective);
barPkN(:,:,n_Measurements_Effective) = barPkk(:,:,end);


for i=(n_Measurements_Effective-1):-1:1

%     barPkp1k = E_Cov(:,:,i+1);
%     barPkk = E_Cov(:,:,i);
%     X_n_n = E_A(:,i);
%     X_n_1_N = ES_A(:,i+1);
    barAk = barPkk(:,:,i) / barPkp1k(:,:,i);
%     X_n_N = X_n_n + A_n * ( X_n_1_N - X_n_1_n );
    X_n_N = E_A(:,i) + (kron(eye(n_Gene),barAk) * ( ES_A(:,i+1) - E_AApr(:,i))); %E_A(:,i-1) is the apriori estimate for time i, since F = I
%     V_n_N = V_n_n + A_n * ( V_n_1_N - V_n_1_n ) * A_n';
    barPkN(:,:,i) = barPkk(:,:,i) + barAk * (barPkN(:,:,i+1) - barPkp1k(:,:,i)) *barAk.';


    % Projection of the space onto a sparse space
    cvx_solver sdpt3
    cvx_begin quiet
    variable X_n_n_Constr(n_Gene.^2)

    minimize ((X_n_n_Constr - X_n_N)'*(X_n_n_Constr - X_n_N) + Alpha.*sum(abs(X_n_n_Constr(:))));
%     subject to norm(X_n_n_Constr,1) < 100

    cvx_end

    X_n_N = X_n_n_Constr;

%     X_n_1_N = X_n_N;

    ES_A(:,i)=X_n_N;

end
	 */


	return 0;
}

CKalmanFilter::~CKalmanFilter()
{
	if(matPkk != NULL)
	{
		for(int i =0;i < nTimePts;i++)
		{
			for(int j = 0; j < nGenes;j++)
			{
				matPkk[i][j].Release();
			}
			delete [] matPkk[i];
		}
		delete [] matPkk;
	}

	if(matPkp1k != NULL)
	{
		for(int i =0;i < nTimePts;i++)
		{
			for(int j = 0; j < nGenes;j++)
			{
				matPkp1k[i][j].Release();
			}
			delete [] matPkp1k[i];
		}
		delete [] matPkp1k;
	}

	if(matPkN != NULL)
	{
		for(int i =0;i < nTimePts;i++)
		{
			for(int j = 0; j < nGenes;j++)
			{
				matPkN[i][j].Release();
			}
			delete [] matPkN[i];
		}
		delete [] matPkN;
	}

	if(pEstA != NULL)
	{
		for(int i =0;i < nTimePts ;i++)
		{
			delete [] pEstA[i];
		}
		delete [] pEstA;
	}

	if(pEstAApr != NULL)
	{
		for(int i =0;i < nTimePts ;i++)
		{
			delete [] pEstAApr[i];
		}
		delete [] pEstAApr;
	}


	if(pEstSmtA != NULL)
	{
		for(int i =0;i < nTimePts ;i++)
		{
			delete [] pEstSmtA[i];
		}
		delete [] pEstSmtA;
	}

	if(matSmthEst != NULL)
	{
		for(int i =0;i < (nTimePts - 1);i++)
		{
			matSmthEst[i].Release();
		}
		delete [] matSmthEst;
	}

	if(matEstA != NULL)
	{
		for(int i =0;i < (nTimePts - 1);i++)
		{
			matEstA[i].Release();
		}
		delete [] matEstA;
	}

	if(matEstAApr != NULL)
	{
		for(int i =0;i < (nTimePts - 1);i++)
		{
			matEstAApr[i].Release();
		}
		delete [] matEstAApr;
	}

	if(pgEstA != NULL)
	{
		for(int i = 0;i < nTimePts;i++)
		{
			pgEstA[i].Destroy();
		}
		delete [] pgEstA;
	}

	if(pgEstAApr != NULL)
	{
		for(int i = 0;i < nTimePts;i++)
		{
			pgEstAApr[i].Destroy();
		}
		delete [] pgEstAApr;
	}

	if(pgEstSmtA != NULL)
	{
		for(int i = 0;i < nTimePts;i++)
		{
			pgEstSmtA[i].Destroy();
		}
		delete [] pgEstSmtA;
	}
}

