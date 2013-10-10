/*
 * CLasso.cpp
 *
 *  Created on: May 4, 2013
 *      Author: root
 */

#include "CLasso.h"
#define DELTALAMBDA 0.1

#include <omp.h>
CLasso::CLasso()
{
	iDim = 0;
	dLambda = 1;
	iPathLen = 1;
	dThresh = 1e-5;
	iMaxIter = 100;
	iVerbose = 1;
}
CLasso::CLasso(CMatrix* mata, CVector* vecb, CVector* vecxhat,double lambda,double rho,double alpha,MKL_INT maxiter,MKL_INT iverbose)
{
	iDim = 0;
	dLambda = 1;
	iPathLen = 1;
	dThresh = 1e-5;
	iMaxIter = maxiter;
	double abstol = 1e-4;
	double reltol = 1e-2;


	CVector vecAtb(mata->cols,1);

	mata->MultiplyV(&vecAtb,vecb,true);
//	printf("vecAtb\n");
//	vecAtb.Print();


	//compute the cholesky factors of the matrix A
	//IN our case it is going to be an eye mat

	CMatrix Lo;//(mata->rows,mata->cols);
	CMatrix Up;//(mata->rows,mata->cols);
	CMatrix InvUp;
	CMatrix InvLo;


	if(mata->rows > mata->cols) // skinny
	{

		/*
		 * L = chol( A'*A + rho*speye(n), 'lower' );
		 */

		CMatrix eye(mata->rows,mata->cols);
		Lo.Create(mata->rows,mata->cols);
		InvLo.Create(mata->rows,mata->rows);
		Up.Create(mata->rows,mata->cols);
		InvUp.Create(mata->rows,mata->cols);

		eye.Eye();
		eye.Scale(rho);

		mata->Copy(&Lo);
		mata->MMMul(&Lo,&Up,1,0,true,false);
		Lo.Add(&Up,&eye);

		LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',Lo.cols,Lo.mat,Lo.cols);
		//zero out the upper part
		for(int i =0;i < Lo.rows;i++)
		{
			for(int j =0;j < Lo.cols;j++)
			{
				if(j <= i )
					continue;
				Lo.mat[i * Lo.cols + j] = 0;
			}
		}
		Lo.Copy(&Up,true);
		Lo.InvMat(&InvLo);
		Up.InvMat(&InvUp);
		Lo.Release();
		Up.Release();
//		printf("Lower Triangular matrix\n");
//		InvLo.Print();
	}
	else //fat
	{
		Lo.Create(mata->rows,mata->rows);
		InvLo.Create(mata->rows,mata->rows);
		Up.Create(mata->rows,mata->rows);
		InvUp.Create(mata->rows,mata->rows);

		CMatrix matAt(mata->cols,mata->rows);

//		printf("Matrix A\n");
//		mata->Print();
		Lo.Eye();
		mata->Copy(&matAt,true);
//		printf("second matrix\n");
//		matAt.Print();
		mata->MMMul(&matAt,&Lo,(1/rho),1,false,false);
//		Lo.Print();
		LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',Lo.cols,Lo.mat,Lo.cols);
		//zero out the upper part
		for(int i =0;i < Lo.rows;i++)
		{
			for(int j =0;j < Lo.cols;j++)
			{
				if(j <= i )
					continue;
				Lo.mat[i * Lo.cols + j] = 0;
			}
		}
//		printf("Lower Tri Mat \n");
//		Lo.Print();
		Lo.InvMat(&InvLo);
//		printf("InvLo\n");
//		InvLo.Print();
		Lo.Copy(&Up,true);
		Up.InvMat(&InvUp);
//		printf("InvUp\n");
//		InvUp.Print();
		Lo.Release();
		Up.Release();
	}

	//invert Lo and Up and save them

	if(iverbose)
	{
		printf("%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "iter","r norm", "eps pri", "s norm", "eps dual", "objective");
	}

	//ADMM Solver
	CVector vecX(mata->cols,1);
	CVector vecZ(mata->cols,1);
	CVector vecU(mata->cols,1);


	CVector histObjVal(maxiter,1);
	CVector histRNorm(maxiter,1);
	CVector histSNorm(maxiter,1);
	CVector histEpsPri(maxiter,1);
	CVector histEpsDual(maxiter,1);


	for(int k =0;k < maxiter;k++)
	{
		CVector vecQ(mata->cols,1);
		CVector vecTmp(mata->cols,1);
		CVector vecTmp2(mata->rows,1);
		CVector vecTmp3(mata->rows,1);
		CVector vecZold(mata->cols,1);
		CVector vecXhat(mata->cols,1);


		vecQ.Copy(&vecAtb); //q = Atb
		//tmp = z - u
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecU,1,-1);
		vecQ.AddI(&vecTmp,1,rho); // q = Atb + rho * (z - u)
//		printf("VecQ\n");
//		vecQ.Print();

		if(mata->rows > mata->cols ) // skinny matrix
		{
			//solve the cholesky factorization
			//x = U \ (L \ q)
			InvLo.MultiplyV(&vecTmp3,&vecQ);
			InvUp.MultiplyV(&vecX,&vecTmp3);
//			printf("Lasso X estimate\n");
//			vecX.Print(true);
		}
		else // fat matrix
		{
			//	x = q / rho - (A' * ( U \ (L \ (A * q))) / rho^2
			//vecTmp2 = A * q
			mata->MultiplyV( &vecTmp2, &vecQ);
//			printf("VecTemp2\n");
//			vecTmp2.Print();
			//vecTmp = U \ ( L \ vecTmp2)
			//vecTmp = inv(U) * ( inv(L) * vecTmp2)

			//vecTmp = inv(L) * vecTmp2
			InvLo.MultiplyV(&vecTmp3,&vecTmp2);
			InvUp.MultiplyV(&vecTmp2,&vecTmp3);
//			LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'L',Lo.cols,1,Lo.mat,Lo.cols,vecTmp2.data,1);
//			printf("VecTemp2\n");
//			vecTmp2.Print();
			//vecX = vecQ
			vecX.Copy(&vecQ);
//			printf("vecX\n");
//			vecX.Print();
			// vecQ = A' * vecTmp2
			mata->MultiplyV( &vecQ,&vecTmp2,true);
//			printf("vecQ\n");
//			vecQ.Print();
			vecX.AddI(&vecQ,1/rho,-1 / (rho * rho));
//			printf("vecX\n");
//			vecX.Print();
		}

//
//		printf("vecX after cholesky soln\n");
//		vecX.Print();
		//zold = z
		vecZold.Copy(&vecZ);
		//x_hat = alpha*x + (1- alpha) * zold
		vecXhat.Copy(&vecX);
		vecXhat.AddI(&vecZold,alpha,(1 - alpha));
//		printf("vecXhat after alpha update\n");
//		vecXhat.Print();
		//z = shrinkage(x_hat + u, lambda/rho);
		//vecTmp = x_hat
		vecTmp.Copy(&vecXhat);
		//tmp = tmp + u
		vecTmp.AddI(&vecU,1,1);
		//this updates z
		Shrinkage(&vecTmp,lambda / rho,&vecZ);
//		printf("vecZ after shrinkage\n");
//		vecZ.Print();


//	    % u-update
//	    u = u + (x_hat - z);
		vecTmp.Copy(&vecXhat); //tmp  = x_hat
		vecTmp.AddI(&vecZ,1,-1);// tmp = tmp - z
		vecU.AddI(&vecTmp);// u = u + tmp

//		printf("vecU after update\n");
//		vecU.Print();

		//find the objective value
		histObjVal.data[k] = Objective(mata,vecb,lambda,&vecX,&vecZ);
		//history.rnom = norm(x - z)
		vecTmp.Copy(&vecX);
		vecTmp.AddI(&vecZ,1,-1);
//		printf("x - z\n");
//		vecTmp.Print();
		histRNorm.data[k] = vecTmp.Norm(2);
		//history.s_norm = norm(-rho * (z - zold))
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecZold,-1 * rho,rho);
//		printf("-rho * ( z - zold)\n");
//		vecTmp.Print();
		histSNorm.data[k] = vecTmp.Norm(2);
		//    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		double normX = vecX.Norm(2);
		vecTmp.Copy(&vecZ);
		vecTmp.ScaleI(-1);
		double normZ = vecTmp.Norm(2);
		if(normX > normZ)
		{
			histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normX;
		}
		else
		{
			histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normZ;
		}

		//    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		vecTmp.Copy(&vecU);
		vecTmp.ScaleI(rho);
		histEpsDual.data[k] = sqrt(mata->cols) * abstol + reltol * vecTmp.Norm(2);

		if(iverbose)
		{
//			fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
//            history.r_norm(k), history.eps_pri(k), ...
//            history.s_norm(k), history.eps_dual(k), history.objval(k));
			printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
					histRNorm.data[k],histEpsPri.data[k],histSNorm.data[k],
					histEpsDual.data[k],histObjVal.data[k]);
		}

//		  if (history.r_norm(k) < history.eps_pri(k) && ...
//		       history.s_norm(k) < history.eps_dual(k))
//		         break;
//		    end
		if(histRNorm.data[k] < histEpsPri.data[k]  &&  histSNorm.data[k] < histEpsDual.data[k])
			break;
	}

	vecxhat->Copy(&vecZ);


}
MKL_INT CLasso::DoLasso(CVector* vecb, CVector* vecxhat,double lambda,double rho,double alpha,MKL_INT maxiter,MKL_INT iverbose)
{

//	iverbose = 1;
	double abstol = 1e-4;
	double reltol = 1e-2;

//	CVector vecAtb(vecb->length,1);

//	mata->MultiplyV(&vecAtb,vecb);
//	printf("vecAtb\n");
//	vecAtb.Print();


	//compute the cholesky factors of the matrix A
	//IN our case it is going to be an eye mat
//
//	CMatrix Lo(mata->rows,mata->cols);
//	CMatrix Up(mata->rows,mata->cols);
//	CMatrix eye(mata->rows,mata->cols);
	CSparseMatrix spmatLo(vecb->length,vecb->length);
	spmatLo.Eye(sqrt(rho + 1)); //scaled eye to rho and cholesky factor
	spmatLo.Inv(true);//inv the diagonal matrix and since this is diagonal so upper diagonal or lower does not make a difference
//
//	eye.Eye();
//	eye.Scale(rho);
//
//	mata->Copy(&Lo);
//	mata->MMMul(&Lo,&Up,1,0,true,false);
//	Lo.Add(&Up,&eye);
//
//	LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',Lo.cols,Lo.mat,Lo.cols);
//	Lo.Copy(&Up,true);
//	printf("Lower Triangular matrix\n");
//	Lo.Print();


	if(iverbose)
	{
		printf("%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "iter","r norm", "eps pri", "s norm", "eps dual", "objective");
	}

	//ADMM Solver
	CVector vecX(vecb->length,1);
	CVector vecZ(vecb->length,1);
	CVector vecU(vecb->length,1);


	CVector histObjVal(maxiter,1);
	CVector histRNorm(maxiter,1);
	CVector histSNorm(maxiter,1);
	CVector histEpsPri(maxiter,1);
	CVector histEpsDual(maxiter,1);


	for(int k =0;k < maxiter;k++)
	{
		CVector vecQ(vecb->length,1);
		CVector vecTmp(vecb->length,1);
		CVector vecZold(vecb->length,1);
		CVector vecXhat(vecb->length,1);


		vecQ.Copy(vecb); //q = Atb
		//tmp = z - u
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecU,1,-1);
		vecQ.AddI(&vecTmp,1,rho); // q = Atb + rho * (z - u)

//		if(mata->rows >= mata->cols )
//		{
//			//solve the cholesky factorization
//			//x = U \ (L \ q)
//			LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'L',mata->cols,1,Lo.mat,mata->cols,vecQ.data,1);
			spmatLo.MulMatVec(&vecQ,&vecTmp); // tmp = L \ q;
			spmatLo.MulMatVec(&vecTmp,&vecX); // vecX = U \ tmp;
//		}
//		else
//			assert(false);//
//
//		printf("vecX after cholesky soln\n");
//		vecX.Print();
		//zold = z
		vecZold.Copy(&vecZ);
		//x_hat = alpha*x + (1- alpha) * zold
		vecXhat.Copy(&vecX);
		vecXhat.AddI(&vecZold,alpha,(1 - alpha));
//		printf("vecXhat after alpha update\n");
//		vecXhat.Print();
		//z = shrinkage(x_hat + u, lambda/rho);
		//vecTmp = x_hat
		vecTmp.Copy(&vecXhat);
		//tmp = tmp + u
		vecTmp.AddI(&vecU,1,1);
		//this updates z
		Shrinkage(&vecTmp,lambda / rho,&vecZ);
//		printf("vecZ after shrinkage\n");
//		vecZ.Print();


//	    % u-update
//	    u = u + (x_hat - z);
		vecTmp.Copy(&vecXhat); //tmp  = x_hat
		vecTmp.AddI(&vecZ,1,-1);// tmp = tmp - z
		vecU.AddI(&vecTmp);// u = u + tmp

//		printf("vecU after update\n");
//		vecU.Print();

		//find the objective value
		histObjVal.data[k] = Objective(NULL,vecb,lambda,&vecX,&vecZ);
		//history.rnom = norm(x - z)
		vecTmp.Copy(&vecX);
		vecTmp.AddI(&vecZ,1,-1);
//		printf("x - z\n");
//		vecTmp.Print();
		histRNorm.data[k] = vecTmp.Norm(2);
		//history.s_norm = norm(-rho * (z - zold))
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecZold,-1 * rho,rho);
//		printf("-rho * ( z - zold)\n");
//		vecTmp.Print();
		histSNorm.data[k] = vecTmp.Norm(2);
		//    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		double normX = vecX.Norm(2);
		vecTmp.Copy(&vecZ);
		vecTmp.ScaleI(-1);
		double normZ = vecTmp.Norm(2);
		if(normX > normZ)
		{
			histEpsPri.data[k] = sqrt(vecb->length) * abstol + reltol * normX;
		}
		else
		{
			histEpsPri.data[k] = sqrt(vecb->length) * abstol + reltol * normZ;
		}

		//    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		vecTmp.Copy(&vecU);
		vecTmp.ScaleI(rho);
		histEpsDual.data[k] = sqrt(vecb->length) * abstol + reltol * vecTmp.Norm(2);

		if(iverbose)
		{
//			fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
//            history.r_norm(k), history.eps_pri(k), ...
//            history.s_norm(k), history.eps_dual(k), history.objval(k));
			printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
					histRNorm.data[k],histEpsPri.data[k],histSNorm.data[k],
					histEpsDual.data[k],histObjVal.data[k]);
		}

//		  if (history.r_norm(k) < history.eps_pri(k) && ...
//		       history.s_norm(k) < history.eps_dual(k))
//		         break;
//		    end
		if(histRNorm.data[k] < histEpsPri.data[k]  &&  histSNorm.data[k] < histEpsDual.data[k])
			break;
	}

	vecxhat->Copy(&vecZ);

	return 0;
}
double CLasso::Objective(CMatrix* matA,CVector* vecB,double lambda,CVector* vecX,CVector* vecZ)
{
	//   p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
	CVector tmp(vecB->length,1);
	if(matA == NULL)
	{
		tmp.Copy(vecX);
	}
	else
	{
		matA->MultiplyV(&tmp,vecX); //tmp = A * x
	}
	tmp.AddI(vecB,1,-1);  // tmp = tmp - b
//	printf("Ax -b\n");
//	tmp.Print();
	double s = 0;
	for(int i =0;i < tmp.length;i++)
	{
		s = s + tmp.data[i] * tmp.data[i];
	}
	s = (s / 2) + (lambda * vecZ->Norm(1));
	return s;
}
MKL_INT CLasso::Shrinkage(CVector* x,double  kappa, CVector* res)
{
	 //the shrinkange function is  z = max( 0, x - kappa ) - max( 0, -x - kappa )

	CVector vecKappa(x->length,x->bColVector);
	CVector tmp1(x->length,x->bColVector);
	CVector tmp2(x->length,x->bColVector);

	vecKappa.Ones();
	vecKappa.ScaleI(kappa);

	tmp1.Copy(x);
	tmp1.AddI(&vecKappa,1,-1);
	//following line assigns res = max(0, x - kappa)
	tmp1.Max(0,res);

	tmp2.Copy(x);
	tmp2.AddI(&vecKappa,-1,-1);
	//the followin does the same for tmp1 = max(0, -x - kappa)
	tmp2.Max(0,&tmp1);

	//finally res = res - tmp1
	res->AddI(&tmp1,1,-1);

	return 0;
}
MKL_INT CLasso::BestLambda(CMatrix* mata, CVector* vecb, CVector* vecxhat,double* dlambda,double rho,double alpha,MKL_INT maxiter,MKL_INT iverbose)
{
	iDim = 0;
	iPathLen = 1;
	dThresh = 1e-5;
	iMaxIter = maxiter;
	double abstol = 1e-4;
	double reltol = 1e-2;

	double delta = 0.1;
	int numLambda = (dlambda[1] - dlambda[0]) / delta;
	double* pMse = new double[numLambda];

	//find the minimum Mse and select that lambda for fine search
	double min = pMse[0];
	int minidx = 0;
	pMse[0] = INFINITY;

	double ltime = omp_get_wtime();


	CVector vecAtb(mata->cols,1);

	mata->MultiplyV(&vecAtb,vecb,true);
//	printf("vecAtb\n");
//	vecAtb.Print();


	//compute the cholesky factors of the matrix A

	CMatrix Lo;//(mata->rows,mata->cols);
	CMatrix Up;//(mata->rows,mata->cols);
	CMatrix InvUp;
	CMatrix InvLo;


	if(mata->rows > mata->cols) // skinny
	{

		/*
		 * L = chol( A'*A + rho*speye(n), 'lower' );
		 */

		CMatrix eye(mata->rows,mata->cols);
		Lo.Create(mata->rows,mata->cols);
		InvLo.Create(mata->rows,mata->rows);
		Up.Create(mata->rows,mata->cols);
		InvUp.Create(mata->rows,mata->cols);

		eye.Eye();
		eye.Scale(rho);

		mata->Copy(&Lo);
		mata->MMMul(&Lo,&Up,1,0,true,false);
		Lo.Add(&Up,&eye);

		LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',Lo.cols,Lo.mat,Lo.cols);
		//zero out the upper part
		for(int i =0;i < Lo.rows;i++)
		{
			for(int j =0;j < Lo.cols;j++)
			{
				if(j <= i )
					continue;
				Lo.mat[i * Lo.cols + j] = 0;
			}
		}
		Lo.Copy(&Up,true);
		Lo.InvMat(&InvLo);
		Up.InvMat(&InvUp);
		Lo.Release();
		Up.Release();
//		printf("Lower Triangular matrix\n");
//		InvLo.Print();
	}
	else //fat
	{
		Lo.Create(mata->rows,mata->rows);
		InvLo.Create(mata->rows,mata->rows);
		Up.Create(mata->rows,mata->rows);
		InvUp.Create(mata->rows,mata->rows);

		CMatrix matAt(mata->cols,mata->rows);

//		printf("Matrix A\n");
//		mata->Print();
		Lo.Eye();
		mata->Copy(&matAt,true);
//		printf("second matrix\n");
//		matAt.Print();
		mata->MMMul(&matAt,&Lo,(1/rho),1,false,false);
//		Lo.Print();
		LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',Lo.cols,Lo.mat,Lo.cols);
		//zero out the upper part
		for(int i =0;i < Lo.rows;i++)
		{
			for(int j =0;j < Lo.cols;j++)
			{
				if(j <= i )
					continue;
				Lo.mat[i * Lo.cols + j] = 0;
			}
		}
//		printf("Lower Tri Mat \n");
//		Lo.Print();
		Lo.InvMat(&InvLo);
//		printf("InvLo\n");
//		InvLo.Print();
		Lo.Copy(&Up,true);
		Up.InvMat(&InvUp);
//		printf("InvUp\n");
//		InvUp.Print();
		Lo.Release();
		Up.Release();
	}

	//invert Lo and Up and save them

	if(iverbose)
	{
		printf("%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n", "iter","r norm", "eps pri", "s norm", "eps dual", "objective");
	}

	//ADMM Solver
	CVector vecX(mata->cols,1);
	CVector vecZ(mata->cols,1);
	CVector vecU(mata->cols,1);


	CVector histObjVal(maxiter,1);
	CVector histRNorm(maxiter,1);
	CVector histSNorm(maxiter,1);
	CVector histEpsPri(maxiter,1);
	CVector histEpsDual(maxiter,1);

	CVector vecyhat(mata->rows,1);


	//coarse lambda search

	for(int l = 0;l < numLambda;l++)
	{
		double lambda = dlambda[0] + (l* DELTALAMBDA);
		for(int k =0;k < maxiter;k++)
		{
			CVector vecQ(mata->cols,1);
			CVector vecTmp(mata->cols,1);
			CVector vecTmp2(mata->rows,1);
			CVector vecTmp3(mata->rows,1);
			CVector vecZold(mata->cols,1);
			CVector vecXhat(mata->cols,1);


			vecQ.Copy(&vecAtb); //q = Atb
			//tmp = z - u
			vecTmp.Copy(&vecZ);
			vecTmp.AddI(&vecU,1,-1);
			vecQ.AddI(&vecTmp,1,rho); // q = Atb + rho * (z - u)
	//		printf("VecQ\n");
	//		vecQ.Print();

			if(mata->rows > mata->cols ) // skinny matrix
			{
				//solve the cholesky factorization
				//x = U \ (L \ q)
				InvLo.MultiplyV(&vecTmp3,&vecQ);
				InvUp.MultiplyV(&vecX,&vecTmp3);
	//			printf("Lasso X estimate\n");
	//			vecX.Print(true);
			}
			else // fat matrix
			{
				//	x = q / rho - (A' * ( U \ (L \ (A * q))) / rho^2
				//vecTmp2 = A * q
				mata->MultiplyV( &vecTmp2, &vecQ);
	//			printf("VecTemp2\n");
	//			vecTmp2.Print();
				//vecTmp = U \ ( L \ vecTmp2)
				//vecTmp = inv(U) * ( inv(L) * vecTmp2)

				//vecTmp = inv(L) * vecTmp2
				InvLo.MultiplyV(&vecTmp3,&vecTmp2);
				InvUp.MultiplyV(&vecTmp2,&vecTmp3);
	//			LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'L',Lo.cols,1,Lo.mat,Lo.cols,vecTmp2.data,1);
	//			printf("VecTemp2\n");
	//			vecTmp2.Print();
				//vecX = vecQ
				vecX.Copy(&vecQ);
	//			printf("vecX\n");
	//			vecX.Print();
				// vecQ = A' * vecTmp2
				mata->MultiplyV( &vecQ,&vecTmp2,true);
	//			printf("vecQ\n");
	//			vecQ.Print();
				vecX.AddI(&vecQ,1/rho,-1 / (rho * rho));
	//			printf("vecX\n");
	//			vecX.Print();
			}

	//
	//		printf("vecX after cholesky soln\n");
	//		vecX.Print();
			//zold = z
			vecZold.Copy(&vecZ);
			//x_hat = alpha*x + (1- alpha) * zold
			vecXhat.Copy(&vecX);
			vecXhat.AddI(&vecZold,alpha,(1 - alpha));
	//		printf("vecXhat after alpha update\n");
	//		vecXhat.Print();
			//z = shrinkage(x_hat + u, lambda/rho);
			//vecTmp = x_hat
			vecTmp.Copy(&vecXhat);
			//tmp = tmp + u
			vecTmp.AddI(&vecU,1,1);
			//this updates z
			Shrinkage(&vecTmp,lambda / rho,&vecZ);
	//		printf("vecZ after shrinkage\n");
	//		vecZ.Print();


	//	    % u-update
	//	    u = u + (x_hat - z);
			vecTmp.Copy(&vecXhat); //tmp  = x_hat
			vecTmp.AddI(&vecZ,1,-1);// tmp = tmp - z
			vecU.AddI(&vecTmp);// u = u + tmp

	//		printf("vecU after update\n");
	//		vecU.Print();

			//find the objective value
			histObjVal.data[k] = Objective(mata,vecb,lambda,&vecX,&vecZ);
			//history.rnom = norm(x - z)
			vecTmp.Copy(&vecX);
			vecTmp.AddI(&vecZ,1,-1);
	//		printf("x - z\n");
	//		vecTmp.Print();
			histRNorm.data[k] = vecTmp.Norm(2);
			//history.s_norm = norm(-rho * (z - zold))
			vecTmp.Copy(&vecZ);
			vecTmp.AddI(&vecZold,-1 * rho,rho);
	//		printf("-rho * ( z - zold)\n");
	//		vecTmp.Print();
			histSNorm.data[k] = vecTmp.Norm(2);
			//    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
			double normX = vecX.Norm(2);
			vecTmp.Copy(&vecZ);
			vecTmp.ScaleI(-1);
			double normZ = vecTmp.Norm(2);
			if(normX > normZ)
			{
				histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normX;
			}
			else
			{
				histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normZ;
			}

			//    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
			vecTmp.Copy(&vecU);
			vecTmp.ScaleI(rho);
			histEpsDual.data[k] = sqrt(mata->cols) * abstol + reltol * vecTmp.Norm(2);

			if(iverbose)
			{
	//			fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
	//            history.r_norm(k), history.eps_pri(k), ...
	//            history.s_norm(k), history.eps_dual(k), history.objval(k));
				printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
						histRNorm.data[k],histEpsPri.data[k],histSNorm.data[k],
						histEpsDual.data[k],histObjVal.data[k]);
			}

	//		  if (history.r_norm(k) < history.eps_pri(k) && ...
	//		       history.s_norm(k) < history.eps_dual(k))
	//		         break;
	//		    end
			if(histRNorm.data[k] < histEpsPri.data[k]  &&  histSNorm.data[k] < histEpsDual.data[k])
				break;
		}

		//calculate the Mse and store it in an array
		//decide if its good enough
//		mat.MultiplyV(&vecYhat,&vec);maxiter
		mata->MultiplyV(&vecyhat,&vecZ);
		vecyhat.AddI(vecb,1,-1);
		pMse[l] =vecyhat.Norm(2);
	}

	//sort the mse and find the most suitable lambda for finer search
	min = pMse[0];
	minidx = 0;
	for(int l = 1; l < numLambda;l++)
	{
		if(pMse[l] < min)
		{
			min = pMse[l];
			minidx = l;
		}
	}
	delete [] pMse;

	double newdelta = 0.001;
	double minlambda = (minidx - 1) * delta + dlambda[0];
	double maxlambda = (minidx + 1) * delta + dlambda[0];
	numLambda = (maxlambda - minlambda) / newdelta;
	pMse = new double[numLambda];
	pMse[0] = INFINITY;

	for(int l = 1; l < numLambda;l++)
	{
		double lambda = l * newdelta + minlambda;
		for(int k =0;k < maxiter;k++)
		{
			CVector vecQ(mata->cols,1);
			CVector vecTmp(mata->cols,1);
			CVector vecTmp2(mata->rows,1);
			CVector vecTmp3(mata->rows,1);
			CVector vecZold(mata->cols,1);
			CVector vecXhat(mata->cols,1);


			vecQ.Copy(&vecAtb); //q = Atb
			//tmp = z - u
			vecTmp.Copy(&vecZ);
			vecTmp.AddI(&vecU,1,-1);
			vecQ.AddI(&vecTmp,1,rho); // q = Atb + rho * (z - u)maxiter
	//		printf("VecQ\n");
	//		vecQ.Print();

			if(mata->rows > mata->cols ) // skinny matrix
			{
				//solve the cholesky factorization
				//x = U \ (L \ q)
				InvLo.MultiplyV(&vecTmp3,&vecQ);
				InvUp.MultiplyV(&vecX,&vecTmp3);
	//			printf("Lasso X estimate\n");
	//			vecX.Print(true);
			}
			else // fat matrix
			{
				//	x = q / rho - (A' * ( U \ (L \ (A * q))) / rho^2
				//vecTmp2 = A * q
				mata->MultiplyV( &vecTmp2, &vecQ);
	//			printf("VecTemp2\n");
	//			vecTmp2.Print();
				//vecTmp = U \ ( L \ vecTmp2)
				//vecTmp = inv(U) * ( inv(L) * vecTmp2)

				//vecTmp = inv(L) * vecTmp2
				InvLo.MultiplyV(&vecTmp3,&vecTmp2);
				InvUp.MultiplyV(&vecTmp2,&vecTmp3);
	//			LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'L',Lo.cols,1,Lo.mat,Lo.cols,vecTmp2.data,1);
	//			printf("VecTemp2\n");
	//			vecTmp2.Print();
				//vecX = vecQ
				vecX.Copy(&vecQ);
	//			printf("vecX\n");
	//			vecX.Print();
				// vecQ = A' * vecTmp2
				mata->MultiplyV( &vecQ,&vecTmp2,true);
	//			printf("vecQ\n");
	//			vecQ.Print();
				vecX.AddI(&vecQ,1/rho,-1 / (rho * rho));
	//			printf("vecX\n");
	//			vecX.Print();
			}

	//
	//		printf("vecX after cholesky soln\n");
	//		vecX.Print();
			//zold = z
			vecZold.Copy(&vecZ);
			//x_hat = alpha*x + (1- alpha) * zold
			vecXhat.Copy(&vecX);
			vecXhat.AddI(&vecZold,alpha,(1 - alpha));
	//		printf("vecXhat after alpha update\n");
	//		vecXhat.Print();
			//z = shrinkage(x_hat + u, lambda/rho);
			//vecTmp = x_hat
			vecTmp.Copy(&vecXhat);
			//tmp = tmp + u
			vecTmp.AddI(&vecU,1,1);
			//this updates z
			Shrinkage(&vecTmp,lambda / rho,&vecZ);
	//		printf("vecZ after shrinkage\n");
	//		vecZ.Print();


	//	    % u-update
	//	    u = u + (x_hat - z);
			vecTmp.Copy(&vecXhat); //tmp  = x_hat
			vecTmp.AddI(&vecZ,1,-1);// tmp = tmp - z
			vecU.AddI(&vecTmp);// u = u + tmp

	//		printf("vecU after update\n");
	//		vecU.Print();

			//find the objective value
			histObjVal.data[k] = Objective(mata,vecb,lambda,&vecX,&vecZ);
			//history.rnom = norm(x - z)
			vecTmp.Copy(&vecX);
			vecTmp.AddI(&vecZ,1,-1);
	//		printf("x - z\n");
	//		vecTmp.Print();
			histRNorm.data[k] = vecTmp.Norm(2);
			//history.s_norm = norm(-rho * (z - zold))
			vecTmp.Copy(&vecZ);
			vecTmp.AddI(&vecZold,-1 * rho,rho);
	//		printf("-rho * ( z - zold)\n");
	//		vecTmp.Print();
			histSNorm.data[k] = vecTmp.Norm(2);
			//    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
			double normX = vecX.Norm(2);
			vecTmp.Copy(&vecZ);
			vecTmp.ScaleI(-1);
			double normZ = vecTmp.Norm(2);
			if(normX > normZ)
			{
				histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normX;
			}
			else
			{
				histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normZ;
			}

			//    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
			vecTmp.Copy(&vecU);
			vecTmp.ScaleI(rho);
			histEpsDual.data[k] = sqrt(mata->cols) * abstol + reltol * vecTmp.Norm(2);

			if(iverbose)
			{
	//			fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
	//            history.r_norm(k), history.eps_pri(k), ...
	//            history.s_norm(k), history.eps_dual(k), history.objval(k));
				printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
						histRNorm.data[k],histEpsPri.data[k],histSNorm.data[k],
						histEpsDual.data[k],histObjVal.data[k]);
			}

	//		  if (history.r_norm(k) < history.eps_pri(k) && ...
	//		       history.s_norm(k) < history.eps_dual(k))
	//		         break;
	//		    end
			if(histRNorm.data[k] < histEpsPri.data[k]  &&  histSNorm.data[k] < histEpsDual.data[k])
				break;
		}

		mata->MultiplyV(&vecyhat,&vecZ);
		vecyhat.AddI(vecb,1,-1);
		pMse[l] =vecyhat.Norm(2);

	}
	min = pMse[0];
	minidx = 0;

	for(int l = 1; l < numLambda;l++)
	{
		if(pMse[l] < min)
		{
			min = pMse[l];
			minidx = l;
		}
	}

	delete [] pMse;


	// perform the lasso one more time to get the final estimate
	double lambda = minidx *newdelta + minlambda;

	for(int k =0;k < maxiter;k++)
	{
		CVector vecQ(mata->cols,1);
		CVector vecTmp(mata->cols,1);
		CVector vecTmp2(mata->rows,1);
		CVector vecTmp3(mata->rows,1);
		CVector vecZold(mata->cols,1);
		CVector vecXhat(mata->cols,1);


		vecQ.Copy(&vecAtb); //q = Atb
		//tmp = z - u
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecU,1,-1);
		vecQ.AddI(&vecTmp,1,rho); // q = Atb + rho * (z - u)
//		printf("VecQ\n");
//		vecQ.Print();

		if(mata->rows > mata->cols ) // skinny matrix
		{
			//solve the cholesky factorization
			//x = U \ (L \ q)
			InvLo.MultiplyV(&vecTmp3,&vecQ);
			InvUp.MultiplyV(&vecX,&vecTmp3);
//			printf("Lasso X estimate\n");
//			vecX.Print(true);
		}
		else // fat matrix
		{
			//	x = q / rho - (A' * ( U \ (L \ (A * q))) / rho^2
			//vecTmp2 = A * q
			mata->MultiplyV( &vecTmp2, &vecQ);
//			printf("VecTemp2\n");
//			vecTmp2.Print();
			//vecTmp = U \ ( L \ vecTmp2)
			//vecTmp = inv(U) * ( inv(L) * vecTmp2)

			//vecTmp = inv(L) * vecTmp2
			InvLo.MultiplyV(&vecTmp3,&vecTmp2);
			InvUp.MultiplyV(&vecTmp2,&vecTmp3);
//			LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'L',Lo.cols,1,Lo.mat,Lo.cols,vecTmp2.data,1);
//			printf("VecTemp2\n");
//			vecTmp2.Print();
			//vecX = vecQ
			vecX.Copy(&vecQ);
//			printf("vecX\n");
//			vecX.Print();
			// vecQ = A' * vecTmp2
			mata->MultiplyV( &vecQ,&vecTmp2,true);
//			printf("vecQ\n");
//			vecQ.Print();
			vecX.AddI(&vecQ,1/rho,-1 / (rho * rho));
//			printf("vecX\n");
//			vecX.Print();
		}

//
//		printf("vecX after cholesky soln\n");
//		vecX.Print();
		//zold = z
		vecZold.Copy(&vecZ);
		//x_hat = alpha*x + (1- alpha) * zold
		vecXhat.Copy(&vecX);
		vecXhat.AddI(&vecZold,alpha,(1 - alpha));
//		printf("vecXhat after alpha update\n");
//		vecXhat.Print();
		//z = shrinkage(x_hat + u, lambda/rho);
		//vecTmp = x_hat
		vecTmp.Copy(&vecXhat);
		//tmp = tmp + u
		vecTmp.AddI(&vecU,1,1);
		//this updates z
		Shrinkage(&vecTmp,lambda / rho,&vecZ);
//		printf("vecZ after shrinkage\n");
//		vecZ.Print();


//	    % u-update
//	    u = u + (x_hat - z);
		vecTmp.Copy(&vecXhat); //tmp  = x_hat
		vecTmp.AddI(&vecZ,1,-1);// tmp = tmp - z
		vecU.AddI(&vecTmp);// u = u + tmp

//		printf("vecU after update\n");
//		vecU.Print();

		//find the objective value
		histObjVal.data[k] = Objective(mata,vecb,lambda,&vecX,&vecZ);
		//history.rnom = norm(x - z)
		vecTmp.Copy(&vecX);
		vecTmp.AddI(&vecZ,1,-1);
//		printf("x - z\n");
//		vecTmp.Print();
		histRNorm.data[k] = vecTmp.Norm(2);
		//history.s_norm = norm(-rho * (z - zold))
		vecTmp.Copy(&vecZ);
		vecTmp.AddI(&vecZold,-1 * rho,rho);
//		printf("-rho * ( z - zold)\n");
//		vecTmp.Print();
		histSNorm.data[k] = vecTmp.Norm(2);
		//    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
		double normX = vecX.Norm(2);
		vecTmp.Copy(&vecZ);
		vecTmp.ScaleI(-1);
		double normZ = vecTmp.Norm(2);
		if(normX > normZ)
		{
			histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normX;
		}
		else
		{
			histEpsPri.data[k] = sqrt(mata->cols) * abstol + reltol * normZ;
		}

		//    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
		vecTmp.Copy(&vecU);
		vecTmp.ScaleI(rho);
		histEpsDual.data[k] = sqrt(mata->cols) * abstol + reltol * vecTmp.Norm(2);

		if(iverbose)
		{
//			fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
//            history.r_norm(k), history.eps_pri(k), ...
//            history.s_norm(k), history.eps_dual(k), history.objval(k));
			printf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
					histRNorm.data[k],histEpsPri.data[k],histSNorm.data[k],
					histEpsDual.data[k],histObjVal.data[k]);
		}

//		  if (history.r_norm(k) < history.eps_pri(k) && ...
//		       history.s_norm(k) < history.eps_dual(k))
//		         break;
//		    end
		if(histRNorm.data[k] < histEpsPri.data[k]  &&  histSNorm.data[k] < histEpsDual.data[k])
			break;
	}

	vecxhat->Copy(&vecZ);
	printf("time: %lf\n",omp_get_wtime() - ltime);

	return 0;
}
CLasso::~CLasso() {

}

