//============================================================================
// Name        : Kalman1.cpp
// Author      : Jehandad
// Version     :
// Copyright   : All rights reserved
// Description : Hello World in C++, Ansi-style
//============================================================================


#include "Kalman1.h"
#include "CRndStream.h"
#include "CTimeVaryingNW.h"
#include "CSparseMatrix.h"
#include "CLasso.h"
#include "CMatrix3D.h"
#include "CKalmanFilter.h"

#include <matio.h>
#include <omp.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/utsname.h>
#include <sys/types.h>

#include <mpi.h>
#include <signal.h>

#include <ios>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define DIR_NAME "/scratch/jxkhan/DataIdx%dGene%dMeas%dTimePts%d"
#define NODE_LOG "/scratch/jxkhan/Node%dLog.txt"

#define STACK_SIZE (1024 * 1024)

#define MINLAMBDA 0
#define MAXLAMBDA 2
#define NUM_REALIZATIONS 1
#define NUM_IDX 1
//#define GENEARR {10}
//#define OBSARR {7}
#define DEALLOC_COV 1
#define NGENE 3218
#define NOBS 	6
#define NTIMEPTS 11

#define DATAGENE 666
#define DATAIDX 777
#define DATAREALZ 888
#define DATADIE 999

#define NODEBUSY 1
#define NODEIDLE 0
#define NUMDATAELEMS 3
#define USEMPI 1

//Define the error reporting macro

struct ChildData
{
	int rank;
	int g;
	int indx;
	int crealizaion;

};

struct NWData
{
	int indx;
	int reliza;
	int pch;
	int splvl;
};


MKL_INT CalcError(CMatrix* pTrue, CMatrix* pEst,CVector* pererr,CVector* dist, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts);
MKL_INT CalcErrorGene(CVector* pTrue, CVector* pEst,CVector* pererr,CVector* dist, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts);
void process_mem_usage(double& vm_usage, double& resident_set);
MKL_INT CalcROCErrorGene(CVector* pTrue, CVector* pEst,CVector* pVecDetEdge,CVector* pvecActEdge, CVector* pvecDetZero,CVector* pvecActZero, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts);
static int  CloneFunc(void *arg);
MKL_INT ReadMat3d(char* filename,char* varname,int dim,CMatrix* retmat);
MKL_INT ReadMat(char* filename,char* varname,CMatrix* retmat);
static int  CreateNWProc(void *arg);
MKL_INT SaveMat3d(char* filename, int nVars, ...);
MKL_INT SaveMat(char* filename,int nVars, ... );


int main(int argc, char *argv[])
{
//printf("Intializing Parallel Realizations\n");

int numprocs = 0, namelen, id = 0,mpirank = 0;
char processor_name[MPI_MAX_PROCESSOR_NAME];

#ifdef USEMPI
setbuf(stdout, NULL);
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
MPI_Get_processor_name(processor_name,&namelen);
printf("MPI: Process %d on %s out of %d\n",mpirank,processor_name,numprocs);

#endif


//MKL_INT arrGene = {10,10,10,10};
//MKL_INT nGene = 50;
//MKL_INT nTimePoints = 3;
//MKL_INT arrObservations[] = {2,5,7,10};
//MKL_INT nObservations = 35;//arrObservations[idx];
double perChange = 0.2;

	{




#ifdef USEMPI
		if(mpirank == 0)
		{
			double start = omp_get_wtime();
			for(int idx = 0; idx < NUM_IDX;idx++)
			{
				MKL_INT nGene = NGENE;
				MKL_INT nTimePoints = NTIMEPTS;
				MKL_INT nObservations = NOBS;
//				double arrSparsityLvl[NUM_IDX] = {0.1, 0.2, 0.3,0.4,0.5};

				for(int curRealization = 0; curRealization < NUM_REALIZATIONS;curRealization++)
				{
					double vm,rss;
					process_mem_usage(vm,rss);
					//The master node generates the synthetic network and writes it to the disk
					printf("[%d]Realization %d Index %d Time: %lf Generating Network \n",mpirank,curRealization,idx, omp_get_wtime()- start);
					printf("[%d]Realization %d Index %d Time: %.0f Mem Usage VM: %f RSS:%f\n",mpirank,curRealization,idx,omp_get_wtime() - start,vm,rss);
					char strDirName[256] = {0};
					char strFileName[256] = {0};
					char strNodeFile[256] = {0};
					char strResultFile[256] = {0};
					sprintf(strDirName,DIR_NAME,idx,nGene,nObservations,nTimePoints);
					umask(0);
					printf("[%d]Realization %d Index %d Time: %lf Creating Dir\n",mpirank,curRealization,idx,omp_get_wtime()- start);
					mkdir(strDirName,S_IRWXU | S_IRWXG | S_IRWXO);
					sprintf(strFileName,"%s/Realization%d.mat",strDirName,curRealization);
					sprintf(strResultFile,"%s/Results%d.mat",strDirName,curRealization);


					//instantiate the process
					char *stack;                    /* Start of stack buffer */
					char *stackTop;                 /* End of stack buffer */
					pid_t pid = 0;
					struct utsname uts;
					stack =(char*)  malloc(STACK_SIZE);
					int st;
					if (stack == NULL)
					{
						printf("Error Allocating Stack\n");
					}
					stackTop = stack + STACK_SIZE;  /* Assume stack grows downward */

					NWData data;
					data.indx = idx;
					data.pch = perChange;
					data.reliza = curRealization;
					data.splvl = 0.18;//arrSparsityLvl[idx];

//					pid = clone(CreateNWProc, stackTop, /*CLONE_NEWUTS | */ SIGCHLD | CLONE_VFORK /*| CLONE_NEWPID*/, &data);
//					if (pid == -1)
//					{
//						printf("Error Creating Clone\n");
//					}
//					do
//					{
//						if(waitpid(pid, &st, WUNTRACED | WCONTINUED) == -1)    /* Wait for child */
//						{
//							printf("Error WaitPid\n");
//						}
//					}while(!WIFEXITED(st) && !WIFSIGNALED(st));
//
//
//					CTimeVaryingNW nw;
//					nw.CreateNW(nGene,nTimePoints,nObservations,perChange,0.18,time(NULL),curRealization);
//
//					SaveMat3d(strFileName,3,
//							(nw.matGeneInteractions),nw.nTimePts,"A_Array",
//							nw.pMatX,nw.nTimePts,"X_Array",
//							nw.pMatY,nw.nTimePts,"Y_Array");
					//also write the same info to the results file

		// TODO the CTimeVaryingNW class should automatically dealloc its variables as soon as the variables go out of scope
					mkl_free_buffers();
					printf("[%d]Realization %d Index %d Time: %lf Network Generated \n",mpirank,curRealization,idx, omp_get_wtime()- start);
					//now we will wait send the gene numbers, realization and index to the target processes
					//for each index default zero
						//for each realization default zero
					CMatrix* matEstA = new CMatrix[nTimePoints + 1];
					CMatrix* matEstSmtA = new CMatrix[nTimePoints + 1];
					//allocate the matrices

					for(int t = 0;t < (nTimePoints + 1); t++)
					{
						matEstA[t].Create(nGene,nGene);
						matEstSmtA[t].Create(nGene,nGene);
					}
					CMatrix  matEstError(1,nTimePoints);
					CMatrix  matActEdge(1,nTimePoints);
					CMatrix matDetEdge(1,nTimePoints);
					CMatrix matActZero(1,nTimePoints);
					CMatrix matDetZero(1,nTimePoints);


					MPI_Status status;
					for(MKL_INT g =0; g < (numprocs - 1);g++)
					{
						MPI_Send(&g,1,MPI_INT,g+1,DATAGENE,MPI_COMM_WORLD);//the new job description
						MPI_Send(&idx,1,MPI_INT,g+1,DATAIDX,MPI_COMM_WORLD);
						MPI_Send(&curRealization,1,MPI_INT,g+1,DATAREALZ,MPI_COMM_WORLD);
					}

					printf("1\n");
					for(MKL_INT g = numprocs - 1 ;g < nGene;g++)
					{
						//the first one to get back will get the next gene
						int newnode = 0;
						MPI_Recv(&newnode,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
						//read in its results and store them in the central file
						if(status.MPI_TAG == DATAGENE)
						{
							CMatrix* pgEstA = new CMatrix[nTimePoints+1];
							CMatrix* pgEstSmtA = new CMatrix[nTimePoints + 1];

							sprintf(strNodeFile,"%s/Node%d.mat",strDirName,status.MPI_SOURCE);

							for(int t = 0;t < (nTimePoints + 1);t++)
							{
								if(ReadMat3d(strNodeFile,"EstKalmanApos",t,&(pgEstA[t])) == -1)
								{
									eprintf("Readmat failed");
									return -1;
								}
								if(ReadMat3d(strNodeFile,"EstKalmanSmt",t,&(pgEstSmtA[t])) == -1)
								{
									eprintf("Readmat failed");
									return -1;
								}

								if(matEstA[t].SetVec(&(pgEstA[t]),newnode,true) == -1)
									eprintf("SetVec Failed");
								if(matEstSmtA[t].SetVec(&(pgEstSmtA[t]),newnode,true) == -1)
									eprintf("SetVec Failed");
							}
							/*
							CMatrix mat;
							CMatrix tmp(1,nTimePoints);
							ReadMat(strNodeFile,"pKalmanErr",&mat);
	//						printf("Old Error prior to gene %d\n",newnode);
	//						matEstError.Print();
							matEstError.Copy(&tmp);
							matEstError.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matActEdge",&mat);
							matActEdge.Copy(&tmp);
							matActEdge.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matDetEdge",&mat);
							matDetEdge.Copy(&tmp);
							matDetEdge.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matActZero",&mat);
							matActZero.Copy(&tmp);
							matActZero.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matDetZero",&mat);
							matDetZero.Copy(&tmp);
							matDetZero.Add(&tmp,&mat);

	//						printf("New error after gene %d\n",newnode);
	//						matEstError.Print();*/
							// now we can remove the file
							remove(strNodeFile);
							for(int t = 0;t < (nTimePoints + 1) ;t++)
							{
								pgEstA[t].Release();
								pgEstSmtA[t].Release();
							}
							delete [] pgEstA;
							delete [] pgEstSmtA;
						}
						if(idx != 0)
						{
							int a =0;
							a = a+1;
						}
						//send it the next gene and other periphernalia
//						printf("[%d] gene = %d \n",mpirank,g);
						MPI_Send(&g,1,MPI_INT,status.MPI_SOURCE,DATAGENE,MPI_COMM_WORLD);
						MPI_Send(&idx,1,MPI_INT,status.MPI_SOURCE,DATAIDX,MPI_COMM_WORLD);
						MPI_Send(&curRealization,1,MPI_INT,status.MPI_SOURCE,DATAREALZ,MPI_COMM_WORLD);
					}
					for(int i =1; i < numprocs;i++)
					{
						//wait for all processes to report completion
						int  newnode = 0;
						MPI_Recv(&newnode,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
						//read in its results and store them in the central file
						if(status.MPI_TAG == DATAGENE)
						{
							CMatrix* pgEstA = new CMatrix[nTimePoints+1];
							CMatrix* pgEstSmtA = new CMatrix[nTimePoints + 1];

							sprintf(strNodeFile,"%s/Node%d.mat",strDirName,status.MPI_SOURCE);

							for(int t = 0;t < (nTimePoints + 1);t++)
							{
								ReadMat3d(strNodeFile,"EstKalmanApos",t,&(pgEstA[t]));
								ReadMat3d(strNodeFile,"EstKalmanSmt",t,&(pgEstSmtA[t]));

								if(matEstA[t].SetVec(&(pgEstA[t]),newnode,true) == -1)
									eprintf("SetVec Failed");
								if(matEstSmtA[t].SetVec(&(pgEstSmtA[t]),newnode,true) == -1)
									eprintf("SetVec Failed");
							}

							for(int t = 0;t < (nTimePoints + 1) ;t++)
							{
								pgEstA[t].Release();
								pgEstSmtA[t].Release();
							}
							delete [] pgEstA;
							delete [] pgEstSmtA;
							/*

							CMatrix mat;
							CMatrix tmp(1,nTimePoints);
							ReadMat(strNodeFile,"pKalmanErr",&mat);
	//						printf("Old Error prior to gene %d\n",newnode);
	//						matEstError.Print();
							matEstError.Copy(&tmp);
							matEstError.Add(&tmp,&mat);
	//						printf("New error after gene %d\n",newnode);
							ReadMat(strNodeFile,"matActEdge",&mat);
							matActEdge.Copy(&tmp);
							matActEdge.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matDetEdge",&mat);
							matDetEdge.Copy(&tmp);
							matDetEdge.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matActZero",&mat);
							matActZero.Copy(&tmp);
							matActZero.Add(&tmp,&mat);

							ReadMat(strNodeFile,"matDetZero",&mat);
							matDetZero.Copy(&tmp);
							matDetZero.Add(&tmp,&mat);
*/
							// now we can remove the file
							remove(strNodeFile);
						}
					}
					matEstError.Scale(1/nGene);
					//now we can save the file for one realization
					SaveMat3d(strResultFile,2,
							matEstA,nTimePoints + 1,"matEstA",
							matEstSmtA,nTimePoints + 1,"matEstSmtA");/*,
							&matEstError,1,"matEstErr",
							&matActEdge,1,"matActEdge",
							&matDetEdge,1,"matDetEdge",
							&matActZero,1,"matActZero",
							&matDetZero,1,"matDetZero");*/
					//release the acquired memory
					for(int t = 0;t < (nTimePoints + 1);t++)
					{
						matEstA[t].Release();
						matEstSmtA[t].Release();
					}
					delete [] matEstA;
					delete [] matEstSmtA;
					matEstA = NULL;
					matEstSmtA = NULL;
				}
			}

			//when all the realization / indices and genes are done
			for(int i = 1;i < numprocs;i++)
			{
				MPI_Send(&i,1,MPI_INT,i,DATADIE,MPI_COMM_WORLD);//tell the nodes the job is done
			}

		}
		else // MPI Slave
		{
#endif
			while(1)
			{
				int curRealization = 0;
				int idx = 0;
				int gene = 2;
#ifdef USEMPI
				gene = -1;
				//wait for and receive the index for parameters, the gene number and the realization number from the MPI framework
				for(int i =0;i < NUMDATAELEMS;i++)
				{
					MPI_Status status;
					int data = 0;
					MPI_Recv(&data,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					switch(status.MPI_TAG)
					{
					case DATAGENE:
						gene = data;
						break;
					case DATAIDX:
						idx = data;
						break;
					case DATAREALZ:
						curRealization = data;
						break;
					case DATADIE:
						printf("MPI Process with rank %d exiting after receving DATADIE\n",mpirank);
						MPI_Finalize();
						return 0;
						break;
					}
				}

				double start = omp_get_wtime();

				char *stack;                    /* Start of stack buffer */
				char *stackTop;                 /* End of stack buffer */
				pid_t pid = 0;
				struct utsname uts;
				stack =(char*)  malloc(STACK_SIZE);
				int status;
				if (stack == NULL)
				{
					printf("Error Allocating Stack\n");
				}
				stackTop = stack + STACK_SIZE;  /* Assume stack grows downward */
				ChildData d;
				d.crealizaion = curRealization;
				d.g = gene;
				d.indx = idx;
				d.rank = mpirank;
//				CloneFunc(&d);
				pid = clone(CloneFunc, stackTop, /*CLONE_NEWUTS | */ SIGCHLD | CLONE_VFORK /*| CLONE_NEWPID*/, &d);
				if (pid == -1)
				{
					printf("Error Creating Clone\n");
				}
				if(wait(&status) == -1)
					printf("Error waiting for Child\n");
				//error checking
				if(WIFSIGNALED(status) || (WEXITSTATUS(status) == -1))
				{
					eprintf("Child returned error");
					MPI_Finalize();
					return 0;
				}
/*				do
				{
					if(wait( &status) == -1)    // Wait for child 
					{
						printf("Error WaitPid\n");
					}
				}while(!WIFEXITED(status));*/

				printf("[%d]Realization %d Index %d:Time Taken = %lf \n\n",mpirank,curRealization,idx,omp_get_wtime() - start); // complete
				MPI_Send(&gene,1,MPI_INT,0,DATAGENE,MPI_COMM_WORLD);//inform the master that the job is done
				mkl_free_buffers();
				double vm,rss;
				process_mem_usage(vm,rss);
				printf("[%d]Realization %d Index %d Time: %.0f Mem Usage VM: %f RSS:%f\n",mpirank,curRealization,idx,omp_get_wtime() - start,vm,rss);

#endif

			}
		}
	}

#ifdef USEMPI
	MPI_Finalize();
#endif
	printf("Simulation Run Completed Successfully\n");
	return 0;
}

MKL_INT CalcError(CMatrix* pTrue, CMatrix* pEst,CVector* pererr,CVector* dist, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts)
{
	CMatrix tmp(numgenes,numgenes);
	for(int z = 0; z < ntimepts;z++)
	{
//		printf("True Interaction Matrix\n");
//		pTrue[z].Print();

//		printf("Estimated Matrix\n");
//		pEst[z].Print();

//		printf("Difference Matrix\n");

		double trueEdgeCount = 0;
		tmp.Add(&(pTrue[z]),&(pEst[z]),1,-1);
//		tmp.Print();
		for(int i =0;i < tmp.rows;i++)
		{
			for(int j =0;j < tmp.cols;j++)
			{
//				printf("%f <= %f \n",fabs(tmp.mat[i * tmp.cols + j]),
//						(alpha * fabs(pTrue[z].mat[i * numgenes + j])));
				if(fabs(tmp.mat[i * tmp.cols + j]) <= (alpha * fabs(pTrue[z].mat[i * numgenes + j])))
				{
					//true edge
					trueEdgeCount++;
				}
			}
		}
		pererr->data[z] = (1 - (trueEdgeCount / (numgenes * numgenes))) * 100;
//		printf("%f \n",(1 - (trueEdgeCount / (numgenes * numgenes))) * 100);
	}

	return 0;
}

MKL_INT CalcROCErrorGene(CVector* pTrue, CVector* pEst,CVector* pVecDetEdge,CVector* pvecActEdge, CVector* pvecDetZero,CVector* pvecActZero, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts)
{
	CVector tmp(numgenes,1);
	for(int z = 0; z < ntimepts;z++)
	{
		double numactualedges = 0;
		double numdetectededges = 0;

		double numactualzero = 0;
		double numdetectedzero = 0;

		tmp.Copy(&(pTrue[z]));
		tmp.AddI(&(pEst[z + 1]),1,-1);
		for(int i =0;i < tmp.length;i++)
		{
			if(fabs(pTrue[z].data[i]) > 0 )
			{
				numactualedges++;
				if(fabs(tmp.data[i]) <= (alpha * fabs(pTrue[z].data[i]))) //an actual edge exists
				{
					//true edge
					numdetectededges++; //an edge that was detected as an edge
				}
			}
			else
			{
				numactualzero++;
				if(fabs(pEst[z+1].data[i] < 0.1))
					numdetectedzero++;
			}
		}
		pvecDetZero->data[z] = numdetectedzero;
		pvecActZero->data[z] = numactualzero;
		pVecDetEdge->data[z] = numdetectededges;
		pvecActEdge->data[z] = numactualedges;
	}
	return 0;
}


MKL_INT CalcErrorGene(CVector* pTrue, CVector* pEst,CVector* pererr,CVector* dist, double alpha,MKL_INT nobs, MKL_INT numgenes,MKL_INT ntimepts)
{
	CVector tmp(numgenes,1);
	for(int z = 0; z < ntimepts;z++)
	{
//		printf("True Interaction Matrix\n");
//		pTrue[z].Print(true);

//		printf("Estimated Matrix\n");
//		pEst[z + 1].Print(true);

//		printf("Difference Matrix\n");

		double trueEdgeCount = 0;
		tmp.Copy(&(pTrue[z]));
		tmp.AddI(&(pEst[z + 1]),1,-1);
//		tmp.Print(true);
		for(int i =0;i < tmp.length;i++)
		{
//				printf("%f <= %f \n",fabs(tmp.mat[i * tmp.cols + j]),
//						(alpha * fabs(pTrue[z].mat[i * numgenes + j])));
			if(fabs(tmp.data[i]) <= (alpha * fabs(pTrue[z].data[i])))
			{
				//true edge
				trueEdgeCount++;
			}
		}
		pererr->data[z] = (1 - (trueEdgeCount / (numgenes))) * 100;
//		printf("%f \n",(1 - (trueEdgeCount / (numgenes))) * 100);
	}

	return 0;
}




void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;

   vm_usage = vm_usage / 1024;
   resident_set = resident_set / 1024;
}

static int CloneFunc(void* arg)
{
	ChildData* cd = (ChildData*)arg;
	setbuf(stdout, NULL);


	MKL_INT nGene = NGENE;
	MKL_INT nTimePoints = NTIMEPTS;
	MKL_INT nObservations = NOBS;//arrObservations[idx];

	double vm,rss;
	double start = omp_get_wtime();
	int gene = cd->g;
	int mpirank = cd->rank;
	int idx = cd->indx;
	int curRealization = cd->crealizaion;

	printf("[%d] New data gene = %d idx = %d curRealization = %d \n",mpirank, gene,idx,curRealization);
	printf("[%d] New params nGene = %d nTimePoints = %d nObservations = %d\n",mpirank, nGene,nTimePoints, nObservations);
	printf("[%d]Realization %d Index %d Time: %.0f \n",mpirank,curRealization,idx,omp_get_wtime() - start);

//			MKL_INT curgene = 0;
	char strDirName[256] = {0};
	char strFileName[256] = {0};
	char strLogFileName[256] = {0};
	char strNodeFile[256] = {0};

	sprintf(strDirName,DIR_NAME,idx,nGene,nObservations,nTimePoints);
//	printf(strDirName);
	sprintf(strFileName,"%s/Realization%d.mat",strDirName,curRealization);
//	printf(strFileName);
	sprintf(strNodeFile,"%s/Node%d.mat",strDirName,mpirank);
//	printf(strNodeFile);
	sprintf(strLogFileName,NODE_LOG,mpirank);
//			FILE fidLog = fopen(strLogFileName,"a");

	//load the X and Y from the data file and set them in the respective variables
	// we need the entire X and just gene'th row of the Y matrix

	printf("[%d]Relaization:%d Index %d Time: %lf Begin Estimate Initial Condition \n",mpirank,curRealization,idx, omp_get_wtime()- start);

	CVector ic;

	double ltime = omp_get_wtime();

	printf("[%d]Realization %d Index %d Time: %lf Calc IC for gene %d \n",mpirank,curRealization,idx, omp_get_wtime()- start,gene);
	double dLambda[] = {MINLAMBDA,MAXLAMBDA};
	CMatrix matH(100,100);
	CVector vec2(nObservations,1);

//	printf("Hello\n");
	char str[] = "Y_Array";
	ReadMat3d(strFileName,str,0,&matH);
//				printf("Y_Array\n");
//				matH.Print();
	matH.GetVec(&vec2,gene,true);
//				printf("Vector\n");
//				vec2.Print(true);
	vec2.bColVector = true;

	ReadMat3d(strFileName,"X_Array",0,&matH);
//				printf("X_Array\n");//matH.Print();
	matH.Transpose();
	//constrain the initial condition using best lasso
	double delta = 0.1;
	int numLambda = (dLambda[1] - dLambda[0]) / delta;
	double dLargestLambda = 0;
	double* gcv = new double[numLambda];

	//find the minimum Mse and select that lambda for fine search
	double min = gcv[0];
	int minidx = 0;
	gcv[0] = INFINITY;


	ic.Create(nGene,true);
	CVector vec(nGene,1);
	//search for the lasso paramter
//#pragma omp parallel for schedule(dynamic)
	for(int k =1; k < numLambda;k++)
	{
//				printf("Lasso for Lambda = %f\n",k * delta + dlambda[0]);
		CVector vecYhat(nObservations,1);
//				printf("mata\n");
//				mat.Print();
//				printf("Observation vector\n");
//				vec2.Print();

		double lambda = k*delta + dLambda[0];
		double pLambda = 0;

		CLasso lasso(&matH,&vec2,&vec,lambda,1,1,100,0);
//				printf("Constrained Soln\n");
//				vec.Print();

		if(vec.IsZero() && dLargestLambda == 0)
		{
			dLargestLambda = lambda;
			gcv[k] = INFINITY;
		}
		else if( dLargestLambda == 0)
		{
			pLambda = 	vec.Cardinality();
			matH.MultiplyV(&vecYhat,&vec);
			vecYhat.AddI(&vec2,1,-1);

			gcv[k] = (vecYhat.Norm(2) / pow(( 1 - (pLambda / nGene) ),2)) / nGene;
		}
		else
			gcv[k] = INFINITY;

	}
	//sort the Mse to find the minimum
	min = gcv[0];
	minidx = 0;
	for(int k = 1; k < numLambda;k++)
	{
		if(gcv[k] < min)
		{
			min = gcv[k];
			minidx = k;
		}
	}
	//fine tune the parameter

	delete [] gcv;

	double newdelta = 0.001;
	double minlambda = (minidx - 1) * delta + dLambda[0];
	double maxlambda = (minidx + 1) * delta + dLambda[0];
	numLambda = (maxlambda - minlambda) / newdelta;
	gcv = new double[numLambda];
	gcv[0] = INFINITY;
	//the second search loop for fine grained search
//			printf("search for a finer lambda\n");
//#pragma omp parallel for schedule(dynamic)
	for(int k =1; k < numLambda;k++)
	{
		CVector vecYhat(nObservations,1);
		double lambda = k*delta + dLambda[0];
		double pLambda = 0;

		CLasso lasso(&matH,&vec2,&vec,k*newdelta + minlambda,1,1,100,0);
//				printf("Constrained Soln\n");
//				vec.Print();
		if(vec.IsZero() && dLargestLambda == 0)
		{
			dLargestLambda = lambda;
			gcv[k] = INFINITY;
		}
		else if( dLargestLambda == 0)
		{
			pLambda = 	vec.Cardinality();
			matH.MultiplyV(&vecYhat,&vec);
			vecYhat.AddI(&vec2,1,-1);

			gcv[k] = (vecYhat.Norm(2) / pow(( 1 - (pLambda / nGene) ),2)) / nGene;
		}
		else
			gcv[k] = INFINITY;
	}
	//find the min msetrueEdgeCount
	min = gcv[0];
	minidx = 0;

	for(int k = 1; k < numLambda;k++)
	{
		if(gcv[k] < min)
		{
			min = gcv[k];
			minidx = k;
		}
	}

	delete [] gcv;
	//restimate the connectivity row using the selected lasso
//			printf("Best fine lambda %f, corresponding Mse %f \n",minidx * newdelta + minlambda,pMse[minidx]);
	//find the optimum estimate (again) and overwrite the unconstrained one
	CLasso lasso(&matH,&vec2,&vec,minidx*newdelta + minlambda,1,1,100,0);
	ic.Copy(&vec);
//			printf("Final estimated vector ic %d Time %lf\n",gene,omp_get_wtime() - ltime);
//			ic.Print(true);

	printf("[%d]Realization %d Index %d Time: %lf Initial Condition Estimate Generated\n",mpirank,curRealization,idx, omp_get_wtime()- start);


	//Run Kalman filter on the network for each time and just one gene


	/*
	 * The Kalman filter object
	 */
	CKalmanFilter filter;
	//extract the observation vector from the file and store it in an array of CVectors
	CVector* pvecY = new CVector[nTimePoints];
	CMatrix* matX = new CMatrix[nTimePoints];
	for(int t = 0;t < nTimePoints;t++)
	{
		CMatrix mat(nObservations,nGene);
		ReadMat3d(strFileName,"Y_Array",t,&mat);
		pvecY[t].Create(nObservations,1);
		mat.GetVec(&(pvecY[t]),gene,true);
		ReadMat3d(strFileName,"X_Array",t,&(matX[t]));

	}
	if(filter.EstimateGene(matX,pvecY,nGene,nTimePoints,nObservations,&ic,dLambda,gene) == -1)
		return -1;

	//release the X annd Y matrices
	for(int t = 0;t < nTimePoints;t++)
	{
		pvecY[t].Destroy();
		matX[t].Release();
	}


	delete [] pvecY;
	delete [] matX;
	ic.Destroy();
	mkl_free_buffers();
	printf("[%d]Realization %d Index %d Time: %lf Gene: %d Kalman Estimate Complete\n",mpirank,curRealization,idx,omp_get_wtime()- start,gene);
	/*
	 * End Kalman Filter
	 */
/*
	CVector vecErr(nTimePoints,1);
	CVector vecDist(nTimePoints,1);




	printf("[%d]Realization %d Index %d Time: %lf Calculating Error\n",mpirank,curRealization,idx,omp_get_wtime()- start);
	//load the interaction matrices
	CVector* pvecA = new CVector[nTimePoints];
	for(int t = 0;t < nTimePoints;t++)
	{
		CMatrix mat(nGene,nGene);
		ReadMat3d(strFileName,"A_Array",t,&mat);
		pvecA[t].Create(nGene,1);
		mat.GetVec(&(pvecA[t]),gene,true);
		pvecA[t].bColVector = true;
	}
	CVector vecActEdge(nTimePoints,1);
	CVector vecDetEdge(nTimePoints,1);
	CVector vecActZero(nTimePoints,1);
	CVector vecDetZero(nTimePoints,1);

	CalcErrorGene(pvecA,filter.pgEstSmtA,&vecErr,&vecDist,0.3,nObservations,nGene,nTimePoints);
	CalcROCErrorGene(pvecA,filter.pgEstSmtA,&vecDetEdge,&vecActEdge,&vecDetZero,&vecActZero,0.3,nObservations,nGene,nTimePoints);

	//release the matrices / vectors
	for(int t = 0;t < nTimePoints;t++)
	{
		pvecA[t].Destroy();
	}
	delete [] pvecA;
	pvecA = NULL;
	*/
/*
	 * Calculate the error metrics and then save the results to the mat file
	 */
	printf("[%d]Realization %d Index %d Time: %lf Begining to save data\n",mpirank,curRealization,idx,omp_get_wtime()- start);
	//convert the CVector arrays to CMatrix arrays and then save them to a mat file
	CMatrix* mEstSmtA = new CMatrix[nTimePoints + 1];
	CMatrix* mEstA	 = new CMatrix[nTimePoints + 1];
	for(int t = 0;t < (nTimePoints + 1);t++)
	{
		mEstSmtA[t].Create(1,nGene);
//					filter.pgEstSmtA[t].Print();
		if(mEstSmtA[t].SetVec(&(filter.pgEstSmtA[t]),0,true) == -1)
			eprintf("SetVec Failed");

		mEstA[t].Create(1,nGene);
		if(mEstA[t].SetVec(&(filter.pgEstA[t]),0,true) == -1)
			eprintf("SetVec Failed");
	}
	CMatrix matErr(1,nTimePoints);
	CMatrix matDist(1,nTimePoints);
	CMatrix matNObs(1,1);
	CMatrix matNGene(1,1);
	CMatrix matNTimePts(1,1);
	CMatrix matActEdge(1,nTimePoints);
	CMatrix matDetEdge(1,nTimePoints);
	CMatrix matDetZero(1,nTimePoints);
	CMatrix matActZero(1,nTimePoints);

/*
	matErr.SetVec(&vecErr,0,true);
//		matErr.Print();
	matDist.SetVec(&vecDist,0,true);
	matActEdge.SetVec(&vecActEdge,0,true);
	matDetEdge.SetVec(&vecDetEdge,0,true);
	matActZero.SetVec(&vecActZero,0,true);
	matDetZero.SetVec(&vecDetZero,0,true);*/

	matNObs.Set(0,0,nObservations);
	matNGene.Set(0,0,nGene);
	matNTimePts.Set(0,0,nTimePoints);

	//check and report if the file already exists

	SaveMat3d(strNodeFile,2,
			mEstA,nTimePoints + 1,"EstKalmanApos",
			mEstSmtA,nTimePoints + 1,"EstKalmanSmt");

	SaveMat(strNodeFile,5,
			&matErr,"pKalmanErr",
			&matDist,"dKalmanErr",
			&matNObs,"nObs",
			&matNGene,"nGenes",
			&matNTimePts,"nTimePts");/*,
			&matActEdge,"matActEdge",
			&matDetEdge,"matDetEdge",
			&matActZero,"matActZero",
			&matDetZero,"matDetZero");*/

	for(int t = 0; t < (nTimePoints + 1);t++)
	{
		mEstSmtA[t].Release();
		mEstA[t].Release();
	}
	delete [] mEstA;
	delete [] mEstSmtA;
	mEstA = NULL;
	mEstSmtA = NULL;

	return 0;
}


MKL_INT ReadMat3d(char* filename,char* varname,int dim,CMatrix* retmat)
{
	mat_t* matfile;
	matvar_t *matvar;
	matvar_t *cellvar;
	size_t dims[2];
	int   start[2]={0,0},stride[2]={1,1},edge[2]={0};

	matfile = Mat_Open(filename,MAT_ACC_RDONLY);
	if(matfile == NULL)
	{
		printf("ReadMat3d read file error file: %s variable %s dim %d\n",filename,varname,dim);
		return -1;
//		assert(false);
	}

	cellvar = Mat_VarRead(matfile,varname);
	if(cellvar == NULL)
	{
		printf("ReadMat3d read cell error file: %s variable %s dim %d\n",filename,varname,dim);
		return -1;
//		assert(false);
	}

	matvar = Mat_VarGetCell(cellvar,dim);
	if(matvar == NULL)
	{
		printf("ReadMat3d read mat error file: %s variable %s dim %d\n",filename,varname,dim);
		return -1;
//		assert(false);
	}
//	Mat_VarPrint(cellvar,1);
//	Mat_VarPrint(matvar,1);

	retmat->Create(matvar->dims[1],matvar->dims[0]);

	memcpy(retmat->mat,matvar->data,sizeof(double) * matvar->dims[0] * matvar->dims[1]);

//
//	edge[0] = matvar->dims[0];
//	edge[1] = matvar->dims[1];
////	Mat_VarGetSize
//	Mat_Rewind(matfile);
//	int err = Mat_VarReadData(matfile,matvar,retmat->mat,start,stride,edge);
////	int err = Mat_VarReadDataLinear(matfile,matvar,retmat->mat,0,1,retmat->rows * retmat->cols);
	retmat->Transpose();
//	retmat->Print();
	Mat_VarFree(matvar);
//	Mat_VarFree(cellvar);

	Mat_Close(matfile);
	return 0;
}

MKL_INT ReadMat(char* filename,char* varname,CMatrix* retmat)
{
	mat_t* matfile;
	matvar_t *matvar;
	matvar_t *cellvar;
	size_t dims[2];
	int   start[2]={0,0},stride[2]={1,1},edge[2]={0};

	matfile = Mat_Open(filename,MAT_ACC_RDONLY);
	if(matfile == NULL)
	{
		printf("ReadMat read file error\n");
		assert(false);
	}

	matvar = Mat_VarRead(matfile,varname);
	if(matvar == NULL)
	{
		printf("ReadMat Error: Variable not found filename %s variable %s \n",filename,varname);
		return -1;
	}

//	Mat_VarPrint(cellvar,1);
//	Mat_VarPrint(matvar,1);

	retmat->Create(matvar->dims[1],matvar->dims[0]);

	memcpy(retmat->mat,matvar->data,sizeof(double) * matvar->dims[0] * matvar->dims[1]);

//
//	edge[0] = matvar->dims[0];
//	edge[1] = matvar->dims[1];
////	Mat_VarGetSize
//	Mat_Rewind(matfile);
//	int err = Mat_VarReadData(matfile,matvar,retmat->mat,start,stride,edge);
////	int err = Mat_VarReadDataLinear(matfile,matvar,retmat->mat,0,1,retmat->rows * retmat->cols);
	retmat->Transpose();
//	retmat->Print();
	Mat_VarFree(matvar);
//	Mat_VarFree(cellvar);

	Mat_Close(matfile);
	return 0;
}


static int  CreateNWProc(void *arg)
{
	NWData* d = (NWData*)arg;

	char strDirName[256] = {0};
	char strFileName[256] = {0};
	int nGene = NGENE;
	int nObservations = NOBS;
	int nTimePoints = NTIMEPTS;

	int idx = d->indx;
	int curRealization = d->reliza;
	int perChange = d->pch;
	int SparsityLvl = d->splvl;

	sprintf(strDirName,DIR_NAME,idx,nGene,nObservations,nTimePoints);
	sprintf(strFileName,"%s/Realization%d.mat",strDirName,curRealization);

	CTimeVaryingNW nw;
	nw.CreateNW(nGene,nTimePoints,nObservations,perChange,SparsityLvl,time(NULL),curRealization);

	SaveMat3d(strFileName,3,
			(nw.matGeneInteractions),nw.nTimePts,"A_Array",
			nw.pMatX,nw.nTimePts,"X_Array",
			nw.pMatY,nw.nTimePts,"Y_Array");

	return 0;
}
MKL_INT SaveMat3d(char* filename, int nVars, ...)
{
	/* triplets would be CMatrix**, depth of the matrix array and the name of the variable */
	mat_t *matfile;
	matvar_t **matvar;
	matvar_t *cellvar;
	size_t dims[2];

	CMatrix* m = NULL;
	int depth = 0;
	char* varname = NULL;

	va_list vl;
	va_start(vl,nVars);
	matfile = Mat_Open(filename,MAT_ACC_RDWR | MAT_FT_MAT73	);

	for(int i =0; i < nVars;i++)
	{
		m = va_arg(vl,CMatrix*);
		depth = va_arg(vl,int);
		varname = va_arg(vl,char*);

		matvar = new matvar_t*[depth+1];


		for(int j =0;j < depth;j++)
		{
			dims[0] = m[j].rows;
			dims[1] = m[j].cols;
			m[j].Transpose();
			matvar[j] = Mat_VarCreate(NULL,MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,m[j].mat,MAT_F_DONT_COPY_DATA);

		}
		matvar[depth] = NULL;


		dims[0] = depth;
		dims[1] = 1;

		cellvar = Mat_VarCreate(varname,MAT_C_CELL,MAT_T_CELL,2,dims,
                matvar,0);
		Mat_VarWrite( matfile, cellvar, MAT_COMPRESSION_ZLIB);
//		Mat_VarPrint(cellvar,1);
//		Mat_VarPrint(matvar[0],1);

		for(int j= 0 ; j < depth;j++)
		{
			m[j].Transpose();
			Mat_VarFree(matvar[j]);
		}
//		Mat_VarFree(cellvar);
		delete [] matvar;

	}
	Mat_Close(matfile);
	return 0;
}
MKL_INT SaveMat(char* filename,int nVars, ... )
{
	/*Pointer to the marix variable and the name by which to save it in pairs */

	mat_t *matfile;
	matvar_t *matvar;
	size_t dims[2];

	va_list vl;
	va_start(vl,nVars);
	matfile = Mat_Open(filename,MAT_ACC_RDWR | MAT_FT_MAT73);

	for(int i =0; i < nVars;i++)
	{

		CMatrix* matrix = NULL;
		char* varname = NULL;
		matrix = va_arg(vl,CMatrix*);
		varname = va_arg(vl,char*);

		dims[0] = matrix->rows;
		dims[1] = matrix->cols;

		matrix->Transpose();


		matvar = Mat_VarCreate(varname,MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,matrix->mat,0);
		Mat_VarWrite( matfile, matvar, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(matvar);


		matrix->Transpose();
	}
	Mat_Close(matfile);
	return 0;
}
