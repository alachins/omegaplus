
/*  
 *  OmegaPlus: A Parallel Tool for Rapid & Scalable Detection of 
 *	       Selective Sweeps in Genome Datasets
 *
 *  Copyright February 2012 by Nikolaos Alachiotis and Pavlos Pavlidis
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an email to
 *  Pavlos Pavlidis (pavlidisp@gmail.com) or
 *  Nikolaos Alachiotis (n.alachiotis@gmail.com)
 *  
 */


#include "OmegaPlus.h"

double gettime(void)
{
	struct timeval ttime;
	gettimeofday(&ttime , NULL);
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

cor_t ABS (cor_t input)
{
	if(input<=0.0f)
		return -input;

	return input;
}

void printHeading (FILE * fp)
{
	fprintf(fp,"\n\n                                      _______________");
	fprintf(fp,"\n\n                                         OmegaPlus"   );
	fprintf(fp,"\n                                      _______________");
	fprintf(fp,"\n\n\n\n OmegaPlus version 3.0.0 released by Nikolaos Alachiotis and Pavlos Pavlidis in December 2014.\n");
}

void printRunInfo (FILE * fp, int argc, char ** argv, int fileFormat, int imputeN, int imputeG, int binary)
{
	fprintf(fp,"\n\n Command:\n\n\t");

	int i;
			
	for(i=0; i<argc; ++i)
		fprintf(fp," %s",argv[i]);
	
	fprintf(fp,"\n\n\n");

	fprintf(fp, "\n Input file format (0:ms, 1:fasta, 2:macs, 3:vcf, 4:sf):\t%d\n", fileFormat);

	int imputeBoth = imputeN + imputeG;

	if(imputeBoth==2)
		fprintf(fp," Gap (-) imputation:\t\t\t\t\t\tON\n Ambiguous character (N) imputation:\t\t\t\tON\n");
	else
	{
		if(imputeG==1)
			fprintf(fp," Gap (-) imputation:\t\t\t\t\t\tON\n Ambiguous character (N) imputation:\t\t\t\tOFF\n");	
		if(imputeN==1)
			fprintf(fp," Gap (-) imputation:\t\t\t\t\t\tOFF\n Ambiguous character (N) imputation:\t\t\t\tON\n");
		if(imputeBoth==0)
			fprintf(fp," Gap (-) imputation:\t\t\t\t\t\tOFF\n Ambiguous character (N) imputation:\t\t\t\tOFF\n");	
	}

	if(borderTol>0)
		fprintf(fp," Omega search strategy:\t\t\t\t\t\t#SNPsLeft - %d < #SNPsRight < #SNPsLeft + %d\n",borderTol, borderTol);
	else
		fprintf(fp," Omega search strategy:\t\t\t\t\t\tExhaustive\n");

	if(fileFormat==FASTA_FORMAT || fileFormat==VCF_FORMAT)
	{
		if(binary==1)
			fprintf(fp," Alignment deduction to binary:\t\t\t\t\tON\n");
		else
			fprintf(fp," Alignment deduction to binary:\t\t\t\t\tOFF\n");
	}

	if(fileFormat!=VCF_FORMAT)
		fprintf(fp,"\n\n");	
}

void introMsg(int argc, char ** argv, FILE * fpInfo, int fileFormat, int imputeN, int imputeG, int binary)
{
	printHeading (stdout);

	printRunInfo (stdout, argc, argv, fileFormat, imputeN, imputeG, binary);
	
	printHeading (fpInfo);

	printRunInfo (fpInfo, argc, argv, fileFormat, imputeN, imputeG, binary);
}

#ifdef _USE_PTHREADS
/*static void pinToCore(int tid)
{
	cpu_set_t cpuset;
         
	CPU_ZERO(&cpuset);    
	CPU_SET(tid, &cpuset);

	if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
	{
		fprintf(stdout, "\n ERROR: Please specify a number of threads that is smaller or equal");
		fprintf(stdout, "\n        to the number of available physical cores (%d).\n\n",tid);
		exit(0);
	}
}*/

#ifdef _USE_PTHREADS_MEMINT
void initializeThreadData(threadData_t * cur, int i, int threads)
{
	cur->threadID=i;
	cur->threadTOTAL=threads;
	cur->threadBARRIER=0;
	cur->threadOPERATION=BUSYWAIT;
	cur->gridPartition=malloc(sizeof(gridPartition_t));
}

void initPartitionInfo(gridPartition_t * gridPartition, int grid, int threads)
{
	int i;
	int divthreads = grid / threads;
	int modthreads = grid % threads;

	gridPartition[0].limitLeft = -1;

	for(i=0;i<threads;i++)
	{
		gridPartition[i].limitRight = gridPartition[i].limitLeft + divthreads + 1;

		if(modthreads!=0)
		{
			gridPartition[i].limitRight++;
			modthreads--;
		}

		if(i<threads-1)
			gridPartition[i+1].limitLeft = gridPartition[i].limitRight - 1;	
	}
	gridPartition[threads-1].limitRight = grid;
}

alignment_struct * copyAlignmentData (alignment_struct * output, alignment_struct * input)
{
	int i,j;
	output = malloc(sizeof(alignment_struct));

	output->states = input->states;
	output->length=input->length;
	output->segsites=input->segsites;

	output->positionsInd = malloc(sizeof(int)*output->segsites);

	for(i=0;i<output->segsites;i++)
		output->positionsInd[i]=input->positionsInd[i];

	output->sequences=input->sequences;

	output->seqtable=input->seqtable;

	output->siteSize=input->siteSize;

	output->compressedSize=input->compressedSize;

	int size= output->states;

	if (output->states==2 || output->states==3)
		size--;		
		
	output->compressedArrays = malloc(sizeof(unsigned int *)*size);

	for(i=0;i<size;i++)
	{
		output->compressedArrays[i] = malloc(sizeof(unsigned int)*output->compressedSize);

		for(j=0;j<output->compressedSize;j++)
			output->compressedArrays[i][j] = input->compressedArrays[i][j];
	}

	return output;
}

void setThreadArgumentsMEMINT (threadData_t * threadData, int tid, alignment_struct * alignment, omega_struct * omega, FILE * fpReport, int matrixSizeMax, gridPartition_t * gridPartition, void * multi_sync_i)
{
	threadData[tid].gridPartition = gridPartition;
	threadData[tid].alignment=alignment;
	threadData[tid].omega=omega;
	threadData[tid].fpReport=fpReport;
	threadData[tid].matrixSizeMax = matrixSizeMax;
#ifdef _USE_PTHREADS_MULTI
	threadData[tid].multi_sync = (multi_sync_t *) multi_sync_i;
#endif
}

#ifdef _USE_PTHREADS_MULTI
int ismyTurn(int partitionIndex,multi_sync_t * multi_sync)
{
	return multi_sync->myTurn[partitionIndex];
}

int getThreadHelp(int partitionIndex,multi_sync_t * multi_sync)
{
	int i;
	int ttid = 0;

	int maxthreads = 2, threadcounter = 0;

	for(i=0;i<multi_sync->size;i++)
	{		
		if(multi_sync->threadAv[i]==1)
		{
			multi_sync->threadAv[i]=0;
			multi_sync->mtid[i]=partitionIndex;
			multi_sync->ttid[i]=ttid++;

			threadcounter++;

			if(threadcounter==maxthreads)
				break;
		}
	}

	for(i=0;i<multi_sync->size;i++)
	{		
		if(multi_sync->mtid[i]==partitionIndex)
		{
			multi_sync->tthreads[i]=ttid+1;
		}
	}

	multi_sync->tthreads[partitionIndex]=ttid+1;
	multi_sync->mtid[partitionIndex]=partitionIndex;	
	
	return ttid;	
}

int getNextIndex(int cur, int total)
{
	int temp = cur + 1;
	
	if(temp==total-1)
		temp = 1;

	return temp;
}

/*void setTurnNext(int partitionIndex,multi_sync_t * multi_sync)
{
	multi_sync->myTurn[partitionIndex]=0;
	
	int i;
	for(i=0;i<multi_sync->size;i++)
		if(multi_sync->mtid[i]==partitionIndex && i!=partitionIndex)
		{
			multi_sync->mtid[i]=-1;
			multi_sync->ttid[i]=-1;
			multi_sync->threadAv[i]=1;
		}

	int next = getNextIndex(partitionIndex, multi_sync->size);

	while(1)
	{
		if(multi_sync->coarsePartDone[next]!=1)
		{
			multi_sync->myTurn[next]=1;
			break;	
		}
		next = getNextIndex(next, multi_sync->size);
	}
}*/

void setTurnNextMaxLoad(int partitionIndex,multi_sync_t * multi_sync)
{
	multi_sync->myTurn[partitionIndex]=0;
	
	int i,max=0, pos=-1;

	for(i=0;i<multi_sync->size;i++)
		if(multi_sync->mtid[i]==partitionIndex && i!=partitionIndex)
		{
			multi_sync->mtid[i]=-1;
			multi_sync->ttid[i]=-1;
			multi_sync->threadAv[i]=1;
		}


	for(i=0;i<multi_sync->size;i++)
		if(multi_sync->load[i]>max)
		{
			max=multi_sync->load[i];
			pos=i;
		}

	if(pos!=-1)
		multi_sync->myTurn[pos]=1;
}

void passMyTurn(int tid,multi_sync_t * multi_sync)
{
	int next = getNextIndex(tid, multi_sync->size);
	int i = multi_sync->size + 1;
	while(i>-1)
	{
		if(multi_sync->coarsePartDone[next]!=1)
		{
			multi_sync->myTurn[tid]=0;
			multi_sync->myTurn[next]=1;
			break;	
		}
		next = getNextIndex(next, multi_sync->size);
		i--;
	}
}

int isLast(multi_sync_t * multi_sync, int tid)
{
	int i;
	int total=0;
	
	for(i=0;i<multi_sync->size;i++)
		total += multi_sync->coarsePartDone[i];

	if(total==multi_sync->size)
		return 1;

	return 0;
}

void setArgsToWorkersMULTI(multi_sync_t * multi_sync, threadArgPWC_t_MULTI * threadArgPWC_multi,int partitionIndex, alignment_struct * alignment, omega_struct * omega, int cvw_i, int firstRowToCompute,float ** myCorrelationMatrix, char * bits_in_16bitsLocal)
{
	int i;

	for(i=0;i<multi_sync->size;i++)
	{
		if(multi_sync->mtid[i]==partitionIndex)
		{
			threadArgPWC_multi[i].alignment = alignment;
			threadArgPWC_multi[i].omega = omega;
			threadArgPWC_multi[i].cvw_i = cvw_i;
			threadArgPWC_multi[i].firstRowToCompute = firstRowToCompute;
			threadArgPWC_multi[i].myCorrelationMatrix = myCorrelationMatrix;
			threadArgPWC_multi[i].bits_in_16bitsLocal = bits_in_16bitsLocal;
		}
	}
}

void terminateThreadsMULTI(multi_sync_t * multi_sync)
{
	int i;
	
	for(i=0;i<multi_sync->size;i++)
		multi_sync->hold[i]=-2;
}

void unleashTheWorkers(multi_sync_t * multi_sync, int partitionIndex, int mode)
{
	int i;

	for(i=0;i<multi_sync->size;i++)
		if(multi_sync->mtid[i]==partitionIndex)
			multi_sync->hold[i]=mode;

}

void synchronizeWorkers(multi_sync_t * multi_sync, int partitionIndex)
{
	int i;
	int status=0;	
	while(1)
	{
		status=0;
		for(i=0;i<multi_sync->size;i++)
			if(multi_sync->mtid[i]==partitionIndex)
				status += multi_sync->hold[i];		

		if(multi_sync->tthreads[partitionIndex]==status)
		{
			for(i=0;i<multi_sync->size;i++)
				if(multi_sync->mtid[i]==partitionIndex)
					multi_sync->hold[i]=-1;		
			
			break;
		}
	
	}

}
#endif

void processPartition (int partitionIndex, gridPartition_t * gridPartition, alignment_struct * alignment, omega_struct * omega, cor_t ** myCorrelationMatrix, void * multi_syncI)
{

	int i, first=0;
	int firstRowToCopy = -1;	    	
	int firstRowToCompute = 1;	
	int firstRowToAdd = 2;

	int lvw_i = gridPartition[partitionIndex].limitLeft, cvw_i;
	int grid = gridPartition[partitionIndex].limitRight;

	char * bits_in_16bitsLocal=NULL;

#ifdef _SHARED
	bits_in_16bitsLocal = bits_in_16bits;
#else
	bits_in_16bitsLocal = malloc(sizeof(char)*(0x1u << 16));

	compute_bits_in_16bitsLocal(bits_in_16bitsLocal);
#endif

#ifdef _USE_PTHREADS_MULTI
	
	int ttid = 0;

	int inithelpLock = 20;
	int helpLock = inithelpLock;

	multi_sync_t * multi_sync = (multi_sync_t *) multi_syncI;
#endif

	for(i=lvw_i+1;i<grid;i++)
	{
		cvw_i=findNextValidOmega(omega, lvw_i, grid);

#ifdef _USE_PTHREADS_MULTI
		if(ttid==0)
		{
			if(validGridP(cvw_i,grid))
			{	
			
				if(first!=0)	
					overlapCorrelationMatrixAdditions (alignment, omega, lvw_i, cvw_i, 
								       &firstRowToCopy, &firstRowToCompute, &firstRowToAdd);

				first=1;

				shiftCorrelationMatrixValues (omega, lvw_i, cvw_i, firstRowToCopy, myCorrelationMatrix);
		
				computeCorrelationMatrixPairwise (alignment, omega, cvw_i, firstRowToCompute, NULL, myCorrelationMatrix, bits_in_16bitsLocal);								

				applyCorrelationMatrixAdditions (omega, cvw_i,firstRowToAdd, myCorrelationMatrix);	
			
				computeOmegas (alignment, omega, cvw_i, NULL, myCorrelationMatrix);

				lvw_i = cvw_i;
			}
		}
		else
		{
			
			if(validGridP(cvw_i,grid))
			{	
			
				if(first!=0)	
					overlapCorrelationMatrixAdditions (alignment, omega, lvw_i, cvw_i, 
								       &firstRowToCopy, &firstRowToCompute, &firstRowToAdd);

				first=1;

				shiftCorrelationMatrixValues (omega, lvw_i, cvw_i, firstRowToCopy, myCorrelationMatrix);
		
				setArgsToWorkersMULTI(multi_sync, threadArgPWC_MULTI,partitionIndex, alignment, omega, cvw_i, firstRowToCompute,myCorrelationMatrix,bits_in_16bitsLocal);

				unleashTheWorkers(multi_sync, partitionIndex,0);

				computeCorrelationsMULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrix, bits_in_16bitsLocal, multi_sync->ttid[partitionIndex], multi_sync->tthreads[partitionIndex]);

				multi_sync->hold[partitionIndex]=1;

				synchronizeWorkers(multi_sync, partitionIndex);

				applyCorrelationMatrixAdditions (omega, cvw_i,firstRowToAdd, myCorrelationMatrix);	
			
				initmaxOmegaValueThreadsMULTI(multi_sync->size);

				unleashTheWorkers(multi_sync, partitionIndex,2);

				computeOmegasMULTI (alignment, omega, cvw_i, myCorrelationMatrix, multi_sync->ttid[partitionIndex], multi_sync->tthreads[partitionIndex], partitionIndex);

				multi_sync->hold[partitionIndex]=1;

				synchronizeWorkers(multi_sync, partitionIndex);

				getmaxOmegaValueThreadsMULTI(omega, cvw_i, multi_sync->size);

				lvw_i = cvw_i;

			}

			helpLock--;
		}

		multi_sync->load[partitionIndex]--;
#else
		if(validGridP(cvw_i,grid))
		{	
	
			if(first!=0)	
				overlapCorrelationMatrixAdditions (alignment, omega, lvw_i, cvw_i, 
							       &firstRowToCopy, &firstRowToCompute, &firstRowToAdd);

			first=1;

			shiftCorrelationMatrixValues (omega, lvw_i, cvw_i, firstRowToCopy, myCorrelationMatrix);

			computeCorrelationMatrixPairwise (alignment, omega, cvw_i, firstRowToCompute, NULL, myCorrelationMatrix, bits_in_16bitsLocal);								

			applyCorrelationMatrixAdditions (omega, cvw_i,firstRowToAdd, myCorrelationMatrix);	
	
			computeOmegas (alignment, omega, cvw_i, NULL, myCorrelationMatrix);

			lvw_i = cvw_i;
		}
#endif
#ifdef _USE_PTHREADS_MULTI	
		if(ttid==0)
		{
			if(ismyTurn(partitionIndex,multi_sync))
			{
				ttid=getThreadHelp(partitionIndex,multi_sync);
				multi_sync->ttid[partitionIndex]=ttid;
			}
		}
		else
		{
			if(helpLock==0)
			{
				helpLock = inithelpLock;

				ttid=0;

				setTurnNextMaxLoad(partitionIndex,multi_sync);

				if(ismyTurn(partitionIndex,multi_sync))
				{
					ttid=getThreadHelp(partitionIndex,multi_sync);
					multi_sync->ttid[partitionIndex]=ttid;
				}
			}
		}
#endif
	}

	//printf("Thread %d at %f\n",partitionIndex,gettime()-mainTime0);
}

void * thread (void * x)
{
	threadData_t * currentThread = (threadData_t *) x;
	
	int i, tid = currentThread->threadID;
	//double time0 = gettime();

	//pinToCore(tid);

	alignment_struct * alignment = NULL;

#ifdef _SHARED
	alignment = currentThread->alignment;
#else
	alignment = copyAlignmentData(alignment, currentThread->alignment);
#endif
	omega_struct * omega = currentThread->omega;

	gridPartition_t * gridPartition = currentThread->gridPartition;

	void * multi_syncT=NULL;
#ifdef _USE_PTHREADS_MULTI

	int cvw_i = -1;
	int firstRowToCompute = -1;
	float ** myCorrelationMatrixN = NULL;
	char * bits_in_16bitsLocal=NULL;

	multi_sync_t * multi_sync = currentThread->multi_sync;

	multi_syncT = multi_sync;
#endif

	int matrixSizeMax = currentThread->matrixSizeMax;

	cor_t ** myCorrelationMatrix = NULL;

	myCorrelationMatrix = createCorrelationMatrix(myCorrelationMatrix,matrixSizeMax);

	processPartition (tid, gridPartition, alignment, omega, myCorrelationMatrix, multi_syncT);

	if(myCorrelationMatrix != NULL)
		for(i=0; i<matrixSizeMax; ++i)
			if(myCorrelationMatrix[i] != NULL)
				free(myCorrelationMatrix[i]);
	
	if(myCorrelationMatrix != NULL)
		free(myCorrelationMatrix);

#ifdef _USE_PTHREADS_MULTI
	multi_sync->coarsePartDone[tid] = 1;	

	if(isLast(multi_sync, tid))
		terminateThreadsMULTI(multi_sync);
	else
	{
		if(ismyTurn(tid,multi_sync)==1)
			passMyTurn(tid,multi_sync);

		multi_sync->threadAv[tid] = 1;

		
		while(1)
		{
			sleep(0);
		
			if(multi_sync->hold[tid]==-2)
				break;

			if(multi_sync->hold[tid]==0)
			{
				alignment = threadArgPWC_MULTI[tid].alignment;
				omega = threadArgPWC_MULTI[tid].omega;
				cvw_i = threadArgPWC_MULTI[tid].cvw_i;
				firstRowToCompute = threadArgPWC_MULTI[tid].firstRowToCompute;
				myCorrelationMatrixN = threadArgPWC_MULTI[tid].myCorrelationMatrix;
				bits_in_16bitsLocal = threadArgPWC_MULTI[tid].bits_in_16bitsLocal;

				computeCorrelationsMULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, multi_sync->ttid[tid], multi_sync->tthreads[tid]);

				multi_sync->hold[tid]=1;
				
				while(multi_sync->hold[tid]==1) sleep(0);
			}

			if(multi_sync->hold[tid]==2)
			{

				alignment = threadArgPWC_MULTI[tid].alignment;
				omega = threadArgPWC_MULTI[tid].omega;
				cvw_i = threadArgPWC_MULTI[tid].cvw_i;
				myCorrelationMatrixN = threadArgPWC_MULTI[tid].myCorrelationMatrix;

				computeOmegasMULTI (alignment, omega, cvw_i, myCorrelationMatrixN, multi_sync->ttid[tid], multi_sync->tthreads[tid], tid);

				multi_sync->hold[tid]=1;
				
				while(multi_sync->hold[tid]==1) sleep(0);
			}
		}
	}
#endif

	//printf("thread %d: %f\n",tid,gettime()-time0);
	return NULL;		
}
#else
void initializeThreadData(threadData_t * cur, int i, int threads)
{
	cur->threadID=i;
	cur->threadTOTAL=threads;
	cur->threadBARRIER=0;
	cur->threadOPERATION=BUSYWAIT;
	cur->threadArgPWC=malloc(sizeof(threadArgPWC_t));
	cur->threadArgCO=malloc(sizeof(threadArgCO_t));
}

void syncThreadsBARRIER(threadData_t * threadData)
{
	int i, threads = threadData[0].threadTOTAL, barrierS=0;

	while(barrierS!=threads)
	{
		barrierS=0;
		for(i=0;i<threads;i++)
			barrierS += threadData[i].threadBARRIER;
	}

	for(i=0;i<threads;i++)
	{
		threadData[i].threadOPERATION=BUSYWAIT;
		threadData[i].threadBARRIER=0;
	}
}

void execFunctionMaster(threadData_t * threadData, int operation)
{
	if(operation==PAIRWISECORRELATION)
		correlationThread(&threadData[0]);

	if(operation==COMPUTEOMEGAS)
		omegasThread(&threadData[0]);
}

void setThreadOperation(threadData_t * threadData, int operation)
{
	int i, threads=threadData[0].threadTOTAL;
	
	for(i=0;i<threads;i++)
		threadData[i].threadOPERATION = operation;	
}

void startThreadOperations(threadData_t * threadData, int operation)
{
	setThreadOperation(threadData, operation);

	execFunctionMaster(threadData,operation);	

	threadData[0].threadBARRIER=1;

	syncThreadsBARRIER(threadData);		
}

void * thread (void * x)
{
	threadData_t * currentThread = (threadData_t *) x;
	
	//int tid = currentThread->threadID;
	
	//pinToCore(tid);

	while(1)
	{
		sleep(0);
		
		if(currentThread->threadOPERATION==EXIT)
			return NULL;
		
		if(currentThread->threadOPERATION==PAIRWISECORRELATION)
		{
			correlationThread(currentThread);

			currentThread->threadBARRIER=1;
			
			while(currentThread->threadBARRIER==1) sleep(0);
		}
		
		if(currentThread->threadOPERATION==COMPUTEOMEGAS)
		{			
			omegasThread(currentThread);

			currentThread->threadBARRIER=1;
			
			while(currentThread->threadBARRIER==1) sleep(0);
		}
	}
	
	return NULL;		
}

void terminateWorkerThreads(pthread_t * workerThreadL, threadData_t * threadData)
{
	int i, threads=threadData[0].threadTOTAL;
	
	for(i=0;i<threads;i++)
		threadData[i].threadOPERATION = EXIT;			

	for(i=1;i<threads;i++)
		pthread_join(workerThreadL[i-1],NULL);
}
#endif
#endif


void initializeAlignmentVariables(alignment_struct * alignment)
{
  
  alignment->positions = NULL;
  alignment->positionsInd = NULL;
  alignment->seqtable=NULL;
  alignment->compressedArrays = NULL;
  alignment->correlationMatrix = NULL;
  alignment->correlationMatrix = NULL;
		

}
	
int main(int argc, char** argv)
{
	mainTime0 = gettime();

	int grid = 0, 
	  alignmentLength=0,
	  minw=0,
	  maxw=0, maxwUSER=-1,
	  minsnps = MINSNPPERWINDOW,
	  nxt_alignment=0,
	  i,
	  noSeparation = 0,
	  generateVCFsamplelist=0,
#if _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
	    lvw_i, // last valid w
	    cvw_i, // current valid w
	    firstRowToCopy = -1, 
	    firstRowToCompute = 1,
	    firstRowToAdd = 2,
#endif
#else
	    lvw_i, // last valid w
	    cvw_i, // current valid w
#ifndef _USE_OPENMP_GENERIC
	    firstRowToCopy = -1, 
	    firstRowToCompute = 1,
	    firstRowToAdd = 2,
#endif
#endif
	    matrixSizeMax = 0,
            imputeN = 0,
	    imputeG = 0,
            binary = 0,
	    alignmentIndex=1,
	    fileFormat=OTHER_FORMAT,
	    threads=0,
	    resultType=-1,
	    ld=-1,
	    filterOut=0;

	    borderTol=-1;
	    int memLimit=0;
#ifdef _USE_OPENMP_GENERIC
	    
	    float memLimit_f=0.0;
	    int genGridList_size=0;
	    int * genGridList=NULL;
	    int genLoad_size=0;
	    int	prev_group_left_index = -1;
	    int	prev_group_right_index = -1;
	    float group_memsize = get_Mem_GroupSize();
	    float totalmem_genoutput=0.0;
	    int leftSNPindex = -1;
	    int rightSNPindex = -1;
	    int first_group_index=-1;
	    int last_group_index=-1;
	    int prev_first_group_index=-1;
	    int prev_last_group_index=-1;
	    int overlap=-1;
	    int overlap_workgroup_map_dim=-1;
            int overlap_workgroup_map_number_of_elements=-1;
	    float *** workgroup_map_ptr=NULL;
	    float ** workgroup_map_ptr_mem=NULL;
	    float *** overlap_workgroup_map_ptr=NULL;
	    float ** overlap_workgroup_map_ptr_mem=NULL;
	    int cur_i, j;
	    int overlap_start, overlap_i, overlap_j;
	    int prev_workgroup_map_dim=-1;
	    int workgroup_map_dim=0;
	    int workgroup_map_number_of_elements=0;
	    int cur_startIndex_tmp, cur_finishIndex_tmp;
	    int prev_startIndex_tmp, prev_finishIndex_tmp;
#endif

	unsigned int seed = 0;

	float maxomegaRealPos = -1.0;
	float maxomegamaxValue = -1.0;
	int maxomegamaxLeftIndex = -1;
	int maxomegamaxRightIndex = -1;

	double time0, time1, totalTimeL=.0, totalTimeG0 = gettime(), totalTimeG1;

  	char** recfile = malloc(sizeof(char*));
  	      *recfile = NULL;
	
	char inputFileName[INFILENAMESIZE],
	     infoFileName[INFILENAMESIZE], 
	     omegaReportFileName[INFILENAMESIZE], 
	     warnFileName[INFILENAMESIZE], 
	     sampleVCFfileName[INFILENAMESIZE];

	sampleVCFfileName[0]='\0';

	void * functionData=NULL;
	functionData = functionData;

	FILE *fpIn=NULL, *fpInfo=NULL, *fpWarnings=NULL, *fpReport=NULL, *fpVCFsamples=NULL;
	
	alignment_struct * alignment;

	omega_struct * omega;

	
   	commandLineParser(argc, argv, inputFileName, &grid, &alignmentLength, &minw, &maxw, recfile, 
			  &minsnps, &imputeN, &imputeG, &binary, &seed, &fileFormat, &threads, &resultType, &ld, &borderTol, &filterOut, &noSeparation, sampleVCFfileName,
			  &generateVCFsamplelist, &memLimit);

	maxwUSER = maxw;

	linkage_disequilibrium = ld;

	if(sampleVCFfileName[0] != '\0')
	{
		if(generateVCFsamplelist==1)
			fpVCFsamples = fopen(sampleVCFfileName,"w");
		else
			fpVCFsamples = fopen(sampleVCFfileName,"r");			
	}


#ifdef _USE_PTHREADS	
	//pinToCore(0); 
	
	workerThreadL=NULL; 

	workerThreadL = malloc (sizeof(pthread_t)*(threads-1));

	threadData_t * threadData = malloc (sizeof(threadData_t)*threads);
	
	functionData = threadData;

	for(i=0;i<threads;i++)
		initializeThreadData(&threadData[i],i,threads);

#ifdef _USE_PTHREADS_MEMINT
	gridPartition_t * gridPartition = malloc(sizeof(gridPartition_t)*threads);

	initPartitionInfo(gridPartition,grid,threads);
#else
	for(i=1;i<threads;i++)
		pthread_create (&workerThreadL[i-1], NULL, thread, (void *) (&threadData[i]));
#endif
#endif

	strcpy(infoFileName,"OmegaPlus_Info.");
	
	strncat(infoFileName,runName,INFILENAMESIZE-strlen(infoFileName));

	strcpy(warnFileName,"OmegaPlus_Warnings.");
	
	strncat(warnFileName,runName,INFILENAMESIZE-strlen(warnFileName));

	strcpy(omegaReportFileName,"OmegaPlus_Report.");

	strncat(omegaReportFileName,runName,INFILENAMESIZE-strlen(omegaReportFileName));


	fpIn = fopen(inputFileName,"r");
	
	fpInfo = fopen(infoFileName,"w");

	if (fileFormat==MS_FORMAT || fileFormat==MACS_FORMAT)
		fpWarnings = fopen(warnFileName, "w");

	fpReport = fopen(omegaReportFileName,"w");	
 
	
	introMsg(argc, argv, fpInfo, fileFormat, imputeN, imputeG, binary);


	alignment = (alignment_struct *)malloc(sizeof(alignment_struct));

	alignment->length = alignmentLength;

	omega = (omega_struct *)malloc(sizeof(omega_struct)*grid);

#ifdef _SHARED
	compute_bits_in_16bits();
#endif
	srand(seed);

	initializeGlobalPointers(alignment);
	

	time0 = gettime();

	nxt_alignment = findFirstAlignment(alignment, fpIn,fpInfo, fileFormat, fpVCFsamples, generateVCFsamplelist, sampleVCFfileName);
				

	while(nxt_alignment==1)
	{
	  
	  initializeAlignmentVariables( alignment );
	  
	  
	  
	  fprintf(stdout," Alignment %d\n",alignmentIndex);
			
	  fprintf(fpInfo," Alignment %d\n",alignmentIndex);
	  
	  if(noSeparation == 0)
	    fprintf(fpReport,"\n//%d\n",alignmentIndex);			
	  else
	    fprintf(fpReport, "\n");
	  
	  if( readAlignment(fpIn,alignment, imputeG, imputeN, binary, fileFormat, fpInfo, filterOut) == 1)
	    {

/*#ifdef _USE_OPENMP_GENERIC
	total_dp_init_time=malloc(sizeof(double)*threads);
	total_dp_update_time=malloc(sizeof(double)*threads);
	total_omega_values_time=malloc(sizeof(double)*threads);
	for(i=0;i<threads;i++)
	{
		total_dp_init_time[i]=0.0;
		total_dp_update_time[i]=0.0;
		total_omega_values_time[i]=0.0;
	}
#endif
*/
		    compressAlignment(alignment);		
		    
		    if (fileFormat==MS_FORMAT || fileFormat == MACS_FORMAT)
		      checkSNIPPositions(fpWarnings, alignment, alignmentIndex);
		    
		    maxw = maxwUSER;		
		    
		    matrixSizeMax = findOmegaBounds(alignment,omega, grid, &maxw, minw, minsnps);
		    
		    if(maxw!=maxwUSER)
		      {
			fprintf(stdout,"\n\n\t\tWARNING: Maximum window (-maxwin) has changed to %d",maxw);
			fprintf(fpInfo,"\n\n\t\tWARNING: Maximum window (-maxwin) has changed to %d",maxw);
		      }

		    fprintf(stdout,"\n\n\t\tProcessing (seconds):\t");		    
		    fprintf(fpInfo,"\n\n\t\tProcessing (seconds):\t");

		    fflush(stdout);
		    
		    
#ifdef _USE_PTHREADS_MEMINT
		    void * multiThreadsArgT=NULL;
#ifdef _USE_PTHREADS_MULTI
		    
		    multi_sync_t * multiThreadsArg = malloc (sizeof(multi_sync_t));
		    
		    multiThreadsArg->size = threads;
		    multiThreadsArg->coarsePartDone = malloc(sizeof(int)*threads);
		    multiThreadsArg->myTurn = malloc(sizeof(int)*threads);
		    multiThreadsArg->threadAv = malloc(sizeof(int)*threads);
		    multiThreadsArg->ttid = malloc(sizeof(int)*threads);
		    multiThreadsArg->mtid = malloc(sizeof(int)*threads);
		    multiThreadsArg->tthreads = malloc(sizeof(int)*threads);
		    multiThreadsArg->hold = malloc(sizeof(int)*threads);
		    multiThreadsArg->load = malloc(sizeof(int)*threads);
		    
		    
		    for(i=0;i<threads;i++)
		      {
			multiThreadsArg->coarsePartDone[i]=0;
			multiThreadsArg->myTurn[i]=0;
			multiThreadsArg->threadAv[i]=0;
			multiThreadsArg->ttid[i]=-1;
			multiThreadsArg->mtid[i]=-1;
			multiThreadsArg->tthreads[i]=0;
			multiThreadsArg->hold[i]=-1;
			multiThreadsArg->load[i]=gridPartition[i].limitRight-gridPartition[i].limitLeft-1;	
		      }
		    
		    multiThreadsArg->myTurn[2]=1;
		    
		    multiThreadsArgT = multiThreadsArg;
		    
		    threadArgPWC_MULTI=malloc(sizeof(threadArgPWC_t_MULTI)*threads);
		    
		    maxOmegaValueThreadsMULTI=malloc(sizeof(float)*threads);
		    maxOmegaLeftIndexThreadsMULTI=malloc(sizeof(int)*threads);
		    maxOmegaRightIndexThreadsMULTI=malloc(sizeof(int)*threads);
		    
		    
		    /*for(i=0;i<threads;i++)
		      {
		      printf("LOAD %d\n",i);
		      
		      printf("%d %d = %d\n",omega[(gridPartition[i].limitLeft+1)].leftIndex,omega[(gridPartition[i].limitRight-1)].rightIndex,omega[(gridPartition[i].limitRight-1)].rightIndex-omega[(gridPartition[i].limitLeft+1)].leftIndex);
		      }*/
		    
#endif
		    
		    for(i=0;i<threads;i++)
		      setThreadArgumentsMEMINT (threadData, i, alignment, omega, fpReport, matrixSizeMax, gridPartition, multiThreadsArgT);
		    
		    for(i=1;i<threads;i++)
		      pthread_create (&workerThreadL[i-1], NULL, thread, (void *) (&threadData[i]));
		    
		    cor_t ** myCorrelationMatrix = NULL;
		    
		    myCorrelationMatrix = createCorrelationMatrix(myCorrelationMatrix,matrixSizeMax);
		    
		    processPartition (0, gridPartition, alignment, omega, myCorrelationMatrix, multiThreadsArgT);
		    
		    
#ifdef _USE_PTHREADS_MULTI
		    multiThreadsArg->coarsePartDone[0] = 1;
		    
		    
		    if(isLast(multiThreadsArg, 0))
		      terminateThreadsMULTI(multiThreadsArg);
		    else
		      {
			if(ismyTurn(0,multiThreadsArg)==1)
			  passMyTurn(0,multiThreadsArg);
			
			multiThreadsArg->threadAv[0] = 1;
			
			while(1)
			  {
			    sleep(0);
			    
			    if(multiThreadsArg->hold[0]==-2)
			      break;
			    
			    if(multiThreadsArg->hold[0]==0)
			      {
				alignment = threadArgPWC_MULTI[0].alignment;
				omega = threadArgPWC_MULTI[0].omega;
				int cvw_i = threadArgPWC_MULTI[0].cvw_i;
				int firstRowToCompute = threadArgPWC_MULTI[0].firstRowToCompute;
				float ** myCorrelationMatrixN = threadArgPWC_MULTI[0].myCorrelationMatrix;
				char * bits_in_16bitsLocal = threadArgPWC_MULTI[0].bits_in_16bitsLocal;
				
				computeCorrelationsMULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, multiThreadsArg->ttid[0], multiThreadsArg->tthreads[0]);
				
				multiThreadsArg->hold[0]=1;
				
				while(multiThreadsArg->hold[0]==1) sleep(0);
			      }
			    
			    if(multiThreadsArg->hold[0]==2)
			      {
				alignment = threadArgPWC_MULTI[0].alignment;
				omega = threadArgPWC_MULTI[0].omega;
				int cvw_i = threadArgPWC_MULTI[0].cvw_i;
				float ** myCorrelationMatrixN = threadArgPWC_MULTI[0].myCorrelationMatrix;
				
				computeOmegasMULTI (alignment, omega, cvw_i, myCorrelationMatrixN, multiThreadsArg->ttid[0], multiThreadsArg->tthreads[0],0);
				
				multiThreadsArg->hold[0]=1;
				
				while(multiThreadsArg->hold[0]==1) sleep(0);
				
			      }
			  }
		      }
		    
#endif
		    
		    for(i=1;i<threads;i++)
		      pthread_join(workerThreadL[i-1],NULL);		
		    
		    for(i=0;i<grid;i++)
		      appendOmegaResultToFile (alignment, omega, i, i, fpReport, resultType);
		    
		    if(myCorrelationMatrix != NULL)
		      for(i=0; i<matrixSizeMax; ++i)
			if(myCorrelationMatrix[i] != NULL)
			  free(myCorrelationMatrix[i]);
		    
		    if(myCorrelationMatrix != NULL)
		      free(myCorrelationMatrix);
#else

#ifdef _USE_OPENMP_GENERIC

	omp_set_num_threads(threads);

	assert(memLimit>1);
	memLimit_f = (float) memLimit;

	genGridList_size = 0;
	genGridList = (int *)malloc(sizeof(int)*grid);

	for(i=0;i<grid;i++)
		genGridList[i]=-1;


	i=0;
	lvw_i=-1;
	cvw_i=findNextValidOmega(omega, lvw_i, grid);

	while(validGridP(cvw_i,grid))
	{

		genGridList[i++]=cvw_i;

		genLoad_size = compute_genLoad(get_SNP_groupID(omega[cvw_i].leftIndex),get_SNP_groupID(omega[cvw_i].rightIndex), genLoad_size, prev_group_left_index, prev_group_right_index);

		prev_group_left_index = get_SNP_groupID(omega[cvw_i].leftIndex);
		prev_group_right_index = get_SNP_groupID(omega[cvw_i].rightIndex);

		totalmem_genoutput = genLoad_size*group_memsize;

		if(totalmem_genoutput>memLimit_f)
		{
			genGridList_size=i;

			leftSNPindex = omega[genGridList[0]].leftIndex;
			rightSNPindex = omega[genGridList[genGridList_size-1]].rightIndex;
				
			first_group_index = get_SNP_groupID(leftSNPindex);
			last_group_index = get_SNP_groupID(rightSNPindex);

			overlap = prev_last_group_index - first_group_index + 1;
				
			if(prev_last_group_index==-1)
				overlap=-1;

			if(overlap>0)
			{
				overlap_workgroup_map_dim = overlap-1;

				overlap_workgroup_map_number_of_elements = (overlap_workgroup_map_dim*overlap_workgroup_map_dim+overlap_workgroup_map_dim)*0.5;

				if(overlap_workgroup_map_ptr!=NULL) { free(overlap_workgroup_map_ptr); overlap_workgroup_map_ptr=NULL; }
				if(overlap_workgroup_map_ptr_mem!=NULL) { free(overlap_workgroup_map_ptr_mem); overlap_workgroup_map_ptr_mem=NULL; }
	
				overlap_workgroup_map_ptr = (float ***)malloc(sizeof(float **)*overlap_workgroup_map_dim);
				overlap_workgroup_map_ptr_mem = (float **)malloc(sizeof(float *)*overlap_workgroup_map_number_of_elements);

				cur_i=0;
				for(i=0;i<overlap_workgroup_map_dim;i++)
				{
					overlap_workgroup_map_ptr[i]=&(overlap_workgroup_map_ptr_mem[cur_i]);
					cur_i+=i+1;
				
					for(j=0;j<i+1;j++)
						overlap_workgroup_map_ptr[i][j]=NULL;
				}

				overlap_start = prev_last_group_index - overlap + 1 - prev_first_group_index;
				overlap_i=overlap_start;
				overlap_j=overlap_start;
				

				for(i=0;i<overlap_workgroup_map_dim;i++)
				{
					for(j=0;j<=i;j++)
					{
						overlap_workgroup_map_ptr[i][j] = workgroup_map_ptr[overlap_i][overlap_j];
						workgroup_map_ptr[overlap_i][overlap_j]=NULL;
						overlap_j++;	
					}
					overlap_j=overlap_start;
					overlap_i++;
				}

				for(i=0;i<prev_workgroup_map_dim;i++)
				{
					for(j=0;j<=i;j++)
					{
						if(workgroup_map_ptr[i][j]!=NULL)
						{
							free(workgroup_map_ptr[i][j]);
							workgroup_map_ptr[i][j]=NULL;
						}
					}
				}

			}

			workgroup_map_dim=last_group_index+1-first_group_index;
			workgroup_map_number_of_elements = (workgroup_map_dim*workgroup_map_dim+workgroup_map_dim)*0.5;

			if(workgroup_map_ptr!=NULL) { free(workgroup_map_ptr); workgroup_map_ptr=NULL; }
			if(workgroup_map_ptr_mem!=NULL) { free(workgroup_map_ptr_mem); workgroup_map_ptr_mem=NULL; }		

			workgroup_map_ptr = (float ***)malloc(sizeof(float **)*workgroup_map_dim);
			workgroup_map_ptr_mem = (float **)malloc(sizeof(float *)*workgroup_map_number_of_elements);

			cur_i=0;
			for(i=0;i<workgroup_map_dim;i++)
			{
				workgroup_map_ptr[i]=&(workgroup_map_ptr_mem[cur_i]);
				cur_i+=i+1;
			
				for(j=0;j<i+1;j++)
					workgroup_map_ptr[i][j]=NULL;
			}

			prev_startIndex_tmp = get_SNP_groupID(omega[genGridList[0]].leftIndex);
			prev_finishIndex_tmp = get_SNP_groupID(omega[genGridList[0]].rightIndex);

			update_workgroup_map_ptr(workgroup_map_ptr,  prev_startIndex_tmp, prev_finishIndex_tmp, first_group_index);

			for(i=1;i<genGridList_size;i++)
			{
				cur_startIndex_tmp = get_SNP_groupID(omega[genGridList[i]].leftIndex);
				cur_finishIndex_tmp = get_SNP_groupID(omega[genGridList[i]].rightIndex);
				update_workgroup_map_partial_ptr(workgroup_map_ptr, cur_startIndex_tmp, cur_finishIndex_tmp, prev_startIndex_tmp, prev_finishIndex_tmp, first_group_index);
				prev_startIndex_tmp = cur_startIndex_tmp;
				prev_finishIndex_tmp = cur_finishIndex_tmp;
			}

			dp_on_tiles_overlap_ptr (first_group_index, last_group_index, workgroup_map_ptr, overlap_workgroup_map_ptr, overlap, first_group_index, alignment,leftSNPindex, rightSNPindex);



			int * threadload = calloc(threads, sizeof(int));
			int * threadloadloc = calloc(threads, sizeof(int));
			int ** threadgrid = malloc(threads*sizeof(int*));
			for(i=0;i<threads;i++)
				threadgrid[i]=malloc(genGridList_size*sizeof(int));

			int thread2addnxt=0;

			//printf("%d %d\n",genGridList_size, threads);
			//assert(genGridList_size>=threads);
			//if(genGridList_size<threads)

			int skip=0;
			int thread_t = threads;

			if(genGridList_size<threads)
			{
				skip =1;
				thread_t = genGridList_size;
			}

			i=0;
			int minload = omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
			int j;
			for(i=0;i<thread_t;i++)
			{
				threadload[i] = omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
				threadgrid[i][threadloadloc[i]]=genGridList[i];
				threadloadloc[i]++;
				if(threadload[i]<=minload)
				{
					minload = threadload[i];
					thread2addnxt = i;
				}
			}

			if(skip==0)
			for(i=threads;i<genGridList_size;i++)
			{
				threadload[thread2addnxt] += omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
				threadgrid[thread2addnxt][threadloadloc[thread2addnxt]]=genGridList[i];
				threadloadloc[thread2addnxt]++;
				minload = threadload[thread2addnxt];

				for(j=0;j<threads;j++)
				{
					if(threadload[j]<minload)
					{
						minload = threadload[j];
						thread2addnxt = j;
					}
				}
			}

			#pragma omp parallel for private (i)
			for(i=0;i<threads;i++)
			{
				int threadload = threadloadloc[i];
				int l;

				for(l=0;l<threadload;l++)
				{
					computeOmegas_generic (omega, threadgrid[i][l], workgroup_map_ptr, first_group_index);
				}
			}

						
			for(i=0;i<threads;i++)
				free(threadgrid[i]);

			free(threadgrid);
			free(threadload);
			free(threadloadloc);

			prev_first_group_index = first_group_index;
			prev_last_group_index = last_group_index;

			for(i=0;i<genGridList_size;i++)
				genGridList[i]=-1;

			i=0;
			prev_group_left_index = -1;
			prev_group_right_index = -1;
			genLoad_size=0;
			prev_workgroup_map_dim = workgroup_map_dim;

		}

		lvw_i = cvw_i;
		cvw_i=findNextValidOmega(omega, lvw_i, grid);

		if(!validGridP(cvw_i,grid) && genGridList[0]!=-1)
		{				

			genGridList_size=i;

			leftSNPindex = omega[genGridList[0]].leftIndex;
			rightSNPindex = omega[genGridList[genGridList_size-1]].rightIndex;
				
			first_group_index = get_SNP_groupID(leftSNPindex);
			last_group_index = get_SNP_groupID(rightSNPindex);

			overlap = prev_last_group_index - first_group_index + 1;
				
			if(prev_last_group_index==-1)
				overlap=-1;


			if(overlap>0)
			{
				overlap_workgroup_map_dim = overlap-1;

				overlap_workgroup_map_number_of_elements = (overlap_workgroup_map_dim*overlap_workgroup_map_dim+overlap_workgroup_map_dim)*0.5;

				if(overlap_workgroup_map_ptr!=NULL) { free(overlap_workgroup_map_ptr); overlap_workgroup_map_ptr=NULL; }
				if(overlap_workgroup_map_ptr_mem!=NULL) { free(overlap_workgroup_map_ptr_mem); overlap_workgroup_map_ptr_mem=NULL; }
	
				overlap_workgroup_map_ptr = (float ***)malloc(sizeof(float **)*overlap_workgroup_map_dim);
				overlap_workgroup_map_ptr_mem = (float **)malloc(sizeof(float *)*overlap_workgroup_map_number_of_elements);

				cur_i=0;
				for(i=0;i<overlap_workgroup_map_dim;i++)
				{
					overlap_workgroup_map_ptr[i]=&(overlap_workgroup_map_ptr_mem[cur_i]);
					cur_i+=i+1;
				
					for(j=0;j<i+1;j++)
						overlap_workgroup_map_ptr[i][j]=NULL;
				}

				overlap_start = prev_last_group_index - overlap + 1 - prev_first_group_index;
				overlap_i=overlap_start;
				overlap_j=overlap_start;
				

				for(i=0;i<overlap_workgroup_map_dim;i++)
				{
					for(j=0;j<=i;j++)
					{
						overlap_workgroup_map_ptr[i][j] = workgroup_map_ptr[overlap_i][overlap_j];
						workgroup_map_ptr[overlap_i][overlap_j]=NULL;
						overlap_j++;	
					}
					overlap_j=overlap_start;
					overlap_i++;
				}

				for(i=0;i<prev_workgroup_map_dim;i++)
				{
					for(j=0;j<=i;j++)
					{
						if(workgroup_map_ptr[i][j]!=NULL)
						{
							free(workgroup_map_ptr[i][j]);
							workgroup_map_ptr[i][j]=NULL;
						}
					}
				}


			}

			workgroup_map_dim=last_group_index+1-first_group_index;
			workgroup_map_number_of_elements = (workgroup_map_dim*workgroup_map_dim+workgroup_map_dim)*0.5;

			if(workgroup_map_ptr!=NULL) { free(workgroup_map_ptr); workgroup_map_ptr=NULL; }
			if(workgroup_map_ptr_mem!=NULL) { free(workgroup_map_ptr_mem); workgroup_map_ptr_mem=NULL; }		

			workgroup_map_ptr = (float ***)malloc(sizeof(float **)*workgroup_map_dim);
			workgroup_map_ptr_mem = (float **)malloc(sizeof(float *)*workgroup_map_number_of_elements);

			cur_i=0;
			for(i=0;i<workgroup_map_dim;i++)
			{
				workgroup_map_ptr[i]=&(workgroup_map_ptr_mem[cur_i]);
				cur_i+=i+1;
			
				for(j=0;j<i+1;j++)
					workgroup_map_ptr[i][j]=NULL;
			}

			prev_startIndex_tmp = get_SNP_groupID(omega[genGridList[0]].leftIndex);
			prev_finishIndex_tmp = get_SNP_groupID(omega[genGridList[0]].rightIndex);

			update_workgroup_map_ptr(workgroup_map_ptr,  prev_startIndex_tmp, prev_finishIndex_tmp, first_group_index);

			for(i=1;i<genGridList_size;i++)
			{
				cur_startIndex_tmp = get_SNP_groupID(omega[genGridList[i]].leftIndex);
				cur_finishIndex_tmp = get_SNP_groupID(omega[genGridList[i]].rightIndex);
				update_workgroup_map_partial_ptr(workgroup_map_ptr, cur_startIndex_tmp, cur_finishIndex_tmp, prev_startIndex_tmp, prev_finishIndex_tmp, first_group_index);
				prev_startIndex_tmp = cur_startIndex_tmp;
				prev_finishIndex_tmp = cur_finishIndex_tmp;
			}

			dp_on_tiles_overlap_ptr (first_group_index, last_group_index, workgroup_map_ptr, overlap_workgroup_map_ptr, overlap, first_group_index, alignment,leftSNPindex, rightSNPindex);



			int * threadload = calloc(threads, sizeof(int));
			int * threadloadloc = calloc(threads, sizeof(int));
			int ** threadgrid = malloc(threads*sizeof(int*));
			for(i=0;i<threads;i++)
				threadgrid[i]=malloc(genGridList_size*sizeof(int));

			int thread2addnxt=0;

			int skip=0;
			int thread_t = threads;

			if(genGridList_size<threads)
			{
				skip =1;
				thread_t = genGridList_size;
			}

			//assert(genGridList_size>=threads);
			i=0;
			int minload = omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
			int j;
			for(i=0;i<thread_t;i++)
			{
				threadload[i] = omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
				threadgrid[i][threadloadloc[i]]=genGridList[i];
				threadloadloc[i]++;
				if(threadload[i]<=minload)
				{
					minload = threadload[i];
					thread2addnxt = i;
				}
			}
	
			if(skip==0)
			for(i=threads;i<genGridList_size;i++)
			{
				threadload[thread2addnxt] += omega[genGridList[i]].rightIndex - omega[genGridList[i]].leftIndex;
				threadgrid[thread2addnxt][threadloadloc[thread2addnxt]]=genGridList[i];
				threadloadloc[thread2addnxt]++;
				minload = threadload[thread2addnxt];

				for(j=0;j<threads;j++)
				{
					if(threadload[j]<minload)
					{
						minload = threadload[j];
						thread2addnxt = j;
					}
				}
			}


			#pragma omp parallel for private (i)
			for(i=0;i<threads;i++)
			{
				int threadload = threadloadloc[i];
				int l;

				for(l=0;l<threadload;l++)
				{
					computeOmegas_generic (omega, threadgrid[i][l], workgroup_map_ptr, first_group_index);
				}
			}

			for(i=0;i<threads;i++)
				free(threadgrid[i]);

			free(threadgrid);
			free(threadload);
			free(threadloadloc);

			leftSNPindex=-1;
			prev_group_left_index = -1;
			prev_group_right_index = -1;


			for(i=0;i<workgroup_map_dim;i++)
			{
				for(j=0;j<=i;j++)
				{
					if(workgroup_map_ptr[i][j]!=NULL)
					{
						free(workgroup_map_ptr[i][j]);
						workgroup_map_ptr[i][j]=NULL;
					}
				}
			}

			if(workgroup_map_ptr!=NULL)
			{
				free(workgroup_map_ptr);
				workgroup_map_ptr=NULL;
			}
			if(workgroup_map_ptr_mem!=NULL)
			{
				free(workgroup_map_ptr_mem);
				workgroup_map_ptr_mem=NULL;
			}				
			i=0;


		}
		
	}


				//for(i=0;i<overlap_workgroup_map_dim;i++)
				//{
				//	for(j=0;j<=i;j++)
				//	{
				//		if(overlap_workgroup_map_ptr[i][j]!=NULL)
				//			assert(0);
				//	}
				//}

	for(i=0;i<grid;i++)
		appendOmegaResultToFile (alignment, omega, i, i, fpReport, resultType);


		
//		free(total_dp_init_time);
//		free(total_dp_update_time);
//		free(total_omega_values_time);
		free(genGridList);
#else
		    
		    alignment->correlationMatrix = createCorrelationMatrix(alignment->correlationMatrix,matrixSizeMax);
		    
		    lvw_i=-1;
		    
		    for(i=0;i<grid;i++)
		      {
			cvw_i=findNextValidOmega(omega, lvw_i, grid);
			
			if(validGridP(cvw_i,grid))
			  {		
			    overlapCorrelationMatrixAdditions (alignment, omega, lvw_i, cvw_i, 
							       &firstRowToCopy, &firstRowToCompute, &firstRowToAdd);
			    
			    shiftCorrelationMatrixValues (omega, lvw_i, cvw_i, firstRowToCopy, alignment->correlationMatrix);
			    
			    computeCorrelationMatrixPairwise (alignment, omega, cvw_i, firstRowToCompute, functionData, NULL,NULL);					
			    
			    applyCorrelationMatrixAdditions (omega, cvw_i,firstRowToAdd,alignment->correlationMatrix);
			    
			    computeOmegas (alignment, omega, cvw_i, functionData,NULL);
			    
			    lvw_i = cvw_i;
			  }
			
			appendOmegaResultToFile (alignment, omega, i, i+1, fpReport, resultType);
		      }			
#endif		    
#endif

	    }
	  else
	    {
		    fprintf(stdout,"\n\n\t\tProcessing:\t\tNOT POSSIBLE\n\n");

		    fprintf(fpInfo,"\n\n\t\tProcessing:\t\tNOT POSSIBLE\n\n");
		    
	    }
	  
	  	time1 = gettime();
		
		totalTimeL = time1-time0; 
		
		fprintf(stdout, "%f\n\n",totalTimeL);

		fprintf(fpInfo, "%f\n\n",totalTimeL);

		maxOmegaResultReport (&maxomegaRealPos, &maxomegamaxValue, &maxomegamaxLeftIndex, &maxomegamaxRightIndex, grid, alignment, omega, fpInfo);

		freeAlignment(alignment, matrixSizeMax);

		time0 = gettime();
		
		nxt_alignment = findNextAlignment(fpIn,fileFormat);
		
		alignmentIndex++;	
	}


#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
	terminateWorkerThreads(workerThreadL,threadData);
#endif
	free(threadData);
#endif
	totalTimeG1 = gettime();

	fprintf(stdout, " Total elapsed time %f seconds.\n\n\n",totalTimeG1-totalTimeG0);

	fprintf(fpInfo, " Total elapsed time %f seconds.\n\n\n",totalTimeG1-totalTimeG0);

	if (fileFormat==MS_FORMAT || fileFormat == MACS_FORMAT) 
		fclose(fpWarnings);

	if(alignment->VCFsample_valid!=NULL)
		free(alignment->VCFsample_valid);

	free(alignment);
	free(omega);
	free(recfile);
	
	if(fpIn!=NULL)
		fclose(fpIn);

	fclose(fpInfo); 

	if(fpReport!=NULL)
		fclose(fpReport);

	if(fpVCFsamples!=NULL)
		fclose(fpVCFsamples);

	return 0;
}
