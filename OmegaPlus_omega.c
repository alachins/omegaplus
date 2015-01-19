
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

inline int min(int a, int b)
{
  if(a < b)
    return a;
  return b;
}

inline int max(int a, int b)
{
  if(a > b)
    return a;
  return b;
}

int moveright(int start, int* intsnps, int snpsize, int pos)
{
	int i = start;

	while(i<snpsize && intsnps[i] <= pos)
		++i;

	return i-1;
}

int moveleft(int start, int* intsnps, int snpsize, int pos)
{
	int i = start;

	while(i>=0 && intsnps[i] >= pos)
		--i;

	return i+1;
}

void boundaries(omega_struct* omstruct, 
		int grid, 
		int* intsnps, 
		int snpsize, 
		int minw, 
		int minsnps, /* the minimum number of snps a window is allowed to have (usually 2)*/
		int maxw, 
		int *changemaxw, 
		int *maxSizeMatrix)
{
	int startsnp = intsnps[0];

	int endsnp = intsnps[snpsize - 1];

	assert(grid>1);

	
	int i, left_boundary, right_boundary, left_min_boundary, right_min_boundary,
	    maxsizematrix = 0, dif;

	float step = (float)(endsnp - startsnp)/(grid-1), omega_position = (float)startsnp;

	if(step < 1)
	{
		fprintf(stdout, "\n\n WARNING: Gridsize is too large (%d) for the region between the first and last SNPs (%d)\n\n", grid, endsnp-startsnp+1);
		
		/* assert(endsnp-startsnp >= grid); */
		
	}

	int checksnp=0;

	for(i=0; i<grid; ++i)
	{
		omstruct[i].omegaRealPos = omega_position;

		left_boundary = omega_position - maxw;

		right_boundary = omega_position + maxw;

		left_min_boundary = omega_position - minw;

		right_min_boundary = omega_position + minw;

		omstruct[i].omegaPos = moveright(checksnp, intsnps, snpsize, omega_position);

		omstruct[i].leftIndex = moveleft(omstruct[i].omegaPos, intsnps, snpsize, left_boundary);

		omstruct[i].rightIndex = moveright(omstruct[i].omegaPos + 1, intsnps, snpsize, right_boundary);

		omstruct[i].leftminIndex = moveleft(omstruct[i].omegaPos, intsnps, snpsize, left_min_boundary);

		omstruct[i].rightminIndex = moveright(omstruct[i].omegaPos, intsnps, snpsize, right_min_boundary);

		omstruct[i].maxValue = 0.0;

		omstruct[i].maxLeftIndex = 0;
		omstruct[i].maxRightIndex = 0;

		/* change the minw if it is too small */
		while(omstruct[i].omegaPos - omstruct[i].leftminIndex + 1  < minsnps)
			--omstruct[i].leftminIndex;

		while( omstruct[i].rightminIndex - omstruct[i].omegaPos  < minsnps)
			++omstruct[i].rightminIndex;

		omstruct[i].valid = 1;

		/* do some checks for the minimum window boundaries */
		if(omstruct[i].leftminIndex < omstruct[i].leftIndex || omstruct[i].rightminIndex > omstruct[i].rightIndex)
			omstruct[i].valid = 0;

		/* do some checks for the maximum window boundaries */
		if(omstruct[i].omegaPos + 1 >= omstruct[i].rightIndex || omstruct[i].omegaPos <= omstruct[i].leftIndex)
			omstruct[i].valid = 0;

		dif = omstruct[i].rightIndex - omstruct[i].leftIndex + 1;

		if(omstruct[i].valid == 1 && dif  > MAXSIZE)
		{
			*changemaxw = 1;
			return;
		}

		if(dif > maxsizematrix)
			maxsizematrix = dif;

		omega_position += step;
	}

	*maxSizeMatrix = maxsizematrix;
}

int findOmegaBounds (alignment_struct * alignment, omega_struct * omega, int grid, int * maxw, int minw, int minsnps)
{
	int changemaxw = 1,
            matrixSizeMax=0,
            maxwc = *maxw;

	while (changemaxw)
	{
		changemaxw = 0;

		boundaries(omega, grid, alignment-> positionsInd, alignment->segsites, 	minw, minsnps, maxwc, &changemaxw, &matrixSizeMax);
	
		if(changemaxw)
			maxwc = (int)(DECREASE * maxwc);
	}

	*maxw = maxwc;

	return matrixSizeMax;	
}

int findNextValidOmega(omega_struct *omega, int lvw_i, int grid)
{
	int i=lvw_i+1;
	
	while(i<grid && omega[i].valid!=1)
		i++;
	
	return i;		
}

int validGridP(int cvw_i, int grid)
{
	if(cvw_i>=0 && cvw_i<grid)
		return 1;
	
	return 0;
}

float computeOmega (float LS, float RS, float TS, int k, int ksel2, int m, int msel2)
{
	float   numerator = (LS + RS) / (ksel2 + msel2);

	float denominator = (TS - LS - RS) / (k*m) + DENOMINATOR_OFFSET;

	float omega =  numerator / denominator;

	return omega;
}

void computeOmegaValues (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, void * threadData)
{
	float LS, RS, TS, tmpW = 0.0, maxW=0.0;

	int i, j, ksel2, msel2, k, m, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex,

	rightMinIndexORIG = rightMinIndex,

	rightMaxIndexORIG = rightMaxIndex;

#ifdef _UNROLL
	int vw = 2;
	int iter, iterations = (rightMaxIndex-rightMinIndex+1) / vw;
	int finaliterations = (rightMaxIndex-rightMinIndex+1) % vw;
	int omegaSNIPIndexPlusOne = omegaSNIPIndex + 1;
	float maxW_0=0.0, maxW_1=0.0;
	int maxLeftIndex_0=0, maxRightIndex_0=0, maxLeftIndex_1=0, maxRightIndex_1=0;
#endif

#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT

	threadData_t * threadDataL = (threadData_t *) threadData;

	int t=0, 	

	tid = threadDataL->threadID,
	
	threads = threadDataL->threadTOTAL;
#endif
#endif	
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{

#ifdef _USE_PTHREADS

#ifndef _USE_PTHREADS_MEMINT
	  
	if(t%threads==tid)
		{
#endif
#endif	
			LS = correlationMatrix[omegaSNIPIndex][i];

			k = omegaSNIPIndex - i + 1;
		
			ksel2 = (k * (k-1)) / 2;

			if(borderTol > 0)
			{
				rightMinIndex = rightMinIndexORIG;

				rightMaxIndex = rightMaxIndexORIG;

			    //fprintf(stderr, "---------------------\nrightMinIndex: %d, rightMaxIndex: %d\n", rightMinIndex, rightMaxIndex);

				int leftSNPs = omegaSNIPIndex - i + 1;
				int equalRightPosition = omegaSNIPIndex + leftSNPs;
  
				rightMinIndex = max(rightMinIndex, equalRightPosition - borderTol);
				rightMaxIndex = min(rightMaxIndex, equalRightPosition + borderTol);
			    
			}
#ifdef _UNROLL
			j = rightMinIndex;

			for(iter=0;iter<iterations;iter++)
			{

				int j_0 = j;
				int j_1 = j_0 + 1;

				j=j+vw;

				float RS_0 = correlationMatrix[j_0][omegaSNIPIndexPlusOne];
				float RS_1 = correlationMatrix[j_1][omegaSNIPIndexPlusOne];

				int m_0 = j_0 - omegaSNIPIndex;
				int m_1 = j_1 - omegaSNIPIndex;
	
				int mmin1_0 = m_0 - 1;
				int mmin1_1 = m_1 - 1;

				int mmin1multm_0 = mmin1_0 * m_0;
				int mmin1multm_1 = mmin1_1 * m_1;

				int msel2_0 = mmin1multm_0 / 2;
				int msel2_1 = mmin1multm_1 / 2;
					
				float TS_0 = correlationMatrix[j_0][i];
				float TS_1 = correlationMatrix[j_1][i];

				int ksel2plusmsel2_0 = ksel2 + msel2_0;
				int ksel2plusmsel2_1 = ksel2 + msel2_1;

				float LSplusRS_0 = LS+RS_0;
				float LSplusRS_1 = LS+RS_1;


				float numerator_0 = LSplusRS_0 / ksel2plusmsel2_0;
				float numerator_1 = LSplusRS_1 / ksel2plusmsel2_1;

				float TSminusLS_0 = TS_0 - LS;
				float TSminusLS_1 = TS_1 - LS;

				float TSminusLSminusRS_0 = TSminusLS_0 - RS_0;
				float TSminusLSminusRS_1 = TSminusLS_1 - RS_1;

				int kmultm_0 = k * m_0;
				int kmultm_1 = k * m_1;

				float denominator_0 = TSminusLSminusRS_0 / kmultm_0;
				float denominator_1 = TSminusLSminusRS_1 / kmultm_1;

				float denominatorplusOFFSET_0 = denominator_0 + DENOMINATOR_OFFSET;
				float denominatorplusOFFSET_1 = denominator_1 + DENOMINATOR_OFFSET;

				float tmpW_0 =  numerator_0 / denominatorplusOFFSET_0;
				float tmpW_1 =  numerator_1 / denominatorplusOFFSET_1;
	
				if(tmpW_0>maxW_0)
				{
					maxW_0 = tmpW_0;
					maxLeftIndex_0 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_0 = j_0 + omega[omegaIndex].leftIndex;
				}				
	
				if(tmpW_1>maxW_1)
				{
					maxW_1 = tmpW_1;
					maxLeftIndex_1 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_1 = j_1 + omega[omegaIndex].leftIndex;
				}			

			}

			if(maxW_0>maxW_1)
			{
				maxW = maxW_0;
				maxLeftIndex = maxLeftIndex_0;
				maxRightIndex = maxRightIndex_0;
			}
			else
			{
				maxW = maxW_1;
				maxLeftIndex = maxLeftIndex_1;
				maxRightIndex = maxRightIndex_1;
			}


			if(finaliterations!=0)			
			{
			
				RS = correlationMatrix[j][omegaSNIPIndexPlusOne];
				m = j - omegaSNIPIndex;	
				int mmin1 = m - 1;
				int mmin1multm = mmin1 * m;
				msel2 = mmin1multm / 2;					
				TS = correlationMatrix[j][i];
				int ksel2plusmsel2 = ksel2 + msel2;
				float LSplusRS = LS+RS;
				float numerator = LSplusRS / ksel2plusmsel2;
				float TSminusLS = TS - LS;
				float TSminusLSminusRS = TSminusLS - RS;
				int kmultm = k * m;
				float denominator = TSminusLSminusRS / kmultm;
				float denominatorplusOFFSET = denominator + DENOMINATOR_OFFSET; 
				tmpW =  numerator / denominatorplusOFFSET;
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}

			}

			maxW_0 = maxW;
			maxLeftIndex_0 = maxLeftIndex;
			maxRightIndex_0 = maxRightIndex;
	
			maxW_1 = maxW;
			maxLeftIndex_1 = maxLeftIndex;
			maxRightIndex_1 = maxRightIndex;		

#else

		
			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				RS = correlationMatrix[j][omegaSNIPIndex+1];

				m = j - omegaSNIPIndex;
	
				msel2 = (m * (m-1)) / 2;
					
				TS = correlationMatrix[j][i];

				tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}
			}

#endif

#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
		}
		t++;
#endif
#endif
	}
#ifdef _USE_PTHREADS
#ifndef _USE_PTHREADS_MEMINT
	threadDataL->threadArgCO->maxValue = maxW;
	threadDataL->threadArgCO->maxLeftIndex = maxLeftIndex;
	threadDataL->threadArgCO->maxRightIndex = maxRightIndex;	
#else
	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
#endif
#else
	omega[omegaIndex].maxValue = maxW;
	omega[omegaIndex].maxLeftIndex  = maxLeftIndex;
	omega[omegaIndex].maxRightIndex = maxRightIndex;
#endif
}

#ifdef _USE_PTHREADS
#ifdef _USE_PTHREADS_MEMINT
void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues (omega, omegaIndex, correlationMatrix, NULL);
}
#ifdef _USE_PTHREADS_MULTI
void initmaxOmegaValueThreadsMULTI(int size)
{
	int j;

	for(j=0;j<size;j++)
		maxOmegaValueThreadsMULTI[j]=-1.0;
}

void getmaxOmegaValueThreadsMULTI(omega_struct * omega, int cvw_i, int size)
{

	int j;

	omega[cvw_i].maxValue = maxOmegaValueThreadsMULTI[0];
	omega[cvw_i].maxLeftIndex  = maxOmegaLeftIndexThreadsMULTI[0];
	omega[cvw_i].maxRightIndex = maxOmegaRightIndexThreadsMULTI[0];

	for(j=1;j<size;j++)
		if(maxOmegaValueThreadsMULTI[j]>omega[cvw_i].maxValue)
		{
			omega[cvw_i].maxValue = maxOmegaValueThreadsMULTI[j];
			omega[cvw_i].maxLeftIndex  = maxOmegaLeftIndexThreadsMULTI[j];
			omega[cvw_i].maxRightIndex = maxOmegaRightIndexThreadsMULTI[j];
		}
}

void computeOmegaValuesMULTI (omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, int ttid, int tthreads, int tid)
{
	float LS, RS, TS, tmpW = 0, maxW=0;

	int i, j, ksel2, msel2, k, m, maxLeftIndex=0, maxRightIndex=0,
	
	omegaSNIPIndex = omega[omegaIndex].omegaPos - omega[omegaIndex].leftIndex,

	leftMinIndex = omega[omegaIndex].leftminIndex - omega[omegaIndex].leftIndex,

	leftMaxIndex = omega[omegaIndex].leftIndex - omega[omegaIndex].leftIndex,
	
	rightMinIndex = omega[omegaIndex].rightminIndex - omega[omegaIndex].leftIndex,

	rightMaxIndex = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex,

	rightMinIndexORIG = rightMinIndex,

	rightMaxIndexORIG = rightMaxIndex;

	int t =0;

#ifdef _UNROLL
	int vw = 2;
	int iter, iterations = (rightMaxIndex-rightMinIndex+1) / vw;
	int finaliterations = (rightMaxIndex-rightMinIndex+1) % vw;
	int omegaSNIPIndexPlusOne = omegaSNIPIndex + 1;
	float maxW_0=0.0, maxW_1=0.0;
	int maxLeftIndex_0=0, maxRightIndex_0=0, maxLeftIndex_1=0, maxRightIndex_1=0;
#endif
	
	for (i=leftMinIndex;i>=leftMaxIndex;i--) // Left Side
	{

		if(t%tthreads==ttid)
		{
			LS = correlationMatrix[omegaSNIPIndex][i];

			k = omegaSNIPIndex - i + 1;
		
			ksel2 = (k * (k-1)) / 2;

			if(borderTol > 0)
			{
				rightMinIndex = rightMinIndexORIG;

				rightMaxIndex = rightMaxIndexORIG;

			    //fprintf(stderr, "---------------------\nrightMinIndex: %d, rightMaxIndex: %d\n", rightMinIndex, rightMaxIndex);

				int leftSNPs = omegaSNIPIndex - i + 1;
				int equalRightPosition = omegaSNIPIndex + leftSNPs;
  
				rightMinIndex = max(rightMinIndex, equalRightPosition - borderTol);
				rightMaxIndex = min(rightMaxIndex, equalRightPosition + borderTol);
			    
			}
		
			#ifdef _UNROLL
			j = rightMinIndex;

			for(iter=0;iter<iterations;iter++)
			{

				int j_0 = j;
				int j_1 = j_0 + 1;

				j=j+vw;

				float RS_0 = correlationMatrix[j_0][omegaSNIPIndexPlusOne];
				float RS_1 = correlationMatrix[j_1][omegaSNIPIndexPlusOne];

				int m_0 = j_0 - omegaSNIPIndex;
				int m_1 = j_1 - omegaSNIPIndex;
	
				int mmin1_0 = m_0 - 1;
				int mmin1_1 = m_1 - 1;

				int mmin1multm_0 = mmin1_0 * m_0;
				int mmin1multm_1 = mmin1_1 * m_1;

				int msel2_0 = mmin1multm_0 / 2;
				int msel2_1 = mmin1multm_1 / 2;
					
				float TS_0 = correlationMatrix[j_0][i];
				float TS_1 = correlationMatrix[j_1][i];

				int ksel2plusmsel2_0 = ksel2 + msel2_0;
				int ksel2plusmsel2_1 = ksel2 + msel2_1;

				float LSplusRS_0 = LS+RS_0;
				float LSplusRS_1 = LS+RS_1;


				float numerator_0 = LSplusRS_0 / ksel2plusmsel2_0;
				float numerator_1 = LSplusRS_1 / ksel2plusmsel2_1;

				float TSminusLS_0 = TS_0 - LS;
				float TSminusLS_1 = TS_1 - LS;

				float TSminusLSminusRS_0 = TSminusLS_0 - RS_0;
				float TSminusLSminusRS_1 = TSminusLS_1 - RS_1;

				int kmultm_0 = k * m_0;
				int kmultm_1 = k * m_1;

				float denominator_0 = TSminusLSminusRS_0 / kmultm_0;
				float denominator_1 = TSminusLSminusRS_1 / kmultm_1;

				float denominatorplusOFFSET_0 = denominator_0 + DENOMINATOR_OFFSET;
				float denominatorplusOFFSET_1 = denominator_1 + DENOMINATOR_OFFSET;

				float tmpW_0 =  numerator_0 / denominatorplusOFFSET_0;
				float tmpW_1 =  numerator_1 / denominatorplusOFFSET_1;
	
				if(tmpW_0>maxW_0)
				{
					maxW_0 = tmpW_0;
					maxLeftIndex_0 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_0 = j_0 + omega[omegaIndex].leftIndex;
				}				
	
				if(tmpW_1>maxW_1)
				{
					maxW_1 = tmpW_1;
					maxLeftIndex_1 = i + omega[omegaIndex].leftIndex;
					maxRightIndex_1 = j_1 + omega[omegaIndex].leftIndex;
				}			

			}

			if(maxW_0>maxW_1)
			{
				maxW = maxW_0;
				maxLeftIndex = maxLeftIndex_0;
				maxRightIndex = maxRightIndex_0;
			}
			else
			{
				maxW = maxW_1;
				maxLeftIndex = maxLeftIndex_1;
				maxRightIndex = maxRightIndex_1;
			}


			if(finaliterations!=0)			
			{
			
				RS = correlationMatrix[j][omegaSNIPIndexPlusOne];
				m = j - omegaSNIPIndex;	
				int mmin1 = m - 1;
				int mmin1multm = mmin1 * m;
				msel2 = mmin1multm / 2;					
				TS = correlationMatrix[j][i];
				int ksel2plusmsel2 = ksel2 + msel2;
				float LSplusRS = LS+RS;
				float numerator = LSplusRS / ksel2plusmsel2;
				float TSminusLS = TS - LS;
				float TSminusLSminusRS = TSminusLS - RS;
				int kmultm = k * m;
				float denominator = TSminusLSminusRS / kmultm;
				float denominatorplusOFFSET = denominator + DENOMINATOR_OFFSET; 
				tmpW =  numerator / denominatorplusOFFSET;
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}

			}

			maxW_0 = maxW;
			maxLeftIndex_0 = maxLeftIndex;
			maxRightIndex_0 = maxRightIndex;
	
			maxW_1 = maxW;
			maxLeftIndex_1 = maxLeftIndex;
			maxRightIndex_1 = maxRightIndex;		

#else

		
			for(j=rightMinIndex;j<=rightMaxIndex;j++) // Right Side
			{
				RS = correlationMatrix[j][omegaSNIPIndex+1];

				m = j - omegaSNIPIndex;
	
				msel2 = (m * (m-1)) / 2;
					
				TS = correlationMatrix[j][i];

				tmpW = computeOmega(LS, RS, TS, k, ksel2, m, msel2);
	
				if(tmpW>maxW)
				{
					maxW = tmpW;
					maxLeftIndex = i + omega[omegaIndex].leftIndex;
					maxRightIndex = j + omega[omegaIndex].leftIndex;
				}
			}

#endif

		}
		t++;
	}

	maxOmegaValueThreadsMULTI[tid]=maxW;
	maxOmegaLeftIndexThreadsMULTI[tid]=maxLeftIndex;
	maxOmegaRightIndexThreadsMULTI[tid]=maxRightIndex;
}

void computeOmegasMULTI (alignment_struct * alignment, omega_struct * omega, int omegaIndex, cor_t ** correlationMatrix, int ttid, int tthreads, int tid)
{
	computeOmegaValuesMULTI (omega, omegaIndex, correlationMatrix, ttid, tthreads, tid);
}

#endif
#else
void computeOmegasThread (alignment_struct * alignment, omega_struct * omega, int omegaIndex, threadData_t * threadData)
{	
	computeOmegaValues (omega, omegaIndex, alignment->correlationMatrix, threadData);
}

void omegasThread(threadData_t * currentThread)
{
	alignment_struct * alignment = currentThread->threadArgCO->alignment;

	omega_struct * omega = currentThread->threadArgCO->omega;

	int omegaIndex = currentThread->threadArgCO->omegaIndex;

        computeOmegasThread (alignment,omega,omegaIndex,currentThread);	
}

void setThreadArgumentsCO(threadData_t * threadData, int tid, alignment_struct * alignment, omega_struct * omega, int omegaIndex)
{
	threadData[tid].threadArgCO->alignment=alignment;
	threadData[tid].threadArgCO->omega=omega;
	threadData[tid].threadArgCO->omegaIndex=omegaIndex;
}

void getAllThreadMaxOmega(threadData_t * threadData, omega_struct * omega, int omegaIndex)
{
	int i, threads = threadData[0].threadTOTAL;

	omega[omegaIndex].maxValue = threadData[0].threadArgCO->maxValue;
	omega[omegaIndex].maxLeftIndex = threadData[0].threadArgCO->maxLeftIndex;
	omega[omegaIndex].maxRightIndex = threadData[0].threadArgCO->maxRightIndex;

	for(i=1;i<threads;i++)
	{
		if(threadData[i].threadArgCO->maxValue>omega[omegaIndex].maxValue)
		{
			omega[omegaIndex].maxValue = threadData[i].threadArgCO->maxValue;
			omega[omegaIndex].maxLeftIndex = threadData[i].threadArgCO->maxLeftIndex;
			omega[omegaIndex].maxRightIndex = threadData[i].threadArgCO->maxRightIndex;
		}
	}
}

void computeOmegaValues_THREADS (alignment_struct * alignment, omega_struct * omega, int omegaIndex, threadData_t * threadData)
{
	int i, threads = threadData[0].threadTOTAL;

	for(i=0;i<threads;i++)
		setThreadArgumentsCO(threadData, i, alignment, omega, omegaIndex);

	startThreadOperations(threadData, COMPUTEOMEGAS);

	getAllThreadMaxOmega(threadData,omega, omegaIndex);
}

void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{

	threadData_t * threadDataL = (threadData_t *) threadData;

	computeOmegaValues_THREADS (alignment, omega, omegaIndex, threadDataL);
}
#endif
#else
void computeOmegas (alignment_struct * alignment, omega_struct * omega, int omegaIndex, void * threadData, cor_t ** correlationMatrix)
{
	computeOmegaValues (omega, omegaIndex, alignment->correlationMatrix, NULL);
}
#endif

void appendOmegaResultToFile (alignment_struct * alignment, omega_struct * omega, int omegaIndex, int gridIndex, FILE * fpOut, int resultType)
{
	if (resultType==RESULTS_ALL)
		if(omega[omegaIndex].valid)
			fprintf(fpOut,"%.4f\t%f\t%d\t%d\t%d\n", omega[omegaIndex].omegaRealPos, omega[omegaIndex].maxValue, alignment->positionsInd[omega[omegaIndex].maxLeftIndex], 
							      alignment->positionsInd[omega[omegaIndex].maxRightIndex], omega[omegaIndex].valid );
		else
			fprintf(fpOut,"%.4f\t%f\t%d\t%d\t%d\n", omega[omegaIndex].omegaRealPos, 0.0, 0, 0, 0);
	else
		if(omega[omegaIndex].valid)
			fprintf(fpOut,"%.4f\t%f\n", omega[omegaIndex].omegaRealPos, omega[omegaIndex].maxValue);
		else
			fprintf(fpOut,"%.4f\t%f\n", omega[omegaIndex].omegaRealPos, 0.0);
}

void maxOmegaResultReport (float * maxomegaRealPos, float * maxomegamaxValue, int * maxomegamaxLeftIndex, int * maxomegamaxRightIndex, int grid, alignment_struct * alignment, omega_struct * omega, FILE * fpInfo)
{
	*maxomegaRealPos = 0.0;
	*maxomegamaxValue = 0.0;
	*maxomegamaxLeftIndex = -1;
	*maxomegamaxRightIndex = -1;

	int i;

	for(i=0;i<grid;i++)
	{
		if(omega[i].maxValue>*maxomegamaxValue)
		{
			*maxomegamaxValue = omega[i].maxValue;
			*maxomegaRealPos = omega[i].omegaRealPos;
			*maxomegamaxLeftIndex = alignment->positionsInd[omega[i].maxLeftIndex];
			*maxomegamaxRightIndex = alignment->positionsInd[omega[i].maxRightIndex];		
		}
	}

	fprintf(stdout, "\t\tMax Omega:\t\t%f\n",*maxomegamaxValue);
	fprintf(stdout, "\t\tLocation:\t\t%.2f\n",*maxomegaRealPos);
	fprintf(stdout, "\t\tLeftmost SNP:\t\t%d\n",*maxomegamaxLeftIndex);
	fprintf(stdout, "\t\tRightmost SNP:\t\t%d\n\n\n",*maxomegamaxRightIndex);

	fprintf(fpInfo, "\t\tMax Omega:\t\t%f\n",*maxomegamaxValue);
	fprintf(fpInfo, "\t\tLocation:\t\t%.2f\n",*maxomegaRealPos);
	fprintf(fpInfo, "\t\tLeftmost SNP:\t\t%d\n",*maxomegamaxLeftIndex);
	fprintf(fpInfo, "\t\tRightmost SNP:\t\t%d\n\n\n",*maxomegamaxRightIndex);
	
}
