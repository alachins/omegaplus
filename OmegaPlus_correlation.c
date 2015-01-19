
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

/* LUT-based population counter */
int iterated_bitcount(unsigned int n)
{
	int count=0;    
    
	while(n)
	{
		count += n & 0x1u ;    
		n >>= 1 ;
	}

	return count;
}

#ifdef _SHARED
void compute_bits_in_16bits(void)
{
	unsigned int i;    
    
	for (i = 0; i < (0x1u<<16); i++)
		bits_in_16bits[i] = iterated_bitcount(i);
}

unsigned int precomputed16_bitcount (unsigned int n)
{
	/* works only for 32-bit unsigned int*/	    
	return bits_in_16bits [n         & 0xffffu]	
	    +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}
#else
void compute_bits_in_16bitsLocal(char * bits_in_16bitsLocal)
{
	unsigned int i;    
    
	for (i = 0; i < (0x1u<<16); i++)
		bits_in_16bitsLocal[i] = iterated_bitcount(i);
}

unsigned int precomputed16_bitcountLocal (char * bits_in_16bitsLocal, unsigned int n)
{
	/* works only for 32-bit unsigned int*/	    
	return bits_in_16bitsLocal [n         & 0xffffu]	
	    +  bits_in_16bitsLocal [(n >> 16) & 0xffffu] ;
}
#endif

cor_t computeCorrelationValueBIN_RSQUARE(int sequences, unsigned int * accumXvec)
{	
	cor_t result;
	cor_t sequencesTotal = (cor_t) sequences;

	if (sequencesTotal==0)
		return 0.;

	cor_t numerator=0.0, 
              denominator=0.0;

	cor_t denomA = (cor_t)accumXvec[0]/sequencesTotal;
	cor_t denomB = (cor_t)accumXvec[1]/sequencesTotal;
	cor_t denomC = (cor_t)accumXvec[2]/sequencesTotal;
	cor_t denomD = (cor_t)accumXvec[3]/sequencesTotal;

	denominator = denomA * denomB * denomC * denomD;

	if (denominator==0.0)
		return 0.0;

	numerator  = ((cor_t)accumXvec[4])/sequencesTotal;
	numerator -= ((cor_t)accumXvec[1]/(cor_t)sequences) * ((cor_t)accumXvec[3]/(cor_t)sequences);
	numerator *= numerator;

	result = (sequences * numerator) /denominator;

	return result;
}

cor_t computeCorrelationValueBIN_DOM(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	if ( ( (m1 >= 0.5) && (m2 >= 0.5) ) || 	( (m1 < 0.5) && (m2 < 0.5 ) ) )
		result =  m11*m22 - m12*m21;
	else if( ( (m1 < 0.5) && (m2 >= 0.5 ) )|| ( (m1 >= 0.5) && (m2 < 0.5) ) )
		result =  m21*m12 - m11*m22;
	else{
		fprintf(stderr, "m1: %e, m2: %e, m11: %e, m12: %e, m21: %e, m22: %e\n", m1, m2, m11, m12, m21, m22);
		assert(0);
	}

	return result;
}

cor_t computeCorrelationValueBIN_JUSTD(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	result = m11 - m1*m2;

	return result;
}

cor_t computeCorrelationValueBIN_DOM2(int sequences, unsigned int * accumXvec)
{
	cor_t result;

	cor_t m1  = (cor_t)accumXvec[1];
	cor_t m11 = (cor_t)accumXvec[4];
	cor_t m12 = m1 - m11;
	cor_t m2  = (cor_t)accumXvec[3];
	cor_t m21 = m2 - m11;
	cor_t m22 = (cor_t)accumXvec[0] - m21;

	if(!(accumXvec[2]-m12==m22))	
	{
		printf("Equality Failed! (%f - %f)\n",accumXvec[2]-m12,m22);
		assert(accumXvec[2]-m12==m22);
	}

	m11/=sequences;
	m22/=sequences;
	m12/=sequences;
	m21/=sequences;
	m1/=sequences;
	m2/=sequences;

	if ( ( (m1 >= 0.5) && (m2 >= 0.5) ) || 	( (m1 < 0.5) && (m2 < 0.5 ) ) )
		result =  m11*m22 - m12*m21;
	else if( ( (m1 < 0.5) && (m2 >= 0.5 ) )|| ( (m1 >= 0.5) && (m2 < 0.5) ) )
		result =  m21*m12 - m11*m22;
	else{
		fprintf(stderr, "m1: %e, m2: %e, m11: %e, m12: %e, m21: %e, m22: %e\n", m1, m2, m11, m12, m21, m22);
		assert(0);
	}

  	cor_t denom = (m1*m2*(1.-m1)*(1.-m2));

  	result = abs(result)/denom;

	return result;
}

cor_t computeCorrelationValueBIN(int sequences, unsigned int * accumXvec)
{
	if(linkage_disequilibrium == RSQUARE)
		return computeCorrelationValueBIN_RSQUARE(sequences, accumXvec);

	if(linkage_disequilibrium == DOM)
		return computeCorrelationValueBIN_DOM(sequences, accumXvec);

	if(linkage_disequilibrium == ABSDOM)
		return ABS(computeCorrelationValueBIN_DOM(sequences, accumXvec));

	if(linkage_disequilibrium == JUSTD)
		return computeCorrelationValueBIN_JUSTD(sequences, accumXvec);

	if(linkage_disequilibrium == ABSD)
		return ABS(computeCorrelationValueBIN_JUSTD(sequences, accumXvec));

	if(linkage_disequilibrium == ABSDOM2)
		return ABS(computeCorrelationValueBIN_DOM2(sequences, accumXvec));

	assert(linkage_disequilibrium==999);

	return 0.0;
}

cor_t computeCorrelationValueDNA(int sequences, cor_t pairwiseCorrelationMatrix[4][4], unsigned int * valid)
{
	int i, j;
	int valid1 = 4, valid2 = 4;

	cor_t correlation = 0.0;

	for(i=0; i<4; ++i)
		if(valid[i] == 0)
			valid1--;

	for(j=0; j<4; ++j)
		if(valid[j+4] == 0)
			valid2--;

	for(i=0; i<4; ++i)
		if(valid[i]!=0)
			for(j=0;j<4;j++)
				if(valid[i]!=0 && valid[j+4]!=0)
					correlation += pairwiseCorrelationMatrix[i][j];

	/* the formula is from: Zaykin et al. 2008, DOI: 10.1534/genetics.108.089409 
	In Zaykin et al. they multiply by (valid1 - 1)(valid2 - 1) * sequences / (valid1 * valid2)
	*/

	if(valid1*valid2!=0)
		correlation *= ( (((valid1 - 1.0)*(valid2 - 1.0))*sequences ) /(valid1 * valid2));

	return correlation;    
}

#ifdef _SHARED
void count01Combs (int total, unsigned int inputL, unsigned int inputR, unsigned int * accumXvec)
{
	int limit=ENTRIES_PER_INT;
	
	if (total<limit)
		limit = total;

	int t1 = precomputed16_bitcount (inputL);
	int t2 = precomputed16_bitcount (inputR);

	unsigned int combAND = inputL & inputR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += limit - t1;
	accumXvec[2] += limit - t2;

	accumXvec[4] +=  precomputed16_bitcount (combAND);	
}

int count01GAPCombs (unsigned int inputL, unsigned int inputR, unsigned int inputLVld, unsigned int inputRVld, unsigned int * accumXvec)
{
	
	unsigned int validComb = inputLVld & inputRVld;	
	unsigned int inputUL = inputL & validComb;
	unsigned int inputUR = inputR & validComb;	

	int validSampleSize = precomputed16_bitcount (validComb);
	
	int t1 = precomputed16_bitcount (inputUL); 
	int t2 = precomputed16_bitcount (inputUR); 

	unsigned int combAND = inputUL&inputUR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += validSampleSize - t1;
	accumXvec[2] += validSampleSize - t2;

	accumXvec[4] +=  precomputed16_bitcount (combAND);

	return validSampleSize;	
}

cor_t computePairwiseCorrelationBIN (alignment_struct * alignment, int s_i, int s_j)
{	
	int i, total;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val;

	total = alignment->sequences;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		count01Combs (total, s_i_val, s_j_val, accumXvec);
	
		total -= ENTRIES_PER_INT;
	}
		
	assert(accumXvec[0]+accumXvec[1]==alignment->sequences);
	assert(accumXvec[2]+accumXvec[3]==alignment->sequences);	

	return computeCorrelationValueBIN(alignment->sequences, accumXvec);	
}

float computePairwiseCorrelationBINGAPS (alignment_struct * alignment, int s_i, int s_j)
{
	int i, sequences=0;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		s_i_vld = alignment->compressedArrays[1][s_i*alignment->siteSize+i];
		s_j_vld = alignment->compressedArrays[1][s_j*alignment->siteSize+i];

		sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec);		
	}
	
	assert(accumXvec[0]+accumXvec[1]==sequences);
	assert(accumXvec[2]+accumXvec[3]==sequences);	

	return computeCorrelationValueBIN(sequences, accumXvec);	
}

float computePairwiseCorrelationDNA (alignment_struct * alignment, int s_i, int s_j)
{
	int i,j,k, total=0;
	unsigned int accumXvec[5],
                     validVec[8];

	cor_t pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val;
   	
	for(i=0;i<8;i++)
		validVec[i]=0;
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<5;k++)
				accumXvec[k]=0;

			total = alignment->sequences;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];

				count01Combs (total, s_i_val, s_j_val, accumXvec);

				total -= ENTRIES_PER_INT;
			}
			
			if(accumXvec[1]==0)
				break;

			validVec[i] = 1;			
			
			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(alignment->sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (alignment->sequences, pairwiseCorrelationMatrix, validVec);	
}

float computePairwiseCorrelationDNAGAPS (alignment_struct * alignment, int s_i, int s_j)
{
	
	int i,j,k,sequences;
	
	unsigned int    accumXvec[5],
                     validVec[8];

	float pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;        

        for(i=0;i<8;i++)
		validVec[i]=0;	
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			sequences=0;

			for(k=0;k<5;k++)
				accumXvec[k]=0;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];

				s_i_vld = alignment->compressedArrays[4][s_i*alignment->siteSize+k];
				s_j_vld = alignment->compressedArrays[4][s_j*alignment->siteSize+k];

				sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec);
			}
				
                        if(accumXvec[1]==0)
				break;

			validVec[i] = 1;
			
			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (sequences, pairwiseCorrelationMatrix, validVec);	
}
#else
void count01Combs (int total, unsigned int inputL, unsigned int inputR, unsigned int * accumXvec, char * bits_in_16bitsLocal)
{
	int limit=ENTRIES_PER_INT;
	
	if (total<limit)
		limit = total;

	int t1 = precomputed16_bitcountLocal (bits_in_16bitsLocal,inputL);
	int t2 = precomputed16_bitcountLocal (bits_in_16bitsLocal,inputR);

	unsigned int combAND = inputL & inputR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += limit - t1;
	accumXvec[2] += limit - t2;

	accumXvec[4] +=  precomputed16_bitcountLocal (bits_in_16bitsLocal,combAND);	
}

int count01GAPCombs (unsigned int inputL, unsigned int inputR, unsigned int inputLVld, unsigned int inputRVld, unsigned int * accumXvec, char * bits_in_16bitsLocal)
{	
	unsigned int validComb = inputLVld & inputRVld;	
	unsigned int inputUL = inputL & validComb;
	unsigned int inputUR = inputR & validComb;	

	int validSampleSize = precomputed16_bitcountLocal (bits_in_16bitsLocal,validComb);
	
	int t1 = precomputed16_bitcountLocal (bits_in_16bitsLocal,inputUL); 
	int t2 = precomputed16_bitcountLocal (bits_in_16bitsLocal,inputUR); 

	unsigned int combAND = inputUL&inputUR;

	accumXvec[1] += t1;  
	accumXvec[3] += t2;

	accumXvec[0] += validSampleSize - t1;
	accumXvec[2] += validSampleSize - t2;

	accumXvec[4] +=  precomputed16_bitcountLocal (bits_in_16bitsLocal,combAND);

	return validSampleSize;	
}

cor_t computePairwiseCorrelationBIN (alignment_struct * alignment, int s_i, int s_j, char * lookuptable)
{	
	int i, total;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val;

	total = alignment->sequences;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		count01Combs (total, s_i_val, s_j_val, accumXvec, lookuptable);
	
		total -= ENTRIES_PER_INT;
	}
		
	assert(accumXvec[0]+accumXvec[1]==alignment->sequences);
	assert(accumXvec[2]+accumXvec[3]==alignment->sequences);	

	return computeCorrelationValueBIN(alignment->sequences, accumXvec);	
}

float computePairwiseCorrelationBINGAPS (alignment_struct * alignment, int s_i, int s_j, char * lookuptable)
{
	int i, total, sequences=0;
	unsigned int accumXvec[5];
	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;

	total = alignment->sequences;

	for(i=0;i<5;i++)
		accumXvec[i]=0;
	
	for(i=0;i<alignment->siteSize;i++)
	{
		s_i_val = alignment->compressedArrays[0][s_i*alignment->siteSize+i];
		s_j_val = alignment->compressedArrays[0][s_j*alignment->siteSize+i];

		s_i_vld = alignment->compressedArrays[1][s_i*alignment->siteSize+i];
		s_j_vld = alignment->compressedArrays[1][s_j*alignment->siteSize+i];

		sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec, lookuptable);		
	}
	
	assert(accumXvec[0]+accumXvec[1]==sequences);
	assert(accumXvec[2]+accumXvec[3]==sequences);	

	return computeCorrelationValueBIN(sequences, accumXvec);	
}

float computePairwiseCorrelationDNA (alignment_struct * alignment, int s_i, int s_j, char * lookuptable)
{
	int i,j,k, total=0;
	unsigned int accumXvec[5],
                     validVec[8];

	cor_t pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val;
   	
	for(i=0;i<8;i++)
		validVec[i]=0;
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<5;k++)
				accumXvec[k]=0;

			total = alignment->sequences;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];

				count01Combs (total, s_i_val, s_j_val, accumXvec, lookuptable);

				total -= ENTRIES_PER_INT;
			}
			
			if(accumXvec[1]==0)
				break;

			validVec[i] = 1;			
			
			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(alignment->sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (alignment->sequences, pairwiseCorrelationMatrix, validVec);	
}

float computePairwiseCorrelationDNAGAPS (alignment_struct * alignment, int s_i, int s_j, char * lookuptable)
{	
	int i,j,k,sequences;
	
	unsigned int    accumXvec[5],
                     validVec[8];

	float pairwiseCorrelationMatrix [4][4]; // A->0 C->1 G->2 T->3

	unsigned int s_i_val, s_j_val, s_i_vld, s_j_vld;        

        for(i=0;i<8;i++)
		validVec[i]=0;	
	
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			sequences=0;

			for(k=0;k<5;k++)
				accumXvec[k]=0;

			for(k=0;k<alignment->siteSize;k++)
			{
				s_i_val = alignment->compressedArrays[i][s_i*alignment->siteSize+k];
				s_j_val = alignment->compressedArrays[j][s_j*alignment->siteSize+k];

				s_i_vld = alignment->compressedArrays[4][s_i*alignment->siteSize+k];
				s_j_vld = alignment->compressedArrays[4][s_j*alignment->siteSize+k];

				sequences += count01GAPCombs (s_i_val, s_j_val, s_i_vld, s_j_vld, accumXvec, lookuptable);
			}
				
                        if(accumXvec[1]==0)
				break;

			validVec[i] = 1;

			if(accumXvec[3]==0)
				continue;

			validVec[j+4] = 1;

			pairwiseCorrelationMatrix[i][j]=computeCorrelationValueBIN(sequences,accumXvec);
		}
	}
	return computeCorrelationValueDNA (sequences, pairwiseCorrelationMatrix, validVec);	
}
#endif

#ifdef _USE_PTHREADS
#ifdef _USE_PTHREADS_MEMINT
#ifdef _USE_PTHREADS_MULTI
void computeCorrelationsBIN_MULTI(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable, int ttid, int tthreads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j, k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%tthreads==ttid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;
			
#ifdef _SHARED
				correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j);
#else
				correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j, lookuptable);
#endif	
			}
			k++;
		}
	}	
}

void computeCorrelationsBINGAPS_MULTI(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable, int ttid, int tthreads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j, k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%tthreads==ttid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
				correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j);
#else
				correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j, lookuptable);
#endif			
			}
			k++;
		}
	}
}

void computeCorrelationsDNA_MULTI(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable, int ttid, int tthreads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j, k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%tthreads==ttid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
				correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j);
#else
				correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j, lookuptable);
#endif
			}
			k++;
		}
	}
}

void computeCorrelationsDNAGAPS_MULTI(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable, int ttid, int tthreads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j, k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%tthreads==ttid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
				correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j);
#else
				correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j, lookuptable);
#endif
			}
			k++;
		}
	}
}

void computeCorrelationsMULTI(alignment_struct * alignment, omega_struct * omega, int cvw_i, int firstRowToCompute, float ** myCorrelationMatrixN, char * bits_in_16bitsLocal, int ttid, int tthreads)
{

	if (firstRowToCompute==-1)
		return;

	switch(alignment->states)
	{	
		case 2: computeCorrelationsBIN_MULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, ttid, tthreads); 
			break;
		case 3: computeCorrelationsBINGAPS_MULTI(alignment,omega,cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, ttid, tthreads); 
			break;
		case 4: computeCorrelationsDNA_MULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, ttid, tthreads);  
			break;
		case 5: computeCorrelationsDNAGAPS_MULTI(alignment, omega, cvw_i, firstRowToCompute, myCorrelationMatrixN, bits_in_16bitsLocal, ttid, tthreads);  
			break;
		default: assert(0);
	}
	
}
#endif
void computeCorrelationsBIN(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;
			
#ifdef _SHARED
			correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j);
#else
			correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j, lookuptable);
#endif	
		}
	}	
}

void computeCorrelationsBINGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
			correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j);
#else
			correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j, lookuptable);
#endif
		}
	}
}

void computeCorrelationsDNA(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
			correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j);
#else
			correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j, lookuptable);
#endif
		}
	}
}

void computeCorrelationsDNAGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, cor_t ** correlationMatrix, char * lookuptable)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

#ifdef _SHARED
			correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j);
#else
			correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j, lookuptable);
#endif
		}
	}
}

void computeCorrelationMatrixPairwise(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable)
{

	if (firstRowIndex==-1)
		return;

	switch(alignment->states)
	{
		case 2: computeCorrelationsBIN(alignment,omega,omegaIndex, firstRowIndex, myCorrelationMatrix, lookuptable); 
			break;
		case 3: computeCorrelationsBINGAPS(alignment,omega,omegaIndex, firstRowIndex, myCorrelationMatrix, lookuptable); 
			break;
		case 4: computeCorrelationsDNA(alignment, omega, omegaIndex, firstRowIndex, myCorrelationMatrix, lookuptable); 
			break;
		case 5: computeCorrelationsDNAGAPS(alignment, omega, omegaIndex, firstRowIndex, myCorrelationMatrix, lookuptable); 
			break;
		default: assert(0);
	}	
}
#else
void computeCorrelationsBIN(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, int tid, int threads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i, j, s_i, s_j, k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%threads==tid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

				alignment->correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j);
			}
			k++;
		}
	}	
}

void computeCorrelationsBINGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, int tid, int threads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j,k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%threads==tid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

				alignment->correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j);
			}
			k++;
		}
	}
}

void computeCorrelationsDNA(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, int tid, int threads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j,k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%threads==tid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

				alignment->correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j);
			}
			k++;
		}
	}
}

void computeCorrelationsDNAGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, int tid, int threads)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j,k=0;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			if(k%threads==tid)
			{
				s_i = i + omega[omegaIndex].leftIndex;
				s_j = j + omega[omegaIndex].leftIndex;

				alignment->correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j);
			}
			k++;
		}
	}
}

void correlationThread(threadData_t * currentThread)
{
	alignment_struct * alignment = currentThread->threadArgPWC->alignment;

	omega_struct * omega = currentThread->threadArgPWC->omega;

	int omegaIndex = currentThread->threadArgPWC->omegaIndex;

	int firstRow=currentThread->threadArgPWC->firstRow;

	switch(alignment->states)
	{	
		case 2: computeCorrelationsBIN(alignment, omega, omegaIndex, firstRow,currentThread->threadID, currentThread->threadTOTAL); 
			break;
		case 3: computeCorrelationsBINGAPS(alignment,omega,omegaIndex, firstRow,currentThread->threadID, currentThread->threadTOTAL); 
			break;
		case 4: computeCorrelationsDNA(alignment, omega, omegaIndex, firstRow,currentThread->threadID, currentThread->threadTOTAL); 
			break;
		case 5: computeCorrelationsDNAGAPS(alignment, omega, omegaIndex, firstRow,currentThread->threadID, currentThread->threadTOTAL); 
			break;
		default: assert(0);
	}		
}

void setThreadArgumentsPWC(threadData_t * threadData, int tid, alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow)
{
	threadData[tid].threadArgPWC->alignment=alignment;
	threadData[tid].threadArgPWC->omega=omega;
	threadData[tid].threadArgPWC->omegaIndex=omegaIndex;
	threadData[tid].threadArgPWC->firstRow=firstRow;
}

void computeCorrelations_THREADS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow, threadData_t * threadData)
{
	int i, threads = threadData[0].threadTOTAL;

	for(i=0;i<threads;i++)
		setThreadArgumentsPWC(threadData, i, alignment, omega, omegaIndex, firstRow);
	
	startThreadOperations(threadData,PAIRWISECORRELATION);
}

void computeCorrelationMatrixPairwise(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable)
{

	if (firstRowIndex==-1)
		return;

	threadData_t * threadDataL = (threadData_t *)threadData;

	computeCorrelations_THREADS(alignment,omega,omegaIndex, firstRowIndex, threadDataL); 

}
#endif
#else
void computeCorrelationsBIN(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationBIN(alignment, s_i, s_j);
		}
	}	
}

void computeCorrelationsBINGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationBINGAPS(alignment, s_i, s_j);
		}
	}
}

void computeCorrelationsDNA(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationDNA(alignment, s_i, s_j);
		}
	}
}

void computeCorrelationsDNAGAPS(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRow)
{
	assert(omega[omegaIndex].rightIndex>=omega[omegaIndex].leftIndex);
	
	if (omega[omegaIndex].rightIndex==omega[omegaIndex].leftIndex)
		return;
	
	int i,j, s_i, s_j;

	int LinesToComputeTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex;

	for(i=firstRow;i<=LinesToComputeTotal;i++)
	{
		for(j=i-1;j>=0;j--)
		{
			s_i = i + omega[omegaIndex].leftIndex;
			s_j = j + omega[omegaIndex].leftIndex;

			alignment->correlationMatrix[i][j] = computePairwiseCorrelationDNAGAPS(alignment, s_i, s_j);
		}
	}
}

void computeCorrelationMatrixPairwise(alignment_struct * alignment, omega_struct * omega, int omegaIndex, int firstRowIndex, void * threadData, cor_t ** myCorrelationMatrix, char * lookuptable)
{

	if (firstRowIndex==-1)
		return;

	switch(alignment->states)
	{
		case 2: computeCorrelationsBIN(alignment,omega,omegaIndex, firstRowIndex); 
			break;
		case 3: computeCorrelationsBINGAPS(alignment,omega,omegaIndex, firstRowIndex); 
			break;
		case 4: computeCorrelationsDNA(alignment, omega, omegaIndex, firstRowIndex); 
			break;
		case 5: computeCorrelationsDNAGAPS(alignment, omega, omegaIndex, firstRowIndex); 
			break;
		default: assert(0);
	}	
}
#endif

cor_t ** createCorrelationMatrix(cor_t ** correlationMatrix, int matrixSize)
{
	int i;

	correlationMatrix = malloc (sizeof(cor_t *)*matrixSize);
	
	for(i=0;i<matrixSize;i++)
	{
		correlationMatrix[i] = malloc (sizeof(cor_t)*(i+1));

		correlationMatrix[i][i] = 0.0;
	}

	return correlationMatrix;
}

void applyCorrelationMatrixAdditions (omega_struct * omega, int omegaIndex, int firstRowIndex, cor_t ** correlationMatrix)
{
	int i,j;

	if (firstRowIndex==-1)
		return;
	
	int LinesToUpdateTotal = omega[omegaIndex].rightIndex - omega[omegaIndex].leftIndex; 

	for(i=firstRowIndex;i<=LinesToUpdateTotal;i++)
	{
		for(j=i-2;j>=0;j--)
		{
			correlationMatrix[i][j] =   correlationMatrix[i][j]  
						  + correlationMatrix[i-1][j]  
						  + correlationMatrix[i][j+1]  
						  - correlationMatrix[i-1][j+1];
		}
	}
}

void overlapCorrelationMatrixAdditions (alignment_struct * alignment, 
                                            omega_struct * omega, 
                                            int lvw_i, 
					    int cvw_i, 
                                            int * firstRowToCopy, 
                                            int * firstRowToCompute, 
                                            int * firstRowToAdd)
{

	*firstRowToCopy    = -1; // Do not shift any matrix values.
	*firstRowToCompute =  1; // Start recomputing correlations from the beginning of the matrix.
	*firstRowToAdd     =  2; // Start recomputing correlation distances from the beginning of the matrix.

	if (lvw_i==-1)
		return;

	int prevIndexA = omega[lvw_i].leftIndex,
	    prevIndexW = omega[lvw_i].rightIndex,
	    curIndexA  = omega[cvw_i].leftIndex,
	    curIndexW  = omega[cvw_i].rightIndex;

	assert(curIndexA>=prevIndexA);
	assert(curIndexW>=prevIndexW);

	int firstRowToCopyL = curIndexA - prevIndexA + 1,
	    prevLastRow     = prevIndexW - prevIndexA;

	if (firstRowToCopyL <= prevLastRow) // if matrices overlap
	{
		if(firstRowToCopyL==1)
			*firstRowToCopy = -1; // Do not shift any matrix values.	
		else
			*firstRowToCopy = firstRowToCopyL;	

		if(curIndexW==prevIndexW) // Same current matrix as the previous one.
		{
			*firstRowToCompute = -1; // Do not recompute anything.
			*firstRowToAdd = -1; // Do not recompute anything.
			return;					
		}
		
		*firstRowToCompute = prevIndexW - curIndexA + 1; // Compute only new rows.
		*firstRowToAdd = prevIndexW - curIndexA + 1; // Compute only new rows.
		return;	
	}
	return;		
}

void shiftCorrelationMatrixValues (omega_struct * omega, int lvw_i, int cvw_i, int firstRowToCopy, cor_t ** correlationMatrix)
{
	if (firstRowToCopy==-1 || lvw_i == -1)
		return;

	int i, j, c_i, c_j;

	int lastRowToCopy = omega[lvw_i].rightIndex - omega[lvw_i].leftIndex;

	int rowsToCopy = lastRowToCopy - firstRowToCopy + 1;

	c_i = firstRowToCopy;

	for(i=1;i<=rowsToCopy;i++)
	{
		c_j = c_i-1;

		for(j=i-1;j>=0;j--)
		{
			correlationMatrix[i][j] = correlationMatrix[c_i][c_j];
			c_j--;			
		}
		c_i++;
	}	
}
