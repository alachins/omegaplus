
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

int bitsPerEntry (alignment_struct * alignment)
{
	if(alignment->states==2) return 1;
	if(alignment->states==3) return 1;
	if(alignment->states==4) return 2;
	if(alignment->states==5) return 3;
	return 1;	
}

int getIntsPerSite(alignment_struct * alignment, int entriesPerInt)
{	
	int div = alignment->sequences / entriesPerInt;
	int mod = alignment->sequences % entriesPerInt;

	int ints = mod==0?div:div+1;

	return ints;
}

unsigned int charToIntMap(char input, int states)
{
	if(states==2)
	{
		if(input==ZERO)	return 0;
		if(input==ONE)	return 1;	
	}

	if(states==3)
	{
		if(input==ZERO)	return 0;
		if(input==ONE)	return 1;
		if(input==GAP)	return 2;		
	}

	if(states==4)
	{
		if(input==AD || input == ad)	return 0;
		if(input==CY || input == cy)	return 1;
		if(input==GU || input == gu)	return 2;
		if(input==TH || input == th)	return 3;
	}

	if(states==5)
	{
		if(input==AD || input == ad)	return 0;
		if(input==CY || input == cy)	return 1;
		if(input==GU || input == gu)	return 2;
		if(input==TH || input == th)	return 3;
		if(input==GAP)	return 4;
	}

	return 0;	
}

void initializeAlignmentCompression (alignment_struct * alignment)
{
	int i, 
            size= alignment->states,
	    div = alignment->sequences / REGISTER_WIDTH,
	    mod = alignment->sequences % REGISTER_WIDTH;
	
	alignment->siteSize = mod==0?div:div+1;

	alignment->compressedSize = alignment->siteSize * alignment->segsites; 

	if (alignment->states==2 || alignment->states==3)
		size--;		
		
	alignment->compressedArrays = malloc(sizeof(unsigned int *)*size);

	for(i=0;i<size;i++)
		alignment->compressedArrays[i] = malloc(sizeof(unsigned int)*alignment->compressedSize);	
}

void mapCharToCodeBIN (unsigned int * code, unsigned int * valid, char in)
{
	switch(in)
	{
		case ZERO: 
			*code  = 0u;
			*valid = 1u;
			break;		
		case ONE:
			*code  = 1u;
			*valid = 1u;
			break;
		case GAP:
			*code  = 2u;
			*valid = 0u;
			break;
		case UN:
			*code  = 2u;
			*valid = 0u;
			break;
		default:
			assert(0);
	}
}

void mapCharToCodeDNA (unsigned int * codeA, unsigned int * codeC, unsigned int * codeG, unsigned int * codeT, unsigned int * valid, char in)
{
  switch(in)
    {
      
    case ad:
    case AD: 
      *codeA = 1u;
      *codeC = 0u;
      *codeG = 0u;
      *codeT = 0u;
      *valid = 1u;
      break;
      
    case cy:
    case CY: 
      *codeA = 0u;
      *codeC = 1u;
      *codeG = 0u;
      *codeT = 0u;
      *valid = 1u;
      break;
      
    case gu:
    case GU: 
      *codeA = 0u;
      *codeC = 0u;
      *codeG = 1u;
      *codeT = 0u;
      *valid = 1u;
      break;
      
    case th:
    case TH: 
      *codeA = 0u;
      *codeC = 0u;
      *codeG = 0u;
      *codeT = 1u;
      *valid = 1u;
      break;		
      
    case GAP:
      *codeA = 0u;
      *codeC = 0u;
      *codeG = 0u;
      *codeT = 0u;
      *valid = 0u;
      break;
    case UN:
      *codeA = 0u;
      *codeC = 0u;
      *codeG = 0u;
      *codeT = 0u;
      *valid = 0u;
      break;
    default:
      assert(0);
    }
}


void compressAlignmentBIN(alignment_struct *alignment)
{
	int i,j,l,m,
            compLimit,
	    compLeft;
	    
	unsigned int tmpEntry, 
		     tmpValid, 
                     compEntry, 
                     compValid; 
	
	initializeAlignmentCompression (alignment);
	
	m=0;
	for(i=0;i<alignment->segsites;i++)
	{		
		l=0;

		compEntry = 0u;
		compValid = 0u;

		compLeft  = alignment->sequences;
		compLimit = REGISTER_WIDTH;

		if(compLeft<compLimit)
			compLimit = compLeft;

		for(j=0;j<alignment->sequences;j++)
		{

			mapCharToCodeBIN (&tmpEntry, &tmpValid, alignment->seqtable[j][i]);		
			
			compEntry = compEntry<<1|tmpEntry;
			compValid = compValid<<1|tmpValid;
		
			l++;
			if(l==compLimit)
			{	
				l=0;
			
				alignment->compressedArrays[0][m]=compEntry; // Vector of Ones
				
				if(alignment->states==3)
					alignment->compressedArrays[1][m]=compValid; // Valid Vector
				
				m++;

				compEntry = 0u;
				compValid = 0u;	
				
				compLeft -= REGISTER_WIDTH;
				compLimit = REGISTER_WIDTH;

				if(compLeft<compLimit)
					compLimit = compLeft;		
			}			
		}
	}
}

void compressAlignmentDNA(alignment_struct *alignment)
{
	int i,j,l,m,
            compLimit,
	    compLeft;
	    
	unsigned int tmpEntryA,
		     tmpEntryC,
		     tmpEntryG,
		     tmpEntryT, 
		     tmpValid, 
                     compEntryA,
		     compEntryC,
		     compEntryG,
		     compEntryT,
		     compValid; 
	
	initializeAlignmentCompression (alignment);
	
	m=0;
	for(i=0;i<alignment->segsites;i++)
	{		
		l=0;

		compEntryA = 0u;
		compEntryC = 0u;
		compEntryG = 0u;
		compEntryT = 0u;
		compValid  = 0u;

		compLeft  = alignment->sequences;
		compLimit = REGISTER_WIDTH;

		if(compLeft<compLimit)
			compLimit = compLeft;

		for(j=0;j<alignment->sequences;j++)
		{

			mapCharToCodeDNA (&tmpEntryA, &tmpEntryC, &tmpEntryG, &tmpEntryT, &tmpValid, alignment->seqtable[j][i]);		
			
			compEntryA = compEntryA<<1|tmpEntryA;
			compEntryC = compEntryC<<1|tmpEntryC;
			compEntryG = compEntryG<<1|tmpEntryG;
			compEntryT = compEntryT<<1|tmpEntryT;

			compValid = compValid<<1|tmpValid;
		
			l++;
			if(l==compLimit)
			{	
				l=0;
			
				alignment->compressedArrays[0][m]=compEntryA; // A Vector
				alignment->compressedArrays[1][m]=compEntryC; // C Vector
				alignment->compressedArrays[2][m]=compEntryG; // G Vector
				alignment->compressedArrays[3][m]=compEntryT; // T Vector
				
				if(alignment->states==5)
					alignment->compressedArrays[4][m]=compValid; // Valid Vector
				
				m++;

				compEntryA = 0u;
				compEntryC = 0u;
				compEntryG = 0u;
				compEntryT = 0u;
				compValid  = 0u;	
				
				compLeft -= REGISTER_WIDTH;
				compLimit = REGISTER_WIDTH;

				if(compLeft<compLimit)
					compLimit = compLeft;		
			}			
		}
	}
}

void compressAlignment(alignment_struct *alignment)
{
	switch(alignment->states)
	{
		case BINARY:
		case BINARY_WITH_GAPS:
			compressAlignmentBIN(alignment);
			break;
		case DNA:
		case DNA_WITH_GAPS:
			compressAlignmentDNA(alignment);
			break;		
		default:
		     	assert(0);	
	}
}
