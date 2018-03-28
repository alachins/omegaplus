
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

void printVersion (FILE * fp)
{
	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.0.0\n\n");
	fprintf(fp,"\tReleased:\t\tMarch 2012\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. Multi-grain parallelism (OmegaPlus-M)\n");
	fprintf(fp,"\t\t\t\t2. Approximated search for the best omega value per position (-b)\n");
	fprintf(fp,"\t\t\t\t3. Default minimum number of SNPs per window set to 5 (MINSNPPERWINDOW)\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.1.0\n\n");
	fprintf(fp,"\tReleased:\t\tJuly 2012\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. Supports VCF input file format\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.0\n\n");
	fprintf(fp,"\tReleased:\t\tOctober 2012\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. no-singletons command line flag to exclude the singletons\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.1\n\n");
	fprintf(fp,"\tReleased:\t\tOctober 2012\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. Minor bug fix in the ms parser\n");
	fprintf(fp,"\n\n");
	
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.2\n\n");
	fprintf(fp,"\tReleased:\t\tJanuary 2013\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. Fixed bug in the VCF parser associated with the handling of missing data\n");
	fprintf(fp,"\n\n");

		
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.4\n\n");
	fprintf(fp,"\tReleased:\t\tDecember 2013\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. add the -noSeparator flag to suppress printing of the // lines that  separate the datasets\n");
	fprintf(fp,"\n\n");

	
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.5\n\n");
	fprintf(fp,"\tReleased:\t\tJanuary 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. The parser of VCF can handle the '.' state in the ALT field\n");
	fprintf(fp,"\n\n");

	
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.7\n\n");
	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. Reads in small letters (FASTA) using the toupper function\n");
	fprintf(fp,"\t\t\t\t2. Reads in multiple frequences in the AF field (INFO)\n");
	fprintf(fp,"\n\n");


	
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.8\n\n");
	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"It can read ^M for the vcf files\n");
	fprintf(fp,"\n\n");

	
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.9\n\n");
	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Accepting . in INFO and PASS fields\n");
	fprintf(fp,"\n\n");


		
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.10\n\n");
	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Distance between grid points can be less than 1\n");
	fprintf(fp,"\n\n");


	
		
	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.2.11\n\n");
	fprintf(fp,"\tReleased:\t\tMay 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Skips empty MS alignments and alignments (MS, FASTA, VCF, MACS) that eventually have no SNPs \n");
	fprintf(fp,"\n\n");


	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t2.3.0\n\n");
	fprintf(fp,"\tReleased:\t\tNovember 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"1. It can process only a subset of the VCF samples in the input file\n");
	fprintf(fp,"\t\t\t\t2. Minor code modifications (unroll and jam)\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.0.0\n\n");
	fprintf(fp,"\tReleased:\t\tDecember 2014\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Faster parallel algorithm in OmegaPlus-G\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.0.1\n\n");
	fprintf(fp,"\tReleased:\t\tFebruary 2017\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Includes \"reports\" option\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.0.2\n\n");
	fprintf(fp,"\tReleased:\t\tMarch 2018\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Includes \"maf\" option\n");
	fprintf(fp,"\n\n");

	fprintf(fp, "\t=====================================\n");

	fprintf(fp,"\n\n");
	fprintf(fp,"\tVersion:\t\t3.0.3\n\n");
	fprintf(fp,"\tReleased:\t\tApril 2018\n\n");
	fprintf(fp,"\tComments:\t\t");
	fprintf(fp,"Includes \"mbs\" option\n");
	fprintf(fp,"\n\n");
	
}

void printHelp (FILE * fp)
{
	fprintf(fp,"\n\n\n");

	fprintf(fp," OmegaPlus | OmegaPlus-F | OmegaPlus-C | OmegaPlus-M\n");
	fprintf(fp,"\t -name runName\n");
	fprintf(fp,"\t -input inputFile\n");
	fprintf(fp,"\t -grid gridNumber\n");
	fprintf(fp,"\t -minwin mininumWindow\n");
	fprintf(fp,"\t -maxwin maximumWindow\n");
	fprintf(fp,"\t[-length alignmentLength]\n");
	fprintf(fp,"\t[-impute N|GAP]\n");
	fprintf(fp,"\t[-seed randomSeed]\n");
	fprintf(fp,"\t[-threads numberOfThreads]\n");	
	fprintf(fp,"\t[-binary]\n");
	fprintf(fp,"\t[-h|-help]\n");
	fprintf(fp,"\t[-all]\n");
	fprintf(fp,"\t[-minsnps minimumNumber]\n");
	fprintf(fp,"\t[-ld ldType]\n");
	fprintf(fp,"\t[-b maxSNPdiff]\n");
	fprintf(fp,"\t[-no-singletons]\n");
	fprintf(fp,"\t[-v|version]\n");
	fprintf(fp,"\t[-noSeparator]\n");
	fprintf(fp,"\t[-sampleList]\n");
	fprintf(fp,"\t[-sampleList_out]\n");
	fprintf(fp,"\t[-reports]\n");
	fprintf(fp,"\t[-maf threshold]\n");
	fprintf(fp,"\t[-mbs]\n");
	fprintf(fp,"\n\n");
	
	fprintf(fp,"\t-name <STRING>\t\tSpecifies a name for the run and the output files.\n\n");
	fprintf(fp,"\t-input <STRING>\t\tSpecifies the name of the input alignment file.\n");
	fprintf(fp,"\t      \t\t\tSupported file formats: MS-like and MaCS-like for binary data and FASTA for DNA data.\n\n");
	fprintf(fp,"\t-grid <INTEGER>\t\tSpecifies the number of omegas to be computed in the alignment. \n\n");
	fprintf(fp,"\t-minwin <INTEGER>\tSpecifies the minimum window to be used for computing linkage disequilibrium values between SNPs. \n\n");
	fprintf(fp,"\t-maxwin <INTEGER>\tSpecifies the maximum window to be used for computing linkage disequilibrium values between SNPs. \n\n");
	fprintf(fp,"\t-length <INTEGER>\tSpecifies the alignment length. Required only for MS-like and MaCS-like input files.\n\n");
	fprintf(fp,"\t-impute <N or GAP>\tEnables the random imputation of the following character (N or GAP) to valid \n");
	fprintf(fp,"\t       \t\t\talphabet characters. To enable the imputation of both N and GAP symbols use -impute twice.\n\n");
	fprintf(fp,"\t-seed <INTEGER>\t\tSpecifies a seed to initialize the random number generator used for the imputation of N \n");
	fprintf(fp,"\t     \t\t\tand GAP symbols as well as the deduction of DNA alignments to binary. Required with -impute\n");
	fprintf(fp,"\t     \t\t\tand -binary.\n\n");
	fprintf(fp,"\t-threads <INTEGER>\tSpecifies the number of threads. Required by OmegaPlus-F, OmegaPlus-C, and OmegaPlus-M.\n\n");
	fprintf(fp,"\t-binary\t\t\tConverts DNA alignments to binary format.\n\n");
	fprintf(fp,"\t-h|-help\t\tDisplays this help message.\n\n");
	fprintf(fp,"\t-all\t\t\tDisplays additional information in the result file.\n\n");
	fprintf(fp,"\t-minsnps <INTEGER>\tSpecifies the minimum number of SNPs per sub-region to calculate omega values.\n\n");
	fprintf(fp, "\t-ld <STRING>\t\tSpecifies the type of linkage disequilibrium measurement.\n");
	fprintf(fp,"\t      \t\t\tSupported types:\n");
	fprintf(fp,"\t      \t\t\tRSQUARE - the r^2 which is the correlation coefficient between a pairs of sites (default)\n");
	fprintf(fp,"\t      \t\t\tD - the D statistic\n");
	fprintf(fp,"\t      \t\t\tABSD - the absolute value of D\n");
	fprintf(fp,"\t      \t\t\tDOM - the D_omega statistic\n");
	fprintf(fp,"\t      \t\t\tABSDOM - the absolute value of D_omega\n");
	fprintf(fp,"\t      \t\t\tABSDOM2 - the absolute value of D_omega normalized by the product of frequencies\n\n");
	fprintf(fp,"\t-b <INTEGER>\t\tEnables approximate search for the best omega value per position:\n");
	fprintf(fp,"\t      \t\t\tthe number of SNPs on the left and right windows does not differ more than <INTEGER>.\n\n");
	fprintf(fp,"\t-no-singletons\t\tExcludes the singletons from the analysis.\n\n");
	fprintf(fp,"\t-v|-version\t\tDisplays version information.\n\n");
	fprintf(fp,"\t-noSeparator\t\tTo suppress printing the // flag that separates datasets. Useful for meta-processing (particularly with R).\n\n");
	fprintf(fp,"\t-sampleList <STRING>\tTo be used with VCF files in order to specify which samples to be included in the analysis.\n\n");
	fprintf(fp,"\t-sampleList_out <STRING>\tTo generate a list of VCF samples in the input VCF file.\n\n");
	fprintf(fp,"\t-reports\t\tTo generate each alignment report in a separate file.\n\n");
	fprintf(fp,"\t-maf <FLOAT>\t\tTo exclude SNPs with minor allele frequency < threshold.\n\n");
	fprintf(fp,"\t-mbs\t\t\tTo specify that the input file is in mbs format (will be processed as if it was ms).\n\n");
	fprintf(fp,"\n\n");
}

int flagMatch(FILE *fp, char flag[], int flaglength, char tmp)
{
	int counter = 0;

	while(counter < flaglength)
	{
		if(tmp != flag[counter])
			break;

		tmp = fgetc(fp);

		++counter;
	}

	return counter;
}

int getFileFormat (FILE * fp)
{
  
  char tmp;

  tmp=fgetc(fp);

  char macsFLAG[] = "COMMAND";

  int flaglength = 7,j;

  char vcfFLAG[] = "##fileformat=VCF";

  int flaglengthVCF = 16;

  int format = OTHER_FORMAT;

  while(tmp!=EOF)
    {
      if(tmp=='/')
	{
	  tmp = fgetc(fp);
	  
	  if(tmp=='/')
	    {
	      format = MS_FORMAT;
	      break;
	    }
	  else
	    tmp = fgetc(fp);				
	}
      else
	{
	  if(tmp=='>')
	    {
	      format = FASTA_FORMAT;
	      break;
	    }
	  else
	    {
	      j = flagMatch(fp, macsFLAG, flaglength, tmp);

	      			    
	      if(j == flaglength)
		{
		  format = MACS_FORMAT;
		  break;
		}
		
	      else
		{

		  fseek(fp, -j - 1, SEEK_CUR);

		  tmp = fgetc( fp);
		  
		  j = flagMatch(fp, vcfFLAG, flaglengthVCF, tmp);

		  if(j == flaglengthVCF)
		    {
		      format = VCF_FORMAT;
		      break;
		    }
		  
		  else
		    tmp = fgetc(fp);
		}
	    }
	}		
    }

  //fprintf(stderr, "\nThe format of the input file is: %d (0: ms, 1: fasta, 2: macs, 3: vcf)\n", format);
  
  return format;
}

void commandLineParser(int argc, char** argv, 
		       char * infile, 
		       int * grid, 
                       int * length, 
                       int * minw, 
                       int * maxw, 
                       char ** recfile, 
                       int * minsnps, 
                       int * imputeN, 
                       int * imputeG, 
                       int * binary, 
                       unsigned int * seed,
		       int * fileFormat,
		       int * threads,
		       int * results,
		       int * ld,
		       int * borderTol,
		       int * filterOut,
		       int *noSeparator,
		       char * samplefile_i,
		       int * generateVCFsamplelist,
		       int *memLimit,
		       int * reports,
		       double * maf,
		       int * fileFormatMBS)
{
  int i, nameSet = 0, fileSet=0, gridSet=0, lengthSet=0, minwSet=0, maxwSet=0, seedSet=0, imputeSet=0, binarySet=0, seedCheck=0;
	char impute, model[100];
	FILE *fp;

#ifdef _USE_PTHREADS
	int threadsSet = 0;
#endif

#ifdef _USE_OPENMP_GENERIC
	int memLimitSet=0;
	int threadsSet = 0;
#endif
	
	strcpy(runName,"x");
	*results = RESULTS_DEFAULT;
	*ld = RSQUARE;

	for(i=1; i<argc; ++i)
	{
	  if(!strcmp(argv[i], "-noSeparator") || !strcmp(argv[i], "-noseparator") )
	    {
	      *noSeparator = 1;
	      continue;
	    }

		if(!strcmp(argv[i], "-name")) 
		{ 
			if (i!=argc-1)
			{
				strcpy(runName,argv[++i]);
				nameSet = 1;
			}
			continue;
		}

    		if(!strcmp(argv[i], "-input")) 
		{ 
			if (i!=argc-1)
			{			
				strcpy(infile,argv[++i]);

				fp=fopen(infile,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",infile);
					exit(0);
				}
				else
				{
					*fileFormat = getFileFormat (fp);

				
					if(*fileFormat==FASTA_FORMAT || *fileFormat==VCF_FORMAT) 
						lengthSet=1;

					if(*fileFormat==VCF_FORMAT)
						seedCheck=1;
								
					fclose(fp);

					fileSet=1;
				}
			}
			continue;
		}

    		if(!strcmp(argv[i], "-grid"))
		{
			if (i!=argc-1)
			{
				*grid = atoi(argv[++i]);

				if(*grid!=0)
					gridSet=1;
			} 
			continue;
		}

		if(!strcmp(argv[i], "-seed"))
		{
			if (i!=argc-1)
			{
				*seed = atoi(argv[++i]);

				if(*seed!=0)
					seedSet=1;
			} 
			continue;
		}

		if(!strcmp(argv[i], "-impute"))
		{
			if (i!=argc-1)
			{
				impute = argv[++i][0];
		
				switch(impute)
				{
					case 'N':
						*imputeN = 1;
						imputeSet=1;
						break;
					case 'G':
						*imputeG = 1;
						imputeSet=1;
						break;
					default:
						assert(0); 
				}
			}
			continue;
		}

		if(!strcmp(argv[i], "-binary"))
		{ 
			*binary = 1;
			binarySet=1;
			continue;
		}

		if(!strcmp(argv[i], "-no-singletons"))
		{ 
			*filterOut = 1;
			continue;
		}

		if(!strcmp(argv[i], "-all"))
		{ 
			*results = RESULTS_ALL; 
			continue;
		}

    		if(!strcmp(argv[i], "-length"))
		{ 
			if (i!=argc-1)
			{
				*length = atoi(argv[++i]);
			
				if(*length!=0)
					lengthSet=1;
			} 
			continue;
		}

		if(!strcmp(argv[i], "-ld"))
		{ 
			if (i!=argc-1)
			{
				strcpy(model,argv[++i]);
			
				if(!strcmp(model, "RSQUARE"))
					*ld = RSQUARE;
				if(!strcmp(model, "DOM"))
					*ld = DOM;
				if(!strcmp(model, "ABSDOM"))
					*ld = ABSDOM;
				if(!strcmp(model, "D"))
					*ld = JUSTD;
				if(!strcmp(model, "ABSD"))
					*ld = ABSD;
				if(!strcmp(model, "ABSDOM2"))
					*ld = ABSDOM2;
				
			} 
			continue;
		}

    		if(!strcmp(argv[i], "-help")||!strcmp(argv[i], "-h"))
		{ 
			printHeading (stdout);

			printHelp (stdout);

			exit(0);
		}

		if(!strcmp(argv[i], "-version")||!strcmp(argv[i], "-v"))
		{ 
			//printHeading (stdout);

			printVersion (stdout);

			exit(0);
		}

    		if(!strcmp(argv[i], "-minwin"))
		{ 
			if (i!=argc-1)
			{
				*minw = atoi(argv[++i]);

				if(*minw!=0)
					minwSet=1;
 			}
			continue; 
		}

    		if(!strcmp(argv[i], "-maxwin"))
		{ 
			if (i!=argc-1)
			{
				*maxw = atoi(argv[++i]);
			
				if(*maxw!=0)
					maxwSet=1;
			}
			continue; 
		}

		if(!strcmp(argv[i], "-minsnps"))
		{ 
			if (i!=argc-1)
			{
				*minsnps = atoi(argv[++i]); 
			}			
			continue; 
		}
		
		if(!strcmp(argv[i], "-b"))
		{
			if(i!=argc-1)
			{
				*borderTol = atoi(argv[++i]);
		      	}
		    	continue;
		}

		/*if(!strcmp(argv[i], "-recmap"))
		{ 
			if (i!=argc-1)
			{
				*recfile = malloc(INFILENAMESIZE*sizeof(char));
      				strcpy(*recfile, argv[++i]);
			}
      			continue; 
    		}*/

		if(!strcmp(argv[i], "-sampleList"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList\n\n");
					exit(0);
				}
		
				strcpy(samplefile_i,argv[++i]);

				fp=fopen(samplefile_i,"r");

				if (fp==NULL)
				{
					fprintf(stderr, "\n ERROR: File %s does not exist.\n\n",samplefile_i);
					exit(0);
				}
				else
				{
					fclose(fp);
				}
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList\n\n");
				exit(0);
			}	

			continue;
		}


		if(!strcmp(argv[i], "-sampleList_out"))
		{ 
			if (i!=argc-1)
			{	
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList_out\n\n");
					exit(0);
			}
		
				strcpy(samplefile_i,argv[++i]);

				*generateVCFsamplelist=1;

				// add check for existing file and ask to overwrite it
			}
			else
			{
				fprintf(stderr, "\n ERROR: File name (string) is missing after argument -sampleList_out\n\n");
				exit(0);
			}	

			continue;
		}

		if(!strcmp(argv[i], "-reports"))
		{
			*reports = 1;
			continue;
		}

		if(!strcmp(argv[i], "-maf"))
		{
		    if (i != argc-1)
			{
				
				if(argv[i+1][0]=='-')
				{
					fprintf(stderr, "\n ERROR: Value is missing after argument -maf\n\n");
					exit(0);
				}
				
				*maf = atof(argv[++i]);

				if(*maf<0.0 || *maf>1.0)
				{
					fprintf(stderr, "\n ERROR: Invalid MAF value (valid: 0.0-1.0)\n\n");
					exit(0);
				}
			} 
			else
			{
				fprintf(stderr, "\n ERROR: Value is missing after argument -maf\n\n");
				exit(0);	
			}
			continue;
		}
	
		if(!strcmp(argv[i], "-mbs"))
		{
			*fileFormatMBS = 1;
			continue;
		}



#ifdef _USE_PTHREADS
		if(!strcmp(argv[i], "-threads"))
		{
			if (i!=argc-1)
			{
				*threads = atoi(argv[++i]);

				if(*threads!=0)
					threadsSet=1;						
			} 
			continue;
		}
#endif

#ifdef _USE_OPENMP_GENERIC

		if(!strcmp(argv[i], "-memLimit"))
		{
			if (i!=argc-1)
			{
				*memLimit = atoi(argv[++i]);

				memLimitSet=1;
			} 
			continue;
		}

		if(!strcmp(argv[i], "-threads"))
		{
			if (i!=argc-1)
			{
				*threads = atoi(argv[++i]);

				if(*threads!=0)
					threadsSet=1;						
			} 
			continue;
		}

#endif

		fprintf(stderr, "\n ERROR: %s is not a valid command line parameter\n\n",argv[i]);
		exit(0);
	}

	if (nameSet==0)
	{
		fprintf(stderr, "\n ERROR: Please specify a name for this run with -name\n\n");
		exit(0);
	}

	if (fileSet==0)
	{
		fprintf(stderr, "\n ERROR: Please specify an alignment with -input\n\n");
		exit(0);
	}

	if (gridSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify the number of omegas with -grid\n\n");
		exit(0);
	}

	if (lengthSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify the alignment length with -length\n\n");
		exit(0);
	}

	if (minwSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify a minimum window size with -minwin\n\n");
		exit(0);
	}

	if (maxwSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify a maximum window size with -maxwin\n\n");
		exit(0);
	}
	
	if ((imputeSet==1 || binarySet==1 || seedCheck==1) && seedSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify a random seed with -seed\n\n");
		exit(0);
	}

	if(*maxw<*minw)
	{
		fprintf(stderr, "\n ERROR: Maximum window size must be greater than minimum window size\n\n");
		exit(0);
	}

#ifdef _USE_PTHREADS
	if (threadsSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify the number of threads with -threads\n\n");
		exit(0);
	}
	else
	{
		if(*threads<=1 && *generateVCFsamplelist==0)
		{
			fprintf(stderr, "\n ERROR: Please specify a number of threads greater than 1\n\n");
			exit(0);
		}
	}
#endif

#ifdef _USE_OPENMP_GENERIC

	if (memLimitSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify a memory limit with -memLimit\n\n");
		exit(0);
	}
	else
	{
		if(*memLimit<=1 && *generateVCFsamplelist==0)
		{
			fprintf(stderr, "\n ERROR: Please specify a memory limit larger than 1 MB\n\n");
			exit(0);
		}
	}

	if (threadsSet==0 && *generateVCFsamplelist==0)
	{
		fprintf(stderr, "\n ERROR: Please specify the number of threads with -threads\n\n");
		exit(0);
	}
	else
	{
		if(*threads<=1 && *generateVCFsamplelist==0)
		{
			fprintf(stderr, "\n ERROR: Please specify a number of threads greater than 1\n\n");
			exit(0);
		}
	}

#endif


}

int isSpace(char ent)
{
	if(ent == 9 || ent == 32) // 9 -> '\t' , 32 -> ' '
		return 1;

	return 0;
}


void ignoreLineSpaces(FILE *fp, char *ent)
{
	while(*ent==' '|| *ent == 9) // horizontal tab
		*ent = fgetc(fp);  
}

int isEndOfLine(char ent)
{
	if(ent == 10 || ent == 13 || ent == 32)
		return 1;

	return 0;
}

int isEndOfLine2(char ent)
{
	if(ent == 10 || ent == 13)
		return 1;

	return 0;
}



void ignoreNewLineChars(FILE *fp, char *ent)
{

  int brk = 0;
  
  while( *ent == 10 || *ent == 13)
    {
      *ent = fgetc(fp);
      
      if( *ent != 10 && *ent != 13)
	{
	  ungetc(*ent, fp);
	  
	  brk = 1;
	  
	  break;
	}
    }
  if(brk == 1)
    *ent = 10;

}


int getNextString(FILE *fp, char ** word, int *readEOL, int *readEOF, int *wordLength)
{

	*readEOL = *readEOF = 0;

	char ent = fgetc(fp);

	int i=0;

	ignoreLineSpaces(fp, &ent);

	if(ent == EOF)
	{
		*readEOF = 1;
		(*word)[0] = '\0';
		return 0;
	}

	if(isEndOfLine2(ent))
	{
		*readEOL = 1;
		(*word)[0] = '\0';
		return 0;
	}

	while(isSpace(ent)==0 && isEndOfLine2(ent)==0 && ent != EOF)
	{

		if(i+1 >= (*wordLength) )
		{
			(*wordLength) = (*wordLength) << 1; 

			(*word) = realloc( (*word), (*wordLength) * sizeof(char) );		
		}

		(*word)[i++] = ent;

		ent = fgetc(fp);
	}

	(*word)[i] = '\0';

	if(isEndOfLine2(ent))
		*readEOL = 1;

	return 1;
}

int getNextString_all_lines(FILE *fp, char ** word, int *readEOL, int *readEOF, int *wordLength)
{
  
  *readEOL = *readEOF = 0;

  char ent = fgetc(fp);

  int i=0;

  while(isSpace(ent) || isEndOfLine2(ent))
  {
	ent = fgetc(fp);
  }


 if(ent == EOF)
    {
      *readEOF = 1;
	    
      (*word)[0] = '\0';
	    
      return 0;
    }


  while(isSpace(ent)==0 && isEndOfLine2(ent)==0 && ent != EOF)
    {

      if( i+1 >= (*wordLength) )
	{
	  (*wordLength) = (*wordLength) << 1; 

	  (*word) = realloc( (*word), (*wordLength) * sizeof(char) );
		    
	}
	
      (*word)[i++] = ent;
      
      ent = fgetc(fp);
    }
  
  (*word)[i] = '\0';
  
  if(isEndOfLine2(ent))
    *readEOL = 1;
  	
  return 1;
}

int findFirstAlignment(alignment_struct * alignment, FILE *fp, FILE *fpInfo, int format, FILE *fpVCFsamples, int generateVCFsamplelist, char * vcf_samples_filename)
{
	char tmp;
	int temp, i;

	if(format==FASTA_FORMAT)
	{
		tmp=fgetc(fp);

		while(tmp!=EOF)
		{
			if(tmp=='>') // first sequence
				return 1;
			else
				tmp = fgetc(fp);
		}
	}
	else if(format == MS_FORMAT)
	{
	    	tmp=fgetc(fp);

		while(tmp!=EOF)
		{
			if(tmp=='/')
			{
				tmp = fgetc(fp);

				if(tmp=='/')
					return 1;
				else
					tmp = fgetc(fp);				
			}
			else
				tmp = fgetc(fp);
		}
	}
	else if(format == MACS_FORMAT)
	{
		char word[1000];
		int sitenumber;

		while( fscanf(fp, "%s", word) )
		{
			if(!strcmp(word, "SITE:") )
			{
				temp = fscanf(fp, "%d", &sitenumber);
				temp = temp;
				if(sitenumber == 0)
					return 1;
			}
		}
	}
	else if(format == VCF_FORMAT)
	{

		char ** string = (char **) malloc (sizeof(char*));
		
		(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);
		
		int status=-1, eol=0, eof=0, maxLength=STRINGLENGTH;

		VCF_header_lines = 0;

		while(1)
		{
			status = getNextString (fp, string, &eol, &eof, &maxLength);

			if(eol==1)
				VCF_header_lines++;
			
			if(status==1)
				if(strcmp((*string),"#CHROM")==0)
				{
					VCF_header_lines++;
					break;
				}

			if(eof==1)
				assert(0);		
		}

		

		
		char * headerFields[VCF_HLENGTH];
		int fieldInd=1;
		int VCFsamples=0;
		
		headerFields[0] = "#CHROM";
		headerFields[1] = "POS";
		headerFields[2] = "ID";
		headerFields[3] = "REF";
		headerFields[4] = "ALT";
		headerFields[5] = "QUAL";
		headerFields[6] = "FILTER";
		headerFields[7] = "INFO";
		headerFields[8] = "FORMAT";

		
		
		while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
		{

			if(fieldInd<VCF_HLENGTH)
				if(strcmp(headerFields[fieldInd],(*string))!=0)
				{
					fprintf(stderr, "\n\n ERROR: VCF header field %s is missing.\n\n",headerFields[fieldInd]);
					assert(0);
				}
			
			fieldInd++;

			if(fieldInd>=VCF_HLENGTH)	
				break;

			if(eol==1 || eof==1)
				assert(0);			

		}

		int sampleList_size=1, sampleList_index=0;
		char ** sampleList = (char **) malloc (sizeof(char *)*sampleList_size);

		while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
		{

			if(sampleList_index==sampleList_size)
			{
				sampleList_size++;
				sampleList = realloc(sampleList, sampleList_size * sizeof(char*));
			}

			sampleList[sampleList_index] = malloc(sizeof(char)*STRINGLENGTH);
			strcpy(sampleList[sampleList_index++], *string);

			VCFsamples++;

			if(eol==1)
				break;

			if(eof==1)
				assert(0);			
		
		}

		assert(VCFsamples!=0);

		if(generateVCFsamplelist==1)
		{
			for(i=0;i<VCFsamples;i++)
			{
				fprintf(fpVCFsamples,"%s\n",sampleList[i]);
			}

			fclose(fpVCFsamples);

			fprintf(stdout, "\n A list of VCF samples has been stored in file %s!\n\n",vcf_samples_filename);
			fprintf(fpInfo, "\n A list of VCF samples has been stored in file %s!\n\n",vcf_samples_filename);

			return 0;
		}

		fprintf(stdout, " Total number of samples in the VCF:\t\t\t\t%d\n",VCFsamples);
		fprintf(fpInfo, " Total number of samples in the VCF:\t\t\t\t%d\n",VCFsamples);

		alignment->VCFsamples = VCFsamples;
		alignment->VCFsample_valid = malloc(sizeof(int)*VCFsamples);

		for(i=0;i<VCFsamples;i++)
			alignment->VCFsample_valid[i]=0; // set all samples to invalid


		if(fpVCFsamples==NULL)
		{
			for(i=0;i<VCFsamples;i++)
				alignment->VCFsample_valid[i]=1; // set all samples to valid


			fprintf(stdout, " Samples excluded from the analysis:\t\t\t\t0\n\n\n");
			fprintf(fpInfo, " Samples excluded from the analysis:\t\t\t\t0\n\n\n");

		}
		else
		{
			int validVCFsamples_matched=0;

			while(getNextString_all_lines (fpVCFsamples, string, &eol, &eof, &maxLength)==1)
			{
				//printf("%s\n",*string);

				for(i=0;i<VCFsamples;i++)
				{
					if(!strcmp(*string, sampleList[i]))
					{
						//printf("Match at %d %s %s\n",i, *string, sampleList[i]);
						alignment->VCFsample_valid[i]=1;
						validVCFsamples_matched++;
						break;
					}
				}

				if(i==VCFsamples)
				{
					fprintf(stdout, " *** NOT FOUND: Sample \"%s\" does not exist in the VCF file ***\n",*string);
					fprintf(fpInfo, " *** NOT FOUND: Sample \"%s\" does not exist in the VCF file ***\n",*string);

				}
			
			}

				fprintf(stdout, " Samples excluded from the analysis:\t\t\t\t%d\n\n\n",VCFsamples-validVCFsamples_matched);
				fprintf(fpInfo, " Samples excluded from the analysis:\t\t\t\t%d\n\n\n",VCFsamples-validVCFsamples_matched);
		}
	
		assert(VCFsamples!=0);

		
		tmp = fgetc(fp);
		
		ungetc(tmp, fp);
		
		ignoreNewLineChars(fp, &tmp);


		getNextString (fp, string, &eol, &eof, &maxLength);

		strncpy(VCF_alignment_name, *string, MAX_CHROM_NAME_VCF);

		assert(strlen(VCF_alignment_name)!=0);

		free(*string);
		free(string);

		for(i=0;i<VCFsamples;i++)
			free(sampleList[i]);

		free(sampleList);
	
		return 1;
	}

	return 0;		
}

int findNextAlignment(FILE *fp, int fileFormat)
{
	char stop,tmp;
	int temp;

	if (fileFormat == MS_FORMAT)
		stop = '/';
	else
		stop = '>';
	
	
	if(fileFormat == MS_FORMAT || fileFormat == FASTA_FORMAT)
	{
		tmp=fgetc(fp);

		while(tmp!=EOF)
		{
			if(tmp==stop)
				return 1;	
			else
	  			tmp = fgetc(fp);
		}	    	    
	}
	
	if(fileFormat == MACS_FORMAT)
	{
		char word[100];
		int sitenumber;

		int nextAl=1;

		while( (nextAl = fscanf(fp, "%s", word)))
		{

			if( nextAl <= 0)
				return 0;

			if(!strcmp(word, "SITE:") )
			{
				temp = fscanf(fp, "%d", &sitenumber);
				temp = temp;

				if(sitenumber == 0)
					return 1;
			}
		}
	}

	if(fileFormat == VCF_FORMAT)
	{
		if(nxtVCFalignment==0)
			return 1;
	}

	return 0;
}

int mapCharToInt(char a)
{
	if(a == ZERO) 	return 0;
	if(a == ONE) 	return 1;
	if(a == GAP || a == UN) return 2;
	if(a == AD) 	return 3;
	if(a == CY) 	return 4;
	if(a == GU) 	return 5;
	if(a == TH) 	return 6;
	return 7;
}



void removeNonPolymorphicSitesBIN(alignment_struct *alignment, FILE * fp, int filterOut, double maf)
{	
	int i,j, dif, rm;
	int setVec[STATESALL];

	int k = 0;

	int * polpositions = malloc(alignment->segsites * sizeof(int));
	
	if(filterOut==0)
	{
		for(i=0;i<alignment->segsites;i++)
		{
			for(j=0;j<STATESALL;j++)
				setVec[j]=0;		
		
			rm=1;

			for(j=0;j<alignment->sequences;j++)
			{
				setVec[mapCharToInt(alignment->seqtable[j][i])] += 1;

				dif  = setVec[0]!=0?1:0;
				dif += setVec[1]!=0?1:0;		

				if(dif==2)
				{
					rm=0;		
					break;
				}
			}
			
			if(rm == 0)
			{
				polpositions[k] = i;
				k++;
			}
		}
	}
	else
	{
		for(i=0;i<alignment->segsites;i++)
		{
			for(j=0;j<STATESALL;j++)
				setVec[j]=0;		
		
			rm=1;

			for(j=0;j<alignment->sequences;j++)
				setVec[mapCharToInt(alignment->seqtable[j][i])] += 1;
	
			dif  = setVec[0]!=0?1:0;
			dif += setVec[1]!=0?1:0;		

			if(dif==2)
			{
				rm=0;

				if(setVec[0] == 1 || setVec[1] == 1)
					rm = 1;
			}

			if(rm == 0)
			{
				polpositions[k] = i;
				k++;
			}
		}		
	}


	double maf_thres = maf;
	if(maf_thres>=0.0 && maf_thres<=1.0)
	{
		for(j=0;j<k;j++)
		{
			int a1c = 0;
			int a2c = 0;
			int aTc = 0;
	
			for(i=0;i<alignment->sequences;i++)
			{
				a1c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 1? 1:0;
				a2c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 0? 1:0;
				aTc += mapCharToInt(alignment->seqtable[i][polpositions[j]]) != 2? 1:0;
			}

			assert(aTc<=alignment->sequences);
			assert(a1c<aTc);
			assert(a2c<aTc);

			double a1f = ((double)a1c)/((double)aTc);
			double a2f = ((double)a2c)/((double)aTc);
			if(a1f<maf_thres || a2f<maf_thres)
			{
				polpositions[j]=-1;
			}
		}

		int k_new = 0;
		for(j=0;j<k;j++)
		{
			if(polpositions[j]!=-1)
				polpositions[k_new++] = polpositions[j];
		}
		assert(k_new<=k);
		k=k_new;
	}
	

	if(k==alignment->segsites)
	{
		free(polpositions);
		return;
	}

	alignment->poltable = malloc(alignment->sequences * sizeof(char*));

	for(i=0; i<alignment->sequences; ++i)
		alignment->poltable[i] = malloc(k * sizeof(char));
	
	for(j=0; j<k; ++j)
		alignment->positionsInd[j] = alignment->positionsInd[polpositions[j]];
	
	for(i=0; i<alignment->sequences; ++i)
	{
		for(j=0; j<k; ++j)
			alignment->poltable[i][j] = alignment->seqtable[i][polpositions[j]];
	
		free(alignment->seqtable[i]);

		alignment->seqtable[i] = alignment->poltable[i];

		alignment->poltable[i] = NULL;
	}

	free(alignment->poltable);

	fprintf(stdout,"\n\t\tDiscarded sites:\t%d",alignment->segsites-k);
	
	fprintf(fp,"\n\t\tDiscarded sites:\t%d",alignment->segsites-k);

	alignment->segsites = k;
	
	free(polpositions);	
}

void removeNonPolymorphicSitesDNA(alignment_struct *alignment, FILE * fp, int filterOut, double maf)
{	
	int i,j, dif, rm;
	int setVec[STATESALL];

	int * polpositions = malloc(alignment->segsites * sizeof(int));
	
	int k = 0;

	if(filterOut==0)
	{
		for(i=0;i<alignment->segsites;i++)
		{
			for(j=0;j<STATESALL;j++)
				setVec[j]=0;
		
			rm=1;

			for(j=0;j<alignment->sequences;j++)
			{
		
				setVec[mapCharToInt(alignment->seqtable[j][i])] += 1;

				dif  = setVec[3]!=0?1:0;
				dif += setVec[4]!=0?1:0;
				dif += setVec[5]!=0?1:0;
				dif += setVec[6]!=0?1:0;

				if(dif>=2)
				{
					rm=0;
					break;
				}
			}		
		
			if(rm == 0)
			{
				polpositions[k] = i;
				k++;
			}
		}
	}
	else
	{
		for(i=0;i<alignment->segsites;i++)
		{
			for(j=0;j<STATESALL;j++)
				setVec[j]=0;
		
			rm=1;

			for(j=0;j<alignment->sequences;j++)		
				setVec[mapCharToInt(alignment->seqtable[j][i])] += 1;

			dif  = setVec[3]!=0?1:0;
			dif += setVec[4]!=0?1:0;
			dif += setVec[5]!=0?1:0;
			dif += setVec[6]!=0?1:0;

			if(dif>=2)
			{
				rm=0;

				if(setVec[3] == 1 || setVec[4] == 1 || setVec[5] == 1 || setVec[6] == 1)
					rm = 1;				
			}
		
			if(rm == 0)
			{
				polpositions[k] = i;
				k++;
			}
		}
	}

	double maf_thres = maf;
	if(maf_thres>=0.0 && maf_thres<=1.0)
	{
		for(j=0;j<k;j++)
		{
			int a1c = 0;
			int a2c = 0;
			int a3c = 0;
			int a4c = 0;
			int aTc = 0;
	
			for(i=0;i<alignment->sequences;i++)
			{
				a1c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 3? 1:0;
				a2c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 4? 1:0;
				a3c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 5? 1:0;
				a4c += mapCharToInt(alignment->seqtable[i][polpositions[j]]) == 6? 1:0;
				aTc += mapCharToInt(alignment->seqtable[i][polpositions[j]]) != 2? 1:0;
			}

			assert(aTc<=alignment->sequences);
			assert(a1c<aTc);
			assert(a2c<aTc);
			assert(a3c<aTc);
			assert(a4c<aTc);

			double a1f = ((double)a1c)/((double)aTc);
			double a2f = ((double)a2c)/((double)aTc);
			double a3f = ((double)a3c)/((double)aTc);
			double a4f = ((double)a4c)/((double)aTc);
			if((a1c!=0 && a1f<maf_thres) || (a2c!=0 && a2f<maf_thres) || (a3c!=0 && a3f<maf_thres) || (a4c!=0 && a4f<maf_thres))
			{
				polpositions[j]=-1;
			}
		}

		int k_new = 0;
		for(j=0;j<k;j++)
		{
			if(polpositions[j]!=-1)
				polpositions[k_new++] = polpositions[j];
		}
		assert(k_new<=k);
		k=k_new;
	}

	if(k==alignment->segsites)
	{
		free(polpositions);
		return;
	}

	alignment->poltable = malloc(alignment->sequences * sizeof(char*));

	for(i=0; i<alignment->sequences; ++i)
		alignment->poltable[i] = malloc(k * sizeof(char));
	
	for(j=0; j<k; ++j)
		alignment->positionsInd[j] = alignment->positionsInd[polpositions[j]];
	
	for(i=0; i<alignment->sequences; ++i)
	{
		for(j=0; j<k; ++j)
			alignment->poltable[i][j] = alignment->seqtable[i][ polpositions[j] ];
	
		free(alignment->seqtable[i]);

		alignment->seqtable[i] = alignment->poltable[i];

		alignment->poltable[i] = NULL;
	}

	free(alignment->poltable);
	
	fprintf(stdout,"\n\t\tDiscarded sites:\t%d",alignment->segsites-k);

	fprintf(fp,"\n\t\tDiscarded sites:\t%d",alignment->segsites-k);

	alignment->segsites = k;

	free(polpositions);	
}

void removeNonPolymorphicSites(alignment_struct *alignment, FILE * fp, int filterOut, double maf)
{	
	if(alignment->states==2 || alignment->states==3)
		removeNonPolymorphicSitesBIN(alignment,fp,filterOut, maf);	

	if(alignment->states==4 || alignment->states==5)
		removeNonPolymorphicSitesDNA(alignment,fp,filterOut, maf);			
}

int getDataType(int * states)
{
	int i,accum=0;

	if(states[STATESALL-1]!=0) // Invalid character was found
		return -1;
	
	for(i=3;i<STATESALL;i++) 
		accum += states[i];
		
	if((states[0]!=0 || states[1]!=0) && accum==0) // BIN
	{
		if(states[2]!=0) // with gaps, Ns or both
			return 3;
		
		return 2;		
	}

	accum=0;

	for(i=0;i<2;i++)
		accum += states[i];
	
	if((states[3]!=0 || states[4]!=0 || states[5]!=0 || states[6]!=0) && accum==0) // DNA
	{
		if(states[2]!=0) // with gaps, Ns or both
			return 5;
		
		return 4;		
	}

	return -1;
}

void countStates(int *numStates, int *gap, int *states)
{
	int i;

	if (*gap==0)
		if (states[2]!=0)
			*gap=1;

	*numStates=0;

	states[2]=0;

	for(i=0;i<STATESALL;i++)
	{
		if(states[i]!=0)
		{
			states[i]=0;
			(*numStates)++;
		}
	}
}

int determineStates(alignment_struct * alignment)
{
	int i,j, 
	    statesNum=0,
	    states[STATESALL];

	for(i=0; i<STATESALL;++i)
		states[i] = 0;
	
	for(i=0;i<alignment->sequences;i++)
		for(j=0;j<alignment->segsites;j++)		
			states[mapCharToInt(alignment->seqtable[i][j])]++;
		
	statesNum = getDataType(states);

	return statesNum;	
}

char getCharBIN (char input, char state0, char state1)
{
	if(input==GAP)
 		return GAP;
	
	if(input == UN)
	  	return UN;

	if(input==state0)
		return ZERO;

	if(input==state1)
		return ONE;

	return 'X';
}

void convertAlignmentBIN (alignment_struct * alignment)
{
	int i,j;
	char state0, state1, tmp;

	for(j=0;j<alignment->segsites;j++)
	{
		state0='x';
		state1='x'; 		

		for(i=0;i<alignment->sequences;i++)
		{
			tmp = alignment->seqtable[i][j];

			if(state0=='x' && (tmp!=GAP && tmp!=UN)) 
				state0=tmp;

			if(state0!='x' && state1=='x' && tmp!=state0 && tmp!=GAP && tmp!=UN)
				state1=tmp;

			alignment->seqtable[i][j] = getCharBIN(tmp,state0,state1);
		}
	}
}

void binaryDeduction (alignment_struct * alignment)
{
	if(alignment->states==2 || alignment->states==3) // BINARY or BINARY with GAPS
		return;

	int i,j,
            withGaps=0, 
            states[STATESALL],
            siteStates,
            maxStates=-1;

	for(i=0; i<STATESALL;++i)
		states[i] = 0;
	
	for(j=0;j<alignment->segsites;j++)
	{
		for(i=0;i<alignment->sequences;i++)
			states[mapCharToInt(alignment->seqtable[i][j])]++;
		
		countStates(&siteStates, &withGaps, states);

		if (siteStates>maxStates)
			maxStates = siteStates;
	}

	assert(maxStates!=-1);

	if (maxStates==2)
	{
		alignment->states = maxStates + withGaps;
		convertAlignmentBIN(alignment);
	}
	
	return;	
}

void freeAlignment(alignment_struct *alignment, int matrixSizeMax)
{
	int i=0;

	if(alignment->positions != NULL)
		free(alignment->positions);

	if(alignment->positionsInd != NULL)
		free(alignment->positionsInd);

	for(i=0; i<alignment->sequences; ++i)
		if(alignment->seqtable[i] != NULL)
			free(alignment->seqtable[i]);

	if(alignment->seqtable!=NULL)
		free(alignment->seqtable);

	if(alignment->compressedArrays != NULL)
		if (alignment->states==2 || alignment->states==3)
			for(i=0;i<alignment->states-1;i++)
				free(alignment->compressedArrays[i]);

	if(alignment->compressedArrays != NULL)
		if (alignment->states==4 || alignment->states==5)
			for(i=0;i<alignment->states;i++)
				free(alignment->compressedArrays[i]);

	if(alignment->compressedArrays != NULL)
		free(alignment->compressedArrays);
		
	if(alignment->correlationMatrix != NULL)
		for(i=0; i<matrixSizeMax; ++i)
			if(alignment->correlationMatrix[i] != NULL)
				free(alignment->correlationMatrix[i]);

	if(alignment->correlationMatrix != NULL)
		free(alignment->correlationMatrix);
}

char getCharacterImputeBIN (int *input)
{
  	char states[2] = {'0', '1'};
  
  	int total = input[0] + input[1];
  
  	int rn = rand() % total + 1;

 	int j = 0;
  
  	int current = input[j];
  
  	while(rn > current)
      		current += input[++j];
  
  	return states[j];
}

char getCharacterImputeDNA(int *input)
{
  	char states[4] = {'A', 'C', 'G', 'T'};
  
  	/* For DNA, the numbers of A, C, G, T are in positions 3,4,5,6 */
  	int total = input[3] + input[4] + input[5] + input[6];
  
  	int rn = rand() % total + 1;

  	int j = 3;

  	int current = input[j];  

  	while(rn > current)
      		current += input[++j];
  
  	return states[j-3];
}

void imputeStatesBINGAPS(alignment_struct * alignment, int imputeG, int imputeN)
{
	int i,j, states[STATESALL];

	char case1='#', case2='#';

	if(imputeG==1)
		case1 = GAP;

	if(imputeN==1)
		case2 = UN;
	
	for(i=0; i<STATESALL;++i)
		states[i] = 0;
	
	for(j=0;j<alignment->segsites;j++)
	{
		for(i=0; i<STATESALL;++i)
			states[i] = 0;
		
		for(i=0;i<alignment->sequences;i++)
			states[mapCharToInt(alignment->seqtable[i][j])]++;

		if(states[2]!=0)
			for(i=0;i<alignment->sequences;i++)
				if(alignment->seqtable[i][j]==case1 || alignment->seqtable[i][j]==case2)
					alignment->seqtable[i][j] = getCharacterImputeBIN(states);				
	}

	assert(alignment->states==3);

	if(imputeG +  imputeN == 2)
	{
		alignment->states = 2;
		return;
	}	
	
	alignment->states = determineStates(alignment);
}

void imputeStatesDNAGAPS(alignment_struct * alignment, int imputeG, int imputeN)
{
	int i,j, states[STATESALL];

	char case1='#', case2='#';

	if(imputeG==1)
		case1 = GAP;
	
	if(imputeN==1)
		case2 = UN;

	for(i=0; i<STATESALL;++i)
		states[i] = 0;
	
	for(j=0;j<alignment->segsites;j++)
	{
		for(i=0; i<STATESALL;++i)
			states[i] = 0;
		
		for(i=0;i<alignment->sequences;i++)
			states[mapCharToInt(alignment->seqtable[i][j])]++;

		if(states[2]!=0)
			for(i=0;i<alignment->sequences;i++)
				if(alignment->seqtable[i][j]==case1 || alignment->seqtable[i][j]==case2)
					alignment->seqtable[i][j] = getCharacterImputeDNA(states);			
	}

	assert(alignment->states==5);

	if(imputeG +  imputeN == 2)
	{
		alignment->states = 4;
		return;
	}
	
	alignment->states = determineStates(alignment);
}

void imputeStates (alignment_struct * alignment, int imputeG, int imputeN)
{
	if (imputeG==0 && imputeN==0)
		return;

	if(alignment->states==3)
		imputeStatesBINGAPS(alignment, imputeG, imputeN);
	
	if(alignment->states==5)
		imputeStatesDNAGAPS(alignment, imputeG, imputeN);
	
	return;
}

void switchValues (int * input, char * sortedStates, int index0, int index1)
{
	int tmp;
        char tmpC;

	if(input[index1]>input[index0])
	{
		tmp = input[index0];
		tmpC = sortedStates[index0];
		input[index0] = input[index1];
		sortedStates[index0] = sortedStates[index1];
		input[index1] = tmp;
   		sortedStates[index1] = tmpC;
	}
}

void sort(int * freqs, char * sortedStates)
{
	switchValues (freqs, sortedStates, 0, 1);
	switchValues (freqs, sortedStates, 2, 3);
	switchValues (freqs, sortedStates, 0, 2);
	switchValues (freqs, sortedStates, 1, 3);
	switchValues (freqs, sortedStates, 1, 2);	
}

int minorToMajor(int* inp, char* outp)
{
 	char major, alpha[4] = {'A', 'C', 'G', 'T'};
   
 	int total = 0, rv, i;
  
	sort(inp,alpha); 
  
	assert(inp[0] >= inp[1] && inp[1] >= inp[2] && inp[2] >= inp[3]);
	
  	if(inp[2] + inp[3] == 0)
    		return -1;

  	total = inp[0] + inp[1];
  
  	for(i=0; i<2; ++i)
    	{
	      rv= rand() % total + 1;
	      major = (rv <= inp[0] ) ? alpha[0] : alpha[1];
	      outp[i] = alpha[ i+2 ];
	      outp[i+2] = major;
    	}

    	return 1;
}

int numberOfStates(alignment_struct *alignment, int segsites)
{
  	int A=0,C=0,G=0,T=0,i,states=0;
  

  	for(i=0; i<alignment->sequences; ++i)
    	{
	      if(alignment->seqtable[i][segsites] == 'A') A++;
	      if(alignment->seqtable[i][segsites] == 'C') C++;
	      if(alignment->seqtable[i][segsites] == 'G') G++;
	      if(alignment->seqtable[i][segsites] == 'T') T++;
    	}

  	if(A > 0) states++;
  	if(C > 0) states++;
  	if(G > 0) states++;
  	if(T > 0) states++;

  	return states;
}

void binaryDeductionUser (alignment_struct * alignment, int binary)
{
	if(binary==0)
		return;

	if(alignment->states==2 || alignment->states==3)
		return;

	int i,j, states[4];

	char replaceMinorStates[4];

	for(i=0; i<4;++i)
		states[i] = 0;
	
	for(j=0;j<alignment->segsites;j++)
	{
		for(i=0; i<4;++i)
			states[i] = 0;

		for(i=0; i<alignment->sequences; i++)
		{
			int map = mapCharToInt(alignment->seqtable[i][j]);

			if(map >=3 && map <=6)
				states[mapCharToInt(alignment->seqtable[i][j])-3]++;		      	
		}
		
		if(minorToMajor(states, replaceMinorStates)==1)
		{
			for(i=0;i<alignment->sequences;i++)
			{
				if(alignment->seqtable[i][j]==replaceMinorStates[0])
					alignment->seqtable[i][j] = replaceMinorStates[2];

				if(alignment->seqtable[i][j]==replaceMinorStates[1])
					alignment->seqtable[i][j] = replaceMinorStates[3];
			}			
		}
	}
}


void skipToNextAlignmentDelimiter(FILE *fp, char d)
{
  char ent = fgetc( fp );
  
  while( ent != d && ent != EOF)
    {
      ent = fgetc( fp );
    }
  
  
}

int skipLine (FILE * fp)
{
	char tmp;

	while( (tmp = fgetc(fp) ) != EOF)
		if(tmp == 10 || tmp == 13)
			break;
		
	if(tmp == EOF)
		return 1;

	return 0;	
}



int readAlignmentMS(FILE *fp, alignment_struct *alignment, int imputeG, int imputeN, int binary, FILE * fpInfo, int filterOut, double maf, int fileFormatMBS)
{
	int i,y, DIM=2;

	char ent, stringtemp[100];

	

	/* get rid of the first line information */
	while( (ent = fgetc(fp) ) != '\n');
	  
	int temp = fscanf(fp,"%s %d", stringtemp, &alignment->segsites);
	temp = temp;

	alignment->segsites = fileFormatMBS==1? alignment->segsites-1:alignment->segsites;

	if(strcmp(stringtemp, "segsites:") != 0 )
	  {
	    fprintf(stderr, "ERROR! \"segsites:\" are expected but %s has been found\n", stringtemp);
	    assert(0);
	  }

	/* ent = fgetc( fp ); */

	/* while (ent != 10 && ent != 13) */
	/*   ent = fgetc( fp ); */

	if( skipLine( fp ) == 1)
	  {
	    fprintf(stderr, "Unexpected end of file after segsites....\n\n");
	    assert(0);
	  }

	alignment->positions = malloc(sizeof(float)*alignment->segsites);
	assert(alignment->positions!=NULL);

	alignment->positionsInd = malloc(sizeof(int)*alignment->segsites);
	assert(alignment->positionsInd!=NULL);

	if( alignment->segsites > 0)
	  {
	    temp = fscanf(fp, "%s", stringtemp);
	    if(strcmp( stringtemp, "positions:") != 0)
	      {
		fprintf(stderr, "ERROR: \"positions:\" is expected\n\n");
		assert(0);
	      }
	  }



	for(i=0;i<alignment->segsites;i++)
		temp = fscanf(fp,"%f",&alignment->positions[i]);		

	if(fileFormatMBS==1)
	{
		for(i=0;i<alignment->segsites;i++)
		{
			alignment->positionsInd[i] = (int)(alignment->positions[i]);
			printf("%d\t", alignment->positionsInd[i]);
		}
	}
	else
	{
		for(i=0;i<alignment->segsites;i++)
			alignment->positionsInd[i] = (int)(alignment->positions[i] * (float)alignment->length);		
	}
	alignment->seqtable = malloc(sizeof(char*)*DIM);
	assert(alignment->seqtable!=NULL);

	if(alignment -> segsites == 0 )
	  skipToNextAlignmentDelimiter(fp, '/');
	
 	y=0;
	/* do this loop only when there are segsites */
	while( alignment->segsites > 0)
	{
		ent = fgetc(fp);
		
		while(ent=='\n' || ent==' ' || ent==13)
			ent = fgetc(fp);		

		if(ent=='/' || ent==EOF)
			break;		
	
		if(y>=DIM)
		{
			DIM += DIM;
			alignment->seqtable = realloc(alignment->seqtable,sizeof(char*)*DIM);
		}
		
		alignment->seqtable[y] = (char*)malloc(sizeof(char)*alignment->segsites);
		
		i=0;
		while(ent!='\n' && ent!=' ' && ent!=13)
		{
			alignment->seqtable[y][i++]=ent;
			ent = fgetc(fp);
		}
		y++;
	}

	alignment->sequences=y;

	assert(alignment->segsites<=alignment->length);
	
	fprintf(stdout,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	fprintf(fpInfo,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	alignment->states = determineStates(alignment);

	if(alignment->states == -1)
	{
	  fprintf(stderr, "WARNING: Alignment seems to be empty (no states found)\n\n");
	  return 0;
	}	

	removeNonPolymorphicSites(alignment, fpInfo, filterOut, maf);

	if( alignment -> segsites < MINSNPS_THRESHOLD)
	  {
	    return 0;
	  }

	imputeStates(alignment, imputeG, imputeN);

	binaryDeductionUser(alignment, binary);

	binaryDeduction(alignment);

	return 1;
}

void ignoreSpaces(FILE * fp, char * ent)
{
	while(*ent==10 || *ent==13 || *ent==32)
		*ent = fgetc(fp);
}

int readAlignmentMACS(FILE *fp, alignment_struct *alignment, int imputeG, int imputeN, int binary, FILE *fpInfo, int filterOut, double maf)
{
	int  i, j, DIM = 1, DIM2 = 2, prevDIM;
	char ent;
	int indCounter = 0, ninds = 0, nsnp = 0;
	int sitevar;
	char siteflag[100];

	char ** mySnpTable = malloc(DIM*sizeof(char*));

	mySnpTable[0] = malloc(DIM2*sizeof(char));

	alignment->positions = malloc(DIM*sizeof(float));

	// read first line 
	int temp = fscanf(fp, "%f", &alignment->positions[0]);
	temp = temp;

	// read first SNP
	while((ent = fgetc(fp)))
	{
		ignoreLineSpaces(fp, &ent);

		if(isEndOfLine(ent))
			break;

		if(indCounter >= DIM2)
		{
			DIM2 += DIM2;
			mySnpTable[0] = realloc(mySnpTable[0], DIM2*sizeof(char));
		}

		mySnpTable[0][indCounter] = ent;

		indCounter++;
	}

	ninds = indCounter;

	nsnp = 1;

	// now rest of the lines
	while(1)
	{
		if(nsnp >= DIM)
		{
			prevDIM = DIM;

			DIM += DIM;

			mySnpTable = realloc(mySnpTable, DIM*sizeof(char*));

			alignment->positions = realloc(alignment->positions, DIM*sizeof(float));

			for(i=prevDIM; i<DIM; ++i)
				mySnpTable[i] = malloc(ninds+1 * sizeof(char));
		}

		temp = fscanf(fp, "%s",siteflag);

		if(!strcmp(siteflag,"SITE:"))
		{
			temp = fscanf(fp, "%d %f %s",&sitevar, &alignment->positions[nsnp], mySnpTable[nsnp]);

			if(!strcmp(siteflag, "SITE:"))
				nsnp++;
		}
		else
			break;
	}

	alignment->seqtable = malloc(ninds * sizeof(char*));

	for(i=0; i<ninds; ++i)
		alignment->seqtable[i] = malloc(nsnp * sizeof(char));

	for(i=0; i<nsnp; ++i)
	{
		for(j=0; j<ninds; ++j)
			alignment->seqtable[j][i] = mySnpTable[i][j];		
	}

	for(i=0; i<DIM; ++i)
		free(mySnpTable[i]);

	free(mySnpTable);

	alignment->sequences = ninds;

	alignment->segsites = nsnp;

	assert(alignment->segsites<=alignment->length);

	fprintf(stdout,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	fprintf(fpInfo,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	alignment->states = determineStates(alignment);

	if(alignment->states==-1)
	{
		fprintf(stderr, "\n\n ERROR: Empty alignment? ( no states found ).\n\n");
		exit(0);
	}

	alignment->positionsInd = malloc(alignment->segsites*sizeof(int));

	for(i=0;i<alignment->segsites;i++)
		alignment->positionsInd[i] = (int)(alignment->positions[i] * (float)alignment->length);

	removeNonPolymorphicSites(alignment, fpInfo, filterOut, maf);

	if( alignment->segsites < MINSNPS_THRESHOLD)
	  {
	    
	    return 0;
	  }

	imputeStates(alignment, imputeG, imputeN);

	binaryDeductionUser(alignment, binary);

	binaryDeduction(alignment);
	
	return 1;
}

int isValidDNACharacter(char input)
{
	if(input==AD) return 1;
	if(input==CY) return 1;
	if(input==GU) return 1;
	if(input==TH) return 1;
	if(input==GAP) return 1;
	if(input==UN) return 1;

	if(input=='X' || input=='K' || 
	   input=='M' || input=='R' || 
	   input=='Y' || input=='S' || 
	   input=='W' || input=='B' || 
	   input=='V' || input=='H' || 
	   input=='D' || input=='.') return 1;

	return 0;
} 

int isValidVCFBase(char input)
{
	if(input=='A' || input=='C' || 
	   input=='G' || input=='T' || 
	   input=='N' || input=='a' || 
	   input=='c' || input=='g' || 
	   input=='t' || input=='n' || input == '.') return 1;

	return 0;
}

int scanStateVector(char * stateVector, char X)
{
	int i;

	for(i=0;i<MAX_STATES_VCF;i++)
		if(stateVector[i]==X)
			return 1;

	return 0;
}

int getStatesNum (char * stateVector)
{
	int i, states=0;
	
	for(i=0;i<MAX_STATES_VCF;i++)
	{
		if(stateVector[i]!='X')
			states++;
		else
			return states;
	}

	return states;
}

int getStatesREF (char * string, char * stateVector, int line)
{
	int i, j, index=0, elen=0, slen=strlen(string);

	char CharToStore = 'X';

	for(j=0;j<MAX_STATES_VCF;j++)
		stateVector[j] = CharToStore;

	for(i=0;i<slen;i++)
	{

		if(string[i]!=',')
		{
		
			if(string[i]=='<')
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
			
			assert(isValidVCFBase(string[i])==1);

			elen++;

			CharToStore = string[i];
		}

		if(string[i]==',' || i == slen-1)
		{
			if(elen==1)
			{
				if(scanStateVector(stateVector, CharToStore)==0)
				{
					stateVector[index++] = CharToStore;
					elen = 0;	
				}
				else
				{
					fprintf(stderr, "\n\n ERROR: Nucleotide %c (field REF) in line %d appears twice.\n\n", CharToStore, line);
					assert(0); 
				}			
			}
			else
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
		}
	}
	return 1;
}

int getStatesALT (char * string, char * stateVector, int line)
{
	int i, j, index=0, elen=0, slen=strlen(string);

	char CharToStore = 'X';

	for(j=0;j<MAX_STATES_VCF;j++)
		if(stateVector[j] == 'X')
		{
			index = j;
			break;
		}

	assert(index!=0);

	for(i=0;i<slen;i++)
	{

		if(string[i]!=',')
		{

			if(string[i]=='<')
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}

			assert(isValidVCFBase(string[i])==1);

			elen++;

			CharToStore = string[i];
		}

		if(string[i]==',' || i == slen-1)
		{
			if(elen==1)
			{
				if(scanStateVector(stateVector, CharToStore)==0)
				{
					stateVector[index++] = CharToStore;
					elen = 0;	
				}
				else
				{
					fprintf(stderr, "\n\n ERROR: Nucleotide %c (field ALT) in line %d appears twice (in REF or ALT).\n\n", CharToStore, line);
					assert(0);
				}			
			}
			else
			{
				for(j=0;j<MAX_STATES_VCF;j++)
					stateVector[j] = 'X';

				return 0;
			}
		}
	}
	return 1;
}

float getValueAF (char * string, int line)
{
  
  int i, j, len = strlen(string), k = 0;
  char AF_s [MAXAFLENGTH+1];
  AF_s[0] = 0;
	
  float AF=-1.0, freqs[1000];

	

  if(strcmp(".",string)==0 || len<4)
    return -1;

  for(i=0;i<len-3;i++)
    {	
      if(string[i]=='A' && string[i+1]=='F' && string[i+2]=='=')
	{
	  if(i==0)
	    break;
	  else
	    if(string[i-1]==';')
	      break;
	}
    }

  if(i==len-3)
    return -1.0;

  i = i+3;

  for(j=0;j<MAXAFLENGTH;j++)
    {
	  
	  
      if( string[i] == ',')
	{
	  assert( AF_s[0] != 0 );
	  
	  freqs[k] = atof( AF_s );
	  
	  j = 0;
	  
	  ++i;
	  
	  ++k;
	}

      if(string[i]==';' || i==len)
	break;
      
      if(string[i]=='.' || (string[i]>=48 && string[i]<=57) || string[i] == 'e' || string[i] == '-' || string[i] == '+')
	{
	  AF_s[j] = string[i];
	  i++;
	  AF_s[j+1] = 0;
	}
      else
	{
	  fprintf(stderr, "\n\n ERROR: Invalid character (%c) in AF of field INFO in line %d.\n\n", string[i], line );
	  assert(0);
	}
    }
	
  if(strlen(AF_s)==0)
    {
      fprintf(stderr, "\n\n ERROR: AF in field INFO in line %d has no value.\n\n", line);
      assert(strlen(AF_s)!=0);
    }

  
  if(k == 0)
    AF = atof(AF_s);
  else
    {

      freqs[k] = atof(AF_s);
      ++k;
      
      AF = 0.;
      
      for( i = 0; i < k; ++i)
	{
	  if( freqs[i] > 0. && freqs[i] < 1.)
	    {
	      AF = freqs[i];
	      break;
	    }
	}
    }

  if(AF < 0.0 || AF > 1.0)
    {
      fprintf(stderr, "\n\n ERROR: AF (%f) in line %d should be >= 0.0 and <= 1.0.\n\n", AF, line);
      assert(AF>=0.0 && AF<=1.0);
    }

  return AF;
}

int checkVTisSNP (char * string)
{
	int i, len = strlen(string);

	if(strcmp(".",string)==0 || len<4)
		return -1;

	for(i=0;i<len-3;i++)
	{
		if(string[i]=='V' && string[i+1]=='T' && string[i+2]=='=')
		{
			if(i==0)
				break;
			else
				if(string[i-1]==';')
					break;
		}	
	}

	if(i==len-3)
		return -1;

	i = i+3;

	assert(i<=len-3);

	if(string[i]=='S' && string[i+1]=='N' && string[i+2]=='P')
		return 1;
	else
		return 0;
}

int getGTpos(char * string, int line)
{
	int i, len = strlen(string), GTposition = 0;

	for(i=0;i<len-1;i++)
	{
		if(string[i]==':')
			GTposition ++;

		if(string[i]=='G' && string[i+1]=='T')
		{
			if(i==0)
			{				
				if(i+2==len)
					return GTposition;

				else
					if(string[i+2]==':')
						return GTposition;
			}
			else
			{
				if(string[i-1]==':')
				{
					if(i+2==len)
						return GTposition;
					else
						if(string[i+2]==':')
							return GTposition;
				}
			}
		}			
	}

	GTposition = -1;

	return GTposition;
}

int getGTfield (char * string, int GTpos)
{
	int i=0, pos=GTpos, len = strlen(string),j=0, counter=0;

	assert(pos!=-1);	

	for(i=0;i<len;i++)
	{
		if(string[i]!=':')
		{	
			if(pos==0)
			{
				string[j++]=string[i];

				if(string[i]=='|' || string[i]=='/')
					counter++;
			}
		}
		else
		{
			pos--;
			
			if(pos==-1)
			{
				string[j]=0;
				break;
			}
		}
	}

	assert(strlen(string) > 0);	

	return counter+1;
}

void dataShuffleKnuth(char * data, int startIndex, int endIndex)
{
	if(startIndex == endIndex)
		return;

	int i, index;
	char tmp;

	for (i = endIndex; i > startIndex; i--)
	{
		index = startIndex + (rand() % (i - startIndex + 1));

		tmp = data[index];
		data[index] = data[i];
		data[i] = tmp;
	}
}

void getGTdata (char * string, char * stateVector, int statesTotal, char * sampleData)
{
	int i, j=0, index=0, start=0, end=0, len = strlen(string);

	for(i=0;i<len;i++)
	{	
		if(string[i]>=48 && string[i]<=57)
		{
			index = string[i]-48;

			assert(index<statesTotal);
			
			sampleData[j++] = stateVector[index];
		}
		else
		{
  			if(string[i]=='.')
			{
				sampleData[j++] = 'N';
			}

			if(string[i]=='/')
			{
				end++;
			}

			if(string[i]=='|')
			{
				dataShuffleKnuth(sampleData, start, end);
				start = j;
				end = j;
			}			
		}
	}

	dataShuffleKnuth(sampleData, start, end);
	
}

void processSampleVCF (char * string, int GTpos, char ** sampleData, int * sampleDataMemSize, char * stateVector, int statesTotal)
{
	
	int i, dataSize = getGTfield (string,GTpos);

	if(dataSize+1>*sampleDataMemSize)
	{
		*sampleDataMemSize = dataSize+1;

		*sampleData = realloc((*sampleData), (*sampleDataMemSize) * sizeof(char));		
	}
	
	for(i=0;i<(*sampleDataMemSize);i++)
		(*sampleData)[i]=0;

	getGTdata (string, stateVector, statesTotal, *sampleData);

	(*sampleData)[dataSize]=0;
}

int readLine_VCF (FILE * fp, char ** string, int lineIndex, alignment_struct * alignment, int * DIM, char ** sampleData, int * sampleDataMemSize, int * SNP_SZ, char ** tmpSNPtable, int * tmpSNPtableLineIndex)
{

	int eol=0, eof=0, maxLength=STRINGLENGTH;
	int elementIndex = -1, lineSkipped=0, isEOF = 1, i,j=0, s;

	char stateVector[MAX_STATES_VCF], tmp;
	int position=-1, statusREF, statesREF=0, statusALT, statesALT=0, VTisSNP=0, GTpos=-1, sampleIndex=-1;
	float AF=0.5;


	if(lineIndex==0)
		elementIndex++;	

	while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
	{
		
		elementIndex++;

		switch(elementIndex)
		{
			case 0: // CHROM
				if(strcmp(VCF_alignment_name,*string)!=0)
				{
					strncpy(VCF_alignment_name, *string, MAX_CHROM_NAME_VCF);

					assert(strlen(VCF_alignment_name)!=0);

					return 0;
				}
				break;

			case 1: // POS
				
				position = atoi(*string);

				break;


			case 2: // ID
	
				break;


			case 3: // REF
				
				statusREF = getStatesREF (*string, stateVector, lineIndex+1+VCF_header_lines);
				statusREF = statusREF;
				statesREF = getStatesNum (stateVector);
			
				if(statesREF==0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					break;
				}
			
				break;


			case 4: // ALT
				
				statusALT = getStatesALT (*string, stateVector, lineIndex+1+VCF_header_lines);
				statusALT = statusALT;
				statesALT = getStatesNum (stateVector);
				
				if(statesALT != 0 && stateVector[ statesALT - 1] == '.')
				  for( s = statesALT - 1; s < statesALT-1 + statesREF; ++s)
				    stateVector[s] = stateVector[s-statesREF];
				
				if(statesALT==0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					break;
				}

				break;


			case 5: // QUAL
				
				break;


			case 6: // FILTER
				
			  if(strcmp("PASS",*string)!=0 && strcmp(".", *string) != 0)
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					break;
				}

				break;


			case 7: // INFO
				
			  if( strcmp( *string, ".") == 0)
			    {
			      AF = .9;
			      break;
			    }
			  
			  AF = getValueAF(*string, lineIndex+1+VCF_header_lines); // returns -1.0 if AF not found
			  
			  if(AF==0.0 || AF==1.0)
			    {
			      isEOF = skipLine(fp);
			      
			      if (isEOF == 1)
				return -1;
			      
			      lineSkipped=1;
			      
			      break;
			    }
			  
			  VTisSNP = checkVTisSNP(*string); // returns -1.0 if VT not found
			  
			  if(VTisSNP==0)
			    {
			      isEOF = skipLine(fp);
			      
			      if (isEOF == 1)
				return -1;
			      
			      lineSkipped=1;
			      
			      break;
			    }
			  
			  break;
			  

			case 8: // FORMAT

				GTpos = getGTpos(*string, lineIndex+1+VCF_header_lines); // returns -1 if GT not found

				if(GTpos==-1)		
				{
					isEOF = skipLine(fp);

					if (isEOF == 1)
						return -1;

					lineSkipped=1;

					break;
				}
			
				break;

			default:

				break;

		}

		if(elementIndex == VCF_HLENGTH-1 || eol==1 || eof==1 || lineSkipped==1)		
			break;
	}

	if (eof== 1 || position == -1)
		return -1;


	if(lineSkipped!=1)
	{

		(*tmpSNPtableLineIndex)++;

		alignment->segsites++;
	
		if(alignment->segsites>=2)
			assert(position>=alignment->positionsInd[alignment->segsites-2]);

		if(alignment->segsites>=*DIM)
		{
			*DIM = *DIM<<1;
			alignment->positionsInd = (int *) realloc(alignment->positionsInd,*DIM*sizeof(int));
		}

		alignment->positionsInd[alignment->segsites-1] = position;

		sampleIndex =0;

		while(getNextString (fp, string, &eol, &eof, &maxLength)==1)
		{

			assert(sampleIndex<alignment->VCFsamples);

			if(alignment->VCFsample_valid[sampleIndex]==1)
			{

				processSampleVCF(*string, GTpos, sampleData, sampleDataMemSize, stateVector, statesALT);

				if(*tmpSNPtableLineIndex==0)
				{
					*SNP_SZ += strlen(*sampleData);

					tmpSNPtable[*tmpSNPtableLineIndex] = (char *) realloc(tmpSNPtable[*tmpSNPtableLineIndex],*SNP_SZ*sizeof(char));
				}
	
				for(i=0;i<strlen(*sampleData);i++)
				{				
		
					if(j>=*SNP_SZ)
					{
						fprintf(stderr, "\n\n ERROR: There are more than %d nucleotides in line %d. Expected %d (according to the first SNP in line %d).\n\n", j, lineIndex+1+VCF_header_lines, *SNP_SZ, VCF_first_SNP );
						assert(j<*SNP_SZ);
					}

					tmpSNPtable[*tmpSNPtableLineIndex][j++] = (*sampleData)[i];
				}
			}
			sampleIndex++;

			if(eol==1 || eof==1)
				break;
		}

		
		if(eol == 1)
		  {
		    tmp = fgetc(fp);
		    ungetc(tmp, fp);
		    ignoreNewLineChars(fp, &tmp);
		  }


		
		if(j!=*SNP_SZ)
		{

			fprintf(stderr, "\n\n ERROR: There are %d nucleotides in line %d. Expected %d (according to the first SNP in line %d).\n\n", j, lineIndex+1+VCF_header_lines, *SNP_SZ, VCF_first_SNP );
			assert(j==*SNP_SZ);
		}
		if (eof== 1)
			return -1;
	}

	return 1;
}



int readAlignmentVCF(FILE *fp, alignment_struct * alignment, int imputeG, int imputeN, int binary, FILE *fpInfo, int filterOut, double maf)
{
	char ** string = (char **) malloc (sizeof(char*));
	(*string) = (char *) malloc(sizeof(char)*STRINGLENGTH);

	int DIM=2;
	alignment->positionsInd = (int *) malloc(DIM*sizeof(int));

	alignment->segsites = 0;

	printf("\t\tChromosome:\t\t%s\n\n", VCF_alignment_name);

	strcpy(VCF_alignment_name_cur, VCF_alignment_name);

	int inAlignment = 1, lineIndex=-1, i, j;

	int sampleDataMemSize = 2;
	char ** sampleData = (char **) malloc(sizeof(char*));
	(*sampleData) = (char *) malloc(sizeof(char)*sampleDataMemSize);

	int SNP_SZ = 0;
	int SNP_NUM = 10;
	char ** tmpSNPtable = (char **) malloc(SNP_NUM*sizeof(char*));
	tmpSNPtable[0] = (char *) malloc(sizeof(char));
	int prevtmpSNPtableLineIndex,tmpSNPtableLineIndex=-1;

	int memCounter = 1;

	while(inAlignment==1)
	{
		lineIndex++;

		prevtmpSNPtableLineIndex = tmpSNPtableLineIndex;

		inAlignment = readLine_VCF (fp, string, lineIndex, alignment, &DIM, sampleData, &sampleDataMemSize, &SNP_SZ, tmpSNPtable, &tmpSNPtableLineIndex);

		if(inAlignment==1)
		{
			if(tmpSNPtableLineIndex+1>=SNP_NUM)
			{
				SNP_NUM += 10;
				tmpSNPtable = realloc(tmpSNPtable, SNP_NUM*sizeof(char*));
			}

			if(prevtmpSNPtableLineIndex!=tmpSNPtableLineIndex)
			{
				if(prevtmpSNPtableLineIndex == -1)
					VCF_first_SNP = lineIndex + VCF_header_lines + 1;

				tmpSNPtable[tmpSNPtableLineIndex+1] = (char *) malloc(sizeof(char)*SNP_SZ);
				memCounter++;
			}
		}
	}

	nxtVCFalignment = inAlignment;

	if(nxtVCFalignment==0)
		VCF_header_lines = lineIndex + VCF_header_lines;

	alignment->sequences = SNP_SZ;

	alignment->seqtable = malloc(alignment->sequences * sizeof(char*));

	alignment->length = alignment->positionsInd[alignment->segsites-1];

	assert(alignment->segsites<=alignment->length);

	for(i=0;i<alignment->sequences;i++)
		alignment->seqtable[i] = malloc(alignment->segsites * sizeof(char));

	for(i=0;i<alignment->segsites;i++)
		for(j=0;j<alignment->sequences;j++)
			alignment->seqtable[j][i] = tmpSNPtable[i][j];		

	for(i=0;i<memCounter;i++)
		free(tmpSNPtable[i]);

	free(tmpSNPtable);

	free(*sampleData);

	free(sampleData);

	free(*string);

	free(string);

	fprintf(stdout,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	fprintf(fpInfo,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	fflush(stdout);

	alignment->states = determineStates(alignment);

	if(alignment->states==-1)
	{
		fprintf(stderr, "\n\n ERROR: Empty alignment file (no states found).\n\n");
		exit(0);
	}

	removeNonPolymorphicSites(alignment, fpInfo, filterOut, maf);

	if(alignment->segsites < MINSNPS_THRESHOLD)
	  {
	    return 0;
	  }

	imputeStates(alignment, imputeG, imputeN);

	binaryDeductionUser(alignment, binary);

	binaryDeduction(alignment);

	return 1;
}

char filterAmb (char input)
{
	if(input=='X' || input=='K' || 
	   input=='M' || input=='R' || 
	   input=='Y' || input=='S' || 
	   input=='W' || input=='B' || 
	   input=='V' || input=='H' || 
	   input=='D' || input=='.') return UN;

	return input;
}

int readFASTASequence (FILE * fp, alignment_struct * alignment, int sequenceIndex)
{
	int i, DIM = 100,
            curLength=0;

	char seqName[SEQNAMESIZE], ent;
	
	char * temp = fgets(seqName, SEQNAMESIZE, fp);
	temp = temp;
	temp=NULL;
	
	if (alignment->segsites==-1)
	{
		alignment->seqtable[sequenceIndex] = malloc(sizeof(char)*DIM);
		
		ent = toupper(getc(fp));	

		ignoreSpaces(fp,&ent);

		while(ent != EOF && ent != '>' && ent != '/' )
		{			
			if(isValidDNACharacter(ent))
			{
				if(curLength < DIM)
					alignment->seqtable[sequenceIndex][curLength++] = filterAmb (ent);
				else
				{
					DIM += DIM;
					alignment->seqtable[sequenceIndex] = realloc(alignment->seqtable[sequenceIndex], DIM*sizeof(char));

					alignment->seqtable[sequenceIndex][curLength++] = filterAmb (ent);
				}
			}
			else
			{
				if(ent >= 33 && ent <= 126)
				{
					fprintf(stderr,"\n ERROR: Invalid character (%c) at position %d of sequence %d.\n\n", ent, curLength+1, sequenceIndex+1);
					exit(0);
				}				
			}

			ent = toupper(getc(fp));			
		}

		if(ent == EOF || ent == '/')
		{
			fprintf(stderr,"\n ERROR: Alignment contains only one sequence.\n\n");
			exit(0);
		}
		
		alignment->segsites = curLength;
	}
	else
	{
		alignment->seqtable[sequenceIndex] = malloc(sizeof(char)* alignment->segsites);
		
		ent = toupper(getc(fp));

		ignoreSpaces(fp,&ent);

		i = 0;
		while(ent != EOF && ent != '>' && ent != '/' )
		{						
			if(isValidDNACharacter(ent))
				alignment->seqtable[sequenceIndex][i++] = filterAmb (ent);
			else
				if(ent >= 33 && ent <= 126)
				{
					fprintf(stderr,"\n ERROR: Invalid character (%c) at position %d of sequence %d.\n\n", ent, i+1, sequenceIndex+1);
					exit(0);				
				}

			ent = toupper(getc(fp));
		}
		
		if(i != alignment->segsites)
		{
			fprintf(stderr, "\n ERROR: Length of sequence %d is %d. Expected %d.\n\n", sequenceIndex+1, i, alignment->segsites);
			exit(0);	
		}	
	}

	if(ent=='>')
		return 1;

	return 0;
}

int readAlignmentFASTA(FILE *fp, alignment_struct *alignment, int imputeG, int imputeN, int binary, FILE * fpInfo, int filterOut, double maf)
{
	int y, DIM=2, nxt_seq = 1, i;	

	alignment->segsites = -1;

	alignment->seqtable = malloc(sizeof(char*)*DIM);
	
 	y=0;
	while(nxt_seq==1)
	{
		if(y>=DIM)
		{
			DIM += DIM;
			alignment->seqtable = realloc(alignment->seqtable,sizeof(char*)*DIM);
		}

		nxt_seq = readFASTASequence (fp, alignment, y);
	
		y++;
	}
	 
	alignment->sequences=y;
	
	alignment->length = alignment->segsites;

	alignment->positionsInd = malloc(sizeof(int)*alignment->segsites);

	for(i=0;i<alignment->segsites;i++)
		alignment->positionsInd[i]=i+1; 
	
	fprintf(stdout,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	fprintf(fpInfo,"\t\tSequences:\t\t%d\n\t\tSites:\t\t\t%d", alignment->sequences,alignment->segsites);

	alignment->states = determineStates(alignment);

	if(alignment->states==-1)
	{
		fprintf(stderr, "\n\n ERROR: Empty alignment file (no states found).\n\n");
		exit(0);
	}	

	removeNonPolymorphicSites(alignment,fpInfo, filterOut, maf);

	if( alignment -> segsites < MINSNPS_THRESHOLD)
	  {
	    return 0;
	  }
	
	imputeStates(alignment, imputeG, imputeN);

	binaryDeductionUser(alignment, binary);

	binaryDeduction(alignment);

	return 1;
}

int readAlignment(FILE *fp, alignment_struct *alignment, int imputeG, int imputeN, int binary, int format, FILE * fpInfo, int filterOut, double maf, int fileFormatMBS)
{
	if(format==MS_FORMAT)
		return readAlignmentMS(fp, alignment, imputeG, imputeN, binary, fpInfo, filterOut, maf, fileFormatMBS);

	if(format==FASTA_FORMAT)
		return readAlignmentFASTA(fp, alignment, imputeG, imputeN, binary, fpInfo, filterOut, maf);
		
	if(format==MACS_FORMAT)
		return readAlignmentMACS(fp, alignment, imputeG, imputeN, binary, fpInfo, filterOut, maf);

	if(format==VCF_FORMAT)
		return readAlignmentVCF(fp, alignment, imputeG, imputeN, binary, fpInfo, filterOut, maf);		

	return 0;
} 

void checkSNIPPositions (FILE* fp, alignment_struct * alignment, int index)
{
	int i,j=0;

	int t1 = alignment->positionsInd[0],
	    t2;
	
	fprintf(fp, "\n// Alignment %d\n\n",index);
	
	for(i=1;i<alignment->segsites;i++)
	{
		t2 = alignment->positionsInd[i];
	
		if (t2==t1)	
			fprintf(fp, " SNPs %d and %d correspond to the same alignment position: %d\n",j,i,alignment->positionsInd[i]);	
		
		t1 = t2;
                
		j=i;		
	}
	fprintf(fp,"\n");	
}

void initializeGlobalPointers(alignment_struct* alignment)
{
	alignment->positions = NULL;
	alignment->positionsInd = NULL;
	alignment->seqtable = NULL;
	alignment->poltable = NULL;
	alignment->compressedArrays=NULL;
	alignment->correlationMatrix=NULL;
	alignment->VCFsample_valid=NULL;	
}
