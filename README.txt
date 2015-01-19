

===========================================================
OmegaPlus: A scalable tool for rapid detection of selective
	   sweeps in whole-genome datasets 
===========================================================


Copyright (C) 2012 Pavlos Pavlidis and Nikolaos Alachiotis 


OmegaPlus implements a linkage disequilibrium-based approach
to compute the omega statistic to detect selective sweeps 
in whole-genome datasets.


Linux Platforms
---------------

Compile:
	make -f Makefile.gcc|Makefile.PTHREADS.FINE.gcc|Makefile.PTHREADS.COARSE.gcc

Execute the sequential version:
	./OmegaPlus -name TEST -input TEST.fasta -minwin 100 -maxwin 1000 -grid 10000

or the fine-grained parallel version:
	./OmegaPlus-F -name TEST -input TEST.fasta -minwin 100 -maxwin 1000 -grid 10000 -threads 4

or the coarse-grained parallel version:
	./OmegaPlus-C -name TEST -input TEST.fasta -minwin 100 -maxwin 1000 -grid 10000 -threads 4

or the multi-grained parallel version:
	./OmegaPlus-M -name TEST -input TEST.fasta -minwin 100 -maxwin 1000 -grid 10000 -threads 4


Windows Platforms
-----------------

The code has been compiled with Microsoft Visual Studio 2010.

A statically linked executable is provided.

Open a terminal and run as follows:

	OmegaPlus.exe -name TEST -input TEST.fasta -minwin 100 -maxwin 1000 -grid 10000

Only the sequential version has been compiled for Windows.


Change Log
----------

(for an updated change log use the -v or -version command line flags)

February 2012:	OmegaPlus v1.0.0

March 2012:	OmegaPlus v1.0.1
		-- Minor bug fixes
		-- Run examples

March 2012:	OmegaPlus v2.0.0
		-- Multi-grain parallelism
		-- Omega approximate search
		-- MINSNPPERWINDOW (minimum number of SNPs per window) set to 5

July 2012:	OmegaPlus v2.1.0
		-- Supports VCF input file format
		-- Minor bug fixes

October 2012:	OmegaPlus v2.2.0
		-- no-singletons command line flag to exclude the singletons

November 2012:	OmegaPlus v2.2.1
		-- Minor bug fix in the ms parser

January 2013:	OmegaPlus v2.2.2
		-- Fixed bug in the VCF parser associated with the handling of missing data

January 2014:	OmegaPlus v2.2.5
		-- Small input changes, it can now read the '.' state in the ALT field



Citations
---------

If you are publishing results based on OmegaPlus please cite:

N. Alachiotis, A. Stamatakis, and P. Pavlidis. OmegaPlus: A Scalable Tool for Rapid
Detection of Selective Sweeps in Whole-Genome Datasets. Bioinformatics, 2012.

Additionally, if you are using the multi-grain parallelized version please cite also:

N. Alachiotis, P. Pavlidis, and A. Stamatakis. Exploiting Multi-grain Parallelism for 
Efficient Selective Sweep Detection. Algorithms and Architectures for Parallel Processing, 56-68, 2012.


Contact
-------

Pavlos Pavlidis - pavlidisp@gmail.com
Nikos Alachiotis - n.alachiotis@gmail.com


This program is free software; you may redistribute it and/or modify its
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.
