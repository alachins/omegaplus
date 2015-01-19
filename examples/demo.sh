sourceCodeDir=..

# compile the code
cd $sourceCodeDir
make clean -f Makefile.gcc
make -f Makefile.gcc
cd examples

# example 1
# analyze the simulated data ms.out (ms-like format). Assume that the length of the region is 1000000
# and calculate the omega statistic at 1000 equidistant positions. 
$sourceCodeDir/OmegaPlus -name EXAMPLE1 -input ms.out -grid 1000 -length 1000000 -minwin 1000 -maxwin 20000


#example 2
# analyze the simulated data ms.out (ms-like format). Assume that the length of the region is 1000000
# and calculate the omega statistic at 1000 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
$sourceCodeDir/OmegaPlus -name EXAMPLE2 -input ms.out -grid 1000 -length 1000000 -minwin 1000 -maxwin 20000 -all



#example 3
# analyze the simulated data ms.withgaps.out (ms-like format). Assume that the length of the region is 1000000
# and calculate the omega statistic at 1000 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
# The data contains a few gaps: Ignore them (This is the default)
$sourceCodeDir/OmegaPlus -name EXAMPLE3 -input ms.withgaps.out -grid 1000 -length 1000000 -minwin 1000 -maxwin 20000 -all



#example 4
# analyze the simulated data ms.withgaps.out (ms-like format). Assume that the length of the region is 1000000
# and calculate the omega statistic at 1000 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
# The data contains a few gaps: replace them probabilistically by using the rest of the states at each site. Note that you must use a random seed.
$sourceCodeDir/OmegaPlus -name EXAMPLE4 -input ms.withgaps.out -grid 1000 -length 1000000 -minwin 1000 -maxwin 20000 -all -impute GAP -seed 2947



#example 5
# analyze fasta data (fasta.FA).
# and calculate the omega statistic at 100 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
$sourceCodeDir/OmegaPlus -name EXAMPLE5 -input fasta.FA -grid 100 -minwin 100 -maxwin 1000 -all



#example 6
# analyze fasta data (fasta.FA).
# and calculate the omega statistic at 100 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
# the data contains gaps (-) and undefined characters (N). Replace them probabilistically with the remaining characters at the site
$sourceCodeDir/OmegaPlus -name EXAMPLE6 -input fasta.withgaps.FA -grid 100 -minwin 100 -maxwin 1000 -all -impute GAP -impute N -seed 12983


#example 7
# analyze fasta data (fasta.FA).
# and calculate the omega statistic at 100 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
# the data contains gaps (-) and undefined characters (N). Replace them probabilistically with the remaining characters at the site
# prior to the computations convert the data to binary format. For large datasets this will accelerate the computations considerably
$sourceCodeDir/OmegaPlus -name EXAMPLE7 -input fasta.withgaps.FA -grid 100 -minwin 100 -maxwin 1000 -all -impute GAP -impute N -seed 12983 -binary


#example 8
# analyze fasta data (fasta.FA).
# and calculate the omega statistic at 100 equidistant positions. 
# Report the coordinates of the windows that maximize the omega statistic at each position
# the data contains gaps (-) and undefined characters (N). Replace them probabilistically with the remaining characters at the site
# prior to the computations convert the data to binary format. For large datasets this will accelerate the computations considerably
# Force both windows to have similar number of SNPs (+/- 10).
$sourceCodeDir/OmegaPlus -name EXAMPLE8 -input fasta.withgaps.FA -grid 100 -minwin 100 -maxwin 1000 -all -impute GAP -impute N -seed 12983 -binary -b 10









