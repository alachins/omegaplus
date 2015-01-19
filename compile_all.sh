rm *.o
make -f Makefile.gcc
rm *.o
make -f Makefile.PTHREADS.FINE.gcc
rm *.o
make -f Makefile.PTHREADS.COARSE.gcc
rm *.o
make -f ./Makefile.PTHREADS.MULTI.gcc
rm *.o;
rm *.o
make -f ./Makefile.OPENMP.GENERIC.gcc
rm *.o;


