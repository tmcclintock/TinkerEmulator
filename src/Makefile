GSL_INC=/afs/slac.stanford.edu/g/ki/software/gsl/1.15/include/
GSL_LIB=/share/sw/free/gsl/2.1/lib/
test:
	gcc -std=c99 test.c -o ./test.out -fPIC -lm -lgsl -lfftw3 -lgsl -lgslcblas -O3 -ffast-math -funroll-loops -std=c99 -L../../../cosmolike_core/class -lclass
