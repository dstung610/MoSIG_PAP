all: triangular_matrix bubble qsort mergesort

mergesort: mergesort.o
			gcc-8 -fopenmp -o mergesort mergesort.o

mergesort.o: mergesort.c
			gcc-8 -c -O0 -fopenmp mergesort.c

qsort: qsort.o
		gcc-8 -fopenmp -o qsort qsort.o

qsort.o: qsort.c
		gcc-8 -c -O0 -fopenmp qsort.c

bubble: bubble.o
	gcc-8 -fopenmp -o bubble bubble.o

bubble.o: bubble.c
	gcc-8 -c -O0 -fopenmp bubble.c

triangular_matrix: triangular_matrix.o
	gcc-8 -fopenmp -o triangular_matrix triangular_matrix.o

triangular_matrix.o: triangular_matrix.c
	gcc-8 -c -O0 -fopenmp triangular_matrix.c

clean:
	rm -f *~ *.o bubble qsort triangular_matrix mergesort
