
.PHONY: all clean run check

all: run check

run: test.c
	gcc -O0 -fopenmp -o omp_strats omp_strats.c
	./omp_strats $(N) $(INP) $(P) $(S)

check:
	python3 format_checker.py $(INP) output_L_$(S)_$(P).txt output_U_$(S)_$(P).txt

clean:
	rm omp_strats

