
.PHONY: all clean run check

all: run check

run: test.c
	bash compile.sh
	bash run.sh $(N) $(INP) $(P) $(S)

check:
	python3 format_checker.py $(INP) output_L_$(S)_$(P).txt output_U_$(S)_$(P).txt

clean:
	rm -f omp_strats mpi_strat

