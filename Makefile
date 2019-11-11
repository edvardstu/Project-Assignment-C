

CC=gcc
FLAGS := -O3
LINKING := -lblas -fopenmp -lgsl -lgslcblas -lm
DEPS = utilities.h
OBJ = main.o utilities.o


.PHONY: run clean

%.o: %.c $(DEPS)
	$(CC) $(FLAGS) -c $< -o $@ 

main: $(OBJ)
	$(CC) $(FLAGS) -o $@ $^ $(LINKING)

run: main
	./main

clean:
	rm -f main *.o
