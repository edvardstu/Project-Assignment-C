

CC=gcc
FLAGS := -O3
LINKING := -lblas -fopenmp -lgsl -lgslcblas -lm
DEPS = utilities.h interactions.h
OBJ = main.o utilities.o interactions.o


.PHONY: run clean

%.o: %.c $(DEPS)
	$(CC) $(FLAGS) -c $< -o $@

main: $(OBJ)
	$(CC) $(FLAGS) -o $@ $^ $(LINKING)

run: main
	./main

clean:
	rm -f main *.o
