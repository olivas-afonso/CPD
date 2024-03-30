CC=mpicc
CFLAGS=-O2 -fopenmp -g

SRCS=$(wildcard *.c)
EXEC=life3d-mpi

all: $(EXEC)

$(EXEC): $(SRCS)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(EXEC)