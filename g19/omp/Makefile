CC=gcc
CFLAGS=-O2 -fopenmp

SRCS=$(wildcard *.c)
EXEC=life3d-omp

all: $(EXEC)

$(EXEC): $(SRCS)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(EXEC)