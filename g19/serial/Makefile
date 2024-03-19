CC=gcc
CFLAGS=-O2 -fopenmp

SRCS=$(wildcard *.c)
EXEC=life3d

all: $(EXEC)

$(EXEC): $(SRCS)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(EXEC)