CC=gcc
CFLAGS=-O2

SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)
EXEC=my_executable

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
