CC=gcc
CFLAGS+= -Wall -fwrapv -O
LDLIBS=-lm 

.PHONEY: all clean

all: out.txt
	cat $<

out.txt: epsilon Makefile main
	./epsilon > $@
	./main >> $@

main: main.o komplex.o 

main.o komplex.o: komplex.h

epsilon: epsilon.o equals.o

epsilon.o: epsilon.c

equal.o: equals.c 

clean:
	$(RM)  epsilon out.txt main main.o komplex.o epsilon.o equals.o
