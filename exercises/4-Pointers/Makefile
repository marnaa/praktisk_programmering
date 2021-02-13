CC=gcc
CFLAGS+= -Wall
LDLIBS=-lm 

out.txt: epsilon Makefile
	./$< > $@
epsilon: epsilon.o equals.o

epsilon.o: epsilon.c

equal.o: equals.c 

clean:
	$(RM)  epsilon out.txt
