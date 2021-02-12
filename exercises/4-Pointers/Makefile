CC=gcc
CFLAGS+= -Wall
LDLIBS=-lm 

out.txt: epsilon Makefile
	./$< > $@
epsilon: epsilon.c equal.c

clean:
	$(RM)  epsilon out.txt
