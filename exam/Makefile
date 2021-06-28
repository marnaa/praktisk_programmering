CFLAGS =  -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: out.matrix_test.txt out.eigval_check.txt

out.matrix_test.txt: main
	@./$< > $@

out.eigval_check.txt: main

main: main.c lanczos.o eigval_routines.o

lanczos.o: lanczos.c

eigval_routines.o: eigval_routines.c

clean:
	$(RM) *.o main out.* *.png
