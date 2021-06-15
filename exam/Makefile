CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: main

main: main.c

clean:
	$(RM) *.o main out.* *.png
