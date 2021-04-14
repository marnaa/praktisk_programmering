CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.stdout.txt

out.stdout.txt: main
	./$< > $@

clean:
	$(RM) *.o main out.* *.png
