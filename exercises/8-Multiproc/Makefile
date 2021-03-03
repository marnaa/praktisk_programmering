CFLAGS = -Wall -Ofast
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp

.PHONEY: all

all: out.txt omp.out.txt

out.txt: main
	./$< > $@


omp.out.txt: ompmain
	./$< > $@

clean:
	$(RM) ompmain main *.txt
