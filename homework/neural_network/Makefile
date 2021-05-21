CFLAGS = -Wall -O -g
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: comparePlot.png diffeqPlot.png out.data.txt

comparePlot.png: comparePlot.gpi out.compare.txt
	gnuplot $<
diffeqPlot.png: diffeqPlot.gpi out.diffeq.txt
	gnuplot $<

out.data.txt: main
	./$< > $@
out.compare.txt: main
	./$<

out.diffeq.txt: main
	./$<

main: main.c ann.o minimization_routines.o

ann.o: ann.c minimization_routines.o


clean:
	$(RM) *.o main out.* *.png
