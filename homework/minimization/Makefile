CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: out.exA.txt exBPlot.png

exBPlot.png: exBPlot.gpi out.exB.txt
	gnuplot $<

out.exA.txt: main
	./$< > $@

out.exB.txt: main
	./$<

clean:
	$(RM) *.o main out.* *.png
