CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: infWellPlot.png timePlot.png

infWellPlot.png: infWellPlot.gpi out.infWell.txt
	gnuplot $<
timePlot.png: timePlot.gpi out.compare.txt
	gnuplot $<

out.compare.txt out.infWell.txt: main
	./$<
clean:
	$(RM) *.o main *.png out.*
