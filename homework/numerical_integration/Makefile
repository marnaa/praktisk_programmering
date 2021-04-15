CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.stdout.txt piComparePlot.png evalPlot.png infCompareValuePlot.png infCompareEvalPlot.png infCompareErrorPlot.png

piComparePlot.png: piComparePlot.gpi out.piCompare.txt
	gnuplot $<
infCompareValuePlot.png: infCompareValuePlot.gpi out.infCompare.txt
	gnuplot $<
infCompareEvalPlot.png: infCompareEvalPlot.gpi out.infCompare.txt
	gnuplot $<
infCompareErrorPlot.png: infCompareErrorPlot.gpi out.infCompare.txt
	gnuplot $<

evalPlot.png: evalPlot.gpi out.piCompare.txt
	gnuplot $<

out.piCompare.txt out.infCompare.txt: main
	./$<

out.stdout.txt: main
	./$< > $@

clean:
	$(RM) *.o main out.* *.png
