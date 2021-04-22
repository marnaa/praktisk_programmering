CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.data.txt  compareErrorPlot.png compareValuePlot.png

compareErrorPlot.png: compareErrorPlot.gpi out.compare.txt
	gnuplot $<
compareValuePlot.png: compareValuePlot.gpi out.compare.txt
	gnuplot $<
out.data.txt: main
	./$< > $@

out.compare.txt: main
	./$<

clean:
	$(RM) *.o main out.* *.png

