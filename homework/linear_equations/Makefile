CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.gstest.txt timeplot.png compareplot.png gslplot.png


out.gstest.txt out.compare.txt out.time.txt : main
	./$<

timeplot.png: qplot.gpi out.time.txt
	gnuplot $<

compareplot.png: compare.gpi out.compare.txt
	gnuplot $<

gslplot.png: gsltime.gpi out.compare.txt
	gnuplot $<


.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png
