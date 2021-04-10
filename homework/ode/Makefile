CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: harmOscillator.png SIRplot.png SIRplotTc3.png threebody.png

threebody.png: body3plot.gpi out.3body.txt
	gnuplot $<

harmOscillator.png: harmplot.gpi out.harmosc.txt
	gnuplot $<

SIRplot.png: SIRplot.gpi out.disease.txt
	gnuplot $<

SIRplotTc3.png: SIRplotTc3.gpi out.disease3.txt
	gnuplot $<

out.harmosc.txt out.disease.txt out.disease3.txt out.3body.txt: main
	./$<


clean:
	$(RM) *.o main out.* *.png
