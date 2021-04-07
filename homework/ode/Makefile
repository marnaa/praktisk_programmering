CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: harmOscillator.png SIRplot.png

harmOscillator.png: harmplot.gpi out.harmosc.txt
	gnuplot $<

SIRplot.png: SIRplot.gpi out.disease.txt
	gnuplot $<

out.harmosc.txt: main
	./$<

out.disease.txt: main
	./$< > $@

clean:
	$(RM) *.o main out.* *.png  
