CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: out.RosenbockAndE.txt psi0Plot.png convPlot.png convBoundPlot.png

psi0Plot.png: psi0Plot.gpi out.psi0.txt
	gnuplot $<


convPlot.png: convPlot.gpi out.conv.txt
	gnuplot $<

convBoundPlot.png: convBoundPlot.gpi out.convBound.txt
	gnuplot $<

out.psi0.txt out.conv.txt: main
	./$<

out.RosenbockAndE.txt: main
	./$< > $@

main: main.c gs_routines.o ode_routines.o

clean:
	$(RM) *.o out.* main *.png
