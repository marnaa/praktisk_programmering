CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: fitplot.png out.ckandHalflife.txt

fitplot.png: fitanddata.gpi out.fitdat.txt out.fitfunc.txt
	gnuplot $<

out.fitdat.txt out.fitfunc.txt out.ckandHalflife.txt: main
	./$<

main: main.c gs_routines.o

clean:
	$(RM) *.o main *.png out.*


