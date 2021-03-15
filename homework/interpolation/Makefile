CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all
all: out.10.txt out.100.txt out.500.txt out.xy.txt intSin10Point.png \
intSin100Point.png intSin500Point.png out.qmaintest.txt out.cmaintest.txt

intSin10Point.png: sin10.gpi out.10.txt out.xy.txt
	gnuplot $<

intSin100Point.png: sin100.gpi out.100.txt out.xy.txt
	gnuplot $<

intSin500Point.png: sin500.gpi out.500.txt out.xy.txt
	gnuplot $<

out.10.txt out.100.txt out.500.txt out.xy.txt: main Makefile
	./$<

out.qmaintest.txt: qmain
	./$< > $@
out.cmaintest.txt: cmain
	./$< > $@

.PHONEY: clean
clean:
	$(RM) main *.o out* qmain *.png cmain
