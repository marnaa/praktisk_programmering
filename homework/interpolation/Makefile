CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all
all: out.10.txt out.100.txt out.500.txt out.xy.txt out.qdata.txt out.qxy.txt intSin10Point.png \
intSin100Point.png intSin500Point.png qplot.png cplot.png out.cdata.txt out.cxy.txt

intSin10Point.png: sin10.gpi out.10.txt out.xy.txt
	gnuplot $<

intSin100Point.png: sin100.gpi out.100.txt out.xy.txt
	gnuplot $<

intSin500Point.png: sin500.gpi out.500.txt out.xy.txt
	gnuplot $<

qplot.png: qplot.gpi out.qdata.txt out.qxy.txt
	gnuplot $<

cplot.png: cplot.gpi out.cdata.txt out.cxy.txt
	gnuplot $<

out.10.txt out.100.txt out.500.txt out.xy.txt: main Makefile
	./$<

out.qdata.txt out.qxy.txt: qmain
	./$<
out.cdata.txt out.cxy.txt: cmain
	./$<

out.cmaintest.txt: cmain
	./$< > $@

.PHONEY: clean
clean:
	$(RM) main *.o out* qmain *.png cmain
