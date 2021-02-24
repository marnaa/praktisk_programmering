CFLGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: all
all: bessel.gnuplot.png erf.gnuplot.png out.valueExerciseA.txt out.erfValues.txt \
	out.besselValues.txt

erf.gnuplot.png: out.erfValues.txt Makefile
	echo '\
	set terminal png; \
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set yrange [-1.1:1.1];\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error Function";\
	plot "$<" using 1:2 with line title "Calculated",\
	"tabulated_data.txt" using 1:2 with points title "Data"\
	' | gnuplot
bessel.gnuplot.png: out.besselValues.txt 
	echo '\
	set terminal png; \
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Bessel_n Function";\
	plot "$<" using 1:2 with line title "Bes0",\
	"$<" using 1:3 with line title "Bes1",\
	"$<" using 1:4 with line title "Bes2",\
	"besselRoots.txt" using 2:1 with points title "Root0",\
	"besselRoots.txt" using 3:1 with points title "Root1",\
	"besselRoots.txt" using 4:1 with points title "Root2"\
	' | gnuplot


out.valueExerciseA.txt out.erfValues.txt out.besselValues.txt: main Makefile
	./$<

.PHONEY: clean
clean:
	$(RM) main out* *.png
