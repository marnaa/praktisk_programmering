CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: erf.gnuplot.png gam.gnuplot.png main

gam.gnuplot.png: data.txt Makefile 
	echo '\
	set terminal png; \
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot "$<" using 1:5 with line title "tgam from math.h",\
	"$<" using 1:6 with line title "gslsfgam",\
	"$<" using 1:7 with line title "Gamma"\
	,"gam_tabulated_data.txt" using 1:2 with points title "tabulated data"\
	' | gnuplot

erf.gnuplot.png: data.txt tabulated_data.txt
	echo '\
	set terminal png; \
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Errorfunction";\
	plot "$<" using 1:2 with line title "erf from math.h",\
	"$<" using 1:3 with line title "gslsferr",\
	"$<" using 1:4 with line title "Homemade erf"\
	,"tabulated_data.txt" using 1:2 with points title "tabulated data"\
	' | gnuplot

data.txt: main
	./$< > $(@)

main: main.o erf.o gam.o

.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt test


