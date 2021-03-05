CFLAGS = -Wall -Werror 
LDLIBS = -lm

.PHONEY: all
all: out.data.txt exp_func.pdf exp_plot.tex

exp_func.pdf: exp_func.tex exp_plot.tex
	pdflatex $<
	pdflatex $<

out.data.txt: main
	./$< > $@

exp_plot.tex: exp.gpi out.data.txt
	gnuplot $<

clean:
	$(RM) ms.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex
