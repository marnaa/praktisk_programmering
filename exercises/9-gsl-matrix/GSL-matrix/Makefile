CFLAGS = -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

out.txt: main
	./$< > $@ 

.PHONEY: clean
clean:
	$(RM) *.o *.txt
