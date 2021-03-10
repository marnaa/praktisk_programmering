CFLAGS = -Wall -Werror -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean
clean:
	$(RM) main *.o 
