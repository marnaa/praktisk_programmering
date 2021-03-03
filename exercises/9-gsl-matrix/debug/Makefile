CFLAGS = -O -w -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += $(shell gsl-config --libs)
.PHONEY: all
all: output; cat $<
output: main ; ./$< > $@
clean: ; $(RM) main output
