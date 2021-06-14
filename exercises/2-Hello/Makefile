CC=gcc
default: out.txt 
	cat out.txt
out.txt: hello
	./hello > out.txt 
hello: hello.o
	$(CC) -o hello hello.o
hello.o: hello.c
	$(CC) -c hello.c 
.PHONEY: clean
clean:                              # clean is a phoney target
	$(RM) hello.o hello out.txt # clean up the directory
.PHONEY: test
test:                      # test target used for debugging
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
