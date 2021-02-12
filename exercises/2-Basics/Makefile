CC = gcc 
CFLAGS = -O -Wall
LDLIBS = -lm 

out.txt: hello
	./hello > out.txt 

hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c -o hello.o 

hello: hello.o
	$(CC) $(LDFLAGS) hello.o -o hello $(LDLIBS)

