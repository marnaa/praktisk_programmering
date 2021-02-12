CC= gcc
CFLAGS= -O
LDLIBS=-lm 
out.txt: math 
	./math > out.txt
math.o: math.c 
	$(CC) $(CFLAGS) -c math.c -o math.o 
math: math.o 
	$(CC) $(CFLAGS) $^ -o  $@  
clean:
	$(RM) *.o out* math 
