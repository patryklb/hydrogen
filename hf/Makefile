CC=gcc
CFLAGS=-I.
DEPS = mlib.h
OBJ = main.o mlib.o
LIBS=-lm -lstdc++ -llapack -lblas 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
clean:
	rm -f $(OBJ)
