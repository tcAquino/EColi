CC = g++
PTH = ../include
PBIN = ../bin
INC = -I$(PTH)
OPT = -O3 -DNDEBUG
CFLAGS = -Wall -std=c++17 $(OPT)

EColi : EColi.o
	$(CC) $(CFLAGS) $(LIB) -o EColi EColi.o
	rm EColi.o
	mv EColi $(PBIN)/EColi

EColi.o : EColi.cpp
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $<

clean :
	rm -f EColi EColi.o
