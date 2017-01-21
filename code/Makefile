CC = g++
CFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic

default: library

library: objects
	ar rcs libtwinpeaks.a *.o

objects: RNG.o
	$(CC) $(CFLAGS) -c RNG.cpp

