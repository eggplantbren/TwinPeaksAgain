CC = g++
CFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic

default: library executable
	rm -rf *.o

library: objects
	ar rcs libtwinpeaks.a *.o

executable: library
	$(CC) $(CFLAGS) -c main.cpp
	$(CC) -o main main.cpp -ltwinpeaks

objects: RNG.o
	$(CC) $(CFLAGS) -c RNG.cpp

