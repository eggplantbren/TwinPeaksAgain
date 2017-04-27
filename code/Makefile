CC = g++
CFLAGS = -std=c++14 -O3 -Wall -Wextra -pedantic

default: executable
	rm -rf *.o

library: objects
	ar rcs libtwinpeaks.a *.o

executable: library
	$(CC) $(CFLAGS) -c main.cpp
	$(CC) $(CFLAGS) -o main main.cpp -ltwinpeaks

objects:
	$(CC) $(CFLAGS) -c RNG.cpp
	$(CC) $(CFLAGS) -c Utils.cpp
	$(CC) $(CFLAGS) -c Examples/Atoms.cpp
	$(CC) $(CFLAGS) -c Examples/SimpleExample.cpp
	$(CC) $(CFLAGS) -c Examples/SpikeSlab.cpp
