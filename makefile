
all: main.cpp
	g++ -pg -Wall main.cpp lattice.h lattice_iterator.h range.h -o main

clean:
	rm main
