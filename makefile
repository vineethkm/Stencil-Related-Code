
all: src/main.cpp
	mkdir build
	g++ -pg -Wall src/main.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/main
	g++ -pg -Wall src/heat_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/heat_equation

heateq: src/heat_equation.cpp
	mkdir build
	g++ -pg -Wall src/heat_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/heat_equation

main: src/heat_equation.cpp
	mkdir build
	g++ -pg -Wall src/main.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/main

clean:
	rm -r build
