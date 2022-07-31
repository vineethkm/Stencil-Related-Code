
all: src/heat_eq.cpp
	mkdir build
	g++ -pg -Wall src/heat_eq.cpp -o build/heat_equation

fluideq: src/fluid_equation.cpp
	mkdir build
	g++ -pg -Wall src/fluid_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/fluid_equation

clean:
	rm -r build
