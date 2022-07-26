
all: src/main.cpp
	mkdir build
	g++ -pg -Wall src/main.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/main
	g++ -pg -Wall src/heat_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/heat_equation
	g++ -pg -Wall src/fluid_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/fluid_equation

heatequation: src/heat_equation.cpp
	mkdir build
	g++ -pg -Wall src/heat_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/heat_equation

heateq: src/heat_eq.cpp
	mkdir build
	g++ -pg -Wall src/heat_eq.cpp -o build/heat_equation


fluideq: src/fluid_equation.cpp
	mkdir build
	g++ -pg -Wall src/fluid_equation.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/fluid_equation

main: src/main.cpp
	mkdir build
	g++ -pg -Wall src/main.cpp src/lattice.h src/lattice_iterator.h src/range.h -o build/main

clean:
	rm -r build
