
all: src/heat_eq.cpp build
	g++ -pg -O -Wall src/heat_eq.cpp -o build/heat_equation

cuda: src/heat_eq.cu build
	nvcc -O3 src/heat_eq.cu -o build/heat_equation

build:
	mkdir build

clean:
	rm -r build
