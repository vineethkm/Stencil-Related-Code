
all: src/heat_eq.cpp
	mkdir build
	g++ -pg -Wall -O3 src/heat_eq.cpp -o build/heat_equation

cuda: src/heat_eq.cu
	mkdir build
	nvcc src/heat_eq.cu -o build/heat_equation

clean:
	rm -r build
