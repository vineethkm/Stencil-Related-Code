#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>


#define gpuerrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort=true)
{
	if(code!= cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n",cudaGetErrorString(code),file,line);
		if(abort) exit(code);
	}
}

const int m = 400;
const int n = 400;
__constant__ double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
__constant__ double stencil2[3] = {1.0,-2.0,1.0};
__constant__ int mGpu = m;
__constant__ int nGpu = n;
__constant__ double dt = 0.01;


void printLattice(double lattice[m][n])
{
    for(size_t i = 0;i<m;i++)
    {
        for(size_t j = 0;j<n;j++)
        {
            std::cout<<lattice[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

__global__ void solveIter(double *lattice, int* order)
{
    //extern __shared__ char mem[];
    //double lat[m][n];
    //lat[0] = (double*)&mem[0];
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    int r = blockIdx.y*blockDim.y + threadIdx.y;
    int i = r*nGpu + c;
    int pdg = (*order)/2;
    double sum = 0;
    //printf("hello world");
    //printf("(%d,%d)",r,c);
    
    //lattice[i] = 1.0;
    
    // Store lattice in shared memroy
    //mem[i] = lattice[i];
    //__syncthreads();

    //lattice[i] = 1.0;
    //__syncthreads();
    //
    
    // Code for after getting errors fixed
    if(r<pdg || r>=m-pdg || c<pdg || c>=n-pdg)
    	return;

    if(*order == 2)
    {
        sum = 0;
	for(int h = -pdg; h<=pdg; h++)
	{
	    sum += stencil2[h+pdg]*(lattice[(r+h)*nGpu + c]);
	    sum += stencil2[h+pdg]*(lattice[(r*nGpu)+c+h]);
	}
    } else if(*order == 4)
    {
        sum = 0;
	for(int h = -pdg; h<=pdg; h++)
	{
	    sum += stencil4[h+pdg]*lattice[(r+h)*nGpu+c];
	    sum += stencil4[h+pdg]*lattice[(r*nGpu)+c+h];
	}
    }
    __syncthreads();
    //printf("%lf",sum);
    lattice[i] = dt*sum;
    return;
    
}


/*
__global__ void solveIterShared(double *lattice, int* order)
{
    extern __shared__ char mem[];
    //double lat[m][n];
    //lat[0] = (double*)&mem[0];
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    int r = blockIdx.y*blockDim.y + threadIdx.y;
    int i = r*nGpu + c;
    int pdg = (*order)/2;
    double sum = 0;
    //printf("hello world");
    
    //lattice[i] = 1.0;
    
    // Store lattice in shared memroy
    mem[i] = lattice[i];
    __syncthreads();

    lattice[i] = 1.0;
    __syncthreads();
    //


        // Code for after getting errors fixed
    if(r<pdg || r>=m-pdg || c<pdg || c>=n-pdg)
    	return;

    if(*order == 2)
    {
        sum = 0;
	for(int h = -pdg; h<=pdg; h++)
	{
	    sum += stencil2[h+pdg]*lat[r+h][c];
	    sum += stencil2[h+pdg]*lat[r][c+h];
	}
	lattice[i] += dt*sum;
    } else if(*order == 4)
    {
        sum = 0;
	for(int h = -pdg; h<=pdg; h++)
	{
	    sum += stencil4[h+pdg]*lat[r+h][c];
	    sum += stencil4[h+pdg]*lat[r][c+h];
	}
	lattice[i] += dt*sum;
    }
    return;
    
    
}
*/

// Solves the Unsteady 2D Heat equation using Gauss-seidel method
void solve(double lattice[m][n],size_t iterations,int order)
{
    //double dt = 0.01;
    // Cuda allocations
    // 1d array of lattice
    double *flatLattice, *flatLatticeGpu;
    flatLattice = (double*)malloc(m*n*sizeof(double));
    memcpy(flatLattice, lattice[0], m*n*sizeof(double)); 

    gpuerrchk(cudaMalloc(&flatLatticeGpu, m*n*sizeof(double)));
    gpuerrchk(cudaMemcpy(flatLatticeGpu, flatLattice, m*n*sizeof(double),cudaMemcpyHostToDevice));
    
    // Constants
    int *orderGpu;
    int *orderCpu = (int*)malloc(sizeof(int));
    *orderCpu = order;
    gpuerrchk(cudaMalloc(&orderGpu, sizeof(int)));
    gpuerrchk(cudaMemcpy(orderGpu,orderCpu, sizeof(int), cudaMemcpyHostToDevice));

    // Create dim3 variables for grid and block
    dim3 block(20,20);
    dim3 grid(m/20,n/20);
    for(size_t iter = 0; iter<iterations; iter++)
    {
	
	solveIter<<<grid,block>>>(flatLatticeGpu, orderGpu);
        // Wait until Kernel is done n*m*sizeof(double)
        gpuerrchk(cudaDeviceSynchronize());
	if(iter%1000 == 0)
	    std::cout<<"Iteration: "<<iter<<std::endl;
    }
	//solveIter<<<grid,block>>>(flatLatticeGpu, orderGpu);
	//char* err = cudaGetErrorString(cudaPeekAtLastError);
	//std::cout<<err<<std::endl;
        //cudaDeviceSynchronize();
    gpuerrchk(cudaMemcpy(flatLattice, flatLatticeGpu, m*n*sizeof(double),cudaMemcpyDeviceToHost));
    memcpy(lattice[0], flatLattice, m*n*sizeof(double));
    for(int i =0;i<200;i++)
	    std::cout<<flatLattice[i]<<" ";
    std::cout<<std::endl;
    cudaFree(flatLatticeGpu);
    //cudaFree(&dtGpu);
    cudaFree(orderGpu);
    //free(flatLatticeGpu);
    //free(dtCpu);
    //free(orderCpu);
}

int main(int argv, char** argc)
{
    int order = atoi(argc[1]);
    int iterations = atoi(argc[2]);
    double lattice1[m][n];

    // Set all values to 0
    memset(lattice1,0,sizeof(lattice1));

    // Setting Boundary Conditions
    // Left and right boundaries are 100
    if(order == 2)
    {
        for(size_t i = 0;i<m;i++)
        {
            lattice1[i][0] = 100.0;
            lattice1[i][n-1] = 100.0;
        }
    } else if(order == 4)
    {
        for(size_t i = 0;i<m;i++)
        {
            lattice1[i][0] = 100.0;
            lattice1[i][1] = 100.0;
            lattice1[i][n-1] = 100.0;
            lattice1[i][n-2] = 100.0;
        }
    }else return 0;

    // Up and Down boundaries are 0

    solve(lattice1,iterations,order);
    std::cout<<"Output:"<<std::endl;

    printLattice(lattice1);

    // Write lattice values to file
    std::ofstream f("test.pgm",std::ios_base::out
                              |std::ios_base::binary
                              |std::ios_base::trunc
                   );

    f<<"P2\n"<<m<<" "<<n<<"\n"<<100<<"\n";
    for(size_t i = 0; i<m; i++)
    {
        for(size_t j = 0; j<n; j++)
        {
            f<<(int)lattice1[i][j]<<" ";
        }
    }
    f<<std::flush;
    return 0;
}
