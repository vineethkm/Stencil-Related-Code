#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <fstream>

__constant__ double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
__constant__ double stencil2[3] = {1.0,-2.0,1.0};
__constant__ size_t mGpu = 200;
__constant__ size_t nGpu = 200;
__constant__ double dt = 0.01;
const size_t m = 200;
const size_t n = 200;

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
    extern __shared__ double mem[];
    double **lat;
    *lat = mem;
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    int r = blockIdx.y*blockDim.y + threadIdx.y;
    int i = r*nGpu + c;
    int pdg = (*order)/2;
    double sum = 0;
    printf("hello world");
    
    lattice[i] = 1.0;
    /*
    // Store lattice in shared memroy
    lat[r][c] = lattice[i];
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
    */
}

// Solves the Unsteady 2D Heat equation using Gauss-seidel method
void solve(double lattice[m][n],size_t iterations,int order)
{
    //double dt = 0.01;
    // Cuda allocations
    // 1d array of lattice
    double *flatLattice, *flatLatticeGpu;
    flatLattice = (double*)malloc(m*n*sizeof(double));
    memcpy(flatLattice, lattice[0], m*n*sizeof(double)); 

    cudaMalloc(&flatLatticeGpu, m*n*sizeof(double));
    cudaMemcpy(flatLatticeGpu, flatLattice, m*n*sizeof(double),cudaMemcpyHostToDevice);
    
    // Constants
    int *orderGpu;
    int *orderCpu = (int*)malloc(sizeof(int));
    *orderCpu = order;
    cudaMalloc(&orderGpu, sizeof(int));
    cudaMemcpy(orderGpu,orderCpu, sizeof(int), cudaMemcpyHostToDevice);

    // Create dim3 variables for grid and block
    dim3 block(20,20);
    dim3 grid(m/20,n/20);
    for(size_t iter = 0; iter<iterations; iter++)
    {
	
	solveIter<<<grid,block,n*m*sizeof(double)>>>(flatLatticeGpu, orderGpu);
        // Wait until Kernel is done
        cudaDeviceSynchronize();
    }
	//solveIter<<<grid,block,n*m*sizeof(double)>>>(flatLatticeGpu, orderGpu);
	//char* err = cudaGetErrorString(cudaPeekAtLastError);
	std::cout<<err<<std::endl;
        cudaDeviceSynchronize();
    cudaMemcpy(flatLattice, flatLatticeGpu, m*n*sizeof(double),cudaMemcpyDeviceToHost);
    memcpy(lattice[0], flatLattice, m*n*sizeof(double));
    for(int i =0;i<200;i++)
	    std::cout<<flatLattice[i]<<" ";
    std::cout<<std::endl;
    cudaFree(flatLatticeGpu);
    //cudaFree(&dtGpu);
    cudaFree(&orderGpu);
    //free(flatLatticeGpu);
    //free(dtCpu);
    //free(orderCpu);
}

int main(int argv, char** argc)
{
    int order = atoi(argc[1]);
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

    solve(lattice1,10000,order);
    std::cout<<"Output:"<<std::endl;

    //printLattice(lattice1);

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
