#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <fstream>

__constant__ double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
__constant__ double stencil2[3] = {1.0,-2.0,1.0};
__constant__ size_t mGpu = 200;
__constant__ size_t nGpu = 200;
double dt = 0.01;
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

__global__ void solveIter(double *lattice, double* dt, int* order)
{
    extern __shared__ double mem[];
    double **lat;
    //lat[0] = mem;
    int c = blockIdx.x*blockDim.x + threadIdx.x;
    int r = blockIdx.y*blockDim.y + threadIdx.y;
    int i = r*n + c;
    /*
    // Store lattice in shared memroy
    lat[r][c] = lattice[i];
    __syncthreads();

    lattice[i] = lat[r][c];
    // Make use of shared memory to act as a cache for threads within the block
    ///*
    int order = 2;

    double sum = 0;
    for(int h = -order/2; h<=order/2; h+=2)
    {
        sum += buffer[x+h][y];
        sum += buffer[x][y+h];
    }
    lattice[x][y] = sum/4.0;
    */
}

// Solves the Unsteady 2D Heat equation using Gauss-seidel method
void solve(double lattice[m][n],size_t iterations,int order)
{
    int d = order/2;
    int pdg = d;
    double stencil[5];
    memset(stencil,0,sizeof(stencil));
    //double dt = 0.01;
    
    for(int i =0;i<5;i++)
        std::cout<<stencil[i]<<" ";
    std::cout<<std::endl;
    std::cout<<"pe"<<std::endl;
    // Cuda allocations
    // 1d array of lattice
    double *flatLattice, *flatLatticeGpu;
    flatLattice = (double*)malloc(m*n*sizeof(double));
    memcpy(flatLattice, lattice[0], m*n*sizeof(double)); 

    std::cout<<"pe"<<std::endl;
    cudaMalloc(&flatLatticeGpu, m*n*sizeof(double));
    cudaMemcpy(flatLatticeGpu, flatLattice, m*n*sizeof(double),cudaMemcpyHostToDevice);
    
    std::cout<<"pe"<<std::endl;

    // Constants
    double *dtGpu;
    double *dtCpu = (double*)malloc(sizeof(double));
    *dtCpu = dt;
    int *orderGpu;
    int *orderCpu = (int*)malloc(sizeof(int));
    *orderCpu = order;
    cudaMalloc(&dtGpu, sizeof(double));
    cudaMalloc(&orderGpu, sizeof(int));
    cudaMemcpy(dtGpu,dtGpu, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(orderGpu,orderCpu, sizeof(int), cudaMemcpyHostToDevice);

    // Create dim3 variables for grid and block
    dim3 block(20,20);
    dim3 grid(m/20,n/20);

    std::cout<<"pe"<<std::endl;
    double sum = 0;
    for(size_t iter = 0; iter<iterations; iter++)
    {
	/*
        for(size_t i = pdg; i<m-pdg; i++)
        {
            for(size_t j = pdg; j<n-pdg; j++)
            {
                sum = 0;
                for(int h = -d; h<=d; h++)
                {
                    sum += stencil[h+d]*lattice[i+h][j];
                    sum += stencil[h+d]*lattice[i][j+h];
                }
                lattice[i][j] += dt*sum;
            }
        }
	*/
    std::cout<<"pe"<<std::endl;
	solveIter<<<grid,block,n*m*sizeof(double)>>>(flatLatticeGpu, dtGpu, orderGpu);
        // Wait until Kernel is done
        cudaDeviceSynchronize();
    }

    cudaMemcpy(flatLattice, flatLatticeGpu, m*n*sizeof(double),cudaMemcpyDeviceToHost);
    memcpy(lattice[0], flatLattice, m*n*sizeof(double));
    cudaFree(flatLatticeGpu);
    cudaFree(&dtGpu);
    cudaFree(&orderGpu);
    free(flatLatticeGpu);
    free(dtCpu);
    free(orderCpu);
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
