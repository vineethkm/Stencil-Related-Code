#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <fstream>

const double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
const double stencil2[3] = {1.0,-2.0,1.0};
size_t m = 0;
size_t n = 0;

// Prints the lattice array to the screen
void printLattice(double** lattice)
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

// Calculates and Returns Z order
uint32_t calcZOrder(uint16_t xPos, uint16_t yPos)
{
    static const uint32_t MASKS[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
    static const uint32_t SHIFTS[] = {1, 2, 4, 8};

    uint32_t x = xPos;  // Interleave lower 16 bits of x and y, so the bits of x
    uint32_t y = yPos;  // are in the even positions and bits from y in the odd;

    x = (x | (x << SHIFTS[3])) & MASKS[3];
    x = (x | (x << SHIFTS[2])) & MASKS[2];
    x = (x | (x << SHIFTS[1])) & MASKS[1];
    x = (x | (x << SHIFTS[0])) & MASKS[0];

    y = (y | (y << SHIFTS[3])) & MASKS[3];
    y = (y | (y << SHIFTS[2])) & MASKS[2];
    y = (y | (y << SHIFTS[1])) & MASKS[1];
    y = (y | (y << SHIFTS[0])) & MASKS[0];

    const uint32_t result = x | (y << 1);
    return result;
}

// Solves the Steady condition for the 2D Heat Equation using Gauss-Seidel Method
/*
void solveStable(double lattice[m][n],size_t iterations,int order)
{
    double sum = 0;
    //double buffer[m+2][n+2]; 
    //memcpy(buffer,lattice,sizeof(lattice));
    for(size_t iter = 0; iter<iterations; iter++)
    {
        for(size_t i = 1; i<=m; i++)
        {
            for(size_t j = 1; j<=n; j++)
            {
                sum = 0;
                for(int h = -order/2; h<=order/2; h+=2)
                {
                    sum += lattice[i+h][j];
                    sum += lattice[i][j+h];
                }
                lattice[i][j] = sum/4.0;
            }
        }
    }

}
*/

// Solves the Unsteady 2D Heat equation using Gauss-seidel method
// Sequential Traversal
void solve(double **lattice,size_t iterations,int order)
{
    int d = order/2;
    int pdg = d;
    double stencil[5];
    memset(stencil,0,sizeof(stencil));
    double dt = 0.01;
    
    if(order == 2)
        memcpy(stencil,stencil2,sizeof(stencil2));
    else if(order == 4)
        memcpy(stencil,stencil4,sizeof(stencil4));
    else
        return;
    
    for(int i =0;i<5;i++)
        std::cout<<stencil[i]<<" ";
    std::cout<<std::endl;

    double sum = 0;

    for(size_t iter = 0; iter<iterations; iter++)
    {
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
    }
}

// Solves the Unsteady 2D Heat equation using Gauss-seidel method
// Z order Traversal
void solveZOrder(double *lattice,size_t iterations,int order)
{
    int d = order/2;
    int pdg = d;
    double stencil[5];
    memset(stencil,0,sizeof(stencil));
    double dt = 0.01;
    
    if(order == 2)
        memcpy(stencil,stencil2,sizeof(stencil2));
    else if(order == 4)
        memcpy(stencil,stencil4,sizeof(stencil4));
    else
        return;
    
    double sum = 0;

    for(size_t iter = 0; iter<iterations; iter++)
    {
        for(size_t i = pdg; i<m-pdg; i++)
        {
            for(size_t j = pdg; j<n-pdg; j++)
            {
                sum = 0;
                size_t n = calcZOrder(j,i);
                /*
                top    = (((z & 0b10101010) − 1) & 0b10101010) | (z & 0b01010101)
                bottom = (((z | 0b01010101) + 1) & 0b10101010) | (z & 0b01010101)
                left   = (((z & 0b01010101) − 1) & 0b01010101) | (z & 0b10101010)
                right  = (((z | 0b10101010) + 1) & 0b01010101) | (z & 0b10101010)
                */
                for(int h = -d; h<=d; h++)
                {
                    sum += stencil[h+d]*lattice[calcZOrder(j+h,i)];
                    sum += stencil[h+d]*lattice[calcZOrder(j,i+h)];
                }
                lattice[n] += dt*sum;
            }
        }
    }
}

int main(int argv, char* argc[])
{
    int order = atoi(argc[1]);
    m = atoi(argc[2]);
    n = m;
    //double lattice1[m][n];
    double **lattice1 = new double*[m];
    for(size_t i = 0;i<m;i++)
    {
        lattice1[i] = new double[n];
        memset(lattice1[i],0,n*sizeof(double));
    }

    // Set all values to 0
    //memset(lattice1,0,m*n*sizeof(double));
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

    // Load lattice to a 1D array using z ordering
    // Allocate memory for 1D array
    double *latticez = new double[m*n];
    //memset(latticez,0,m*n*sizeof(double));
    for(size_t i =0;i<m;i++)
    {
        for(size_t j=0;j<n;j++)
        {
            latticez[calcZOrder(j,i)] = lattice1[i][j];
        }
    }
    
    solve(lattice1,10000,order);
    solveZOrder(latticez,10000,order);
    std::cout<<"Output:"<<std::endl;
    
    /*
    for(size_t i =0;i<m;i++)
    {
        for(size_t j=0;j<n;j++)
        {
            lattice1[i][j] = latticez[calcZOrder(j,i)];
        }
    }
    */
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
