#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <fstream>

const double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
const double stencil2[3] = {1.0,-2.0,1.0};
const size_t m = 200;
const size_t n = 200;

// Prints the lattice array to the screen
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
void solve(double lattice[m][n],size_t iterations,int order)
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
    //double buffer[m+2][n+2]; 
    //memcpy(buffer,lattice,sizeof(lattice));

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

int main(int argv, char* argc[])
{
    std::cout<<"Hello World"<<std::endl;
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
