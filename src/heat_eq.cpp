#include <iostream>
#include <cstring>
#include <fstream>

//const double stencil4[5] = {-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0};
//const double stencil2[3] = {1.0,-2.0,1.0};
const size_t m = 200;
const size_t n = 200;

void printLattice(double lattice[m+2][n+2])
{
    for(size_t i = 0;i<m+2;i++)
    {
        for(size_t j = 0;j<n+2;j++)
        {
            std::cout<<lattice[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void solveStable(double lattice[m+2][n+2],size_t iterations,int order)
{
    int d = order/2;
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

void solve(double lattice[m+2][n+2],size_t iterations,int order)
{
    int d = order/2;
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

int main()
{
    std::cout<<"Hello World"<<std::endl;
    double lattice1[m+2][n+2];
    double lattice2[m+2][n+2];
    std::cout<<sizeof(lattice1)<<std::endl;

    // Set all values to 0
    memset(lattice1,0,sizeof(lattice1));

    // Setting Boundary Conditions
    // Left and right boundaries are 100
    for(size_t i = 0;i<m+2;i++)
    {
        lattice1[i][0]= 100.0;
        lattice1[i][n+1] = 100.0;
    }

    // Up and Down boundaries are 0

    std::cout<<sizeof(lattice1)<<std::endl;
    // Copy lattice1 to lattice2
    memcpy(lattice2,lattice1,sizeof(lattice1));

    std::cout<<sizeof(lattice1)<<std::endl;
    //double** l = reinterpret_cast<double**> lattice1;
    //printLattice(lattice1);

    //solveStable(lattice1,10000,2);
    //std::cout<<"Output:"<<std::endl;

    //printLattice(lattice1);

    // Write lattice values to file
    std::ofstream f("test.pgm",std::ios_base::out
                              |std::ios_base::binary
                              |std::ios_base::trunc
                   );

    f<<"P2\n"<<m+2<<" "<<n+2<<"\n"<<100<<"\n";
    for(size_t i = 0; i<m+2; i++)
    {
        for(size_t j = 0; j<n+2; j++)
        {
            f<<(int)lattice1[i][j]<<" ";
        }
    }
    f<<std::flush;
    return 0;
}
