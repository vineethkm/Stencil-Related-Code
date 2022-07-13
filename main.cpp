#include<iostream>
#include<array>
#include "lattice.h"

using namespace lattice;
int main()
{
    int n = 10;
    auto itr = lattice_iterator<3>({{0,0,0}},{{n,n,n}});
    lattice_iterator<3> itr_end;
    for(;itr != itr_end; itr++)
    {
        std::cout<<(*itr)[0]<<(*itr)[1]<<(*itr)[2]<<std::endl;
    }

    return 0;
}
