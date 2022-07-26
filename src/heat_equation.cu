#include<iostream>
#include "lattice.h"
#include<array>
#include<algorithm>
#include<vector>
#include<fstream>

using namespace lattice;

template <int O>
double stencil(const int i) {
    const std::array<double,O+1> coeff = {{1.0,-2.0,1.0}};
    return coeff[i+O/2];
}

template <>
double stencil<4>(const int i) {
    const std::array<double,4+1> coeff = {{-1.0/12.0,4.0/3.0,-5.0/2.0,4.0/3.0,-1.0/12.0}};
    return coeff[i+4/2];
}

__global__ void solve(auto domain, std::vector<double> values1,std::vector<double> values0,)
{
    for(const auto& index: domain) {
            const size_t base_index = calculate(index);
            values1[base_index] = values0[base_index];
            for (unsigned int d = 0; d < D; d++) {
                int_d stencil_index = index;
                for (int j = -order/2; j <= order/2; j++) {
                    stencil_index[d] = j+index[d];
                    const double coeff = stencil<order>(j);
                    values1[base_index] += r*coeff*values0[calculate(stencil_index)];
                }
            }
        }
}

int main()
{
    std::cout<<"This is code to simulate 2d heat equation"<<std::endl;
    const unsigned int D = 2;
    typedef std::array<int,D> int_d;
    const size_t m = 50;
    const size_t n = 50;
    const int order = 4;
    const double h = 1.0/m;
    const double dt = 0.9*std::pow(h,2)/(-stencil<order>(0)*D);
    const double r = dt/std::pow(h,2);
    const double t = 1.0;
    const int timesteps = t/dt;
    std::cout << "timesteps = "<<timesteps<<" r = "<<r<<std::endl;

    const int_d min = {{0,0}};
    const int_d max = {{n+order,m+order}};
    const int_d min_domain = {{order/2,order/2}};
    const int_d max_domain = {{n+order/2,m+order/2}};
    const int_d max_left = {{order/2,m+order}};
    const int_d min_right = {{n+order/2,0}};

#if __cplusplus > 201402L
    auto all = make_iterator_range(
                        lattice_iterator<D>(min,max),
                        false);
    auto domain = make_iterator_range(
                        lattice_iterator<D>(min_domain,max_domain),
                        false);
    auto left_boundary = make_iterator_range(
                        lattice_iterator<D>(min,max_left),
                        false);
    auto right_boundary = make_iterator_range(
                        lattice_iterator<D>(min_right,max),
                        false);
#else
    auto all = make_iterator_range(
                        lattice_iterator<D>(min,max),
                        lattice_iterator<D>());
    auto domain = make_iterator_range(
                        lattice_iterator<D>(min_domain,max_domain),
                        lattice_iterator<D>());
    auto left_boundary = make_iterator_range(
                        lattice_iterator<D>(min,max_left),
                        lattice_iterator<D>());
    auto right_boundary = make_iterator_range(
                        lattice_iterator<D>(min_right,max),
                        lattice_iterator<D>());
#endif

    std::vector<double> values0(all.size(),0.0);
    std::vector<double> values1(all.size(),0.0);

    auto calculate = [&](const int_d& index) {
        int ret = 0;
        unsigned int multiplier = 1.0;
        for (size_t i = 0; i<D; ++i) {
            if (i > 0) {
                multiplier *= max[i-1];
            }
            ret += multiplier*index[i];
        }
        return ret;
    };


    for(const auto& index: left_boundary) {
        values0[calculate(index)] = 1.0;
        values1[calculate(index)] = 1.0;
    }
    for(const auto& index: right_boundary) {
        values0[calculate(index)] = 1.0;
        values1[calculate(index)] = 1.0;
    }

    for (int i = 0; i < timesteps; ++i) {
        for(const auto& index: domain) {
            const size_t base_index = calculate(index);
            values1[base_index] = values0[base_index];
            for (unsigned int d = 0; d < D; d++) {
                int_d stencil_index = index;
                for (int j = -order/2; j <= order/2; j++) {
                    stencil_index[d] = j+index[d];
                    const double coeff = stencil<order>(j);
                    values1[base_index] += r*coeff*values0[calculate(stencil_index)];
                }
            }
        }
        values1.swap(values0);
    }



    /*
     * Output an image
     */

    std::ofstream f("test.pgm",std::ios_base::out
                              |std::ios_base::binary
                              |std::ios_base::trunc
                   );

    int maxColorValue = 255;

    std::vector<unsigned char> image(all.size());
    auto minmax = std::minmax_element(values0.begin(),values0.end());
    std::transform(values0.begin(),values0.end(),image.begin(),[&](const double value) {
            uint8_t val = maxColorValue*(value-*std::get<0>(minmax))
                    /(*std::get<1>(minmax)-*std::get<0>(minmax));
            return val;

        });

    f << "P5\n" << (max[1]-min[1]) << " " << (max[0]-min[0]) << "\n" << maxColorValue << "\n";

    for(int i=0; i<max[0]; ++i) {
        const size_t index = i*(max[1]-min[1]);
        f.write( reinterpret_cast<const char*>(&image[index]), (max[1]-min[1]) );
    }



    f << std::flush;
    return 0;
}