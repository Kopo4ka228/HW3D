#include "geom.hpp"
#include <iostream>
#include <fstream>
#include <random>

int main()
{
    const int N = 10000;
    std::ofstream out;
    out.open("input.txt");
    out << N << std::endl;

    for(int i=0; i < N; i++)
    {
        out << rand() <<' '<< rand() << ' ' << rand() << std::endl;
    }
    out.close();
}