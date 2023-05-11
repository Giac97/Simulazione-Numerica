#include "genetic.h"
#include <iostream>

int main(int argc, char const *argv[])
{

    Chromosome ch(34, 0);
    std::cout <<"L^2 = "<< ch.length2() << std::endl;
    std::cout <<"L^1 = "<< ch.length() << std::endl;


    return 0;
}

