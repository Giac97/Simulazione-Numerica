#include "genetic.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <ostream>

City::City(int id, double x, double y)
{
    initCity(id, x, y);
}

void City::initCity(int id, double x, double y)
{
    m_id = id;
    m_x = x;
    m_y = y;
}

double City::dist2(City other)
{
    double dx = this->m_x - other.m_x;
    double dy = this->m_y - other.m_y;

    double d2 = dx * dx + dy * dy;

    return d2;

}

double City::dist(City other)
{
    double d2 = dist2(other);
    return sqrt(d2);
}

double City::getX()
{
    return m_x;
}

Chromosome::Chromosome(int ncities, int method)
{
    initChromosome(ncities, method);
}

void Chromosome::initChromosome(int ncities, int method)
{

    m_ncities = ncities;

    double radius = 10.;

    int seed[4];
    Random rnd;
    int p1, p2;
    std::ifstream  Primes, Seed;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();

    for (int i = 0; i < ncities; i++)
    {
        if (method == 0)
        {
            double Theta = rnd.Rannyu(0, 2. * M_PI);
            double x = radius * cos(Theta);
            double y = radius * sin(Theta);

            City newCity = City(i, x, y);
            m_cities.push_back(newCity);
        }

        if (method >= 1)
        {
            double x = rnd.Rannyu(0., radius);
            double y = rnd.Rannyu(0., radius);

            City newCity = City(i, x, y);
            m_cities.push_back(newCity);
        }

    }
}

void Chromosome::printX()
{
    for(int i = 0; i < m_ncities; i++)
    {
        std::cout << m_cities.at(i).getX() << std::endl;
    }
    std::cout << "N_cities = " << m_ncities << std::endl;
    std::cout << "Cpacity : " << m_cities.capacity() << std::endl;
}

double Chromosome::length2()
{
    double length2 = 0;
    for (int i = 0; i < m_ncities - 2; i++)
    {
        double dist2 = m_cities.at(i).dist2(m_cities.at(i + 1));
        length2 += dist2;
    }
    length2 += m_cities.at(m_ncities - 1).dist2(m_cities.at(0));
    return length2;
}

double Chromosome::length()
{
    double length = 0;
    for (int i = 0; i < m_ncities - 2; i++)
    {
        double dist = m_cities.at(i).dist(m_cities.at(i + 1));
        length += dist;
    }
    length += m_cities.at(m_ncities - 1).dist(m_cities.at(0));
    return length;
}