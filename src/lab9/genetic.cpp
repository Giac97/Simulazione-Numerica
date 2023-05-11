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

/**
 * Returns the squared distance between the city and another city using the formula d2 = dx^2 + dy^2.
 * Mainly used in the evaluation of the cost function L^(2).
 * 
 * @param other Other city from which the distance squared is to be calculated
 * 
 * @returns dist2 The squared distance between the two cities
*/
double City::dist2(City other)
{
    double dx = this->m_x - other.m_x;
    double dy = this->m_y - other.m_y;

    double d2 = dx * dx + dy * dy;

    return d2;

}

/**
 * Returns the distance between the city and another city using the formula d = sqrt(dx^2 + dy^2).
 * Mainly used in the evaluation of the cost function L^(1).
 * 
 * @param other Other city from which the distance is to be calculated
 * 
 * @returns dist The squared distance between the two cities
*/
double City::dist(City other)
{
    double d2 = dist2(other);
    return sqrt(d2);
}

/**
 * Utility function used to get the value of the private variable storing the x coordinate of the city
 * 
 * @returns The x coordinate of the city
*/
double City::getX()
{
    return m_x;
}

Cities::Cities(int ncities, int method)
{
    initCities(ncities, method);
}

void Cities::initCities(int ncities, int method)
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

void Cities::printX()
{
    for(int i = 0; i < m_ncities; i++)
    {
        std::cout << m_cities.at(i).getX() << std::endl;
    }
    std::cout << "N_cities = " << m_ncities << std::endl;
    std::cout << "Cpacity : " << m_cities.capacity() << std::endl;
}

double Cities::length2()
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

/**
 * @brief Returns the L^(1) cost function, defined as the sum of the distance between each city in order plus the distance between the last and the first city
 * 
 * @return double The L^(1) cost function
 */
double Cities::length()
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

/**
 * @brief Returns the vector of cities to be used by other functions and classes
 * 
 * @return std::vector<City> The vector of cities
 */
std::vector<City> Cities::getCities()
{
    return m_cities;
}

