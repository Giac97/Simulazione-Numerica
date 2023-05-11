#pragma once
#include <vector>

class City
{
    private:
        int m_id = 0;
        double m_x = 0.;
        double m_y = 0.;
    
    public:
        City(int id, double x, double y);

        void initCity(int id, double x, double y);

        double getX();

        double dist2(City other);
        double dist(City other);
        
};



class Chromosome
{
    private:
        int m_ncities;
        std::vector<City> m_cities;
    public:
        Chromosome(int ncities, int method);
        void initChromosome(int ncities, int method);
        void printX();

        double length2();
        double length();
};