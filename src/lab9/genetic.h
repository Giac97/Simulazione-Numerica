#pragma once
#include <vector>

/**
 * @brief Data structure for a city containing the main informations about a single city.
 * 
 * The values contained in this class are the x and y coordinates of the city and the id.
 * The id is an integer value used to identify the city for a better storage and presentation espacially when dealing with 
 * multiple cities
*/
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


/**
 * @brief A list of cities with functionalities to compute the cost functions
 * 
 * The Cities is defined as a list of ncities, it contains functionalities which calculate the cost functions L^(1) and L^(2)
 * 
*/
class Cities
{
    private:
        int m_ncities;
        std::vector<City> m_cities;
    public:
        Cities(int ncities, int method);
        void initCities(int ncities, int method);
        void printX();

        double length2();
        double length();

        std::vector<City> getCities();
};

class Population
{
    private:
        int m_n;
        Cities m_cities;
        
    
    public:
        Population(int n_cities, int n, int method);


};