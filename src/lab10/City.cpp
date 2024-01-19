// City.cpp

#include "City.h"
#include <cmath>
#include <vector>

City::City(double x, double y) : x_(x), y_(y) {}

/**
 * Calculates the Euclidean distance between this city and another city.
 *
 * @param other The other city to calculate the distance to.
 * @return The Euclidean distance between the two cities.
 */
double City::calculateDistance(const City& other) const {
    double xDistance = std::abs(x_ - other.x_);
    double yDistance = std::abs(y_ - other.y_);
    return std::sqrt((xDistance * xDistance) + (yDistance * yDistance));
}

/**
 * Overloaded operator to allow printing a City object to the output stream.
 *
 * @param os The output stream.
 * @param city The City object to be printed.
 * @return The modified output stream.
 */
std::ostream& operator<<(std::ostream& os, const City& city) {
    os << city.x_ << " " << city.y_;
    return os;
}

/**
 * Calculates the total distance of a given path.
 *
 * @param path The path represented as a vector of city indices.
 * @param cities The vector of City objects representing all the cities.
 * @return The total distance of the path.
 */
double calculatePathDistance(const std::vector<int>& path, const std::vector<City>& cities) {
    double distance = 0.0;
    for (int i = 0; i < path.size() - 1; ++i) {
        int cityIndex1 = path[i];
        int cityIndex2 = path[i + 1];
        distance += cities[cityIndex1].calculateDistance(cities[cityIndex2]);
    }
    int lastCityIndex = path.back();
    distance += cities[lastCityIndex].calculateDistance(cities[path[0]]);  // Connects the last city to the first city
    return distance;
}
