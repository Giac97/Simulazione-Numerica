#pragma once

#include <iostream>
#include <vector>
class City {
public:
    City(double x = 0, double y = 0);

    double calculateDistance(const City& other) const;

    friend std::ostream& operator<<(std::ostream& os, const City& city);

private:
    double x_;
    double y_;
};


double calculatePathDistance(const std::vector<int>& path, const std::vector<City>& cities);
