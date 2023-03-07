#include "utils.h"
#include <cmath>
#include<armadillo>
using namespace arma;

bool intersect(double d, double l, double r, double theta)
{
    bool cross = false;
    double y_comp = l / 2 * sin(theta);

    if (r + y_comp >= d or r - y_comp <= 0)
        cross = true;

    return cross;

}
