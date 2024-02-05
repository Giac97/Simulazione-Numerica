#include "utils.h"
#include "random.h"
#include <cmath>
#include<armadillo>
using namespace arma;

/**
 * @brief Finds whether a needle dropped at a position r of length l crosses lines separated by a distance d
 *
 * @param d The distance between lines on the floor
 * @param l the length of the needle
 * @param r the point where the needle is dropped
 *
 * @return true if the needle is find to cross the line, else false
 */
bool intersect(double d, double l, double r, double theta)
{
    bool cross = false;
    double y_comp = l / 2 * sin(theta);

    if (r + y_comp >= d or r - y_comp <= 0)
        cross = true;

    return cross;
}

/**
 * Same as the intersect function but doesn't require PI, uses instead a uniform random y component 
 */

bool intersect_no_pi(double d, double l, double r,double x, double y){
    bool cross = false;
    double norm = sqrt(x * x + y * y);     
    double tip = r + l * y / norm;
    if (tip >= d or tip <= 0)
        cross = true;
    
    return cross;

}
