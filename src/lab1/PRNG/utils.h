#pragma once
#include<iostream>
#include<fstream>
#include<armadillo>
using namespace arma;
double error(rowvec AV1, rowvec AV2,double n);
bool intersect(double d, double l, double r, double theta);
