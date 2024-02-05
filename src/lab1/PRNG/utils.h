#pragma once
#include<iostream>
#include<fstream>
#include<armadillo>
#include"random.h"
using namespace arma;
bool intersect(double d, double l, double r, double theta);
//bool intersect_no_pi(double d, double l, double r, int idx, Random rnd);
bool intersect_no_pi(double d, double l, double r, double x, double y);
