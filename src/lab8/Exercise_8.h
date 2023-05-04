#ifndef Exercise_8_h
#define Exercise_8_h

#include "random.h"
int seed[4];
Random rnd;


int Npoints;
double stepMax;

// Functions
double PsiTrial(double x, double sigma, double mu);
double Potential(double x);
bool Accept(double xOld, double xNew, double sigma, double mu);
void Init(void);
#endif 