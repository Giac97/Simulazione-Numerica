#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Exercise_8.h"
#include <armadillo>
using namespace arma;
int main()
{

    Init();

    double xOld, xNew;
    stepMax = 1.0;
    int accepted = 0;
    double sigma = 0.5;
    double mu = 1.5;
    Npoints = 100000;
    colvec points(Npoints, fill::zeros);
    for (int i = 0; i < Npoints - 1; i++)
    {

        double dx = rnd.Rannyu(-stepMax, stepMax);
        
        xOld = points(i);
        xNew = xOld + dx;

        if (Accept(xOld, xNew, sigma, mu))  
        {
            points(i + 1) = xNew;
            accepted++;
        }
            
        else
            points(i + 1) = xOld;

    }
    points.save("points.dat", raw_ascii);
    std::cout << (double)accepted / (double)Npoints << std::endl;
    return 0;   
}


void Init(void)
{
    int p1, p2;
    std::ifstream  Primes, Seed;
    Primes.open("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    Seed.close();
}


/**
Returns the value of the trial wavefunction at x with variational parameters mu, sigma
*/
double PsiTrial(double x, double sigma, double mu)
{
    return exp(-( x - mu ) * ( x - mu ) / (2. * sigma * sigma)) + exp(-( x + mu ) *( x + mu ) / (2. * sigma * sigma));
}

/**
Returns the value of the potential at x
*/
double Potential(double x)
{
    return ( x * x - 2.5 ) * x * x;
}

bool Accept(double xOld, double xNew, double sigma, double mu)
{
    bool accMove = false;
    double psiOld = PsiTrial(xOld, sigma, mu);
    double psiNew = PsiTrial(xNew, sigma, mu);

    double psiOldSq = psiOld * psiOld;
    double psiNewSq = psiNew * psiNew;

    double acc = std::min(1., psiNewSq / psiOldSq);

    if (acc == 1.)
        accMove = true;
    else
    {
        double r = rnd.Rannyu();
        if (r < acc)
            accMove = true;
    }
    return accMove;

}
