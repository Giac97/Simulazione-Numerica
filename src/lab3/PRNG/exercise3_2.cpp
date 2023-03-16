/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <armadillo>

using namespace std;
using namespace arma;
/**
 * This function evaluates statistical uncertainties when using blocking given the average and the square average.
 *
 * @param AV1 Block average
 * @param AV2 Square of the block average
 * @param n Number of blocks
 *
 * @return error
 */

double error(rowvec AV1, rowvec AV2, int n)
{
    if (n == 0)
        return 0.;
    else
    {
        return sqrt((AV2(n) - AV1(n) * AV1(n)) / n);
    }
}


double priceStep(double sOld, double dt, double mu, double sigma, double Z)
{
    double exponent = (mu - 0.5 * sigma * sigma) * dt + sigma * Z * sqrt(dt);
    return sOld * exp(exponent);
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
   //TOtal number of throws
   int M = 10000;

   //Number of blocks
   int N = 100;

   //Number of throws per block
   int L = M / N;
    
   // Initial price
   double S0 = 100.;

   // Delivery time
   double T = 1.;

   //Strike price
   double K = 100.;
   
   // Risk free interest rate
   double r = 0.1;

   // Voltility
   double sigma = 0.25;

    
   // Number steps
   int nSteps = 100;

   // Delta t
   double dt = T / (double) nSteps;


   rowvec ST(M, fill::zeros); 

   rowvec callPrice(M, fill::zeros);

   rowvec putPrice(M, fill::zeros);
   for (int i = 0; i < M; i++)
   {
        colvec S_t(nSteps, fill::zeros);
        S_t(0) = S0;
        for (int j = 0; j < nSteps - 1; j++)
        {
            double Z = rnd.Gauss(0. , 1.);
            S_t(j+1) = priceStep(S_t(j), dt, r, sigma, Z);
        }    
        
        if (i == 100)
            S_t.save("price_sim.dat", raw_ascii);
       ST(i) = S_t(nSteps - 1);
       callPrice(i) = exp(-r * T) * max(0., ST(i) - K);
       putPrice(i)  = exp(-r * T) * max(0., K - ST(i));
   }

   
    
   rowvec callAv(N, fill::zeros);
   rowvec putAv(N, fill::zeros);
    
   rowvec callAv2(N, fill::zeros);
   rowvec putAv2(N, fill::zeros);

   for (int i = 0; i < N; i++)
   {
        double sumCall = 0;
        double sumPut = 0;
        for (int j = 0; j < L; j++)
        {
            int k = j + i * L;
            sumCall += callPrice(k);
            sumPut += putPrice(k);
        }
        callAv(i) = sumCall / L;
        callAv2(i) = callAv(i) * callAv(i);

        putAv(i) = sumPut / L;
        putAv2(i) = putAv(i) * putAv(i);
   }
    
   rowvec callAvProg(N, fill::zeros);
   rowvec callAvProg2(N, fill::zeros);
   rowvec callErr(N, fill::zeros);

   rowvec putAvProg(N, fill::zeros);
   rowvec putAvProg2(N, fill::zeros);
   rowvec putErr(N, fill::zeros);

   for (int i = 0; i < N; i++)
   {
        for (int j = 0; j < i; j++)
        {
            callAvProg(i) += callAv(j);
            putAvProg(i) += putAv(j);

            callAvProg2(i) += callAv2(j);
            putAvProg2(i) += putAv2(j);
        }
        callAvProg(i) /= float(i + 1);
        callAvProg2(i) /= float(i + 1);
        callErr(i) = error(callAvProg, callAvProg2, i);

        putAvProg(i) /= float(i + 1);
        putAvProg2(i) /= float(i + 1);
        putErr(i) = error(putAvProg, putAvProg2, i);
   }    

   std::cout << "Call price at T: " << callAvProg(N - 1) << " +/- " << callErr(N - 1) << std::endl;
   std::cout << "Put price at T: " << putAvProg(N - 1) << " +/- " << putErr(N - 1) << std::endl;


   colvec callCol = callAvProg.st();
   colvec callErrCol = callErr.st();
   
   colvec putCol = putAvProg.st();
   colvec putErrCol = putErr.st();
   vec x = regspace(0, N - 1); 
   
   mat out(N,5, fill::zeros);
   out.col(0) = x;
   out.col(1) = callCol;
   out.col(2) = callErrCol;    
   out.col(3) = putCol;
   out.col(4) = putErrCol; 
   
   out.save("discrete_pricing.dat", raw_ascii);
   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
