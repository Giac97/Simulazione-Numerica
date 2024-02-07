
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include<armadillo>

using namespace arma;
using namespace std;
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
    
   int M = 100;
   int L = 1e4;
   ofstream chiOut;
   chiOut.open("chiout.dat");
   for (int i = 0; i < M; i++)
   {
       double chi = 0.0;
       rowvec n_i = rowvec(M, fill::zeros);

       for (int j = 0; j < L; j++)
       {
           double r = rnd.Rannyu();
           
           for (int k = 0; k < M; k++)
           {
               if (r >= 1.0 / M * k && r < 1.0 / M * (k + 1))
                   n_i(k)++;
           }


       }

       for ( int k = 0; k < M; k++)
           chi += (n_i(k) - L / M) * (n_i(k) - L / M) / (L / M);

       chiOut << chi << std::endl;
   }

   return 0;
}

   

