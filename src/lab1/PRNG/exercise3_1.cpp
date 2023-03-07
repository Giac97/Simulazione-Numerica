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
#include "utils.h"
#include "random.h"

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
   
   //Number of blocks:
   int N = 100;

   //Number of throws:
   int M = 10000;
    
   //Number of throws per block
   int B = M / N;
   //Details of the system; d: distance between lines, l: length of needle 
   double d = 1.0;
   double l = 0.4;
    
   //Initiating variables for r: random point [0, d), theta: random angle [0: pi), cross: whether needle crosses upper or lower line
   double r;
   double theta;
   bool cross;
    
   //initiating empty row vector to fill with block estimation of the value of pi
   rowvec pi(N, fill::zeros);
   rowvec pi2(N, fill::zeros);
   for (int i = 0; i < N; i++)
   {
       double N_hit = 0.;
       for (int j = 0; j < B; j++)
       {
            r = rnd.Rannyu(0., d);
            theta = rnd.Rannyu(0., M_PI);
            cross = intersect(d, l, r, theta);
            if (cross == true)
                N_hit += 1.;
       }

       double pi_est = 2. * l * (double)N / (N_hit * d);
        
       pi(i) = pi_est;
       pi2(i) = pi(i) * pi(i);
   }

    rowvec pi_prog(N, fill::zeros);
	rowvec pi2_prog(N, fill::zeros);
	rowvec err_prog(N, fill::zeros);
    
   for (int i = 0; i < N; i++)
   {
        for (int j = 0; j <= i; j++)
        {
            pi_prog(i) += pi(j);
            pi2_prog(i) += pi2(j);
        }

        pi_prog(i) /= float (i+1);
        pi2_prog(i) /= float (i+1);
        err_prog(i) = error(pi_prog, pi2_prog, i);
   }

	vec x = regspace(0, N - 1);
   for (int i = 0; i < N; i++)
    std::cout << x(i) << "\t"<< pi_prog(i) << "\t" << err_prog(i) << std::endl;
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
