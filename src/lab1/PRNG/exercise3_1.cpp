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

double error(double AV1, double AV2, int n)
{
    if (n == 0)
        return 0.;
    else
    {
        return sqrt((AV2 - AV1 * AV1) /(double) n);
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
   int N = 500;

   //Number of throws:
   int M = 1000000;
    
   //Number of throws per block
   int B = M / N;
   //Details of the system; d: distance between lines, l: length of needle 
   double d = 5.0;
   double l = 3.0;
   double pi = 0.0;
   double pi2 = 0.0;
   double r;
   
   bool cross;
   
   for (int i = 0; i < N; i++)
   {
       double N_hit = 0.;
       for (int j = 0; j < B; j++)
       {
            r = rnd.Rannyu(0., d);
            double theta = rnd.Rannyu(0., M_PI);
            double x = 0.0;
            double y = 0.0;
            double norm = 2.0;
            while (norm > 1.0)
            {
                x = rnd.Rannyu(-1.0,1.0);
                y = rnd.Rannyu(-1.0,1.0);
                norm = sqrt(x * x + y * y);
            }
            cross = intersect_no_pi(d, l, r, x, y);
            if (cross == true)
                N_hit += 1.;
       }

       double pi_est = 2. * l * (double)B / (N_hit * d);
       
       pi += pi_est;
       pi2 += pi_est * pi_est;

       std::cout << i << "\t" << pi / double(i + 1) << "\t" << error(pi / double(i+1), pi2 / double(i+1), i) << std::endl;




   }

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
