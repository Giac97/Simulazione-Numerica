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
#include "utils.h"
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



double integrand(double x)
{
    return M_PI / 2. * cos( M_PI * x / 2.);
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

    
   int N_walkers = 10000;

   //TOtal number of steps
   int M = 100;

   //Number of blocks
   int N = 100;

   //Number of throws per block
   int L = M / N_walkers;
    
    
   mat r(M, 3, fill::zeros);

   rowvec I_mc(N, fill::zeros);
   rowvec I_mc2(N, fill::zeros);

   for (int i = 0; i<M-1; i++)
   {
       int dx =(int) rnd.Rannyu(0.,6);
       
       rowvec mv(3, fill::zeros);
       switch (dx) {
           case 0:
                mv = {1, 0, 0};
                break;
            case 1:
                mv = {-1, 0, 0};
                break;
           case 2:
                mv = {0, 1, 0};
                break;
            case 3:
                mv = {0, -1, 0};
                break;
            case 4:
                mv = {0, 0, 1};
                break;
            case 5:
                mv = {0, 0, -1};
                break;
        }
       r.row(i+1) = r.row(i) + mv;
    }

    r.save("walk.dat", raw_ascii);
    
   for (int i = 0; i < N; i++)
   {
       double SUM = 0.;
       for (int j = 0; j<L; j++)
       {
           int k = j + i * L;
           SUM += r(k);
       }
       I_mc(i) = SUM / L;
       I_mc2(i) = I_mc(i) * I_mc(i);
   }

/**    
    
   rowvec I_prog(N, fill::zeros);
   rowvec I_prog2(N, fill::zeros);
   rowvec I_err(N, fill::zeros); 

   for (int i = 0; i<N; i++)
   {
       for (int j = 0; j <= i; j++)
       {
           I_prog(i) += I_mc(j);
           I_prog2(i) += I_mc2(j);
       }
       I_prog(i) /= float(i + 1);
       I_prog2(i) /= float(i + 1);
       I_err(i) = error(I_prog, I_prog2, i);
   }

   colvec Icol = I_prog.st();
   colvec errcol = I_err.st();
   vec x = regspace(0, N - 1); 
   
   mat out(N,3, fill::zeros);
   out.col(0) = x;
   out.col(1) = Icol;
   out.col(2) = errcol;

	out.save("I_mc_uniform.dat", raw_ascii);
*/
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
