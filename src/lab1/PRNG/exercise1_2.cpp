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
   
   

	int M = 100000; //Number of throws
	int N = 100;		//Number of Blocks
	int L = M / N; 	//Number of throws per block
	
	rowvec r(M, fill::zeros);
	rowvec AV1 (N, fill::zeros);
	rowvec AV2 (N, fill::zeros);
	for (int i = 0; i<M; i++)
	{
		r(i) = rnd.Rannyu();
	}
	
	for (int i = 0; i<= N - 1; i++)
	{
		double SUM = 0.;
		for (int j = 0; j < L; j++)
		{
			int k = j + i * L;
			SUM += (r(k) - 0.5) * (r(k) - 0.5);
		}
		AV1(i) = SUM / L;
		AV2(i) = AV1(i) * AV1(i);
	}
	
	rowvec sum_prog(N, fill::zeros);
	rowvec sum2_prog(N, fill::zeros);
	rowvec err_prog(N, fill::zeros);
	vec x = regspace(0, N - 1);
	
	for (int i = 0; i < N ; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			sum_prog(i) += AV1(j);
			sum2_prog(i) += AV2(j);
		}
		sum_prog(i) /= float (i + 1);
		sum2_prog(i) /= float (i + 1);
		err_prog(i) = error(sum_prog, sum2_prog, i);
	}
	
	for (int i = 0; i < N ; i++)
		std::cout << x(i) << "\t" << sum_prog(i) - 1./12.  << "\t" << err_prog(i) << std::endl;

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
