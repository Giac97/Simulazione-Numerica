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
#include <armadillo>

using namespace std;
using namespace arma;

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
   
   int M = 10000;
    
   int N = 100;
    
   mat expSample(M, N, fill::zeros);
   mat lorentzSample(M, N, fill::zeros);
   mat uniformSample(M, N, fill::zeros); 

   for (int i = 0; i < M; i++)
   {
       for (int j = 0; j < N; j++)
       {
           expSample(i, j) = rnd.Exp(1.);
           lorentzSample(i, j) = rnd.Lorentz(1., 0.);
           uniformSample(i, j) = rnd.Rannyu();
       }
   }

   rowvec sumExp = sum(expSample, 0);
   sumExp.save("sumexp.dat", raw_ascii);

    rowvec sumUnif = sum(uniformSample, 0);
    sumUnif.save("sumunif.dat", raw_ascii);
    
    rowvec sumLorentz = sum(lorentzSample, 0);
    sumLorentz.save("sumlorentz.dat", raw_ascii);


   
   expSample.save("exp.dat", raw_ascii);
    
   lorentzSample.save("lorentz.dat", raw_ascii);
    
    



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
