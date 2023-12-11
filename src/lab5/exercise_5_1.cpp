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
#include <cmath>

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

/**
 *Returns the square of the 1s orbital wavefunction given a position in space as an armadillo rowve
 *
 * @param pt Position (x, y, z) in space
 *
 * @return psi^2 (1s)
 */
double psi2_1_0_0(rowvec pt)
{
    rowvec pt2 = square(pt);
    double r2 = sum(pt2);
    double r = sqrt(r2);
    
    double psi = 1. / sqrt(M_PI) * exp(-r);
    return psi * psi;
}


/**
 *Returns the square of the 2p orbital wavefunction given a position in space as an armadillo rowve
 *
 * @param pt Position (x, y, z) in space
 *
 * @return psi^2 (2p)
 */
double psi2_2_1_0(rowvec pt)
{
    rowvec pt2 = square(pt);
    double r2 = sum(pt2);
    double r = sqrt(r2);
    double theta = acos(pt(2) / r);
    double psi = 1. / 8. * sqrt(2. / M_PI) * r * exp(-r / 2.) * cos(theta);
    return psi * psi;
}

/**
 *Returns the square of the 3d orbital wavefunction given a position in space as an armadillo rowve (test)
 *
 * @param pt Position (x, y, z) in space
 *
 * @return psi^2 (3d)
 */
double psi2_3d_z2(rowvec pt)
{
    rowvec pt2 = square(pt);
    double r2 = sum(pt2);
    double r = sqrt(r2);
    double rho = 2. / 3. * r;
    double theta = acos(pt(2) / r);

    double Radial = 1. / (9. )* sqrt(30.) * rho * rho * pow(M_E, -rho * 0.5);
    double Angular = sqrt(5. / 4.) * (3 * pt2(2) - 1.) * 1. / sqrt(4. * M_PI);
    double psi = Radial * Angular;
    double psi2 = psi * psi;
    return psi2;

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
  
    
   int M = 500000;
   int relax = 1000;
   mat points(M + relax, 3, fill::zeros);
   points(0,0) = 4.;
   points(0,1) = 4.;
   
   //Read the step from the command line
   double delta = atof(argv[1]);
   //Set the number of rejected moves to 0
   int rej = 0;
   for (int i = 0; i < M - 1 + relax; i++)
   {
       //Random displacemnts dx, dy, dz
       double dx = rnd.Rannyu(-delta, delta);
       double dy = rnd.Rannyu(-delta, delta);
       double dz = rnd.Rannyu(-delta, delta);
                
       rowvec trial_move = {dx, dy, dz};
    
       rowvec trial = points.row(i) + trial_move;
       // Compute the square of wavefunction at old and new position, change the function to evaluate different orbital
       double psi2Old = psi2_1_0_0(points.row(i));
       double psi2New = psi2_1_0_0(trial);
       double acc = min(1., psi2New / psi2Old);

       if (acc == 1)
       {
           points.row(i + 1) = trial;
       }
       else
       {
           double r = rnd.Rannyu();
           if (r <= acc)
           {
               points.row(i + 1) = trial;
           }
                
           else
           { 
               rej++;

               points.row(i + 1) = points.row(i);
           }
       }
   }
   

   points.save("psi2_1s.dat", raw_ascii);
    
   double rejRate = double(rej) / double(M);
   std::cout << "***************************************" << std::endl;
   std::cout << "" << std::endl;
   std::cout << "Rejection rate: " << rejRate * 100 << "%" << std::endl;
   std::cout << "" << std::endl;	
   std::cout << "***************************************" << std::endl;
	int N = 100;
	int L = M / N;
	
	rowvec RAV1 (N, fill::zeros);
	rowvec RAV2 (N, fill::zeros);

	for (int i = 0; i < N ; i++)
	{
		double sum = 0;
		for (int j = 0; j < L; j++)
		{
			int k = j + i * L + relax;
			rowvec point2_k = square(points.row(k));
            double r2_k = point2_k(0) + point2_k(1) + point2_k(2);
            double r_k = sqrt(r2_k);
			sum += r_k;
		}
		RAV1(i) = sum / L;
		RAV2(i) = RAV1(i) * RAV1(i);
	}
	rowvec sum_prog(N, fill::zeros);
	rowvec sum2_prog(N, fill::zeros);
	rowvec err_prog(N, fill::zeros);
	vec x = regspace(0, N - 1);
	
	for (int i = 0; i < N ; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			sum_prog(i) += RAV1(j);
			sum2_prog(i) += RAV2(j);
		}
		sum_prog(i) /= float (i + 1);
		sum2_prog(i) /= float (i + 1);
		err_prog(i) = error(sum_prog, sum2_prog, i);
	}
    std::ofstream rAvg;
    rAvg.open("rAvg_1s.dat");
	for (int i = 0; i < N ; i++)
		rAvg << x(i) << "\t" << sum_prog(i)  << "\t" << err_prog(i) << std::endl;	
	
   rnd.SaveSeed();
   rAvg.close();
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
