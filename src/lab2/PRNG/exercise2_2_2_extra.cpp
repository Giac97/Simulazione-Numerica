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

double distance2(rowvec x)
{
    return x(0) * x(0) + x(1) * x(1) + x(2) * x(2);
}

int main(int argc, char *argv[])
{

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open())
    {
        Primes >> p1 >> p2;
    }
    else
        cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    ifstream input("seed.in");
    string property;
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else
        cerr << "PROBLEM: Unable to open seed.in" << endl;

    int N_walkers = 10000;

    // Total number of steps
    int M = 100;

    // Number of blocks
    int N = 100;

    // Creating matrix containg for each walker the squared distance from the origing at each step (0 -> M)
    mat d2_walkers(M, N_walkers, fill::zeros);

    // 
    for (int i = 0; i < N_walkers - 1; i++)
    {

        mat r(M, 3, fill::zeros);

        colvec d2(M, fill::zeros);

        // Executing M random steps for a single walker
        for (int j = 0; j < M - 1; j++)
        {
            int dx = (int)rnd.Rannyu(0., 6);
            double theta = rnd.Rannyu(0, M_PI);
            double phi = rnd.Rannyu(0, 2.0 * M_PI);
            // Row vector containg the displacement vector to move from step j to j+1
            rowvec mv(3, fill::zeros);
            double dist = rnd.Cauchy(1.0, 1.0);
            mv = {dist*sin(theta) * cos(phi), dist*sin(theta) * sin(phi), dist*cos(theta)};
            r.row(j + 1) = r.row(j) + mv;
            d2(j + 1) = distance2(r.row(j + 1));
        }
        if (i == 2)
        {
            r.save("levy_continue_points.dat", raw_ascii);   
        }

        d2_walkers.col(i) = d2;
    }

    // Number of blocks
    int Nb = 100;

    // Number of walkers per block
    int L = N_walkers / Nb;

    // Matrices where values of the averaged distance squared and its square will be stored at each step
    mat Avg1(M, Nb, fill::zeros);
    mat Avg2(M, Nb, fill::zeros);

    
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < Nb; j++)
        {
            double SUM1 = 0.;
            for (int k = 0; k < L; k++)
            {
                int w = k + j * L;
                SUM1 += d2_walkers(i, w);
            }
            Avg1(i, j) = SUM1 / L;
            Avg2(i, j) = Avg1(i, j) * Avg1(i, j);
        }
    }

    mat Avg_prog(M, Nb, fill::zeros);
    mat Avg_prog2(M, Nb, fill::zeros);
    mat Avg_err(M, Nb, fill::zeros);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < Nb; j++)
        {
            for (int k = 0; k < j; k++)
            {
                Avg_prog(i, j) += Avg1(i, k);
                Avg_prog2(i, j) += Avg2(i, k);
            }
            Avg_prog(i, j) /= float(j + 1);
            Avg_prog2(i, j) /= float(j + 1);
            Avg_prog(i, j) = sqrt(Avg_prog(i, j));
            Avg_prog2(i, j) = sqrt(Avg_prog2(i, j));
            Avg_err(i, j) = error(Avg_prog.row(i), Avg_prog2.row(i), j);
        }
    }

    colvec Avg_final = Avg_prog.col(Nb - 1);
    colvec Err_final = Avg_err.col(Nb - 1);

    mat out_rw(M, 2);
    out_rw.col(0) = Avg_final;
    out_rw.col(1) = Err_final;

    out_rw.save("levy_continue_distance.dat", raw_ascii);

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
