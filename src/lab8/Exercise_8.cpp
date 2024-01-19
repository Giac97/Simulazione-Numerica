#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <string>

#include <iomanip>
#include "Exercise_8.h"

#include <vector>



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

    std::cout << std::endl;
    std::cout << "#=================ilBeca==================#" << std::endl;
    std::cout << "#                       ____              #" << std::endl;
    std::cout << "#                     ,'  , `.  ,----..   #" << std::endl;
    std::cout << "#        ,---.     ,-+-,.' _ | /   /   \\  #" << std::endl;
    std::cout << "#       /__./|  ,-+-. ;   , |||   :     : #" << std::endl;
    std::cout << "#  ,---.;  ; | ,--.'|'   |  ;|.   |  ;. / #" << std::endl;
    std::cout << "# /___/ \\  | ||   |  ,', |  ':.   ; /--`  #" << std::endl;
    std::cout << "# \\   ;  \\ ' ||   | /  | |  ||;   | ;     #" << std::endl;
    std::cout << "#  \\   \\  \\: |'   | :  | :  |,|   : |     #" << std::endl;
    std::cout << "#   ;   \\  ' .;   . |  ; |--' .   | '___  #" << std::endl;
    std::cout << "#    \\   \\   '|   : |  | ,    '   ; : .'| #" << std::endl;
    std::cout << "#     \\   `  ;|   : '  |/     '   | '/  : # " << std::endl;
    std::cout << "#      :   \\ |;   | |`-'      |   :    /  #" << std::endl;
    std::cout << "#       '---\" |   ;/           \\   \\ .'   #" << std::endl;
    std::cout << "#             '---'             `---`     #" << std::endl;
    std::cout << "#                                         #" << std::endl;
    std::cout << "#=========================================#" << std::endl;
    std::cout << std::endl << std::endl;
}


/**
Returns the value of the trial wavefunction at x with variational parameters mu, sigma
*/
double waveFunction(double x, double sigma, double mu)
{
    return exp(-( x - mu ) * ( x - mu ) / (2. * sigma * sigma)) + exp(-( x + mu ) *( x + mu ) / (2. * sigma * sigma));
}

double waveFunction2(double x, double sigma, double mu)
{
    return waveFunction( x,  sigma,  mu) * waveFunction( x,  sigma,  mu);
}

/**
Returns the value of the potential at x
*/
double Potential(double x)
{
    return 0.5 * x * x;
    //return ( x * x - 2.5 ) * x * x;
}

/**
 * @brief Uses the Montecarlo algorithm to accept or reject a move
 * 
 * @param xOld Original position
 * @param xNew Proposed position
 * @param sigma Value of the sigma parmeter of the wavefunction
 * @param mu Value of the mu parmeter of the wavefunction
 * @return true accept the move
 * @return false reject the move
 */
bool acceptMove(double xNew)
{
    bool accMove = false;
    double psiOld = waveFunction2(x0, sigma, mu);
    double psiNew = waveFunction2(xNew, sigma, mu);



    double acc = std::min(1., psiNew / psiOld);

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

/**
 * @brief Updates the values of the coordinate using the Metropolis algorithm, proposes a move randomly between -stepMax and +stepMax and accepts it or rejects it
 * IF the move is accepted the position is updated and the accepted counter is increased by 1
 * 
 */
void metroMove()
{
    double delta =  rnd.Rannyu(-stepMax, stepMax);
    double xTrial = x0 + delta;
    double p = waveFunction2(xTrial, sigma, mu) / waveFunction2(x0, sigma, mu);
    double acc = std::min(1.0, p);
    
    
    if (acceptMove(xTrial))
    {
        x0 = xTrial;
        accepted++;
    }
    attempted++;

}

/**
 * @brief Performs a nEquil number of Metropolis steps to find a better starting coordinate.
 * 
 * @param nEquil 
 */
void Equilibrate(int nEquil)
{
    int accepted = 0;
    int attempted = 0;
    for (int i = 0; i < nEquil; i++)
    {
        metroMove();
    }


}

/**
 * @brief Computes the Hamiltonian applied to the wavefunction given the postion and the variational parameters
 * 
 * @param mu Variational parameter, a double
 * @param sigma Variational parameter, a double
 * @param x Coordinate, a double
 * @return H|Psi> 
 */
double evalHamiltonian(double mu, double sigma, double x){

    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;
    double x2 = x * x;

    double psi = waveFunction(x, sigma, mu);
    double V = Potential(x);



    double psiA = exp(-( x - mu ) * ( x - mu ) / (2. * sigma * sigma));
    double psiB = exp(-( x + mu ) * ( x + mu ) / (2. * sigma * sigma));

    double A = psiA *((x - mu) * (x - mu) /sigma4 - (1. / sigma2)) ;
    double B = psiB *((x + mu) * (x + mu) /sigma4 - (1. / sigma2)) ;
    
    double kinPsi = -0.5 * (A + B);
    return kinPsi / psi + V;


}


int main()
{

    //Initialize the random number generator and print the ascii art 
    Init();

    //read the parameters and variables from the input file
    readInput("input.in");

    //Perform a few steps to find a good starting point
    Equilibrate(50);

    //Reset values of Energy and the attempted and accepted move counters
    Energy = 0;
    attempted = 0;
    accepted = 0;

    //If the simAnneal flag is set to false a sampling of the wavefunction and energy is performed at fixed variational parameters
    if (!simAnneal)
        sampleEnergy();
    //If it is set to True starts a simulated annealing run
    else
    {
        annealOutput.open("annealing.out");
        for (int s = 0; s < 100; s++) 
        {
            for (int i = 0; i < tSteps; i++)
            {
                oldEnergy = progrAvgEnergy / nBlocks;

                variationalAttempt();

                step++;

            }
            annealOutput << std::setw(15) << beta << std::setw(15) <<  oldEnergy << std::setw(15) << errEnergy <<std::setw(15) << sigma << std::setw(15) << mu << std::endl;
            std::cout << "=============================================" << std::endl;
            std::cout << "Variational step n: " << s << " beta = " << beta << std::endl;
            std::cout << "Energy = " << oldEnergy << " +/- " << errEnergy << std::endl;
            std::cout << "sigma = " << sigma << " mu =  " << mu << std::endl;
            std::cout << std::endl;

            beta += deltaBeta;

        }

        
        sampleOutput.open("sigmamu.out");
        sampleOutput << "step" << std::setw(15) << "sigma" << std::setw(15) <<  "mu" << std::endl << std::endl;
        double progrAvgSigma = 0;
        double progrAvgSigma2 = 0;
        double progrAvgMu = 0;
        double progrAvgMu2 = 0;
        int finalSteps = 5000;
        for (int i = 0; i < finalSteps; i++)
        {
            oldEnergy = progrAvgEnergy / nBlocks;
            variationalAttempt();
            progrAvgSigma += sigma;
            progrAvgMu += mu;
            progrAvgSigma2 += sigma * sigma;
            progrAvgMu2 += mu * mu;
            sampleOutput << i << std::setw(15) << sigma << std::setw(15) << mu << std::endl;

        }

        sampleOutput.close();
        progrAvgSigma /= finalSteps;
        progrAvgMu /= finalSteps;
        progrAvgMu2 /= finalSteps;
        progrAvgSigma2 /= finalSteps;

        errMu = Error(progrAvgMu, progrAvgMu2, finalSteps);
        errSigma = Error(progrAvgSigma, progrAvgSigma2, finalSteps);

        



           

        //set sim anneal to false to sample energy and wavefunction close to minimum
        simAnneal = false;
        std::cout << "After finding the variational ground state starting a sampling of the wavefunction and energy with the final mu nad sigma values" << std::endl;
        nPoints *= 10;
        sampleEnergy();

        std::cout << "\n #=================================# \n" << std::endl;
        std::cout << "Finished exploration of values" << std::endl;
        std::cout << "sigma = " << std::setw(25) << progrAvgSigma << " +/- " << errSigma << std::endl;
        std::cout << "mu = " << std::setw(25) << progrAvgMu << " +/- " << errMu << std::endl;
        std::cout << "! Final energy = "<< std::setw(25) << oldEnergy << " +/- " << errEnergy << std::endl;
    }
    annealOutput.close();

    return 0;




    

}

/**
 * @brief Reads in the parameters from the input file
 * 
 * @param fileName the name of the input file
 */
void readInput(std::string fileName)
{
    std::ifstream inputFile;

    inputFile.open(fileName);

    inputFile >> simAnneal;
    inputFile >> nBlocks;
    inputFile >> nPoints;
    inputFile >> stepMax;
    inputFile >> x0;
    inputFile >> mu;
    inputFile >> sigma;
    inputFile >> beta;
    inputFile >> deltaBeta;
    std::cout << "VMC For simple 1d system" << std::endl;
    std::cout << "Reading input file...." << std::endl;
    std::cout << std::endl;
    std::cout << "Number of blocks: " << nBlocks << std::endl;
    std::cout << "Number of points: " << nPoints << std::endl;
    std::cout << "Initial value of sigma: " << sigma << std::endl;
    std::cout << "Initial value of mu: " << mu << std::endl;
    std::cout << "Initial inverse temperature: " << beta << std::endl;
    std::cout << "#######################################" << std::endl;

    inputFile.close();


}

/**
 * @brief Samples the wavefunction and computes <H> and its accuracy at a given value of the variational parameters sigma, mu
 * If not doing a simulated annealing run, the block number, acceptance rate and energy is outputted to terminal every block and the sampled postion, energy, and error estimate is written to two files
 * one for the wavefuntion (position.out) one for the energy (energy.out)
 * 
 */
void sampleEnergy()
{
    //opening the output files for writing results:
    outEnergy.open("energy.out");
    outPosition.open("position.out");

    //reinitializing the values of the progressiive average and errors to 0
    progrAvgEnergy = 0.0;
    progrAvgEnergySq = 0.0;

    //loop over the number of blocks:
    for (int i = 0; i < nBlocks; i++)
    {
        Energy = 0.0;
        attempted  = 0;
        accepted = 0;

        //looping over the steps inside the block
        for (int j = 0; j < stepsPerBlock; j++)
        {
            metroMove();
            Energy += evalHamiltonian(mu, sigma, x0);

            //saving position (only if not performing a sim anneal run)
            if (!simAnneal)
                outPosition << x0 << std::endl;
        }
        //If not performing sim anneal, display the acceptance rate
        if (!simAnneal)
        {
            std::cout << "Block number: " << i << std::endl;
            std::cout << "Acceptance rate = " << (double) accepted / attempted * 100 << " %" << std::endl;
            std::cout << "Energy average in block: E = " << Energy / stepsPerBlock << std::endl;
            std::cout << "#-#-#-#-#-#-#-#-#-#-#-#" << std::endl;
        }

        //Performing average over the block values
        Energy /= stepsPerBlock;
        progrAvgEnergy += Energy;
        progrAvgEnergySq += Energy * Energy;

        //Write average energy and error to file (again if not in a sim anneal run)
        if (!simAnneal)
        {
            if (outEnergy.is_open())
            {
                if (i == 0)
                    outEnergy << i << std::setw(15)  << Energy << std::setw(15) << progrAvgEnergy << std::setw(15) << 0 << std::endl;
                else
                    outEnergy<< i << std::setw(15) << Energy << std::setw(15) << progrAvgEnergy / (i + 1) << std::setw(15) << Error(progrAvgEnergy / (i + 1), progrAvgEnergySq / (i + 1), i) << std::endl;
            }
            else
                std::cerr << "ERROR! Couldn't open energy output file!" << std::endl;
        }

    }
    outEnergy.close();
    outPosition.close();

}

/**
 * @brief Estimate of the error
 * 
 * @param avg1 Average value
 * @param avg2 Squared average value
 * @param n Block number
 * @return error estimate 
 */
double Error(double avg1, double avg2, int n) 
{
   if(n==0){
      return 0;
   }
   else{
      return sqrt((avg2-avg1*avg1)/n);
   }
}

/**
 * @brief Attempts, using the Metropolis algorith, a variation of the variational parameters and updates them if the move is accepted
 * 
 */
void variationalAttempt()
{
    // propose a variation of the parameters proportional to the temperature (higher temperatures would give a bigger possible variation)
    deltaMu = rnd.Rannyu(-dMuMax, dMuMax) / beta;
    deltaSigma = rnd.Rannyu(-dSigmaMax, dSigmaMax) / beta;

    mu += deltaMu;
    sigma += deltaSigma;

    sampleEnergy();
    attempted++;
    newEnergy = progrAvgEnergy / nBlocks;
    //std::cout << newEnergy << std::endl;
    errEnergy = Error(progrAvgEnergy/nBlocks,progrAvgEnergySq/nBlocks, (nBlocks-1));
    deltaEnergy = newEnergy - oldEnergy;
    double acc = std::min(1.0, exp(-deltaEnergy * beta));
    
    if (acc == 1.0)
    {
        oldEnergy = newEnergy;
        accepted++;
    }
    else if(rnd.Rannyu() < acc)
    {
        oldEnergy = newEnergy;
        accepted++;
    }
    else
    {
        mu -= deltaMu;
        sigma -= deltaSigma;
    }
}
