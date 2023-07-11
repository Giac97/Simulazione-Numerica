#ifndef Exercise_8_h
#define Exercise_8_h

#include "random.h"
#include <string>

int seed[4];
Random rnd;

double sigmaNew, muNew;

int nPoints = 100000;
int nBlocks = 100;
int stepsPerBlock = (int)(nPoints / nBlocks);
int tSteps = 50;
int accepted;
int attempted;


bool simAnneal = false;


double beta = 1.0;
double deltaBeta = 0.05;
double mu;
double sigma;
double deltaMu;
double deltaSigma;
double dMuMax = 1;
double dSigmaMax = 1;
double errMu;
double errSigma;
double x0;

double Energy;
double stepMax = 2.0;

double oldEnergy;
double newEnergy;
double deltaEnergy;
double errEnergy;

int step = 0;

double progrAvgEnergy;
double progrAvgEnergySq;

std::ofstream outEnergy;
std::ofstream outPosition;
std::ofstream annealOutput;
std::ofstream sampleOutput;

// Functions
double waveFunction(double x, double sigma, double mu);
double waveFunction2(double x, double sigma, double mu);
void readInput(std::string fileName);
double Potential(double x);
bool Accept(double xOld, double xNew, double sigma, double mu);
void Init(void);
void sampleEnergy();
void variationalAttempt();
double Error(double avg1, double avg2, int n);
#endif 