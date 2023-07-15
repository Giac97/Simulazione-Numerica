#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "City.h"
#include "Genetic.h"
int seed[4];

void Init(void)
{
    int p1, p2;
    std::ifstream  Primes, Seed;
    Primes.open("Primes");
    Primes >> p1 >> p2;
    Primes.close();

    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    Seed.close();
}

int main() {
    Init();
    std::cout << "Starting the Genetic algorithm for the TSP" << std::endl;
    
    readInput("input.in");

    std::cout << std::endl;
    std::cout << "This will take a few minutes...\n\n" << std::endl;
    
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Generate random cities
    
    std::cout << "Cities generation started \n\n" << std::endl;
    std::vector<City> cities = cityGeneration(CITY_GEN);
    std::cout << "Cities generated \n\n" << std::endl;
    Population population(POPULATION_SIZE, NUM_CITIES);
    population.initialize();
    std::cout << "Initial population generated \n\n" << std::endl;
    std::ofstream outputFile("distance.out");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening file: distance.out" << std::endl;
        return 1;
    }

    std::ofstream bestPathsFile("best_paths.txt");
    if (!bestPathsFile.is_open()) {
        std::cerr << "Error opening file: best_paths.txt" << std::endl;
        return 1;
    }

    for (int i = 0; i < NUM_GENERATIONS; ++i) {
        population.evaluateFitness(cities);
        population.evolve(cities, outputFile, bestPathsFile, i + 1);
    }

    outputFile.close();
    bestPathsFile.close();

    std::cout << "Best paths are written to the file: best_paths.txt" << std::endl;

    return 0;
}
