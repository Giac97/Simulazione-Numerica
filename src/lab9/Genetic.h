// GeneticAlgorithm.h
#pragma once

#include <vector>
#include <fstream>
#include <limits>
#include <string>
#include <cmath>
#include "City.h"
#include "random.h"

extern Random rnd;
// Number of cities
extern int NUM_CITIES;

// Number of individuals in each generation
extern int POPULATION_SIZE;

// Number of generations
extern int NUM_GENERATIONS;

// Mutation rate (probability of swapping two cities)
extern double MUTATION_RATE;

//Method of city genertion/reading, 1 on a circle, 2 inside square
extern int CITY_GEN;

//NUmber of steps between writing path, 0 to not output anything, 
extern int EVERY;


class Population {
public:
    Population(int populationSize, int numCities);
    void initialize();
    void evaluateFitness(const std::vector<City>& cities);
    std::vector<int> getBestPath() const;
    void evolve(const std::vector<City>& cities, std::ofstream& outputFile, std::ofstream& bestPathsFile, int generation);

private:
    int populationSize_;
    int numCities_;
    std::vector<std::vector<int>> population_;
    std::vector<double> fitness_;

    std::vector<int> generateRandomPermutation(int numCities);
    std::vector<int> crossover(const std::vector<int>& parent1, const std::vector<int>& parent2);
    void mutate(std::vector<int>& path);
    std::vector<int> tournamentSelection(const std::vector<City>& cities);
};

void readInput(std::string fName);
std::vector<City> cityGeneration(int Method);
