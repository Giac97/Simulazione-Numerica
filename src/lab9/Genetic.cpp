// GeneticAlgorithm.cpp

#include "Genetic.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ostream>
#include "random.h"
#include "City.h"

Random rnd;
int NUM_CITIES;
int NUM_GENERATIONS;
double MUTATION_RATE;
int CITY_GEN;
int EVERY;
int POPULATION_SIZE = 500;
void readInput(std::string fName)
{


    std::ifstream inputFile;

    inputFile.open(fName);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: input.in" << std::endl;
        
    }
    inputFile >> NUM_CITIES;
    inputFile >> NUM_GENERATIONS;
    inputFile >> MUTATION_RATE;
    inputFile >> CITY_GEN;
    inputFile >> EVERY;

    std::cout << std::endl;
    std::cout << "<==============================>" << std::endl;
    std::cout << "Number of cities: " << NUM_CITIES << std::endl;
    std::cout << "Number of generations: " << NUM_GENERATIONS << std::endl;
    std::cout << "Mutation rate: " << MUTATION_RATE * 100 << " %" << std::endl;
    std::cout << "<==============================>\n" << std::endl;



    inputFile.close();

}
Population::Population(int populationSize, int numCities) : populationSize_(populationSize), numCities_(numCities) {
    population_.resize(populationSize_);
    fitness_.resize(populationSize_);
}

void Population::initialize() {
    for (int i = 0; i < populationSize_; ++i) {
        population_[i] = generateRandomPermutation(numCities_);
    }
}

void Population::evaluateFitness(const std::vector<City>& cities) {
    for (int i = 0; i < populationSize_; ++i) {
        fitness_[i] = calculatePathDistance(population_[i], cities);
    }
}

std::vector<int> Population::getBestPath() const {
    auto minFitness = std::min_element(fitness_.begin(), fitness_.end());
    int index = std::distance(fitness_.begin(), minFitness);
    return population_[index];
}

void Population::evolve(const std::vector<City>& cities, std::ofstream& outputFile, std::ofstream& bestPathsFile, int generation) {
    std::vector<std::vector<int>> newPopulation(populationSize_);

    for (int i = 0; i < populationSize_; ++i) {
        if (rnd.Rannyu() > 0.5){
        std::vector<int> parent1 = tournamentSelection(cities);
        std::vector<int> parent2 = tournamentSelection(cities);
        std::vector<int> offspring = crossover(parent1, parent2);
        mutate(offspring);
        newPopulation[i] = offspring;
        }
        else
            newPopulation[i] = population_[i];
    }

    population_ = newPopulation;

    evaluateFitness(cities);

    // Output the best distance to the "distance.out" file
    double bestDistance = calculatePathDistance(getBestPath(), cities);
    outputFile << bestDistance << std::endl;

    // Output the best path to the "best_paths.txt" file
    if (generation == 1 || generation % EVERY == 0)
    {
        bestPathsFile << "Generation " << generation << std::endl;
        bestPathsFile << "Best distance: " << bestDistance << std::endl;
        bestPathsFile << "Best path:" << std::endl;
        std::vector<int> bestPath = getBestPath();
        for (int i = 0; i < bestPath.size(); ++i) {
            int cityIndex = bestPath[i];
            bestPathsFile << cities[cityIndex] << std::endl;
        }
        bestPathsFile << std::endl;
    }
    // Calculate percent completed
    double percentCompleted = (static_cast<double>(generation + 1) / NUM_GENERATIONS) * 100.0;

    // Update the progress bar

    //left side of progress bar
    std::string progressBar = "[";

    //Number of sections of the progress bar
    int numCompleted = static_cast<int>(percentCompleted / 2.0);
    //fill up to number of completed sections with "="
    for (int i = 0; i < numCompleted; ++i) {
        progressBar += "=";
    }
    //fill the rest with "space"
    for (int i = numCompleted; i < 50; ++i) {
        progressBar += " ";
    }

    //right side of progress bar
    progressBar += "]";

    // Print the progress bar
    std::cout << "\rProgress: " << progressBar << "  " << percentCompleted << " %";
    std::cout.flush();
}


std::vector<int> Population::generateRandomPermutation(int numCities) {
    std::vector<int> permutation(numCities);
    for (int i = 0; i < numCities; ++i) {
        permutation[i] = i;
    }
    std::random_shuffle(permutation.begin() + 1, permutation.end() - 1);  // Exclude the first and last city
    return permutation;
}

std::vector<int> Population::crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size = parent1.size();
    std::vector<int> offspring(size, -1);

    // Copy the first city from parent1 to offspring
    offspring[0] = parent1[0];
    offspring[offspring.size() - 1] = parent1[offspring.size() - 1];
    // Perform Ordered Crossover (OX)
    int startPos = static_cast<int>( rnd.Rannyu(0, size));
    int endPos = static_cast<int>(rnd.Rannyu(0, size));

    if (startPos > endPos) {
        std::swap(startPos, endPos);
    }

    // Copy the selected segment from parent1 to offspring
    for (int i = startPos; i <= endPos; ++i) {
        offspring[i] = parent1[i];
    }

    // Fill the remaining positions in offspring using parent2
    int parent2Index = 0;
    for (int i = 0; i < size; ++i) {
        if (offspring[i] == -1) {
            while (std::find(offspring.begin(), offspring.end(), parent2[parent2Index]) != offspring.end()) {
                ++parent2Index;
            }
            offspring[i] = parent2[parent2Index];
            ++parent2Index;
        }
    }
    offspring[offspring.size() - 1] = parent1[offspring.size() - 1];
    return offspring;
}

void Population::mutate(std::vector<int>& path) {
    for (int i = 1; i < path.size() - 1; ++i) {  // Exclude the first city
        if (rnd.Rannyu() < MUTATION_RATE) {
            int j = rand() % (path.size() - 1) + 1;  // Exclude the first city
            std::swap(path[i], path[j]);
        }
    }
}

std::vector<int> Population::tournamentSelection(const std::vector<City>& cities) {
    int tournamentSize = 5;
    int randomIndex;
    double bestDistance = std::numeric_limits<double>::max();
    std::vector<int> bestPath;

    for (int i = 0; i < tournamentSize; ++i) {
        randomIndex = rand() % populationSize_;
        std::vector<int> currentPath = population_[randomIndex];
        double currentDistance = calculatePathDistance(currentPath, cities);

        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            bestPath = currentPath;
        }
    }

    return bestPath;
}


std::vector<City> cityGeneration(int Method)
{
    std::vector<City> cities;
    if (Method == 1)
    {
        for (int i = 0; i < NUM_CITIES; ++i) {
            double th = rnd.Rannyu(0, 2 * 3.1415926);
            double R = 10;
            double x = R * cos(th);
            double y = R * sin(th);
            cities.push_back({ x, y });
        }
        cities[cities.size() - 1] = cities[0];
    }
   
    else if (Method == 2)
    {
        int L = 10;
        for (int i = 0; i < NUM_CITIES; ++i) {

            double x = rnd.Rannyu(-L, L);
            double y = rnd.Rannyu(-L, L);
            cities.push_back({ x, y });
        }
        cities[cities.size() - 1] = cities[0];

    }
    return cities;
}


