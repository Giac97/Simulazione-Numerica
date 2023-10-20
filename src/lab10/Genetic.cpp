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
bool COMMUNICATE;
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
    inputFile >> COMMUNICATE;

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

/**
 * @brief obtains the best path accoridng to the distance
 * 
 * @return std::vector<int> The ordered list of city ids
 */
std::vector<int> Population::getBestPath() const {
    auto minFitness = std::min_element(fitness_.begin(), fitness_.end());
    int index = std::distance(fitness_.begin(), minFitness);
    return population_[index];
}

/**
 * @brief evolves a single step of the genetic algorithms performing crossovers and mutations
 * 
 * 
 * @param cities 
 * @param generation 
 */
void Population::evolve(const std::vector<City>& cities, int generation) {
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
    //outputFile << bestDistance << std::endl;

    // Output the best path to the "best_paths.txt" file
    // if (generation == 1 || generation % EVERY == 0)
    // {
    //     bestPathsFile << "Generation " << generation << std::endl;
    //     bestPathsFile << "Best distance: " << bestDistance << std::endl;
    //     bestPathsFile << "Best path:" << std::endl;
    //     std::vector<int> bestPath = getBestPath();
    //     for (int i = 0; i < bestPath.size(); ++i) {
    //         int cityIndex = bestPath[i];
    //         bestPathsFile << cities[cityIndex] << std::endl;
    //     }
    //     bestPathsFile << std::endl;
    // }
    // // Calculate percent completed
    // double percentCompleted = (static_cast<double>(generation + 1) / NUM_GENERATIONS) * 100.0;

    // // Update the progress bar

    // //left side of progress bar
    // std::string progressBar = "[";

    // //Number of sections of the progress bar
    // int numCompleted = static_cast<int>(percentCompleted / 2.0);
    // //fill up to number of completed sections with "="
    // for (int i = 0; i < numCompleted; ++i) {
    //     progressBar += "=";
    // }
    // //fill the rest with "space"
    // for (int i = numCompleted; i < 50; ++i) {
    //     progressBar += " ";
    // }

    //right side of progress bar
    //progressBar += "]";

    // Print the progress bar
    //std::cout << "\rProgress: " << progressBar << "  " << percentCompleted << " %";
    //std::cout.flush();
}


std::vector<int> Population::generateRandomPermutation(int numCities) {
    std::vector<int> permutation(numCities);
    for (int i = 0; i < numCities; ++i) {
        permutation[i] = i;
    }
    std::random_shuffle(permutation.begin() + 1, permutation.end() - 1);  // Exclude the first and last city
    return permutation;
}

/**
 * @brief Performs an ordered crossovers between to parents selected using the tournament
 * selection, picks two ordered points at random of the two paths and exchanges them ensuring 
 * no repetitions
 * 
 * @param parent1 
 * @param parent2 
 * @return std::vector<int> The offspring from the crossover of the two parents
 */
std::vector<int> Population::crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size = parent1.size();
    std::vector<int> offspring(size, -1);

    // Copy the first city from parent1 to offspring
    offspring[0] = parent1[0];
    offspring[offspring.size() - 1] = parent1[offspring.size() - 1];
    // Perform Ordered Crossover (OX)
    int startPos = static_cast<int>( rnd.Rannyu(0, size-1));
    int endPos = static_cast<int>(rnd.Rannyu(0, size-1));

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

/**
 * @brief Performs, if the condition is verified accoridng to a predefined mutation rate
 * one of the possible mutations to a path
 * 
 * @param path Vector of integers representing city ids
 */
void Population::mutate(std::vector<int>& path) {

    if (rnd.Rannyu() < MUTATION_RATE)
    {
        double r = rnd.Rannyu();

        //perform a standard swap of two cities at index i, j
        if (r < 0.5)
        {
            int j = rand() % (path.size() - 1) + 1;
            int i = rand() % (path.size() - 1) + 1;
            std::swap(path[i], path[j]);
        }

        //invert between index i and j
        if (r < 1 && r > 0.5)
        {
            int j = rand() % (path.size() - 1) + 1;
            int i = rand() % (path.size() - 1) + 1;
            if (j < i)
                std::swap(i, j);
            
            int len = j - 1;

            std::reverse(path.begin() + 1, path.begin() + j);

        }


    }

}

/**
 * @brief Performs a tournament selection to select a random path among different paths
 * 
 * 
 * @param cities 
 * @return std::vector<int> 
 */
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

/**
 * @brief Changes the best path of a population to a given path, used for communication
 * between processors of the best path
 * 
 * @param newBest Vector of city ids to be substituted to the best path
 */
void Population::setBest(std::vector<int> newBest)
{
    auto minFitness = std::min_element(fitness_.begin(), fitness_.end());
    int index = std::distance(fitness_.begin(), minFitness);
    population_[index] = newBest;
}

/**
 * @brief Writes the best distance of a population for a given processor.
 * The distance is appended to a file, so if restarting the distance files should be deleted
 * 
 * 
 * @param rank The processor id
 * @param cities The list of cities ordered according to the best path
 */
void Population::writeBestLengthMPI(int rank, const std::vector<City> &cities)
{
    std::string fname;
    if (CITY_GEN == 4)
        fname = "./out_ita/distance_"+std::to_string(rank)+".out";
    else
        fname = "./out/distance_"+std::to_string(rank)+".out";
    double bestDistance = calculatePathDistance(getBestPath(), cities);
    std::ofstream distanceMP;
    distanceMP.open(fname, std::ios_base::app);
    distanceMP << bestDistance << std::endl;
    distanceMP.close();
}

/**
 * @brief Writes the best path for a given generation to a text file for each processor.
 * The path is appended to the already existing file, so if a new simulation is run the best practice is
 * to remove the path files already present
 * 
 * @param rank The processor id
 * @param cities The list of cities in the path order
 * @param generation The generation number
 */
void Population::writeBestPathMPI(int rank,const std::vector<City> &cities, int generation)
{
    std::ofstream pathMP;
    std::string fname;
    if (CITY_GEN == 4)
        fname = "./out_ita/path/path_"+std::to_string(rank)+".out";
    else
        fname = "./out/path/path_"+std::to_string(rank)+".out";
    double bestDistance = calculatePathDistance(getBestPath(), cities);
    pathMP.open(fname, std::ios_base::app);
    pathMP << "Generation " << generation << std::endl;
    pathMP << "Best distance: " << bestDistance << std::endl;
    pathMP << "Best path:" << std::endl;
    std::vector<int> bestPath = getBestPath();
    for (int i = 0; i < bestPath.size(); ++i) {
        int cityIndex = bestPath[i];
        pathMP << cities[cityIndex] << std::endl;
    }
    pathMP << std::endl;
    
}
/**
 * @brief Method that generates the list of cities according to the specification in the input file 
 * 
 * 
 * @param Method an integer specifying the method used to generate cities: 1 -> Generates cities randomly distributed on a circle. 2 -> Generate cities randomly distributed inside a square. 3 -> Read the coordinates of the cities from a txt file
 *  
 * @return std::vector<City> 
 */
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

    else if (Method == 3)
    {
        std::ifstream cityFile;
        cityFile.open("American_capitals.dat");
        double x, y;
        std::string state, name;
        cityFile >> state >> state >> state >> state; //reads the first line
        for (int i = 0; i < NUM_CITIES; i++)
        {
            
            cityFile >> state >> name >> x >> y;
            std::cout << state << " " << name<< " " << x<< " " << y  << std::endl;
            cities.push_back({x, y});
        }
        cities[cities.size() - 1] = cities[0];

    }

    else if (Method == 4)
    {
        std::ifstream cityFile;
        cityFile.open("capoluoghi.txt");
        double x, y;
        for (int i = 0; i < NUM_CITIES; i++)
        {
            
            cityFile >> x >> y;
            //std::cout << state << " " << name<< " " << x<< " " << y  << std::endl;
            cities.push_back({x, y});
        }
        cities[cities.size() - 1] = cities[0];
    }
    return cities;
}

