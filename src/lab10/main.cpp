#include <mpi.h>
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

int main(int argc, char* argv[]) {
    Init();
    std::cout << "Starting the Genetic algorithm for the TSP" << std::endl;
    
    readInput("input.in");

    //Initializing MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;


    std::cout << std::endl;
    std::cout << "This will take a few minutes...\n\n" << std::endl;
    
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // Generate cities from input
    
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

    bool Migrate = false;
    int migrateEvery = 20;
    int counter = 0;
    int sender, receiver;
    for (int i = 0; i < NUM_GENERATIONS; ++i) {
        if (counter == migrateEvery)
            Migrate = true;
        if (!Migrate)
        {
            population.evaluateFitness(cities);
            population.evolve(cities, i + 1);
            counter++;
        }
        else if (!COMMUNICATE)
        {
            population.evaluateFitness(cities);
            population.evolve(cities, i + 1);
            counter++;
        }
        else if (Migrate && COMMUNICATE == true)
        {
            if (rank == 0)
            {
                sender = (int)rnd.Rannyu(0, size);
                receiver = (int)rnd.Rannyu(0, size);
                while (sender == receiver)
                {
                    receiver = (int)rnd.Rannyu(0, size);
                }
                std::cout << std::endl << "Migration (giver, receiver): " << ", " << sender << ", " << receiver << std::endl;
            }
            MPI_Bcast(&sender, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
            MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
            int migrator[NUM_CITIES];
            if (rank == sender)
            {
                std::vector<int> select = population.getBestPath();
                
                for (int j = 0; j < NUM_CITIES; j++)
                {
                    migrator[j] = select[j];
                }
                MPI_Send(migrator, NUM_CITIES, MPI_INTEGER, receiver, 1, MPI_COMM_WORLD);
            }

            if (rank == receiver)
            {
                MPI_Recv(migrator, NUM_CITIES, MPI_INTEGER, sender, 1, MPI_COMM_WORLD, &stat);
                std::vector<int> newBest;
                for (int j = 0; j < NUM_CITIES; j++)
                {
                    newBest.push_back(migrator[j]);
                }
                population.setBest(newBest);
            }
            counter = 0;
            Migrate = false;
            population.evaluateFitness(cities);
            population.evolve(cities, i + 1);

        }
        population.writeBestLengthMPI(rank, cities);
        if (i == 1 || i % EVERY == 0)
            population.writeBestPathMPI(rank, cities, i);
    }

    outputFile.close();
    bestPathsFile.close();

    std::cout << "Best paths are written to the file: best_paths.txt" << std::endl;
    MPI_Finalize();
    return 0;
}
