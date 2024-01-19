# Genetic Algorithm
The present code evaluates the best path connecting a number of cities using a genetic algorithm.

## Input File
The options for the simulation are given in an input file in which  the following parameters are specified in the following order:

Number of cities: an integer giving the number of cities, must be greater than 2
Number of generations: How many generations to simulate
Mutation rate: probability of randomly mutating a path
City generation: How cities are generated: 
- 1 on a circle
- 2 inside a square
- 3 American capitals
- 4 Italian capoluoghi
For methods 3 and 4 the files "American_capitals.dat" and "capoluoghi.txt" are required to be in the folder, the correct number of cities also has to be specified manually, else the code will only read the first few cities or gnerate an error if the number given is larger than the number of cities
Output the best path every N generations

## Outputs
The code will output two files: best_paths.txt containing the best path of the generation every N steps as assigned in the input file and distance.out containing the shortest distance of each generation
