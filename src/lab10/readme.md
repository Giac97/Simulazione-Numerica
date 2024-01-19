# Genetic Algorithm
The present code evaluates the best path connecting a number of cities using a genetic algorithm.

## Input File
The options for the simulation are given in an input file in which  the following parameters are specified in the following order:

Number of cities: an integer giving the number of cities, must be greater than 2
Number of generations: How many generations to simulate
Mutation rate: probability of randomly mutating a path
City generation: How cities are generated: 1 on a circle, 2 inside a square
Output the best path every N generations

## Outputs
The code will output two files: best_paths.txt containing the best path of the generation every N steps as assigned in the input file and distance.out containing the shortest distance of each generation
