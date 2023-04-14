# Ising Model 1d
## Edits to the input file
A functionality has been added to the input file to allow for different options when generating the initial configuration as a "restart" parameter to be added at the end of the input file. Depending on the value given to this parameter the simulation will use as the starting configuration either a completely random configuration, the final configuration from the last simulation run (usefule when simulating a temperature sweep) or a perfectly ordered configuration with every spin "up".
This parameter is an integer number:
- 0: A random configuration is generated
- 1: The configuration is read from the final configuration of the last run
- >=2: An ordered configuration is generated

## Notes on the folder structure
The main code folder contains all the ouputs from each simulation, the code, the Makefile and some bash scripts to automate simulations, the outputs folder contains two folders, one for the Gibbs and one from the Metropolis simulation. Each will contain individual files containg the value of the internal energy, susceptibility, magnetization and specific heat, together with the uncertainty, at the final block for each temperature.
Note that, as with the other codes, the notebook is stored in the notebook folder, together with a copy of the formatted outputs.
