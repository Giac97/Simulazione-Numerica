# Exercise 5
## Compile and run
To compile run "./compile_5-1.sh" this will create an executable named "ex5". The code has to be run providing after the executable name a float value for the montecarlo move step size, for example:

> ex5 0.1

Will run the code with a maximum step size of 0.1

## Wavefunction
The code can sample positions for the 1s, 2p and 3d orbitals of the hydrogen atom, however, in order to change the orbital it is necessary to change the code and re-compile it.
For ease of use two codes are created with their respictive compialtion scripts for the 1s and 2p orbitals:
> exercise_5-1.cpp 

contains the code for the 1s orbital, it produces two files:

> psi2_1s.dat  
> rAvg_1s.dat

the first containing the sampled points, the second the block averaged expectation value of r.

The file
> exercise_5-1_2p.cpp

contains the code for the 2p orbital, produces two files:

> psi2_2p.dat  
> rAvg_2p.dat

the first containing the sampled points, the second the block averaged expectation value of r.
The code will output the sampled points and the average radius together with the block average and the error
