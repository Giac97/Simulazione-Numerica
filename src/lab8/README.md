# VMC
The code requires an input file, in order the required parameters are:

- simAnneal: either a 0 or a 1, determines whether to perform the simulated annealing in order to find the values of sigma or mu (if set to 1) or just sample the wavefunction for the given values of $\sigma$ and $\mu$
- nBlocks: The number of blocks for the blocking averages
- nPoints: Number of points per block
- stepMax: Maximum step in coordinate space
- x0: Initial coordinates of the walker
- mu: The value of $\mu$ (initial value if doing the annealing) 
- sigma: The value of $\sigma$ (initial value if doing the annealing) 
- beta: Initial value of $\beta = T^{-1}$
- deltaBeta: Increase in the inverse temperature per block

