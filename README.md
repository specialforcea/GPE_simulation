# GPE_simulation
A Matlab library for numerically solve **Gross-Pitaevskii Equation** in 1D, 2D and 3D with *Time-splitting spectral method*, 
evolution operator is apporximated to first order of t.
Main functionalities:
1. Solve for GPE ground state under user specified potential.
2. Under user specified potential, evolve initial state for time t and return evolved state.

## User guide
scale_parameters.m is the file for user to specify the settings of simulation. Spatial size, mesh number, external potential, 
coupling strength, scattering length can all be changed here. With customized settings:
1. specify stop_crit in Groundstate.m and run it, it computes groundstate with imaginary time evolution method.
2. with specified initial state, evolve_time t, function *dynamic* return the state evolved under GPE for time t.
