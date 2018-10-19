# GPE_simulation
A Matlab library of numerical functions for **Gross-Pitaevskii Equation** in 1D and 2D. 

Library mainly consists of 3 parts, 1d_scalar, 1d_spinor and 2d_spinor. 3d_spinor code is there but never really tested. Speckle analysis and simulation is an application of 1d_spinor, in which more specific functions for analysis of laser speckle potential is made and tested, GPE evolves under speckle potential with and without spin orbit coupling.

## Basic functions in each part:

1. Solve for GPE ground state under user specified potential with **Imaginary time evolution method**.(Default potential is optical dipole trap with frequency specified as Dip_Freq in scale_parameters.m).
2. Under user specified potential (Default is optical lattice), evolve initial state for time t and return final state. Numerical algorithm is **time splitting spectra method**. 

## Guide and Explanation of basic functions
1. scale_parameters.m is the first file to edit, it specified parameters and conditions of simulation. For example, atom mass is Rb_mass for Rb87, N is atom number, Dip_freq is dipole trap frequency. More explanations are in the file.
2. Once simulation condition is set, you can call Groundstate.m and get phi_0 as ground state. In this file, using Thomas Fermi function as initial state, it evolves with imaginary time, you can specify Deltat(each evolution time step) and Stop_crit, which is the stopping criteria of evolution, when the difference of wavefunction before and after one evolution of Deltat is smaller than stop_crit, indicating the convergence of ground state,  the evolution stops. 
3. To evolve a wavefunction under some condition for some time t, dynamic.m is the function to go to. phi_2= dynamic(phi,t_evo,Deltat,c0,Nx,trap,V,k_scale,f,deltax,deltaf,L,xmin,xmax), takes phi, evolve it for t_evo, each time step is Deltat. other parameters are explained in scale_parameters.m. To customize potential under which the phi evolves, go to dynamic.m in line 15, change potential to the one you want, default is dipole trap plus optical lattice, these two can be controlled with parameters trap and V. This function returns final state.
4. There are other ancillary functions you may not need, just ignore them. Any questions please email: specialforceatp@gmail.com.
