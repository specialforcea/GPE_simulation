%%% This file is an example of how to use GPELab

%% Ground state of Gross-Pitaevskii equation with optical potential (see W. Bao, I-L. Chern, F.Y. Lim, 
%% 'Efficient and spectrally accurate numerical methods for computing ground and first excited states in Bose-Einstein condensates', Journal of Computational Physics, Vol. 219, No. 2, pp. 836-854 (2006))


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 5e-2;
Stop_time = [];
Stop_crit = {'Energy',1e-12};
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -16;
xmax = 16;
Nx = 2^10+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 0.5;
Beta = 250;
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D, @(x) x.^2/2+25*sin(pi*x/4).^2);
Physics1D = Nonlinearity_Var1d(Method, Physics1D);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_Choice);

%% Setting informations and outputs
Outputs = OutputsINI_Var1d(Method);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);