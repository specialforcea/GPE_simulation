%%% This file is an example of how to use GPELab (FFT version)

%% Ground state of a Gross-Pitaevskii equation with quadratic potential and cubic nonlinearity in 1D


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'Energy',1e-8};
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -20;
xmax = 20;
Nx = 2^10+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 0.5;
Beta = 100;
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D);
Physics1D = Nonlinearity_Var1d(Method, Physics1D);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_Choice);

%% Setting informations and outputs
save = 0;
Outputs = OutputsINI_Var1d(Method, save);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching computation
%-----------------------------------------------------------

[Phi_1, Outputs] = GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);


%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 1;
Stop_crit = [];
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);

%% Setting the physical problem

Brownian = Brownian_Process2d(Method);

Delta = 0.5;
Beta = 100;
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D);
Physics1D = StochasticPotential_Var1d(Method, Physics1D, @(W,X) (1/2)*X.^2.*W, [], @(t,X) Brownian(t));
Physics1D = Nonlinearity_Var1d(Method, Physics1D);

%% Setting informations and outputs
Solution_save = 0;
Outputs_iterations = 5;
Output_function{1} = @(Phi,X,FFTX) sqrt(Geometry1D.dx)*sqrt(sum(abs(Phi).^2));
Output_name{1} = 'L2-Norm ';
Outputs = OutputsINI_Var1d(Method,Outputs_iterations,Solution_save,Output_function,Output_name);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var1d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching computation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab1d(Phi_1,Method,Geometry1D,Physics1D,Outputs,[],Print);
