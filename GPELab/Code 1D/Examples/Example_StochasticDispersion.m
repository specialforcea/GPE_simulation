%%% This file is an example of how to use GPELab (FFT version)

%% Ground state of a Gross-Pitaevskii equation with quadratic potential and cubic nonlinearity in 1D


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 5;
Stop_crit = [];
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -25;
xmax = 25;
Nx = 2^12+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 1;
Beta = 1;
Brownian = Brownian_Process2d(Method);
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = StochasticDispersion_Var1d(Method, Physics1D, @(W,FFTX) FFTX.^2.*W, [], @(t,X) Brownian(t));
%Physics1D = Dispersion_Var1d(Method, Physics1D);
%V = ones(1,Nx-1);
%Physics1D = StochasticPotential_Var1d(Method, Physics1D, @(W,X) V*X.*W, [], @(t,X) Brownian(t));
Physics1D = Nonlinearity_Var1d(Method, Physics1D, @(Phi,X) -abs(Phi).^2);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var1d(Method, Geometry1D, Physics1D, InitialData_Choice);

%% Setting informations and outputs
Solution_save = 1;
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

[Phi_1, Outputs] = GPELab1d(Phi_0,Method,Geometry1D,Physics1D,Outputs,[],Print);

Draw_Timesolution1d(Outputs,Method,Geometry1D,Figure_Var1d)