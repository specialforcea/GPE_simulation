%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-2;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-10};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -15;
xmax = 15;
ymin = -15;
ymax = 15;
Nx = 2^6+1;
Ny = 2^6+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Beta = 100;
Omega = 0.6;
Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Potential_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);
Physics2D = Gradientx_Var2d(Method, Physics2D);
Physics2D = Gradienty_Var2d(Method, Physics2D);

%% Setting the initial data
InitialData_Choice = 2;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);

%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 50;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);