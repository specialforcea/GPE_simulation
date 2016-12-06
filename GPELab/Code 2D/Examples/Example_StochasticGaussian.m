%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM

clear all

%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-7};
Max_iter = 5.5e4;
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit,Max_iter);
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^6+1;
Ny = 2^6+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Beta = 1000;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Potential_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);

%% Setting the initial data
InitialData_Choice = 2;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);

%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 500;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi_1, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);

%% Setting the method
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-3;
Stop_time = 2;
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time);

%% Adding random gaussian potential
X_0 = 0;
Y_0 = 0;
d = 4;
V_0 = 2;
Brownian = Brownian_Process2d(Method);
Physics2D = StochasticPotential_Var2d(Method, Physics2D, @(W,X,Y) V_0*exp(-((X-X_0).^2+(Y-Y_0).^2)/2*d^2).*W, [], @(t,X,Y) Brownian(t));

%% Setting outputs
Save_solution = 1;
Save_evo = 10;
Outputs = OutputsINI_Var2d(Method,Save_evo,Save_solution);
Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);


%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_1,Method,Geometry2D,Physics2D,Outputs,[],Print);