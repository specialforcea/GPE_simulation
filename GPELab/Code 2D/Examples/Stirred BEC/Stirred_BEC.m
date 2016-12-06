%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-3;
Stop_time = 10;
Stop_crit = {'MaxNorm',1e-7};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;
Nx = 2^9+1;
Ny = 2^9+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Beta = 15000;
Omega = 0;
mu = 0.74;
X_m = @(t) 5*cos(mu*t)*(1-sin(mu*t));
Y_m = @(t) 5*cos(mu*t)*sin(mu*t);
V0 = 100;
d = 0.3;
TimePotential = @(t,X,Y) V0*exp(-((X-X_m(t)).^2+(Y-Y_m(t)).^2)/d^2);
Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Potential_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);
Physics2D = TimePotential_Var2d(Method, Physics2D, TimePotential);

%% Setting the initial data
load Phi_ini

%% Setting informations and outputs
Save_solution = 1;
Save_iter = 10;
Outputs = OutputsINI_Var2d(Method,Save_iter,Save_solution);
Printing = 1;
Evo = 50;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi,Method,Geometry2D,Physics2D,Outputs,[],Print);