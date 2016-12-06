%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Relaxation';
Deltat = 1e-4;
Stop_time = 5;
Stop_crit = [];
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^8+1;
Ny = 2^8+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 1;
Beta = 10;
Physics2D = Physics2D_Var2d(Method, Delta, Beta);
Physics2D = Dispersion_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D);

%% Setting the initial data
Random_phase = Stationary_Gaussian_Field2d(Geometry2D,@(X,Y)exp(-(X.^2+Y.^2)/2));
Random_phase = Random_phase/max(max(Random_phase));
Phi_0{1} = exp(-2i*pi*Random_phase);

%% Setting informations and outputs
Outputs = OutputsINI_Var2d(Method);
Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);