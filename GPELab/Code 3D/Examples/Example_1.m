%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 1;
Type = 'BESP';
Deltat = 1e-3;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-6};
Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
zmin = -10;
zmax = 10;
Nx = 2^6+1;
Ny = 2^6+1;
Nz = 2^6+1;
Geometry3D = Geometry3D_Var3d(xmin,xmax,ymin,ymax,zmin,zmax,Nx,Ny,Nz);

%% Setting the physical problem
Delta = 0.5;
Beta = 500;
Omega = 0.7;
Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
Physics3D = Dispersion_Var3d(Method, Physics3D);
Physics3D = Potential_Var3d(Method, Physics3D);
Physics3D = Nonlinearity_Var3d(Method, Physics3D);
Physics3D = Gradientx_Var3d(Method, Physics3D);
Physics3D = Gradienty_Var3d(Method, Physics3D);

%% Setting the initial data
InitialData_Choice = 2;
Phi_0 = InitialData_Var3d(Method, Geometry3D, Physics3D, InitialData_Choice);

%% Setting informations and outputs
Outputs = OutputsINI_Var3d(Method);
Printing = 1;
Evo = 15;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab3d(Phi_0,Method,Geometry3D,Physics3D,Outputs,[],Print);