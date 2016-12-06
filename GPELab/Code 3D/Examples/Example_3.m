%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP';
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-4};
Method = Method_Var3d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -8;
xmax = 8;
ymin = -8;
ymax = 8;
zmin = -8;
zmax = 8;
Nx = 2^7+1;
Ny = 2^7+1;
Nz = 2^7+1;
Geometry3D = Geometry3D_Var3d(xmin,xmax,ymin,ymax,zmin,zmax,Nx,Ny,Nz);

%% Setting the physical problem
Delta = 0.5;
Beta = 1000;
Beta_coupled = [1,0.5;0.5,1];
Omega = 0;
Kappa = 2;
RashbaDispersion{1,1} = @(FFTX,FFTY,FFTZ) Delta*(FFTX.^2+FFTY.^2+FFTZ.^2);
RashbaDispersion{1,2} = @(FFTX,FFTY,FFTZ) Kappa*FFTX + 1i*Kappa*FFTY;
RashbaDispersion{2,1} = @(FFTX,FFTY,FFTZ) Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{2,2} = @(FFTX,FFTY,FFTZ) Delta*(FFTX.^2+FFTY.^2+FFTZ.^2);
Physics3D = Physics3D_Var3d(Method, Delta, Beta, Omega);
Physics3D = Dispersion_Var3d(Method, Physics3D,RashbaDispersion);
Physics3D = Potential_Var3d(Method, Physics3D);
Physics3D = Nonlinearity_Var3d(Method, Physics3D,Coupled_Cubic3d(Beta_coupled),[],Coupled_Cubic_energy3d(Beta_coupled));
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