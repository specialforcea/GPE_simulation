%%% This file is an example of how to use GPELab (FFT version)

%% GROUND STATE COMPUTATION WITH A ROTATING TERM AND COUPLED NONLINEARITIES


%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP';
Deltat = 5e-2;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-6};
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
Method.Precond = 'FLaplace';
xmin = -4;
xmax = 4;
ymin = -4;
ymax = 4;
Nx = 2^8+1;
Ny = 2^8+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);

%% Setting the physical problem
Delta = 0.5;
Beta = 10;
Beta_coupled = [1,2;2,1];
Omega = 0;
Kappa = 3;
Physics2D = Physics2D_Var2d(Method, Delta, Beta, Omega);
RashbaDispersion{1,1} = @(FFTX,FFTY) Delta*(FFTX.^2+FFTY.^2);
RashbaDispersion{1,2} = @(FFTX,FFTY) Kappa*FFTX - 1i*Kappa*FFTY;
RashbaDispersion{2,1} = @(FFTX,FFTY) Kappa*FFTX + 1i*Kappa*FFTY;
RashbaDispersion{2,2} = @(FFTX,FFTY) Delta*(FFTX.^2+FFTY.^2);
Physics2D = Dispersion_Var2d(Method, Physics2D,RashbaDispersion);
Physics2D = Potential_Var2d(Method, Physics2D);
Physics2D = Nonlinearity_Var2d(Method, Physics2D,Coupled_Cubic2d(Beta_coupled),...
[],Coupled_Cubic_energy2d(Beta_coupled));
Physics2D = Gradientx_Var2d(Method, Physics2D);
Physics2D = Gradienty_Var2d(Method, Physics2D);

%% Setting the initial data
InitialData_Choice = 1;
Phi_0 = InitialData_Var2d(Method, Geometry2D, Physics2D, InitialData_Choice);

%% Setting informations and outputs
Save = 0;
Outputs = OutputsINI_Var2d(Method, Save);
Printing = 1;
Evo = 100;
Draw = 1;
Print = Print_Var2d(Printing,Evo,Draw);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);