%%% This file is an example of how to use GPELab

%% Ground state of a system of Gross-Pitaevskii equations (see W. Bao and Y. Cai, 
%% 'Ground States of Two-component Bose-Einstein Condensates with an Internal Atmoic Josephson Junction', East Asian Journal on Applied Mathmatics, Vol. 1, No. 1, pp. 49-81)


clear all;
%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP';
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'Energy',1e-12};
Max_Iter = 1e6;
Preconditioner = 'FLaplace';
Method = Method_Var1d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit, Max_Iter, Preconditioner);
xmin = -16;
xmax = 16;
Nx = 2^10+1;
Geometry1D = Geometry1D_Var1d(xmin,xmax,Nx);

%% Setting the physical problem
Delta = 0.5;
Rabi_frequency = -1;
Detuning_constant = 0;
Beta = 600;
Beta11 = 1;
Beta12 = 0.94;
Beta22 = 0.97;
Physics1D = Physics1D_Var1d(Method, Delta, Beta);
Physics1D = Dispersion_Var1d(Method, Physics1D);
Physics1D = Potential_Var1d(Method, Physics1D, quadratic_potential_Josephson1d(Detuning_constant, Rabi_frequency));
Physics1D = Nonlinearity_Var1d(Method, Physics1D,Cubic_Josephson1d(Beta11,Beta12,Beta22),[],Cubic_Josephson_energy1d(Beta11,Beta12,Beta22));

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