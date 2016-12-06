
%%% This file is an example of how to use GPELab (FFT version)

%-----------------------------------------------------------
% Setting the data
%-----------------------------------------------------------

%% Setting the method and geometry

Computation = 'Dynamic';
Ncomponents = 1;
Type = 'Splitting';
Deltat = 1e-2;
Stop_time = 1;
Stop_crit = [];
Method = Method_Var3d(Computation,Ncomponents, Type, Deltat, Stop_time,Stop_crit);
xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;
zmin = -2;
zmax = 2;
Nx = 2^7+1;
Ny = 2^7+1;
Nz = 2^7+1;
Geometry3D = Geometry3D_Var3d(xmin,xmax,ymin,ymax,zmin,zmax, Nx, Ny, Nz);

%% Setting the physical problem
Delta = 1;
Beta = 1e-3;
Physics3D = Physics3D_Var3d(Method, Delta, Beta);
Physics3D = Dispersion_Var3d(Method, Physics3D);
Physics3D = Nonlinearity_Var3d(Method, Physics3D);

%% Setting the initial data
A = 0.1;
d = 0.6;
Random_Phase = Stationary_Gaussian_Field3d(Geometry3D,...
@(X,Y,Z) A*exp(-(X.^2 + Y.^2 + Z.^2)/(2*d^2)));
Phi_1{1} = exp(-2i*pi*Random_Phase);

%% Setting informations and outputs
Solution_save = 0;
Outputs_iterations = 10;
Output_function{1} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
Geometry3D.dz*sum(sum(sum(ifftn(FFTX.*fftn(Phi)).*conj(Phi))));
Output_function{2} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
Geometry3D.dz*sum(sum(sum(ifftn(FFTY.*fftn(Phi)).*conj(Phi))));
Output_function{3} = @(Phi,X,Y,Z,FFTX,FFTY,FFTZ) Geometry3D.dx*Geometry3D.dy*...
Geometry3D.dz*sum(sum(sum(ifftn(FFTZ.*fftn(Phi)).*conj(Phi))));
Output_name{1} = 'BEC Momentum X';
Output_name{2} = 'BEC Momentum Y';
Output_name{3} = 'BEC Momentum Z';
Outputs = OutputsINI_Var3d(Method,Outputs_iterations,Solution_save,...
Output_function,Output_name);

Printing = 1;
Evo = 10;
Draw = 1;
Print = Print_Var3d(Printing,Evo,Draw);

View = 3;
Isovalue = 1e-2;
Aspect = 1;
Figure = Figure_Var3d(View,Isovalue,Aspect);

%-----------------------------------------------------------
% Launching simulation
%-----------------------------------------------------------

[Phi,Outputs]= GPELab3d(Phi_1,Method,Geometry3D,Physics3D,Outputs,[],Print,Figure);